#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:49:56 2024

Add coastlines to spilhaus projection by extracting them from cartopy 
and adding breaks where curves spill over edges of projection


@author: crowan
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cartopy.io.shapereader as shapereader
from shapely import Polygon, LineString
from scipy.spatial import distance_matrix

from spilhaus import from_lonlat_to_spilhaus_xy, line_lonlat_to_spilhaus_xy
    
def line_lonlat_to_pretty_spilhaus_xy(coords, prettypoly, limit=11707219.26):  
    """
    Wrapper to line_lonlat_to_spilhaus_xy function, which rotates and filters/splits
    segments that touch edge of map so that they conist of segments that fall within the polygon that fully matches 
    projection space to the boundaries of the ocean basin.

    Parameters
    ----------
    coords : nummpy array of longitude (range -180 to 180) - latitude (range -90 to 90) pairs
    prettypoly : shapely Polygon definiing plot limits, in Spilhaus x-y pairs
    limit: defines 'edge region' of projection space: 99% of absolute outer value of 11825474

    Returns
    -------
    List of nmupy arrays of line segments that fall within the defined polygon: [Spilhaus x (easting), Spillhuas y (northing)]

    """
    
    segments=line_lonlat_to_spilhaus_xy(coords)
    
    result=[]
    if len(segments)>1:
        for segment in segments: 
        # see if any part of segment falls within confines of "pretty" projection space and adds it to returned result
            check=fit_check(segment, prettypoly) 
            if check is not None: 
                result.append(check)     
            # then rotates segment back to border of projection space it crossed and does the same thing.
            rotated_x=segment[-1][1]
            rotated_y=segment[-1][0]
            dx=np.flip(segment[:,0]-segment[-1][0],0)
            dy=np.flip(segment[:,1]-segment[-1][1],0)
            if segment[-1][0]> limit: # E side - N side rotation
                 rotseg=np.column_stack((rotated_x+dy,rotated_y-dx))
            elif segment[-1][1]> limit: # N side - E side rotation
                 rotseg=np.column_stack((rotated_x-dy,rotated_y+dx))
            elif segment[-1][0] < - limit: # W side - S side rotation
                rotseg=np.column_stack((rotated_x+dy,rotated_y-dx))
            elif segment[-1][1] < - limit: # S side - W side rotation
                rotseg=np.column_stack((rotated_x-dy,rotated_y+dx))
            check=fit_check(rotseg, prettypoly)
            if check is not None: 
                result.append(check)
    else: # don't need to bother with the rotation bit for a single segment - not sure if fit check strictly needed either
        check=fit_check(segments[0], prettypoly) 
        if check is not None: 
            result.append(check)  
    return result
    
def fit_check(coords, polygon):
    """
    Checks if linestring defined by coords lies within the polygon, and returns
    the part of it that does.


    Parameters
    ----------
    coords : nummpy array of Spilhaus easting (x) and northing (y) pairs
    polygon : shapely Polygon definiing plot limits, also in Spilhaus x-y pairs

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    if len(coords)>1:
        if LineString(coords).within(polygon): # whole string within Polygon
            return coords
        elif LineString(coords).crosses(polygon): # part of string within Polygon
            # find bits which fall within polygon 
            inside_bits=LineString(coords).intersection(polygon)
            if inside_bits.geom_type== 'LineString':
                x,y=inside_bits.coords.xy
                return np.column_stack(([a for a in x],[b for b in y]))
            elif inside_bits.geom_type== 'MultiLineString': 
                # not sure if this should really happen if polygon is correctly adjusted...
                # joining segments should just trace edge of polygon and simplifies things downstream.
                result=[]
                for subgeom in inside_bits.geoms:
                    x,y=subgeom.xy
                    result.append(np.column_stack(([a for a in x],[b for b in y])))
                return np.concatenate(result)
            else: return None
        else: return None
    else: # we have a single point
        return coords
    

# prettypolygon.txt contains vertices for polygon tracing out mask for prettified map (in Spilhaus projection coordinates)
polydata=pd.read_csv('prettypolygon.txt', sep='\t')        
prettypoly=Polygon([(x,y) for x,y in zip(polydata['x'], polydata['y'])])

# read coastline data from cartopy - 50m and 10m also available?
coast = shapereader.natural_earth(resolution='50m',
                                  category='physical',
                                  name='coastline')

coastlines = shapereader.Reader(coast).geometries()
coast_latlons=[]

for geom in coastlines:
    if geom.geom_type == 'LineString':
        coast_latlons.append(geom.xy)
    elif geom.geom_type == 'MultiLineString':
        for subgeom in geom.geoms:
            coast_latlons.append(subgeom.xy)

extreme=11825474
limit=0.95*extreme
threshold=3e6

# fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)
# ax.plot(polydata['x'], polydata['y'], color='grey', lw=0.5)

result=[]
for coastline in coast_latlons:
    lon, lat=np.array(coastline[0].tolist()),np.array(coastline[1].tolist())
    spil_x, spil_y=from_lonlat_to_spilhaus_xy(lon,lat)
    # calculates change in xy values for adjacent points in linestring
    diff=pd.DataFrame(np.column_stack(([a-b for a,b in zip(spil_x[1:],spil_x[:-1])], 
                                  [a-b for a,b in zip(spil_y[1:],spil_y[:-1])])),
                  columns=['dx','dy'])
    # finds spillover points where xy changes a lot between adjacent points 
    breaks=diff.index[((diff['dx']>threshold) | (diff['dx']<-threshold)) | ((diff['dy']>threshold) | (diff['dy']<-threshold))]
    # if breaks exist split linestring data at breakpoints
    if len(breaks)>0:
        broken_segments=[]
        p=0
        for b in breaks:
            broken_segments.append(np.column_stack(([x for x in spil_x[p:b+1]],[y for y in spil_y[p:b+1]])))
            p=b+1
        broken_segments.append(np.column_stack(([x for x in spil_x[p:]],[y for y in spil_y[p:]]))) #remember the last segment!
        # where the segment spills over the limits of the coordinate projection, to modify to make boundaries truly fit ocean basins
        # need to rotate the segments into the appropriate position "outside" the projection limits, then find which sections
        # actually fit within "pretty" boundaries.
        for item in broken_segments[:-1]: #last sections don't neccessarily end at barrier!
            # ax.plot(item[:,0], item[:,1], color='green', lw=0.2)
            check=fit_check(item, prettypoly)
            if check is not None: 
                result.append(['Bordering', check])          
            rotated_x=item[-1][1]
            rotated_y=item[-1][0]
            dx=np.flip(item[:,0]-item[-1][0],0)
            dy=np.flip(item[:,1]-item[-1][1],0)
            # checks last 
            if item[-1][0]> limit: # E side - N side rotation
                rotseg=np.column_stack((rotated_x+dy,rotated_y-dx))
            elif item[-1][1]> limit: # N side - E side rotation
                rotseg=np.column_stack((rotated_x-dy,rotated_y+dx))
            elif item[-1][0] < - limit: # W side - S side rotation
                rotseg=np.column_stack((rotated_x+dy,rotated_y-dx))
            elif item[-1][1] < - limit: # S side - W side rotation
                rotseg=np.column_stack((rotated_x-dy,rotated_y+dx)) 
            else: # not quite sure why there are these ones..
                if item[0][0]> limit: # E side - N side rotation
                    rotseg=np.column_stack((rotated_x+dy,rotated_y-dx))
                elif item[0][1]> limit: # N side - E side rotation
                    rotseg=np.column_stack((rotated_x-dy,rotated_y+dx))
                elif item[0][0] < - limit: # W side - S side rotation
                    rotseg=np.column_stack((rotated_x+dy,rotated_y-dx))
                elif item[0][1] < - limit: # S side - W side rotation
                    rotseg=np.column_stack((rotated_x-dy,rotated_y+dx))            
            # ax.plot(rotseg[:,0], rotseg[:,1], color='orange', lw=0.2)
            check=fit_check(rotseg, prettypoly)
            if check is not None: 
                result.append(['Bordering', check])
        # ax.plot(broken_segments[-1][:,0], broken_segments[-1][:,1], color='blue', lw=0.2)
        check=fit_check(broken_segments[-1], prettypoly)
        if check is not None: 
            result.append(['Bordering', check])
    else: # if no breaks, just use whole linestring
       # ax.plot(spil_x,spil_y, color='red', lw=0.2)
       coords=np.column_stack(([x for x in spil_x],[y for y in spil_y]))
       if fit_check(coords, prettypoly) is not None:
           # there are some handfall of segments (mostly small islands) that lie outside polygon.
           # technically these should be rotated too but currently just excluding
           if LineString(coords).is_closed:
               result.append(['Enclosed',coords])
           else:
               result.append(['Open',coords])

coastresult=pd.DataFrame(result, columns=['Type','Coords'])

# All coast strings are now within the pretty polygon, but several types:
# 1. "Closed" are closed linestrings that fall entirely within the pretty polygon
# 2. "Enclosed" are unclosed linestrings that have been split by crossing the border of the polygon
# 3. "Open" are unclosed linestrings that fall entirely within the pretty polygon: Includes a two small islands split into two segements

# This gets us the segments we want to join into a single bordering coastline segment
bordering_segments=coastresult[coastresult.Type == "Bordering"]['Coords'].to_list()
bordering_segments.append(coastresult[coastresult.Type == "Open"].iloc[5]['Coords'])

# Bit of brute forcing in here, but it generates a continuous outer coastline without obvious out of sequence jumping
merged_line = bordering_segments[0][::-1]  # Start with the first segment
remaining_lines = bordering_segments[1:]  # segments to join to first: last 2 are tiny and screw things up
# This works. 
i=0
while remaining_lines:
#while i<23:
    # Calculate distances between ends of merged line and ends of each remaining segment
    distances = np.array([np.concatenate(distance_matrix(
        merged_line[[0, -1]], [line[0], line[-1]])) for line in remaining_lines])
    
    # Find the minimum distance and its flattened index
    min_distance_index = np.argmin(distances)
    
    # Determine the nearest line and the attachment points 
    nearest_line_index = min_distance_index // 4 # Integer division to get the row
    attachment_points = min_distance_index % 4  # Modulo operation to get the column (0, 1, 2, or 3)
    
    print(nearest_line_index, attachment_points, len(remaining_lines[nearest_line_index]))
    
    # Attach the nearest line
    if attachment_points == 0: # start of merged -> start of nearest segment 
        merged_line = np.concatenate([merged_line[::-1], remaining_lines[nearest_line_index]])
    elif attachment_points == 1:  # start of merged -> end of nearest segment
        merged_line = np.concatenate([remaining_lines[nearest_line_index], merged_line])
    elif attachment_points == 2: # end of merged -> start of nearest segment
        merged_line = np.concatenate([merged_line, remaining_lines[nearest_line_index]])
    elif attachment_points == 3:  # end of merged -> start of nearest segment
        merged_line = np.concatenate([merged_line, remaining_lines[nearest_line_index][::-1]]) 

    del remaining_lines[nearest_line_index]
    i=i+1


# Drop last point which does something weird, and close polygon.
merged_line=np.concatenate([merged_line[:-1],[merged_line[0]]])

# Combine merged borderling coastline with closed polygons within the ocean basins
final_types, final_coords=['Bordering'],[merged_line]
final_types.extend(coastresult[coastresult.Type == "Enclosed"]['Type'].to_list())
final_coords.extend(coastresult[coastresult.Type == "Enclosed"]['Coords'].to_list())



# Add a couple of islands that are split in two for some reason.
final_types.extend(['Enclosed','Enclosed'])
final_coords.extend([np.concatenate([[coastresult[coastresult.Type == "Open"].iloc[8]['Coords'][-1]],
                                   coastresult[coastresult.Type == "Open"].iloc[5]['Coords'],
                                   coastresult[coastresult.Type == "Open"].iloc[8]['Coords']]),
                    np.concatenate([[coastresult[coastresult.Type == "Open"].iloc[7]['Coords'][-1]],
                                    coastresult[coastresult.Type == "Open"].iloc[6]['Coords'],
                                    coastresult[coastresult.Type == "Open"].iloc[7]['Coords']])
                    ])


# Write to a file in gmt format
outfile = open('Spilhaus_Coastlines_110m.txt','w')
for type, coords in zip(final_types, final_coords):
    outfile.write("> %s\n" % (type))
    for xy in coords:
        outfile.write("%s\t%s\n" % (xy[0], xy[1]))
outfile.close()


fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)
ax.plot(polydata['x'], polydata['y'], color='grey', lw=0.5)
for i,row in coastresult.iterrows():
     ax.plot(row.Coords[:,0], row.Coords[:,1], color='grey', lw=0.5)
ax.plot(merged_line[:,0], merged_line[:,1], color='red', lw=0.5)



fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)
ax.plot(polydata['x'], polydata['y'], color='grey', lw=0.5)
for coastline in coast_latlons:
    lon, lat=np.array(coastline[0].tolist()),np.array(coastline[1].tolist())
    spil_x, spil_y=from_lonlat_to_spilhaus_xy(lon,lat)
    ax.plot(spil_x, spil_y, color='red', lw=0.5)


test=np.concatenate([[coastresult[coastresult.Type == "Open"].iloc[7]['Coords'][-1]],coastresult[coastresult.Type == "Open"].iloc[6]['Coords'],coastresult[coastresult.Type == "Open"].iloc[7]['Coords']])

fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)
ax.plot(polydata['x'], polydata['y'], color='grey', lw=0.5)
for row in bordering_segments[23:]:
     ax.plot(row[:,0], row[:,1], lw=0.5)



fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)
ax.plot(polydata['x'], polydata['y'], color='grey', lw=0.5)
for i,row in coastresult[coastresult.Type == "Enclosed"].iterrows():
    ax.plot(row.Coords[:,0], row.Coords[:,1], color='blue', lw=0.5)
    
 
fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)
ax.plot(polydata['x'][:5], polydata['y'][:5], color='grey', lw=0.5)
for i,row in coastresult[coastresult.Type == "Bordering"].iterrows():
    ax.plot(row.Coords[:,0], row.Coords[:,1], color='red', lw=0.5)
row=coastresult[coastresult.Type == "Open"].iloc[5]
ax.plot(row.Coords[:,0], row.Coords[:,1], color='blue', lw=0.5)
        
fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)
ax.plot(polydata['x'], polydata['y'], color='grey', lw=0.5)
for i,row in coastresult[coastresult.Type == "Open"].iterrows():
    ax.plot(row.Coords[:,0], row.Coords[:,1], color='red', lw=0.5)