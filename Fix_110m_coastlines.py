#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 21:39:26 2024

Code used to generate Spilhaus_Coastlines_100m.txt
Uses cartopy coastlines - can't just extract and plot as you go because:
    
    1. Line segments spill over edges of projection, creating a mess - need to break at projection boundaries
    2. Boundaries of raw Sphilhaus projection need to be adjusted to create a true '1 ocean basin view'
    
This code performs the neccessary adjustments and exports all coastlines as closed polygons in segmented gmt format. 
The first segment in the file is the 'Bordering' that defines the edges of the world oceans. 
Everything else is 'Enclosed', i.e. plots as a landmass within the ocean basin.

@author: crowan
"""

import numpy as np
import pandas as pd
import cartopy.io.shapereader as shapereader
from shapely import Polygon, LineString
from scipy.spatial import distance_matrix
from spilhaus import from_lonlat_to_spilhaus_xy

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
            x,y=LineString(coords).intersection(polygon).coords.xy
            return np.column_stack(([a for a in x],[b for b in y]))
        else: return None
    else: # we have a single point
        return coords
    

# prettypolygon.txt contains vertices for polygon tracing out mask for prettified map (in Spilhaus projection coordinates)
polydata=pd.read_csv('prettypolygon.txt', sep='\t')        
prettypoly=Polygon([(x,y) for x,y in zip(polydata['x'], polydata['y'])])

# read coastline data from cartopy - 50m and 10m also available?
coast = shapereader.natural_earth(resolution='110m',
                                  category='physical',
                                  name='coastline')

coastlines = shapereader.Reader(coast).geometries()
coast_latlons=[geom.xy for geom in coastlines]

extreme=11825474
limit=0.95*extreme # defining when point close to edge of projection
threshold=3e6 # detection value for strings that spill over to other edges of projection

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
            #ax.plot(item[:,0], item[:,1], color='green', lw=0.2)
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
            #ax.plot(rotseg[:,0], rotseg[:,1], color='orange', lw=0.2)
            check=fit_check(rotseg, prettypoly)
            if check is not None: 
                result.append(['Bordering', check])
        #ax.plot(broken_segments[-1][:,0], broken_segments[-1][:,1], color='blue', lw=0.2)
        check=fit_check(broken_segments[-1], prettypoly)
        if check is not None: 
            result.append(['Bordering', check])
    else: # if no breaks, just use whole linestring
       #ax.plot(spil_x,spil_y, color='red', lw=0.2)
       coords=np.column_stack(([x for x in spil_x],[y for y in spil_y]))
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
bordering_segments=pd.concat([coastresult[coastresult.Type == "Open"][:5],coastresult[coastresult.Type == "Bordering"]])['Coords']

# Bit of brute forcing in here, but it generates a continuous outer coastline without obvious out of sequence jumping
merged_line = bordering_segments.iloc[0][::-1]  # Start with the first segment
remaining_lines = bordering_segments.iloc[1:-2].tolist()  # segments to join to first: last 2 are tiny and screw things up
# This works. 
i=0
#while remaining_lines:
while i<23:
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

# A bit of jiggery-pokery to fix the one place the join is bad and close the bounding coast polygon
merged_line=np.concatenate([[merged_line[-1]],merged_line[:14][::-1],merged_line[14:]])

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
