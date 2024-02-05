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
from shapely import Polygon

from spilhaus import from_lonlat_to_spilhaus_xy

# prettypolygon.txt contains vertices for polygon tracing out mask for prettified map (in Spilhaus projection coordinates)
polydata=pd.read_csv('prettypolygon.txt', sep='\t')        
prettypoly=Polygon([(x,y) for x,y in zip(polydata['x'], polydata['y'])])

# read coastline data from cartopy - 50m and 10m also available?
coast = shapereader.natural_earth(resolution='110m',
                                  category='physical',
                                  name='coastline')

coastlines = shapereader.Reader(coast).geometries()
coast_latlons=[geom.xy for geom in coastlines]

# unit value in converted coordinates which triggers a break (based on visual inspection of histograms of diff.dx and diff.dy)
threshold=3e6

extreme=11825474
limit=0.99*extreme


fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)

for coastline in coast_latlons:
    lon, lat=np.array(coastline[0].tolist()),np.array(coastline[1].tolist())
    spil_x, spil_y=from_lonlat_to_spilhaus_xy(lon,lat)
    # calculates change in xy values for adjacent points in linestring
    diff=pd.DataFrame(np.column_stack(([a-b for a,b in zip(spil_x[1:],spil_x[:-1])], 
                                  [a-b for a,b in zip(spil_y[1:],spil_y[:-1])])),
                 columns=['dx','dy'])
    # finds spillover points where xy changes a lot between adjacent points 
    breaks=diff.index[((diff['dx']>threshold) | (diff['dx']<-threshold)) | ((diff['dy']>threshold) | (diff['dy']<-threshold))]
    # if breaks exist split linestring data at breakpoints and plot individually
    if len(breaks)>0:
       p=0
       for b in breaks:
           ax.plot(spil_x[p:b+1],spil_y[p:b+1], color='blue')
           p=b+1
       ax.plot(spil_x[p:],spil_y[p:], color='blue') #remember the last segment!
    else: # if no breaks, just use whole linestring
       ax.plot(spil_x,spil_y, color='red')  


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
    # if breaks exist split linestring data at breakpoints and plot individually
    if len(breaks)>0:
        broken_segments=[]
        p=0
        for b in breaks:
            broken_segments.append(np.column_stack(([x for x in spil_x[p:b+1]],[y for y in spil_y[p:b+1]])))
            p=b+1
        broken_segments.append(np.column_stack(([x for x in spil_x[p:]],[y for y in spil_y[p:]]))) #remember the last segment!
        result.append(broken_segments)
    


fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)
ax.plot(polydata['x'], polydata['y'], color='grey')

for thing in result:
    for item in thing[0:-1]: #last sections don't neccessarily end at barrier!
        ax.plot(item[:,0], item[:,1], color='green')
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
        ax.plot(rotseg[:,0], rotseg[:,1], color='orange')
    ax.plot(thing[-1][:,0], thing[-1][:,1], color='blue') # Don't care about flipping for final segment
    
  




polydata=pd.read_csv('prettypolygon2.txt', sep='\t')        
prettypoly=Polygon([(x,y) for x,y in zip(polydata['x'], polydata['y'])])

fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)
ax.plot(polydata['x'], polydata['y'], color='grey')

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
    # if breaks exist split linestring data at breakpoints and plot individually
    if len(breaks)>0:
        broken_segments=[]
        p=0
        for b in breaks:
            broken_segments.append(np.column_stack(([x for x in spil_x[p:b+1]],[y for y in spil_y[p:b+1]])))
            p=b+1
        broken_segments.append(np.column_stack(([x for x in spil_x[p:]],[y for y in spil_y[p:]]))) #remember the last segment!
        result.append(broken_segments)
        for item in broken_segments[:-1]: #last sections don't neccessarily end at barrier!
            ax.plot(item[:,0], item[:,1], color='green', lw=0.2)
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
            ax.plot(rotseg[:,0], rotseg[:,1], color='orange', lw=0.2)
        ax.plot(broken_segments[-1][:,0], broken_segments[-1][:,1], color='blue', lw=0.2)
    else: # if no breaks, just use whole linestring
       ax.plot(spil_x,spil_y, color='red', lw=0.2)  




