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

# read coastline data from cartopy
coast = shapereader.natural_earth(resolution='110m',
                                  category='physical',
                                  name='coastline')

coastlines = shapereader.Reader(coast).geometries()
coast_latlons=[geom.xy for geom in coastlines]

# unit value in converted coordinates which triggers a break (based on visual inspection of histograms of diff.dx and diff.dy)
threshold=3e6

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

fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)
ax.plot(polydata['x'], polydata['y'], color='green')


# result=[]
# for coastline in coast_latlons:
#     lon, lat=np.array(coastline[0].tolist()),np.array(coastline[1].tolist())
#     spil_x, spil_y=from_lonlat_to_spilhaus_xy(lon,lat)
#     # calculates change in xy values for adjacent points in linestring
#     diff=pd.DataFrame(np.column_stack(([a-b for a,b in zip(spil_x[1:],spil_x[:-1])], 
#                                   [a-b for a,b in zip(spil_y[1:],spil_y[:-1])])),
#                  columns=['dx','dy'])
#     # finds spillover points where xy changes a lot between adjacent points 
#     breaks=diff.index[((diff['dx']>threshold) | (diff['dx']<-threshold)) | ((diff['dy']>threshold) | (diff['dy']<-threshold))]
#     # if breaks exist split linestring data at breakpoints and plot individually
#     if len(breaks)>0:
#         broken_segments=[]
#         p=0
#         for b in breaks:
#             broken_segments.append(np.column_stack(([x for x in spil_x[p:b+1]],[y for y in spil_y[p:b+1]])))
#             p=b+1
#         np.column_stack(([x for x in spil_x[p:]],[y for y in spil_y[p:]])) #remember the last segment!
#         result.append(broken_segments)
    
# extreme=11825474
# limit=0.99*extreme

for item in result[1:2]:
    xmaxes, xmins = [], []
    ymaxes, ymins = [], []
    alt1=item[0]
    first_flipped=np.flip(item[0],0)
    dx=first_flipped[:,0]-item[1][0][1]
    dy=first_flipped[:,1]-item[1][0][0]
    if np.max(alt1[:,0])> limit: # N side - E side rotation
        to_add=np.column_stack((item[1][0][0]+dy,item[1][0][1]-dx))
    elif np.max(alt1[:,1])> limit: # E side - N side rotation
        to_add=np.column_stack((item[1][0][0]-dy,item[1][0][1]+dx))
    elif np.min(alt1[:,0]) < - limit: # S side - W side rotation
        to_add=np.column_stack((item[1][0][0]+dy,item[1][0][1]-dx))
    elif np.min(alt1[:,1]) < - limit: # W side - S side rotation
        to_add=np.column_stack((item[1][0][0]-dy,item[1][0][1]+dx))   
    alt2=np.flip(to_add,0) # need to back calculate the first segment!
    lastpoint=item[0][-1]
    flag=1
    for thing in item[1:]:
        dx=thing[:,0]-lastpoint[1]
        dy=thing[:,1]-lastpoint[0]
        
        # note transformation is with respect to previous segment
        if np.max(thing[:,0])> limit: # N side - E side rotation
            to_add=np.column_stack((lastpoint[0]+dy,lastpoint[1]-dx))
        elif np.max(thing[:,1])> limit: # E side - N side rotation
            to_add=np.column_stack((lastpoint[0]-dy,lastpoint[1]+dx))
        elif np.min(thing[:,0]) < - limit: # S side - W side rotation
            to_add=np.column_stack((lastpoint[0]+dy,lastpoint[1]-dx))
        elif np.min(thing[:,1]) < - limit: # W side - S side rotation
            to_add=np.column_stack((lastpoint[0]-dy,lastpoint[1]+dx))      
        lastpoint=thing[-1]
        if flag==1: #concatenate rotated with alt1 and non-rotated with alt2
            alt1=np.concatenate((alt1,to_add))
            alt2=np.concatenate((alt2,thing))
            flag=-1
        else: #concatenate non-rotated with alt1 and rotated with alt2
            alt1=np.concatenate((alt1,thing))
            alt2=np.concatenate((alt2,to_add))
            flag=1
    ax.plot(alt1[:,0], alt1[:,1], color='green')
    ax.plot(alt2[:,0], alt2[:,1], color='orange')






