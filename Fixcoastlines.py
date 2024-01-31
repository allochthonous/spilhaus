#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 14:49:56 2024

Add coastlines to spilhaus projection by extracting them from cartopy 
and adding breaks where curves spill over edges of projection


@author: crowan
"""

import matplotlib.pylot as plt
import numpy as np
import pandas as pd
import cartopy.io.shapereader as shapereader

from spilhaus import from_lonlat_to_spilhaus_xy

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



         
