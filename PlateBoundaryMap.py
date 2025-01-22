#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 11:14:14 2024

@author: crowan

Plotting plate boundaries on a basic Spilhaus coastline map

"""

import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from spilhaus import from_lonlat_to_spilhaus_xy, line_lonlat_to_pretty_spilhaus_xy
from shapely import Polygon


def plot_segs_file(file, ax, colour='black', order=1):
    f=open(file, 'r')
    ridgeseg=[]
    flag=0
    for line in f:
        if re.match('^>', line):
            if flag==1: 
                seglat=[]
                seglong=[]
                for data in ridgeseg:
                    data=data.rstrip('\n')
                    seglong.append(float(data.split("\t")[0]))
                    seglat.append(float(data.split("\t")[1]))
                spil_x, spil_y=from_lonlat_to_spilhaus_xy(np.array(seglong),np.array(seglat))
                ax.plot(spil_x, spil_y, color=colour, lw=0.8, zorder=order)
            ridgeseg = []
            flag=1
        else:
            ridgeseg.append(line)
    seglat=[]
    seglong=[]
    for data in ridgeseg:
        data=data.rstrip('\n')
        seglong.append(float(data.split("\t")[0]))
        seglat.append(float(data.split("\t")[1]))
    spil_x, spil_y=from_lonlat_to_spilhaus_xy(np.array(seglong),np.array(seglat))
    ax.plot(spil_x, spil_y, color=colour, lw=0.8, zorder=order)
    f.close()

def extract_segs_file(file):
    output=[]
    f=open(file, 'r')
    ridgeseg=[]
    flag=0
    for line in f:
        if re.match('^>', line):
            if flag==1: 
                seglat=[]
                seglong=[]
                for data in ridgeseg:
                    data=data.rstrip('\n')
                    seglong.append(float(data.split("\t")[0]))
                    seglat.append(float(data.split("\t")[1]))
                output.append(np.column_stack(([x for x in seglong],[y for y in seglat])))
            ridgeseg = []
            flag=1
        else:
            ridgeseg.append(line)
    f.close()
    seglat=[]
    seglong=[]
    for data in ridgeseg:
        data=data.rstrip('\n')
        seglong.append(float(data.split("\t")[0]))
        seglat.append(float(data.split("\t")[1]))
    output.append(np.column_stack(([x for x in seglong],[y for y in seglat])))    
    return output

def extract_coast_file(file):
    # only different from above because also extracts segtype
    # output is final list of lat-long coordinates, segtype the list of segment types listed in the break lines
    output=[]
    segtype=[]
    f=open(file, 'r')
    coastseg=[]
    flag=0
    for line in f: #looping through each line in file
        if re.match('^>', line): # looking for breaklines - where found extract segtype
            segtype.append(line.split(" ")[1].rstrip('\n'))
            if flag==1: # if not the first run have last section data to append to output
                seglat=[]
                seglong=[]
                for data in coastseg:
                    data=data.rstrip('\n')
                    seglong.append(float(data.split("\t")[0]))
                    seglat.append(float(data.split("\t")[1]))
                output.append(np.column_stack(([x for x in seglong],[y for y in seglat])))
            coastseg = []
            flag=1
        else:
            coastseg.append(line) # where no break just appends coordinates to coastseg
    f.close()
    seglat=[]
    seglong=[]
    for data in coastseg:
        data=data.rstrip('\n')
        seglong.append(float(data.split("\t")[0]))
        seglat.append(float(data.split("\t")[1]))
    output.append(np.column_stack(([x for x in seglong],[y for y in seglat])))
    return pd.DataFrame([[a,b] for a,b in zip(segtype, output)], columns=['Type','Coords'])


ridgecolour=(0.7372549019607844, 0.7411764705882353, 0.8627450980392157)


# prettypolygon.txt contains vertices for polygon tracing out mask for prettified map (in Spilhaus projection coordinates)
polydata=pd.read_csv('prettypolygon.txt', sep='\t')        
prettypoly=Polygon([(x,y) for x,y in zip(polydata['x'], polydata['y'])])


    
coastresult=extract_coast_file("./coastlines/Spilhaus_Coastlines_110m.txt")

ridgesegs=extract_segs_file("./Misc/Ridgesegs.txt")
subsegs=extract_segs_file("./Misc/Subsegs.txt")
othersegs=extract_segs_file("./Misc/othersegs.txt")

fig, ax = plt.subplots(1, 1, figsize=(16,16), dpi=300)
ax.plot(polydata['x'], polydata['y'], color='grey', lw=0.5)
for i,row in coastresult.iterrows():
    ax.plot(row.Coords[:,0], row.Coords[:,1], color='grey', lw=0.5)
for ridgeseg in ridgesegs:
    converted=line_lonlat_to_pretty_spilhaus_xy(ridgeseg,prettypoly)
    for coords in converted:
        ax.plot(coords[:,0], coords[:,1], color='red', lw=0.8, zorder=1)
for subseg in subsegs:
    converted=line_lonlat_to_pretty_spilhaus_xy(subseg,prettypoly)
    for coords in converted:
        ax.plot(coords[:,0], coords[:,1], color='purple', lw=0.8, zorder=1)
for otherseg in othersegs:
     converted=line_lonlat_to_pretty_spilhaus_xy(otherseg,prettypoly)
     for coords in converted:
         ax.plot(coords[:,0], coords[:,1], color='orange', lw=0.8, zorder=1)       
