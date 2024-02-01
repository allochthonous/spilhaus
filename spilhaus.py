from typing import Tuple
import numpy.typing as npt
import numpy as np
import pandas as pd
from pyproj import Proj, Transformer
from numpy import sin, cos, tan, arcsin, arctan, arctan2, pi


def from_lonlat_to_spilhaus_xy(
        longitude: npt.NDArray,
        latitude: npt.NDArray
) -> Tuple[npt.NDArray, npt.NDArray]:
    """
    Converts longitude and latitude (degrees N and E) into Spilhaus coordinates

    :param longitude: -180 to 180
    :param latitude: -90 to 90
    :return: Spilhaus x (easting) and y (northing)
    """

    # constants (https://github.com/OSGeo/PROJ/issues/1851)
    e = np.sqrt(0.00669438)
    lat_center_deg = -49.56371678
    lon_center_deg = 66.94970198
    azimuth_deg = 40.17823482
    
    # parameters derived from constants
    lat_center_rad = lat_center_deg * pi / 180 
    lon_center_rad = lon_center_deg * pi / 180 
    azimuth_rad = azimuth_deg * pi / 180
    conformal_lat_center = -pi / 2 + 2 * arctan(
        tan(pi/4 + lat_center_rad/2) *
        ((1 - e * sin(lat_center_rad)) / (1 + e * sin(lat_center_rad))) ** (e / 2)
    )
    alpha = -arcsin(cos(conformal_lat_center) * cos(azimuth_rad))
    lambda_0 = lon_center_rad + arctan2(tan(azimuth_rad), -sin(conformal_lat_center))
    beta = pi + arctan2(-sin(azimuth_rad), -tan(conformal_lat_center))
    
    # coordinates in radians
    lon = longitude * pi / 180
    lat = latitude * pi / 180
    
    # conformal latitude, in radians
    lat_c = -pi / 2 + 2 * arctan(
        tan(pi/4 + lat/2) * ((1 - e * sin(lat)) / (1 + e * sin(lat))) ** (e / 2)
    )
    
    # transformed lat and lon, in degrees
    lat_s = 180 / pi * arcsin(sin(alpha) * sin(lat_c) - cos(alpha) * cos(lat_c) * cos(lon - lambda_0))
    lon_s = 180 / pi * (
        beta + arctan2(
            cos(lat_c) * sin(lon - lambda_0), 
            (sin(alpha) * cos(lat_c) * cos(lon - lambda_0) + cos(alpha) * sin(lat_c))
        )
    )
    
    # projects transformed coordinates onto plane (Adams World in a Square II)
    p = Proj(proj='adams_ws2')
    adams_x, adams_y = p(lon_s, lat_s)
    spilhaus_x = -(adams_x + adams_y) / np.sqrt(2)
    spilhaus_y = (adams_x - adams_y) / np.sqrt(2)
    
    return spilhaus_x, spilhaus_y
    
    
def make_spilhaus_xy_gridpoints(
        spilhaus_res: int = 1000
) -> pd.DataFrame:
    """
    Creates a data frame of Spilhaus coordinates

    :param spilhaus_res: grid resolution
    :return: df with columns `x` and `y` (nrow = spilhaus_res ** 2, ncol = 2)
    """
    
    # regular grid of points in Spilhaus map
    extreme = 11825474
    m = np.linspace(start=-extreme, stop=extreme, num=spilhaus_res)
    gr = np.array(np.meshgrid(m, m)).reshape(2, spilhaus_res ** 2).T
    spilhaus_df = pd.DataFrame({
        'x': gr[:, 0],
        'y': gr[:, 1]
    })
    return spilhaus_df
    
    
def prettify_spilhaus_df(
    spilhaus_df: pd.DataFrame,
    z_lower_bound: float = -np.inf,
    z_upper_bound: float = np.inf
) -> pd.DataFrame:
    """
    Prettifies a Spilhaus data frame

    :param spilhaus_df: raw Spilhaus data frame
    :param z_lower_bound: lower cutoff for valid z values
    :param z_upper_bound: upper cutoff for valid z values
    :return: prettier Spilhaus data frame
    """
    
    # mask for "out of bounds" values (used to filter z at the end)
    spilhaus_df['l'] = (
        (spilhaus_df['z'] <= z_lower_bound) | 
        (spilhaus_df['z'] >= z_upper_bound))
    
    spilhaus_x = np.array(spilhaus_df['x'])
    spilhaus_y = np.array(spilhaus_df['y'])
    spilhaus_z = np.array(spilhaus_df['z'])
    spilhaus_l = np.array(spilhaus_df['l'])
    
    extreme = 11825474 # limit of Spilhaus coordinate system
    spilhaus_res = int(np.sqrt(len(spilhaus_x)))

    # augmented grid points (creates larger spatial grid to tesselate data)
    # order:
    # - original grid
    # - repetition to left (data grid rotated 90º CCW)
    # - repetition to right (data grid rotated 90º CCW)
    # - repetition below (data grid rotation 90º CW)
    # - repetition above (data grid rotated 90º CW)
    
    aug_x = np.concatenate([
        spilhaus_x, 
        spilhaus_x - 2 * extreme,
        spilhaus_x + 2 * extreme,
        spilhaus_x,
        spilhaus_x
    ])
    aug_y = np.concatenate([
        spilhaus_y,
        spilhaus_y,
        spilhaus_y,
        spilhaus_y - 2 * extreme,
        spilhaus_y + 2 * extreme,
    ])
    
    # tesselate grid data, with appropriate rotations 
    aug_z = np.concatenate([
        spilhaus_z, 
        np.fliplr(spilhaus_z.reshape(spilhaus_res, spilhaus_res).T).reshape(-1,),
        np.flip(spilhaus_z.reshape(spilhaus_res, spilhaus_res).T, axis=1).reshape(-1,),
        np.flip(spilhaus_z.reshape(spilhaus_res, spilhaus_res).T, axis=0).reshape(-1,),
        np.flip(spilhaus_z.reshape(spilhaus_res, spilhaus_res).T, axis=0).reshape(-1,),
    ])
    aug_l = np.concatenate([
        spilhaus_l, 
        np.fliplr(spilhaus_l.reshape(spilhaus_res, spilhaus_res).T).reshape(-1,),
        np.flip(spilhaus_l.reshape(spilhaus_res, spilhaus_res).T, axis=1).reshape(-1,),
        np.flip(spilhaus_l.reshape(spilhaus_res, spilhaus_res).T, axis=0).reshape(-1,),
        np.flip(spilhaus_l.reshape(spilhaus_res, spilhaus_res).T, axis=0).reshape(-1,),
    ])
    
    # This is the bit that trims around the 'best' part of the tesselated grid to produce the "pretty" projection
    cutpoint = 1.1 * extreme
    keep = ~(
        aug_l 
        | (aug_x < -cutpoint)
        | (aug_x > cutpoint)
        | (aug_y < -cutpoint)
        | (aug_y > cutpoint) # after defining the absolute extremes, so more trimming here.
        | (aug_y > 1.089e7 - 0.176 * aug_x)
        | (aug_y > 1.6e7 + 0.8333 * aug_x)
        | (aug_x < -0.984e7 - 0.565 * aug_y)
        | (aug_y < -1.378e7 + 0.46 * aug_x)
        | (aug_x > 1.274e7 + 0.172 * aug_y)
        | (aug_y > 1e7 - 0.5 * aug_x)
        | (aug_y > 2.3e7 + aug_x)
        | ((aug_y < 0.29e7) & (aug_x < -1.114e7))
        | ((aug_y < 0.39e7) & (aug_x < -1.17e7))
        | ((aug_y < -1.21e7) & (aug_x > 0.295e7))
        | ((aug_y < -1.2e7) & (aug_x > 0.312e7))
        | ((aug_y < -1.16e7) & (aug_x > 0.4e7))
        | ((aug_y < -1.11e7) & (aug_x > 0.45e7))
    )

    def prettify_axis(u: npt.NDArray) -> npt.NDArray:
        unique_u = np.unique(u)
        n = len(unique_u)
        res_u = np.median(unique_u[1:n] - unique_u[0:(n - 1)])
        return ((u - np.min(u)) / res_u).astype(int)

    pretty_spilhaus_df = pd.DataFrame({
        'x': prettify_axis(aug_x[keep]),
        'y': prettify_axis(aug_y[keep]),
        'z': aug_z[keep]
    }).drop_duplicates(subset=['x', 'y'])
    
    return pretty_spilhaus_df
    

def from_spilhaus_xy_to_lonlat(
        spilhaus_x,
        spilhaus_y
) -> Tuple[npt.NDArray, npt.NDArray]:
    """
    Converts Spilhaus coordinates into longitude and latitude (degrees N and E)

    :param spilhaus_x: Spilhaus easting
    :param spilhaus_y: Spilhaus northing
    :return: longitude (-180 to 180) and latitude (-90 to 90)
    """
    
    # constants
    e = np.sqrt(0.00669438)
    lat_center_deg = -49.56371678
    lon_center_deg = 66.94970198
    azimuth_deg = 40.17823482
    
    # parameters derived from constants
    lat_center_rad = lat_center_deg * pi / 180 
    lon_center_rad = lon_center_deg * pi / 180 
    azimuth_rad = azimuth_deg * pi / 180
    conformal_lat_center = -pi / 2 + 2 * arctan(
        tan(pi/4 + lat_center_rad/2) *
        ((1 - e * sin(lat_center_rad)) / (1 + e * sin(lat_center_rad))) ** (e / 2)
    )
    alpha = -arcsin(cos(conformal_lat_center) * cos(azimuth_rad))
    lambda_0 = lon_center_rad + arctan2(tan(azimuth_rad), -sin(conformal_lat_center))
    beta = pi + arctan2(-sin(azimuth_rad), -tan(conformal_lat_center))
    
    # take spilhaus coordinates and compute transformed coords, in degrees
    itransformer = Transformer.from_crs({"proj":'adams_ws2'}, 4326, always_xy=True)
    adams_x = (spilhaus_y - spilhaus_x) * np.sqrt(2) / 2
    adams_y = - (spilhaus_y + spilhaus_x) * np.sqrt(2) / 2
    lon_s, lat_s = itransformer.transform(adams_x, adams_y)
    
    #transformed coords in radians
    lon_s_rad = lon_s * pi / 180
    lat_s_rad = lat_s * pi / 180
    
    # conformal latitude
    lat_c = arcsin(sin(alpha) * sin(lat_s_rad) + cos(alpha) * cos(lat_s_rad) * cos(lon_s_rad - beta))
    
    # longitude, in radians
    lon = lambda_0 + arctan2(
        cos(lat_s_rad) * sin(lon_s_rad - beta),
        sin(alpha) * cos(lat_s_rad) * cos(lon_s_rad - beta) - cos(alpha) * sin(lat_s_rad)
    )
    
    # latitude (iterative formula from https://mathworld.wolfram.com/ConformalLatitude.html)
    lat = lat_c
    for i in range(10):
        lat = -0.5 * pi + 2 * arctan(
            tan(pi / 4 + lat_c / 2) *
            ((1 + e * sin(lat)) / (1 - e * sin(lat))) ** (e / 2)
        )
        
    # coordinates in degrees
    longitude = ((lon * 180 / pi + 180) % 360) - 180
    latitude = lat * 180 / pi
    
    return longitude, latitude
