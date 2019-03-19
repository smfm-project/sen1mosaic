#!/usr/bin/env python

import argparse
import datetime
import glob
import numpy as np
import os
import pandas
import time

import pdb

def connectToAPI(username, password):
    '''
    Connect to the SciHub API with sentinelsat. Sign up at https://scihub.copernicus.eu/.
    
    Args:
        username: Scihub username. 
        password: Scihub password.        
    '''
    
    import sentinelsat
    
    # Let API be accessed by other functions
    global scihub_api
    
    # Connect to Sentinel API
    scihub_api = sentinelsat.SentinelAPI(username, password, 'https://scihub.copernicus.eu/dhus')
    

def _buildWkt(search_area):
    """
    Function to build a well-known text polygon from a list of extents.
    
    Args:
        search_area: A list in the format [minlon, minlat, maxlon, maxlat]
    
    Returns:
        A wkt POLYGON string.
    """
    
    # Get strings from search_area floats
    lonmin, latmin, lonmax, latmax = [str(float(i)) for i in search_area]
    
    wkt = 'POLYGON((%s %s,%s %s,%s %s,%s %s,%s %s))'%(lonmin,latmin,lonmax,latmin,lonmax,latmax,lonmin,latmax,lonmin,latmin)
    
    return wkt


def search(search_area, start = '20140403', end = datetime.datetime.today().strftime('%Y%m%d'), direction= '*'):
    """search(search_area, start = '20140403', end = datetime.datetime.today().strftime('%Y%m%d'), direction= '*')
    
    Searches for Sentinel-1 GRD IW images that meet conditions of date range and extent.
    
    Args:
        search_area: A list in the format [minlon, minlat, maxlon, maxlat]
        start: Start date for search in format YYYYMMDD. Start date may not precede 20140403, the launch date of Sentinel1-. Defaults to 20140403.
        end: End date for search in format YYYYMMDD. Defaults to today's date.
    
    Returns:
        A pandas dataframe with details of scenes matching conditions.
    """
    
    import sentinelsat

    # Test that we're connected to SciHub
    assert 'scihub_api' in globals(), "The global variable scihub_api doesn't exist. You should run connectToAPI(username, password) before searching the data archive."
    
    # Set up start and end dates
    startdate = sentinelsat.format_query_date(start)
    enddate = sentinelsat.format_query_date(end)
    
    # Build a POLYGON wkt
    search_polygon = _buildWkt(search_area)
    
    # Search data, filtering by options.
    products = scihub_api.query(search_polygon, beginposition = (startdate,enddate),
                         platformname = 'Sentinel-1', producttype = 'GRD', orbitdirection = direction, sensoroperationalmode = 'IW')

    # convert to Pandas DataFrame, which can be searched modified before commiting to download()
    products_df = scihub_api.to_dataframe(products)
    
    if products_df.empty: raise IOError("No products found. Check your search terms.")
    
    return products_df


def removeDuplicates(products_df, data_dir = os.getcwd()):
    '''removeDuplicates(products_df, data_dir = os.getcwd())
    
    Remove images from search results that have already been downloaded
    
    Args:
        products_df: Pandas dataframe from search() function.
        data_dir: Directory containing Sentinel-1 data. Defaults to current working directory.
    
    Returns:
        A dataframe with duplicate files removed.
    '''
    
    # Add trailing '/' to data directory
    data_dir = data_dir.rstrip('/') + '/'
    
    # Build a list of files that are already downloaded
    drop_files = []
    for filename in products_df['filename']:
        if len(glob.glob(data_dir + filename.split('.')[0] + '*')) > 0:
            drop_files.append(filename)
    
    # And drop them from the pandas dataframe
    products_df = products_df[~products_df['filename'].isin(drop_files)]
    
    return products_df


def download(products_df, output_dir = os.getcwd()):
    ''' download(products_df, output_dir = os.getcwd())
    
    Downloads all images from a dataframe produced by sentinelsat.
    
    Args:
        products_df: Pandas dataframe from search() function.
        output_dir: Optionally specify an output directory. Defaults to the present working directory.
    '''
    
    assert os.path.isdir(output_dir), "Output directory doesn't exist."
    
    if products_df.empty == True:
        raise IOError("No products found. Check your search terms.")
        
    else:
        
        print 'Downloading %s product(s)'%str(len(products_df))
        # Download selected products
        scihub_api.download_all(products_df['uuid'], output_dir)


def main(username, password, search_area, start, end, output_dir = os.getcwd(), direction = '*'):
    """main(username, password, search_area, start, end, output_dir = os.getcwd(), direction = '*')
    
    Download Sentinel-1 data from the Copernicus Open Access Hub, specifying a particular tile, date ranges and degrees of cloud cover. This is the function that is initiated from the command line.
    
    Args:
        username: Scihub username. Sign up at https://scihub.copernicus.eu/.
        password: Scihub password.
        search_area: A list in the format [minlon, minlat, maxlon, maxlat] defining the search area.
        start: Start date for search in format YYYYMMDD.
        end: End date for search in format YYYYMMDD.
        output_dir: Optionally specify an output directory. Defaults to the current working directory.
        direction: Orbital direction (either ASCENDING or DESCENDING). If not specified, both will be considered.
    """
    
    # Connect to API
    connectToAPI(username, password)
        
    # Search for files, return a data frame containing details of matching Sentinel-2 images
    products = search(search_area, start = start, end = end, direction = direction)
    
    # Identify images that already exist at the download location and remove them from products dataframe
    products = removeDuplicates(products, data_dir = output_dir)
    
    # Download products
    download(products, output_dir = output_dir)
    
    

if __name__ == '__main__':

    # Set up command line parser
    parser = argparse.ArgumentParser(description = "Download Sentinel-1 IW GRD data from the Copernicus Open Access Hub, specifying a latitude/longitude range, date ranges and ascending/descending orbits. Files that are already present in the destination directory won't be re-downloaded.")

    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    # Required arguments
    required.add_argument('-u', '--user', type = str, required = True, help = "Scihub username")
    required.add_argument('-p', '--password', type = str, metavar = 'PASS', required = True, help = "Scihub password")
    required.add_argument('-a', '--search_area', nargs = 4, metavar = ('LONMIN', 'LATMIN', 'LONMAX', 'LATMAX'), type = float, required = True, help = "Extent of search area, in format <lonmin latmin lonmax latmax>.")
    required.add_argument('-s', '--start', type = str, metavar = 'YYYYMMDD', required = True, help = "Start date for search in format YYYYMMDD.")
    required.add_argument('-e', '--end', type = str, metavar = 'YYYYMMDD', required = True, help = "End date for search in format YYYYMMDD.")    
    
    # Optional arguments
    optional.add_argument('-o', '--output_dir', type = str, metavar = 'PATH', default = os.getcwd(), help = "Specify an output directory. Defaults to the current working directory.")
    optional.add_argument('-d', '--direction', type = str, metavar = 'DIR', default = '*', help = "Specify <ASCENDING> or <DESCENDING> if only a single orbital direction should be downloaded. Defaults to downloading both.")
    
    # Get arguments from command line
    args = parser.parse_args()
    
    # Run through entire processing sequence
    main(args.user, args.password, args.search_area, args.start, args.end, output_dir = args.output_dir, direction = args.direction)
