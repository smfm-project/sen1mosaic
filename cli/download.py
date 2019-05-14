#!/usr/bin/env python

import argparse
import datetime
import os

import sen1mosaic.download

import pdb

##############################################################
### Command line interface for downloading Sentinel-1 data ###
##############################################################


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
    sen1mosaic.download.connectToAPI(username, password)
        
    # Search for files, return a data frame containing details of matching Sentinel-2 images
    products = sen1mosaic.download.search(search_area, start = start, end = end, direction = direction)
    
    # Identify images that already exist at the download location and remove them from products dataframe
    products = sen1mosaic.download.removeDuplicates(products, data_dir = output_dir)
    
    # Download products
    sen1mosaic.download.download(products, output_dir = output_dir)
    
    

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
