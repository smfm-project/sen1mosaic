import argparse
import datetime
import os
import numpy as np

import sen1mosaic.IO
import sen1mosaic.mosaic

import sen2mosaic.core

import pdb

def main(source_files, extent_dest, EPSG_dest, output_res = 20, pol = 'both', start = '20140101', end = datetime.datetime.today().strftime('%Y%m%d'), output_dir = os.getcwd(), output_name = 'S1_output', verbose = False):
    """main(source_files, extent_dest, EPSG_dest, output_res = 20, pol = 'both', start = '20140101', end = datetime.datetime.today().strftime('%Y%m%d'), output_dir = os.getcwd(), output_name = 'S1_output', verbose = False)
    
    Function to run through the entire chain for converting output of sen2Three into custom mosaics. This is the function that is initiated from the command line.
    
    Args:
        source_files: A list of level 3A input files.
        extent_dest: List desciribing corner coordinate points in destination CRS [xmin, ymin, xmax, ymax].
        EPSG_dest: EPSG code of destination coordinate reference system. Must be a UTM projection. See: https://www.epsg-registry.org/ for codes.
        output_res: Optionally specify a value for pixel size in m. Defaults to 20 m.
        pol: Optionally specify a polarisation ('VV' or 'VH'), or process 'both'. Defaults to 'both'.
        output_dir: Optionally specify an output directory. Defaults to current working directory.
        output_name: Optionally specify a string to precede output file names. Defaults to 'S1_output'
        verbose: Set True to print progress.

    """
       
    assert len(extent_dest) == 4, "Output extent must be specified in the format [xmin, ymin, xmax, ymax]"
    assert len(source_files) >= 1, "No source files in specified location."
    assert pol in ['VV', 'VH', 'both'], "Polarisation must be set to 'VV', 'VH' or 'both'. You input %s."%str(pol)
    
    # Test that output directory is writeable
    output_dir = os.path.abspath(os.path.expanduser(output_dir))
    assert os.path.exists(output_dir), "Output directory (%s) does not exist."%output_dir
    assert os.access(output_dir, os.W_OK), "Output directory (%s) does not have write permission. Try setting a different output directory"%output_dir
            
    # Convert band and res list to numpy arrays for indexing
    if pol == 'both':
        pol_list = np.array(['VV', 'VH'])
    else:
        pol_list = np.array([pol])
    
    # Get metadata for output dictionary
    md_dest = sen2mosaic.core.Metadata(extent_dest, output_res, EPSG_dest)
        
    # Keep track of output filenames
    filenames = []
    
    # Process images for each polarisation
    for pol in pol_list:
                    
        # Build composite image for list of input scenes
        filename = sen1mosaic.mosaic.buildComposite(source_files, pol, md_dest, start = start, end = end, output_name = output_name, output_dir = output_dir, verbose = verbose)       
        
        # Keep track of output filenames
        filenames.append(filename)
    
    # Skip if there are no 'VH' images (or no images at all).
    if pol_list.tolist() == ['VV','VH'] and len(filenames) > 1:
        
        # Build VRT output files for straightforward visualisation
        if verbose: print('Building .VRT images for visualisation.')
        
        # Build a VV/VH image
        filename_VVVH = sen1mosaic.mosaic.buildVVVH(filenames[pol_list.tolist().index('VV')]%'mean', filenames[pol_list.tolist().index('VH')]%'mean', md_dest, output_dir = output_dir, output_name = output_name)
        
        # Build a false colour composite image
        # False colour image (VV, VH, VV/VH)
        sen1mosaic.mosaic.buildVRT(filenames[pol_list.tolist().index('VV')]%'mean', filenames[pol_list.tolist().index('VH')]%'mean', filename_VVVH,'%s/%s_VVmean_VHmean_VVVH_R%s.vrt'%(output_dir, output_name, str(output_res)))
        
        # Build an alternative false colour composite image
        sen1mosaic.mosaic.buildVRT(filenames[pol_list.tolist().index('VV')]%'min', filenames[pol_list.tolist().index('VH')]%'min', filenames[pol_list.tolist().index('VV')]%'stdev','%s/%s_VVmin_VHmin_VVstdev_R%s.vrt'%(output_dir, output_name, str(output_res)))
    
    if verbose: print('Processing complete!')


if __name__ == "__main__":
    
    # Set up command line parser
    parser = argparse.ArgumentParser(description = "Collate preprocessed Sentinel-1 data into mosaicked tiles. This script mosaics Sentinel-1 data into a customisable grid square, based on specified UTM coordinate bounds. Files are output as GeoTiffs of mean, min, max, and standard deviation of each available backscatter.")

    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # Required arguments
    required.add_argument('-te', '--target_extent', nargs = 4, metavar = ('XMIN', 'YMIN', 'XMAX', 'YMAX'), type = float, help = "Extent of output image tile, in format <xmin, ymin, xmax, ymax>.")
    required.add_argument('-e', '--epsg', type=int, help="EPSG code for output image tile CRS. This must be UTM. Find the EPSG code of your output CRS as https://www.epsg-registry.org/.")

    # Optional arguments
    optional.add_argument('infiles', metavar = 'S1_FILES', type = str, default = [os.getcwd()], nargs = '*', help = 'Input files from preprocess.py. Specify a valid S1 input file (.dim), multiple files through wildcards, or a directory. Defaults to processing all S1 files in current working directory.')
    optional.add_argument('-st', '--start', type = str, default = '20140101', help = "Start date for tiles to include in format YYYYMMDD. Defaults to processing all dates.")
    optional.add_argument('-en', '--end', type = str, default = datetime.datetime.today().strftime('%Y%m%d'), help = "End date for tiles to include in format YYYYMMDD. Defaults to processing all dates.")
    optional.add_argument('-r', '--resolution', type = int, metavar = 'RES', default = 20, help=  "Output resolution in metres. Defaults to 20 m.")
    optional.add_argument('-o', '--output_dir', type = str, metavar = 'PATH', default = os.getcwd(), help = "Output directory. If nothing specified, downloads will output to the present working directory, given a standard filename.")
    optional.add_argument('-n', '--output_name', type=str, metavar = 'NAME', default = 'S1_output', help="Optionally specify a string to precede output filename.")
    optional.add_argument('-p', '--pol', type=str, metavar = 'POL', default = 'both', help="Specify a single polarisation ('VV' or 'VH') or 'both'. Defaults to processing both.")
    optional.add_argument('-v', '--verbose', action = 'store_true', help = "Print script progress.")
    
    # Get arguments
    args = parser.parse_args()
    
    # Convert arguments to absolute paths    
    infiles = sorted([os.path.abspath(i) for i in args.infiles])   
    output_dir = os.path.abspath(args.output_dir)
    
    # Extract all eligible input files (.dim, or directory containing .dim)
    infiles = sen1mosaic.IO.prepInfiles(infiles, image_type = 'post')
    
    main(infiles, args.target_extent, args.epsg, start = args.start, end = args.end, output_res = args.resolution, pol = args.pol, output_dir = output_dir, output_name = args.output_name, verbose = args.verbose)
    
    
