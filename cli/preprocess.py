#!/usr/bin/env python

import argparse
import datetime
import functools
import numpy as np
import os

import sen1mosaic.core
import sen1mosaic.IO
import sen1mosaic.preprocess

import sen2mosaic.multiprocess

import pdb

####################################################################
### Command line interface for preprocessing Sentinel-1 GRD data ###
####################################################################



def main(infiles, output_dir = os.getcwd(), temp_dir = os.getcwd(), multilook = 2, output_name = 'processed', speckle_filter = False, short_chain = False, noorbit = False, output_units = 'decibels', gpt = '~/snap/bin/gpt', verbose = False):
    """main(infiles, output_dir = os.getcwd(), temp_dir = os.getcwd(), multilook = 2, output_name = 'processed', speckle_filter = False, short_chain = False, noorbit = False, output_units = 'decibels', gpt = '~/snap/bin/gpt', verbose = False)
    
    Preprocess Sentinel-1 GRD IW data from the Copernicus Open Access Data Hub. This functon takes a list of Sentinel-1 input files, and uses the SNAP graph processing tool to generate radiometric/terrain corrected images.
    
    Args:
        infiles: A list of input Sentinel-1 .zip files to process together. Where > 1 file, the files should all be from one overpass. See function splitFiles() for assistance with this.
        output_dir: Directory for output .dim/.data files. Defaults to current working directory.
        temp_dir: Directory to output temporary files. Defaults to current working directory.
        multilook: Multilook integer. Defaults to 2.
        output_name: Name to put in output files for identification. Defaults to 'processed'.
        speckle_filter: Set True to include a Refined Lee speckle filter.
        short_chain: Set True to run a shorter processing chain that omits some optional preprocessing steps at the expense of output quality.
        noorbit: Set True to skip the downloading of a precise Sentinel-1 orbit file.
        output_units: Units to output data, either in decibels or natural units.
        gpt: Path to SNAP graph processing tool. Defaults to ~/snap/bin/gpt.
        verbose: Set True to print progress.
    
    Returns:
        A boolean; True where processing completed successfully and False where something went wrong.
    """
        
    # Process input files    
    output_file = sen1mosaic.preprocess.processFiles(infiles, output_dir = output_dir, temp_dir = temp_dir, multilook = multilook, output_name = output_name, speckle_filter = speckle_filter, short_chain = short_chain, noorbit = noorbit, output_units = output_units, gpt = gpt, verbose = verbose)
    
    # Test that output file has been generated correctly.
    if sen1mosaic.preprocess.testCompletion(output_file, output_dir = output_dir) == False:
        for infile in infiles:
            print('WARNING: %s does not appear to have completed processing successfully.'%infile)
    
    else:
        for infile in infiles:
            print('File %s processed successfully'%infile)
    
    return sen1mosaic.preprocess.testCompletion(output_file, output_dir = output_dir)




if __name__ == '__main__':
    """
    A scipt to pre-process Sentinel-1 IW GRD for mosaicking purposes.
    """
    
    # Set up command line parser
    parser = argparse.ArgumentParser(description = 'Pre-process Sentinel-1 IW GRD data from the Copernicus Open Access Hub to radiometric/terrain corrected images.')

    parser._action_groups.pop()
    positional = parser.add_argument_group('Positional arguments')
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    
    # Positional arguments
    positional.add_argument('infiles', metavar = 'S1_FILES', type = str, default = [os.getcwd()], nargs='*', help='Input files. Specify a valid S1 input file (.zip), multiple files through wildcards, or a directory. Defaults to processing all S1 files in current working directory.')

    # Required arguments
    
    # Optional arguments
    
    optional.add_argument('-o', '--output_dir', metavar = 'PATH', type = str, default = os.getcwd(), help = "Output directory for processed files. Defaults to current working directory.")
    optional.add_argument('-n', '--output_name', metavar = 'STR', type = str, default = 'processed', help = "String to be included in output filenames for identification. Defaults to 'processed'.")
    optional.add_argument('-t', '--temp_dir', metavar = 'PATH', type = str, default = os.getcwd(), help = "Output directory for intermediate files. Defaults to current working directory.")
    optional.add_argument('-ms', '--max_scenes', metavar = 'N', type = int, default = 3, help = "Maximum number of scenes from an overpass to reconstitute and process together. Higher values result in fewer output files with fewer artefacts at scene boundaries, but require more RAM. Defaults to 3 scenes.")
    optional.add_argument('-m', '--multilook', metavar = 'N', type = int, default = 2, help = "Multilooking reduces image noise by degrading output resolution from ~10 x 10 m by a factor. Defaults to 2 (~20 x 20 m output).")
    optional.add_argument('-f', '--speckle_filter', action = 'store_true', help = "Apply a speckle filter (Refined Lee) to output images.")
    optional.add_argument('-u', '--output_units', metavar='UNITS', type=str, default='decibels', help="Output units, set to either decibels (default) or natural.")
    optional.add_argument('-s', '--short', action = 'store_true', help = "Perform a more rapid processing chain, ommitting some nonessential preprocessing steps.")
    optional.add_argument('-no', '--noorbit', action = 'store_true', help = "Skip downloading of a precise orbit file.")
    optional.add_argument('-g', '--gpt', metavar = 'PATH', type = str, default = '~/snap/bin/gpt', help='Path to graph processing tool. Defaults to ~/snap/bin/gpt.')
    optional.add_argument('-st', '--start', type = str, default = '20140101', help = "Start date for tiles to include in format YYYYMMDD. Defaults to processing all dates.")
    optional.add_argument('-en', '--end', type = str, default = datetime.datetime.today().strftime('%Y%m%d'), help = "End date for tiles to include in format YYYYMMDD. Defaults to processing all dates.")
    optional.add_argument('-v', '--verbose', action = 'store_true', help = "Print script progress.")
    optional.add_argument('-p', '--processes', type = int, metavar = 'N', default = 1, help = "Specify a maximum number of tiles to process in parallel. Note: more processes will require more resources. Defaults to 1.")

    #optional.add_argument('-ov', '--overlap', action = 'store_true', help = "Overlap scenes by one, which can be used to corret for artefacts at scene cut points. This requires more storage, and longer ocessing time")
    
    # Parse command line arguments    
    args = parser.parse_args()   
    
    # Extract all eligible input files (.zip, or directory containing .zip)
    infiles = sen1mosaic.IO.prepInfiles(args.infiles, image_type = 'pre')
    
    # Slim down files to those within date range
    infiles = sen1mosaic.preprocess.reduceFilesToTimePeriod(infiles, args.start, args.end)
    
    assert len(infiles) > 0, "No valid input files detected."
    
    # Convert arguments to absolute paths    
    infiles = np.array(sorted([os.path.abspath(i) for i in infiles])) # Also sort, and convert to an array.
    
    output_dir = os.path.abspath(args.output_dir)
    temp_dir = os.path.abspath(args.temp_dir)   
    
    # Determine which images should be processed together as one contiguous overpass
    infiles_split = sen1mosaic.preprocess.splitFiles(infiles, args.max_scenes, overlap = False)
    
    # Keep things simple if using one process
    if args.processes == 1:
        
        for input_files in infiles_split:
        
            # Execute module
            main(input_files, output_dir = args.output_dir, temp_dir = args.temp_dir, multilook = args.multilook, output_name = args.output_name, speckle_filter = args.speckle_filter, short_chain = args.short, noorbit = args.noorbit, output_units = args.output_units, gpt = args.gpt, verbose = args.verbose)
    
    else:
        
        # Or use multiprocessing
        main_partial = functools.partial(main, output_dir = args.output_dir, temp_dir = args.temp_dir, multilook = args.multilook, output_name = args.output_name, speckle_filter = args.speckle_filter, short_chain = args.short, noorbit = args.noorbit, output_units = args.output_units, gpt = args.gpt, verbose = args.verbose)
                    
        sen2mosaic.multiprocess.runWorkers(main_partial, args.processes, infiles_split)

   
