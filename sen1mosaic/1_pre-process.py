
import argparse
import pdb
import numpy as np
import os
import tempfile
import time
import sys
import uuid


def preprocess(infile, outfile, xmlfile = os.path.join(os.path.dirname(__file__), '../cfg/0_calibrate.xml')):
    
    os.system('~/snap/bin/gpt %s -x -Pinputfile=%s -Poutputfile=%s'\
      %(xmlfile,infile,outfile)) # -c 16384M #-c 32768M -q 16


def stitching(infiles, outfile, xmlfile = os.path.join(os.path.dirname(__file__), '../cfg/1_stitch_scenes.xml')):
    
    os.system('~/snap/bin/gpt %s -x -Pinputfiles=%s -Poutputfile=%s'\
      %(xmlfile,infiles,outfile)) # -c 16384M # -c 32768M -q 16


def correction(infile, outfile, xmlfile = os.path.join(os.path.dirname(__file__), '../cfg/2_terrain_correction_short.xml')):
    
    os.system('~/snap/bin/gpt %s -x -Pinputfile=%s -Poutputfile=%s'\
      %(xmlfile,infile,outfile)) # -c 16384M # -c 32768M -q 16


def _getMetaData(infile):
    '''
    Takes a Sentinel-1 filename as input, and returns metadata based on its filename.
    
    Args:
        infile: A string indicating a Sentinel-1 .SAFE or .zip filename
    
    Returns:
        A dictionary containing data on starttime, endtime, and satellite overpass.
    '''
    
    md = {}
    
    md['date'] = infile.split('/')[-1].split('_')[4].split('T')[0]
    md['starttime'] = infile.split('/')[-1].split('_')[4].split('T')[-1]
    md['endtime'] = infile.split('/')[-1].split('_')[5].split('T')[-1]
    md['orbit'] = infile.split('/')[-1].split('_')[6]
    md['datatake'] = infile.split('/')[-1].split('_')[7]
    
    return md


def getContiguousImages(infiles):
    '''
    Use Sentinel-1 metadata from filenames to determine which files should be processesed as part of a single pass.
    
    Args:
        infiles: An array of Sentinel-1 .zip or .SAFE files
    
    Returns:
        An array giving a unique integer to each contiguous group 
    '''
    
    # Sort input files alphabetically
    infiles = np.sort(infiles)
       
    # Initialise group counter and list
    this_group = 1
    group = np.zeros_like(infiles, dtype=np.int)
    
    for n,infile in enumerate(infiles):
        
        # Get start time of current image
        starttime = _getMetaData(infile)['starttime']
                
        # For the first loop, there's no previous file to compare to so skip this test
        if n > 0:
                # If two images are contigous, the endtime of the previous image will equal starttime of current image
            if starttime != endtime:
                # If two images are not contigous, form a new group code
                this_group += 1
            
        # Save end time of this image for the next loop
        endtime = _getMetaData(infile)['endtime']
        
        group[n] = (this_group)
    
    return group


def splitFiles(infiles, max_scenes):
    '''
    Split a 1-d numpy array into overlapping segments of size n. The purpose of this function is to prevent very long chains of Sentinel-1 data being processed together and crashing SNAP. Based on solution in https://stackoverflow.com/questions/36586897/splitting-a-python-list-into-a-list-of-overlapping-chunks.
    
    Args:
        infiles: A numpy array of Sentinel-1 file paths.
        max_scenes: Number of files per segment.
    
    Returns:
        A list of arrays split into segments of size max_scenes.
    '''
    
    overlap = 1  # overlap size
    
    infiles_split = [infiles[i:i + max_scenes - overlap + 1] for i in xrange(0,infiles.shape[0], max_scenes - overlap)]
    
    return infiles_split


def processFiles(infiles, output_dir = os.getcwd(), temp_dir = os.getcwd(), remove = True):
    '''
    A function to pre-process one or more Sentinel-1 IW GRD images in preparation for mosaicking with the SNAP Graph Processing Tool. Images are processed in three steps: 1) Calibration, 2) Reassembly into a single image (if >1 image used from an overpass), and 3) Geometric correction.
    
    Args:
        infiles: An array of Sentinel-1 .zip or .SAFE files
        output_dir: The directory to output the final file. Defaults to current working directory.
        temp_dir: A directory to output intermediate processing steps. These files are deleted after processing.
        remove: A boolean (True/False) value that determines whether intemediate files are deleted. Defaults to True.
    
    Returns:
        The path to the output file.
    '''
    
    # Ensure that output and temporary output directories have a '/'.
    output_dir = output_dir.rstrip('/') + '/'
    temp_dir = temp_dir.rstrip('/') + '/'
    
    # Build an output filename
    md_start = _getMetaData(infiles[0])
    md_end = _getMetaData(infiles[-1])
    
    output_file = output_dir + 'S1_processed_%s_%s_%s_%s_%s'%(md_start['date'],md_start['starttime'],md_end['endtime'],md_start['orbit'], md_start['datatake'])
    
    # Step 1: Run calibration SNAP processing chain
    preprocess_files = []
    
    for infile in infiles:
               
        # Determine a temporary output filename (which must be preceded by original filename. See: http://forum.step.esa.int/t/sliceassembly-op-after-eapphasecorrection-op/1959/5)
        outfile = temp_dir + infile.split('/')[-1][:-4] + '_cal'
            
        # Keep a record of which files have already been processed for each pass
        preprocess_files.append(outfile) 
        
        # Execute Graph Processing Tool
        preprocess(infile, outfile)
    
    # Step 2: Where more than one image, they need to be reassmbled into a single image    
    if len(preprocess_files) > 1: 
                            
        # Format input files to a string separated by commas
        infiles_formatted = ".dim,".join(preprocess_files) + ".dim"
            
        # Select graph that first reassembles multiple images
        outfile = preprocess_files[0][:-4] + '_stitch.dim'
        
        # Execute Graph Processing Tool
        stitching(infiles_formatted, outfile)
    
    # Step 3: Perform geometric correction
    correction(outfile, output_file)
    
    # Tidy up by deleting temporary intermediate files
    if remove:
        for this_file in preprocess_files:
            os.system('rm %s.dim'%this_file)
            os.system('rm -r %s.data'%this_file)
                
        if len(preprocess_files)>1:
            os.system('rm %s'%outfile)
            os.system('rm -r %s.data'%outfile[:-4])
    
    return output_file


def main(infiles, output_dir = os.getcwd(), temp_dir = os.getcwd(), max_scenes = 3):
    '''
    '''
    
    # Convert arguments to absolute paths    
    infiles = np.array(sorted([os.path.abspath(i) for i in infiles])) # Also sort, and convert to an array.
    output_dir = os.path.abspath(output_dir)
    temp_dir = os.path.abspath(temp_dir)
    
    # Determine which images should be processed together as one contiguous overpass
    group = getContiguousImages(infiles)
    
    # Process one group at a time
    for this_group in np.unique(group):
        
        infiles_split = splitFiles(infiles[group == this_group], max_scenes)
        
        for input_files in infiles_split:
        
            output_file = processFiles(input_files, output_dir = output_dir, temp_dir = temp_dir)
    


if __name__ == '__main__':
    """
    A scipt to pre-process Sentinel-1 IW GRD for mosaicking purposes.
    """
    
    # Set up command line parser
    parser = argparse.ArgumentParser(description = 'Pre-process Sentinel-1 IW GRD data from the Copernicus Open Access Hub for mosaicking purposes.')

    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    # Required arguments
    required.add_argument('infiles', metavar='N', type=str, nargs='+', help='Input files. Either specify a valid S1 input file, or multiple files through wildcards.')
    
    # Optional arguments
    optional.add_argument('-o', '--output_dir', type=str, default = os.getcwd(), help="Optionally specify an output directory or file. If nothing specified, we'll apply a standard filename and output to the present working directory.")
    optional.add_argument('-t', '--temp_dir', type=str, default = os.getcwd(), help="Optionally specify a temporary output directory. If nothing specified, we'll output intermediate files to the present working directory.")
    optional.add_argument('-m', '--max_scenes', type=str, default = 3, help="Optionally specify a maximum number of scenes to stitch in one go. If nothing specified, we'll set this to a default of 3 scenes.")

    # Parse command line arguments    
    args = parser.parse_args()
    
    # Execute module
    main(args.infiles, output_dir = args.output_dir, temp_dir = args.temp_dir, max_scenes = args.max_scenes)
        

        

