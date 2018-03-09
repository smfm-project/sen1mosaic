#!/usr/bin/env python

import argparse
import functools
import glob
import multiprocessing
import numpy as np
import os
from osgeo import gdal
import signal
import subprocess
import sys
import time

import pdb


### Functions to enable command line interface with multiprocessing

def _do_work(job_queue, counter=None):
    """
    Processes jobs from  the multiprocessing queue until all jobs are finished
    Adapted from: https://github.com/ikreymer/cdx-index-client
    
    Args:
        job_queue: multiprocessing.Queue() object
        counter: multiprocessing.Value() object
    """
    
    import Queue
        
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    
    while not job_queue.empty():
        try:
            job = job_queue.get_nowait()
            
            main_partial(job)

            num_done = 0
            with counter.get_lock():
                counter.value += 1
                num_done = counter.value
                
        except Queue.Empty:
            pass

        except KeyboardInterrupt:
            break

        except Exception:
            if not job:
                raise


def _run_workers(n_processes, jobs):
    """
    This script is a queuing system that respects KeyboardInterrupt.
    Adapted from: https://github.com/ikreymer/cdx-index-client
    Which in turn was adapted from: http://bryceboe.com/2012/02/14/python-multiprocessing-pool-and-keyboardinterrupt-revisited/
    
    Args:
        n_processes: Number of parallel processes
        jobs: List of input tiles for sen2cor
    """
    
    import psutil 
    
    # Queue up all jobs
    job_queue = multiprocessing.Queue()
    counter = multiprocessing.Value('i', 0)
    
    for job in jobs:
        job_queue.put(job)
    
    workers = []
    
    for i in xrange(0, n_processes):
        
        tmp = multiprocessing.Process(target=_do_work, args=(job_queue, counter))
        tmp.daemon = True
        tmp.start()
        workers.append(tmp)

    try:
        
        for worker in workers:
            worker.join()
            
    except KeyboardInterrupt:
        for worker in workers:
            print 'Keyboard interrupt (ctrl-c) detected. Exiting all processes.'
            # This is an impolite way to kill sen2cor, but it otherwise does not listen.
            parent = psutil.Process(worker.pid)
            children = parent.children(recursive=True)
            parent.send_signal(signal.SIGKILL)
            for process in children:
                process.send_signal(signal.SIGKILL)
            worker.terminate()
            worker.join()
            
        raise



### Primary functions

def _prepInfiles(infiles):
    """
    Function to identify valid input files for processing chain
    
    Args:
        infiles: A list of input files, directories, or tiles for Sentinel-1 inputs.
    Returns:
        A list of all Sentinel-1 IW GRD files in infiles.
    """
    
    # Get absolute path, stripped of symbolic links
    infiles = [os.path.abspath(os.path.realpath(infile)) for infile in infiles]
    
    # List to collate 
    infiles_reduced = []
    
    for infile in infiles:
         
        # Where infile is a directory:
        infiles_reduced.extend(glob.glob('%s/*.zip'%infile))
        
        # Where infile is a .zip file
        if infile.split('/')[-1].split('.')[-1] == 'zip': infiles_reduced.extend([infile])
    
    # Strip repeats (in case)
    infiles_reduced = list(set(infiles_reduced))
    
    # Reduce input files to only IW GRD files
    infiles_reduced = [infile for infile in infiles_reduced if ('_IW_GRDH_' in infile.split('/')[-1])]
    
    return infiles_reduced


def preprocessGraph(infile, outfile, short_chain = False, verbose = False):
    """
    Step 1: Preprocess input data.
    """
    
    if short_chain:
        xmlfile = os.path.join(os.path.dirname(__file__), '../cfg/1_calibrate_short.xml')
    else:
        xmlfile = os.path.join(os.path.dirname(__file__), '../cfg/1_calibrate.xml')
       
    # Execute chain
    os.system('~/snap/bin/gpt %s -x -Pinputfile=%s -Poutputfile=%s'%(xmlfile, infile, outfile))    


def multilookGraph(infile, outfile, multilook, single = True):
    """
    Step 2: Multilook and stitch scenes together.
    """
    
    # Pick processing chain for 1 single image vs a strip of several joined images
    if single:
        xmlfile = os.path.join(os.path.dirname(__file__), '../cfg/2_multilook_single.xml')
    else:
        xmlfile = os.path.join(os.path.dirname(__file__), '../cfg/2_multilook.xml')
    
    # Execute chain
    os.system('~/snap/bin/gpt %s -x -Pinputfiles=%s -Poutputfile=%s -Pmultilook=%s'\
      %(xmlfile, infile, outfile, multilook))


def correctionGraph(infile, outfile, extent, speckle_filter = False, short_chain = False):
    """
    Step 3: Terrain correction and filtering.
    """
    
    if speckle_filter and short_chain:
        xmlfile = os.path.join(os.path.dirname(__file__), '../cfg/3_terrain_correction_filter_short.xml')
    elif speckle_filter and not short_chain:
        xmlfile = os.path.join(os.path.dirname(__file__), '../cfg/3_terrain_correction_filter.xml')
    elif not speckle_filter and short_chain:
        xmlfile = os.path.join(os.path.dirname(__file__), '../cfg/3_terrain_correction_short.xml')
    else:
        xmlfile = os.path.join(os.path.dirname(__file__), '../cfg/3_terrain_correction.xml')
    
    os.system('~/snap/bin/gpt %s -x -Pinputfile=%s -Poutputfile=%s -Pextent=%s'\
      %(xmlfile,infile,outfile, extent))


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


def getExtent(infile, buffer_size = 1000, multilook = 2):
    '''
    Occasional border artifacts are left in Sentinel-1 data in the range direction. We remove pixels from each edge of the image to catch these. To perform this operation, we must get the extent of the image. This does waste data, but must remain until SNAP/Sentinel-1 data formats are consistent.
    
    # See also: http://forum.step.esa.int/t/grd-border-noise-removal-over-ocean-areas/1582/13
    
    Args:
        infile: /path/to/the/Sentinel-1.dim file
        buffer_size: Number of pixels to remove from range direction.
    
    Returns:
        A string with the new extent to use for this file.
    '''
    
    assert infile[-4:] == '.dim', "The input to getExtent() must be a .dim file. You input %s."%str(infile)
    assert type(buffer_size) == int, "buffer_size must be an integer. You input %s."%str(buffer_size)
    assert type(multilook) == int, "multilook must be an integer. You input %s."%str(multilook)
        
    filename = sorted(glob.glob(infile[:-4] + '.data/*.img'))[0]
    
    # Reduce buffer_size in line with degree of multilooking
    buffer_size = int(round(float(buffer_size) / float(multilook)))
    
    ds = gdal.Open(filename,0)
    
    # Multiply the second buffer size by two to account for the extent removed by the first buffer_size
    extent = buffer_size, 0, ds.RasterXSize - (buffer_size * 2), ds.RasterYSize
    
    return ','.join([str(i) for i in extent])


def processFiles(infiles, output_dir = os.getcwd(), temp_dir = os.getcwd(), multilook = 2, speckle_filter = False, short_chain = False, verbose = False):
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
    
    output_file = output_dir + 'S1_L2_%s_%s_%s_%s_%s'%(md_start['date'],md_start['starttime'],md_end['endtime'],md_start['orbit'], md_start['datatake'])
        
    # Step 1: Run calibration SNAP processing chain
    preprocess_files = []
    
    for infile in infiles:
               
        # Determine a temporary output filename (which must be preceded by original filename. See: http://forum.step.esa.int/t/sliceassembly-op-after-eapphasecorrection-op/1959/5)
        outfile = temp_dir + infile.split('/')[-1][:-4] + '_cal'
            
        # Keep a record of which files have already been processed for each pass
        preprocess_files.append(outfile) 
        
        if verbose: print 'Pre-processing %s'%infile
        
        # Execute Graph Processing Tool
        preprocessGraph(infile, outfile, short_chain = short_chain)
        
        outfile += '.dim' # preprocess should not have .dim following file, but correction requires it so add it here


    # Step 2: Perform multilooking
    
    # Where more than one image, they can be be reassembled into a single image. This removes artefacts from boundaries.
    if len(preprocess_files) > 1: 
                            
        # Format input files to a string separated by commas
        infiles_formatted = ".dim,".join(preprocess_files) + ".dim"
            
        # Select graph that first reassembles multiple images
        outfile = preprocess_files[-1] + '_mtl_%st.dim'%str(len(preprocess_files))
        
        single = False
        
        if verbose: print 'Multilooking %s'%infiles_formatted

    # And for case where only one file is input
    else:
        infiles_formatted = preprocess_files[0] + '.dim'
        
        outfile = preprocess_files[0] + '_mtl_1t.dim'
        
        single = True
        
        if verbose: print 'Multilooking %s'%infiles_formatted[0]
    
    # Execute Graph Processing Tool
    multilookGraph(infiles_formatted, outfile, multilook, single = single)
    
    # Step 3: Perform geometric correction
    
    # Get file extent
    extent = getExtent(outfile, multilook = multilook)
    
    if verbose: print 'Geometrically correcting %s'%outfile # outfile = latest file
    
    # Execute Graph Processing Tool
    correctionGraph(outfile, output_file, extent, speckle_filter = speckle_filter, short_chain = short_chain)
    
    # Tidy up by deleting temporary intermediate files
    for this_file in preprocess_files:
        if verbose: print 'Removing %s'%this_file
        os.system('rm %s.dim'%this_file)
        os.system('rm -r %s.data'%this_file)
            
    if verbose: print 'Removing %s'%outfile[:-4]
    os.system('rm %s'%outfile)
    os.system('rm -r %s.data'%outfile[:-4])
    
    return output_file


def testCompletion(L1_files, output_dir = os.getcwd()):
    """
    Function to test whether the processing chains have generated output data
    """
    
    print 'Not testing output file generation. Yet.'
    
    return True
    
    
    
def removeL1(L1_files):
    """
    Function to remove L1 input_files once processing has completed
    """
    
    print 'Not deleting input file. Yet.'



def splitFiles(infiles, max_scenes, overlap = False):
    '''
    Split a 1-d numpy array into overlapping segments of size n. The purpose of this function is to prevent very long chains of Sentinel-1 data being processed together and crashing SNAP. Based on solution in https://stackoverflow.com/questions/36586897/splitting-a-python-list-into-a-list-of-overlapping-chunks.
    
    Args:
        infiles: A numpy array of Sentinel-1 file paths.
        max_scenes: Number of files per segment.
        overlap: Set to True to re-process overlapping scnes to avoid ugly overlaps
    
    Returns:
        A list of arrays split into segments of size max_scenes.
    '''
    
    assert max_scenes > 0, "max_scenes must be > 0, else there won't be any scenes to process."
    if overlap: assert max_scenes > 1, "max_scenes must be > 1 if overlap is specified, else there won't be any scenes to overlap."

    # Get a unique group number. A unique number is assigned to each overpass.
    groups = getContiguousImages(infiles)
    
    # Get overlapping strips of Sentinel-1 data, where overlap == True
    if max_scenes > 1:
        
        infiles_split = []
    
        for group in np.unique(groups):
            
            these_infiles = infiles[groups == group]
            
            if overlap:
                n = max_scenes - 1
            else:
                n = max_scenes
            
            these_infiles_split = [these_infiles[i:i+max_scenes] for i in xrange(0,these_infiles.shape[0], n)]
            
            # This catches case where only one overlapping file is included
            if overlap:
                for i in these_infiles_split:
                    if len(these_infiles_split) > 1 and len(these_infiles_split[-1]) == 1:
                        these_infiles_split = these_infiles_split[:-1]
            
            infiles_split.extend(these_infiles_split)
    
    else:
        
        infiles_split = infiles
        
    return infiles_split


def main(infiles, output_dir = os.getcwd(), temp_dir = os.getcwd(), multilook = 2, speckle_filter = False, short_chain = False, remove = False, verbose = False):
    '''
    '''
    
    output_file = processFiles(infiles, output_dir = output_dir, temp_dir = temp_dir, multilook = multilook, speckle_filter = speckle_filter, short_chain = short_chain, verbose = verbose)
                        
    ##TODO: Test that output file has been generated correctly.
    if testCompletion(infiles, output_dir = output_dir) == False:
        print 'WARNING: %s did not complete processing.'%infile
        
    # TODO: Build a removal function
    if remove: removeL1(infiles_split)




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
    
    # Optional arguments
    optional.add_argument('infiles', metavar='N', type=str, default = [os.getcwd()], nargs='*', help='Input files. Specify a valid S1 input file (.zip), multiple files through wildcards, or a directory. Defaults to current working directory.')
    optional.add_argument('-o', '--output_dir', type=str, default = os.getcwd(), help = "Optionally specify an output directory or file. If nothing specified, we'll apply a standard filename and output to the present working directory.")
    optional.add_argument('-t', '--temp_dir', type=str, default = os.getcwd(), help = "Optionally specify a temporary output directory. If nothing specified, we'll output intermediate files to the present working directory.")
    optional.add_argument('-m', '--max_scenes', type=int, default = 3, help = "Optionally specify a maximum number of scenes to stitch in one go. If nothing specified, we'll set this to a default of 3 scenes.")
    optional.add_argument('-l', '--multilook', type=int, default = 2, help = "Multilook integer. The resolution of the output If not specified, we'll use 2.")
    optional.add_argument('-f', '--speckle_filter', action = 'store_true', help = "Apply a speckle filter (Refined Lee) to Sentinel-1 images. Defaults to False.")
    optional.add_argument('-s', '--short', action = 'store_true', help = "Perform a more rapid processing chain, ommitting some processing steps. Useful for testing. Defaults to False.")
    optional.add_argument('-r', '--remove', action = 'store_true', help = "Optionally delete input scenes after processing complete.")
    optional.add_argument('-v', '--verbose', action = 'store_true', help = "Set script to print progress.")
    optional.add_argument('-p', '--n_processes', type = int, metavar = 'N', default = 1, help = "Specify a maximum number of tiles to process in paralell. Bear in mind that more processes will require more memory. Defaults to 1.")
    optional.add_argument('-ov', '--overlap', action = 'store_true', help = "Overlap scenes by one, which reduces the impact of inter-scene boundaries at the expense of processing time")

    # Parse command line arguments    
    args = parser.parse_args()   
    
    # Extract all eligible input files (.zip, or directory containing .zip)
    infiles = _prepInfiles(args.infiles)
    
    # Convert arguments to absolute paths    
    infiles = np.array(sorted([os.path.abspath(i) for i in infiles])) # Also sort, and convert to an array.
    output_dir = os.path.abspath(args.output_dir)
    temp_dir = os.path.abspath(args.temp_dir)
    
    # Determine which images should be processed together as one contiguous overpass
    infiles_split = splitFiles(infiles, args.max_scenes, overlap = args.overlap)
    
    # Keep things simple if using one process
    if args.n_processes == 1:
            
        for input_files in infiles_split:
        
            # Execute module
            main(args.infiles, output_dir = args.output_dir, temp_dir = args.temp_dir, multilook = args.multilook, speckle_filter = args.speckle_filter, short_chain = args.short, remove = args.remove, verbose = args.verbose)
    
    else:
        
        main_partial = functools.partial(main, output_dir = args.output_dir, temp_dir = args.temp_dir, multilook = args.multilook, speckle_filter = args.speckle_filter, short_chain = args.short, remove = args.remove, verbose = args.verbose)
                    
        _run_workers(args.n_processes, infiles_split)

        

