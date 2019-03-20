#!/usr/bin/env python

import argparse
import datetime
import glob
import math
import numpy as np
import os
from scipy import ndimage
import subprocess

import utilities
import sen2mosaic.utilities

import pdb


## Primary functions

def _createOutputArray(md, dtype = np.uint16):
    '''
    Create an output array from metadata dictionary.
    
    Args:
        md: A metadata object from utiltities.Metadata()
    
    Returns:
        A numpy array sized to match the specification of the metadata dictionary.
    '''
    
    output_array = np.zeros((md.nrows, md.ncols), dtype = dtype)
    
    return output_array


def _createGdalDataset(md, data_out = None, filename = '', driver = 'MEM', dtype = 3, options = []):
    '''
    Function to create an empty gdal dataset with georefence info from metadata dictionary.

    Args:
        md: A metadata object
        data_out: Optionally specify an array of data to include in the gdal dataset.
        filename: Optionally specify an output filename, if image will be written to disk.
        driver: GDAL driver type (e.g. 'MEM', 'GTiff'). By default this function creates an array in memory, but set driver = 'GTiff' to make a GeoTiff. If writing a file to disk, the argument filename must be specified.
        dtype: Output data type. Default data type is a 16-bit unsigned integer (gdal.GDT_Int16, 3), but this can be specified using GDAL standards.
        options: A list containing other GDAL options (e.g. for compression, use [compress'LZW'].

    Returns:
        A GDAL dataset.
    '''
    
    from osgeo import gdal
    
    gdal_driver = gdal.GetDriverByName(driver)
    ds = gdal_driver.Create(filename, md.ncols, md.nrows, 1, dtype, options = options)
    ds.SetGeoTransform(md.geo_t)
    ds.SetProjection(md.proj.ExportToWkt())
    
    # If a data array specified, add it to the gdal dataset
    if type(data_out).__module__ == np.__name__:
        ds.GetRasterBand(1).WriteArray(data_out)
    
    # If a filename is specified, write the array to disk.
    if filename != '':
        ds = None
    
    return ds


def _reprojectImage(ds_source, ds_dest, md_source, md_dest):
    '''
    Reprojects a source image to match the coordinates of a destination GDAL dataset.
    
    Args:
        ds_source: A gdal dataset from _createGdalDataset() containing data to be repojected.
        ds_dest: A gdal dataset from _createGdalDataset(), with destination coordinate reference system and extent.
        md_source: A metadata dictionary created by buildMetadataDictionary() representing the source image.
        md_dest: A metadata dictionary created by buildMetadataDictionary() representing the destination image.
    
    Returns:
        A GDAL array with resampled data
    '''
    
    from osgeo import gdal
    
    proj_source = md_source.proj.ExportToWkt()
    proj_dest = md_dest.proj.ExportToWkt()
    
    # Reproject source into dest project coordinates
    gdal.ReprojectImage(ds_source, ds_dest, proj_source, proj_dest, gdal.GRA_NearestNeighbour)
            
    ds_resampled = ds_dest.GetRasterBand(1).ReadAsArray()
    
    return ds_resampled


def _updateDataArray(data_out, data_resampled, action = 'sum'):
    '''
    Function to update contents of output array based on image_n array.
    
    Args:
        data_out: A numpy array representing the band data to be output.
        data_resampled: A numpy array containing resampled band data to be added to data_out.
        action: 'sum', 'min', or 'max', which respectively adds data_resampled to data_out, replaces pixels in data_out with data_resampled where data_resampled < data_out, and replaces pixels in data_out with data_resampled where data_resampled > data_out.
    
    Returns:
        The data_out array with pixels from data_resampled added.
        
    '''
    
    assert action in ['sum', 'min', 'max'], "Variable 'action' must be set to 'sum', 'min' or 'max'. It was set to %s."%str(action)
    
    # Add usable data to data_out array   
    if action == 'sum':
        mask = data_resampled != 0
        data_out[mask] += data_resampled[mask]
    elif action == 'max':
        mask = np.logical_and(np.logical_or(data_resampled > data_out, data_out == 0), data_resampled != 0)
        data_out[mask] = data_resampled[mask]
    elif action == 'min':
        mask = np.logical_or(data_resampled < data_out, data_out == 0)
        data_out[mask] = data_resampled[mask]

    return data_out


def loadPolarisation(scene, pol, md_dest):
    '''
    Funciton to load and reproject a Sentinel-2 band array.
    
    Args:
        scene: A Sentinel-1 scene of class utilities.LoadScene().
        pol: The name of a polarisation to load (e.g. 'VV')
        md_dest: An object of class utilities.Metadata() to reproject image to.
    
    Returns:
        A numpy array of resampled data
    '''
     
    # Write array to a gdal dataset
    ds_source = _createGdalDataset(scene.metadata, dtype = 6, data_out = scene.getBand(pol))                

    # Create an empty gdal dataset for destination
    ds_dest = _createGdalDataset(md_dest, dtype = 6)
            
    # Reproject source to destination projection and extent
    data_resampled = _reprojectImage(ds_source, ds_dest, scene.metadata, md_dest)
        
    return data_resampled



def generateDataArray(scenes, pol, md_dest, output_dir = os.getcwd(), output_name = 'S1_output', verbose = False):
    """generateDataArray(source_files, pol, md_dest, output_dir = os.getcwd(), output_name = 'S1_output', verbose = False)
    
    Function which generates an output GeoTiff file from list of pre-processed S1 source files for a specified output polarisation and extent.

    Args:
        scenes: A list of pre-processed S1 input files.
        pol: Polarisation to process ('VV' or 'VH')
        md_dest: Dictionary from buildMetaDataDictionary() containing output projection details.
        output_dir: Directory to write output files. Defaults to current working directory.
        output_name: Optionally specify a string to prepend to output files. Defaults to 'S1_output'.
        
    Returns:
        A string with the filename pattern. Returns 'NODATA' where not valid input images.
    """
    
    # Create array to contain output array for this band. Together these arrays are used to calculate the mean, min, max and standard deviation of input images.
    data_num = _createOutputArray(md_dest, dtype = np.int16) # To track number of images for calculating mean
    data_sum = _createOutputArray(md_dest, dtype = np.float32) # To track sum of input images
    data_var = _createOutputArray(md_dest, dtype = np.float32) # To track sum of variance of input images
    data_min = _createOutputArray(md_dest, dtype = np.float32) # To track sum of max value of input images
    data_max = _createOutputArray(md_dest, dtype = np.float32) # To track sum of min value of input images
    
    data_date = _createOutputArray(md_dest, dtype = np.float32) # Data from each image
                   
    # For each source file
    for n, scene in enumerate(scenes):
        
        if verbose: print '    Adding pixels from %s'%scene.filename.split('/')[-1]
        
        
        # Update output arrays if we're finished with the previous date. Skip on first iteration as there's no data yet.
        if n != 0:
            if scene.datetime.date() != last_date:
                data_num = _updateDataArray(data_num, (data_date != 0) * 1, action = 'sum')
                data_sum = _updateDataArray(data_sum, data_date, action = 'sum')
                data_var = _updateDataArray(data_var, data_date ** 2, action = 'sum')
                data_min = _updateDataArray(data_min, data_date, action = 'min')
                data_max = _updateDataArray(data_max, data_date, action = 'max')
                
        # Update date for next loop
        last_date = scene.datetime.date()
        
        # Load data        
        data_resampled = loadPolarisation(scene, pol, md_dest)
        
        # Update array for this date (allowing only 1 measurement per date to be included in sum)
        data_date = _updateDataArray(data_date, data_resampled, action = 'min')
        
        # Tidy up
        ds_source = None
        ds_dest = None
   
    # Update output arrays on final loop
    data_num = _updateDataArray(data_num, (data_date != 0) * 1, action = 'sum')
    data_sum = _updateDataArray(data_sum, data_date, action = 'sum')
    data_var = _updateDataArray(data_var, data_date ** 2, action = 'sum')
    data_min = _updateDataArray(data_min, data_date, action = 'min')
    data_max = _updateDataArray(data_max, data_date, action = 'max')
    
    # Get rid of zeros in cases of no data
    data_num[data_num==0] = 1
    
    # Calculate mean of input data
    data_mean = data_sum / data_num.astype(np.float32)
    
    # Calculate std of input data (See: https://stackoverflow.com/questions/5543651/computing-standard-deviation-in-a-stream). Where standard deviation undefined (< 2 samples), set to 0.
    data_std = np.zeros_like(data_mean)
    data_std[data_num > 1] = ((data_num * data_var - np.abs(data_sum) * data_sum)[data_num > 1]) / ((np.abs(data_sum) * (data_num - 1))[data_num > 1])
    data_std = np.sqrt(data_std)
    data_std[data_num < 2] = 0.
    
    if verbose: print 'Outputting polarisation %s'%pol
    
    # Generate default output filename
    filename = '%s/%s_%s_%s_R%sm.tif'%(output_dir, output_name, '%s', pol, str(md_dest.res))
    
    # Output files (mean, stdev, max, min)
    ds_out = _createGdalDataset(md_dest, data_out = data_mean, filename = filename%'mean', driver='GTiff', dtype = 6, options = ['COMPRESS=LZW'])
    ds_out = _createGdalDataset(md_dest, data_out = data_std, filename = filename%'stdev', driver='GTiff', dtype = 6, options = ['COMPRESS=LZW'])    
    ds_out = _createGdalDataset(md_dest, data_out = data_max, filename = filename%'max', driver='GTiff', dtype = 6, options = ['COMPRESS=LZW'])
    ds_out = _createGdalDataset(md_dest, data_out = data_min, filename = filename%'min', driver='GTiff', dtype = 6, options = ['COMPRESS=LZW'])

    return filename


def buildVVVH(VV_file, VH_file, output_dir = os.getcwd(), output_name = 'S1_output'):
    """buildVVVH(VV_file, VH_file, output_dir = os.getcwd(), output_name = 'S1_output')
    
    Function to build a VV/VH array and output to GeoTiff.
    
    Args:
        VV_file: Path to a VV output Geotiff
        VH_file: Path to a VH output Geotiff
        output_dir: Directory to write output files. Defaults to current working directory.
        output_name: Optionally specify a string to prepend to output files. Defaults to 'S1_output'.
    
    Returns:
        Path to VV/VH GeoTiff
    """
    
    from osgeo import gdal
    
    # Load datasets
    ds_VV = gdal.Open(VV_file,0)
    ds_VH = gdal.Open(VH_file,0)
    
    data_VV = ds_VV.ReadAsArray()
    data_VH = ds_VH.ReadAsArray()
    
    mask = np.logical_or(data_VV == 0., data_VH == 0.)
    
    data_VV[data_VV >= 0] = -0.00001
    data_VH[data_VH >= 0] = -0.00001

    VV_VH = data_VV / data_VH
    
    VV_VH[mask] = 0.
        
    # Output to GeoTiff
    res = str(int(round(ds_VV.GetGeoTransform()[1])))
    filename = '%s/%s_%s_VVVH_R%sm.tif'%(output_dir, output_name, 'mean', res)
    
    gdal_driver = gdal.GetDriverByName('GTiff')
    ds = gdal_driver.Create(filename, ds_VV.RasterXSize, ds_VV.RasterYSize, 1, 6, options = ['COMPRESS=LZW'])
    ds.SetGeoTransform(ds_VV.GetGeoTransform())
    ds.SetProjection(ds_VV.GetProjection())
        
    ds.GetRasterBand(1).WriteArray(VV_VH)
    
    # And write
    ds = None
        
    return filename

    
def buildVRT(red_band, green_band, blue_band, output_path):
    """
    Builds a three band RGB vrt for image visualisation. Outputs a .VRT file.
    
    Args:
        red_band: Filename to add to red band
        green_band: Filename to add to green band
        blue_band: Filename to add to blue band
        output_name: Path to output file
    """
    
    # Remove trailing / from output directory name if present
    output_path = output_path.rstrip('/')
    
    # Ensure output name is a VRT
    if output_path[-4:] != '.vrt':
        output_path += '.vrt'
    
    command = ['gdalbuildvrt', '-separate', '-overwrite']
    command += [output_path, red_band, green_band, blue_band]
    
    subprocess.call(command)


    

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
    md_dest = sen2mosaic.utilities.Metadata(extent_dest, output_res, EPSG_dest)    
    
    # Keep track of output filenames
    filenames = []
    
    # Process images for each polarisation
    for pol in pol_list:
                
        if verbose: print 'Doing polarisation %s'%pol
        
        # Load metadata for all Sentinel-1 datasets
        scenes = [utilities.LoadScene(source_file) for source_file in source_files]
        
        # Sort input scenes
        scenes = utilities.sortScenes(scenes)
        
        # Reduce the pool of scenes to only those that overlap with output tile
        scenes_tile = utilities.getSourceFilesInTile(scenes, md_dest, pol = pol, start = start, end = end, verbose = verbose)
        
        # It's only worth processing a tile if at least one input image is inside tile
        if len(scenes_tile) == 0:
            print "    No data inside specified tile for polarisation %s. Skipping."%pol
            continue
        
        # Combine pixels into output images for each band
        filename = generateDataArray(scenes_tile, pol, md_dest, output_dir = output_dir, output_name = output_name, verbose = verbose)
        
        # Keep track of output filenames
        filenames.append(filename)
    
    # Skip if there are no 'VH' images (or no images at all).
    if pol_list.tolist() == ['VV','VH'] and len(filenames) > 1:
        
        # Build VRT output files for straightforward visualisation
        if verbose: print 'Building .VRT images for visualisation.'
        
        # Build a VV/VH image
        filename_VVVH = buildVVVH(filenames[pol_list.tolist().index('VV')]%'mean', filenames[pol_list.tolist().index('VH')]%'mean', output_dir = output_dir, output_name = output_name)
        
        # Build a false colour composite image
        # False colour image (VV, VH, VV/VH)
        buildVRT(filenames[pol_list.tolist().index('VV')]%'mean', filenames[pol_list.tolist().index('VH')]%'mean', filename_VVVH,'%s/%s_VVmean_VHmean_VVVH_R%s.vrt'%(output_dir, output_name, str(output_res)))
        
        # Build an alternative false colour composite image
        buildVRT(filenames[pol_list.tolist().index('VV')]%'min', filenames[pol_list.tolist().index('VH')]%'min', filenames[pol_list.tolist().index('VV')]%'stdev','%s/%s_VVmin_VHmin_VVstdev_R%s.vrt'%(output_dir, output_name, str(output_res)))
    
    if verbose: print 'Processing complete!'


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
    infiles = utilities.prepInfiles(infiles)
    
    main(infiles, args.target_extent, args.epsg, start = args.start, end = args.end, output_res = args.resolution, pol = args.pol, output_dir = output_dir, output_name = args.output_name, verbose = args.verbose)
    
    
