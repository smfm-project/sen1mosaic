#!/usr/bin/env python

import argparse
import glob
import math
import numpy as np
import os
from scipy import ndimage
import subprocess


import pdb



### Functions for command line interface.

def _prepInfiles(infiles):
    """
    Function to identify valid input files for processing chain
    
    Args:
        infiles: A list of input files, directories, or tiles for Sentinel-1 inputs.
    Returns:
        A list of all Sentinel-1 .dim files in infiles.
    """
    
    # Get absolute path, stripped of symbolic links
    infiles = [os.path.abspath(os.path.realpath(infile)) for infile in infiles]
    
    # List to collate 
    infiles_reduced = []
    
    for infile in infiles:
         
        # Where infile is a directory:
        infiles_reduced.extend(glob.glob('%s/*.dim'%infile))
        
        # Where infile is a .dim file
        if infile.split('/')[-1].split('.')[-1] == 'dim': infiles_reduced.extend([infile])
    
    # Strip repeats (in case)
    infiles_reduced = list(set(infiles_reduced))
    
    # Reduce input files to only Sentinel-1 processed files from preprocess.py
    infiles_reduced = [infile for infile in infiles_reduced if ('S1_' in infile.split('/')[-1])]
    
    return infiles_reduced


## Primary functions

def _createOutputArray(md, dtype = np.uint16):
    '''
    Create an output array from metadata dictionary.
    
    Args:
        md: A metadata dictionary created by buildMetadataDictionary().
    
    Returns:
        A numpy array sized to match the specification of the metadata dictionary.
    '''
    
    output_array = np.zeros((md['nrows'], md['ncols']), dtype = dtype)
    
    return output_array


def _sortSourceFiles(source_files):
    '''
    When building a large mosaic, it's beest for input tiles to be in a consistent order to avoid strange overlap artefacts. This function sorts a list of files in alphabetical order by their filename.
    
    Args:
        source_files: A list of processed Sentinel-1 files
    
    Returns:
        A list of source_files sorted by date.
    '''
    
    source_files.sort()
    
    return source_files


def _loadSourceFile(S1_file, pol):
    '''
    Loads a Sentinel-1 pre-processed file of a given polarisation into a numpy array.
    
    Args:
        S1_file: /path/to/a/ pre-processed S1 .dim file
        pol: Polarisation (either 'VV' or 'VH'
    
    Returns:
        A numpy array.
    '''
    
    from osgeo import gdal
    
    assert pol in ['VV', 'VH'], "Polarisation must be either 'VV' or 'VH'."
    assert S1_file[-4:] == '.dim', "Input file must be a .dim file."
    
    # Remove trailing slash from input filename, if it exists
    S1_file = S1_file.rstrip('/')
    
    # Identify source file following the standardised file pattern
    image_path = glob.glob(S1_file[:-4] + '.data/*_%s*.img'%pol)[0]
       
    # Load the image 
    data = gdal.Open(image_path).ReadAsArray()
    
    return data



def _createGdalDataset(md, data_out = None, filename = '', driver = 'MEM', dtype = 3, options = []):
    '''
    Function to create an empty gdal dataset with georefence info from metadata dictionary.

    Args:
        md: A metadata dictionary created by buildMetadataDictionary().
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
    ds = gdal_driver.Create(filename, md['ncols'], md['nrows'], 1, dtype, options = options)
    ds.SetGeoTransform(md['geo_t'])
    ds.SetProjection(md['proj'].ExportToWkt())
    
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
    
    proj_source = md_source['proj'].ExportToWkt()
    proj_dest = md_dest['proj'].ExportToWkt()
    
    # Reproject source into dest project coordinates
    gdal.ReprojectImage(ds_source, ds_dest, proj_source, proj_dest, gdal.GRA_NearestNeighbour)
            
    ds_resampled = ds_dest.GetRasterBand(1).ReadAsArray()
    
    return ds_resampled


def _testOutsideTile(md_source, md_dest):
    '''
    Function that uses metadata dictionaries from buildMetadatadisctionary() metadata to test whether any part of a source data falls inside destination tile.
    
    Args:
        md_source: A metadata dictionary created by buildMetadataDictionary() representing the source image.
        md_dest: A metadata dictionary created by buildMetadataDictionary() representing the destination image.
        
    Returns:
        A boolean (True/False) value.
    '''
    
    from osgeo import osr
            
    # Set up function to translate coordinates from source to destination
    tx = osr.CoordinateTransformation(md_source['proj'], md_dest['proj'])
    
    # And translate the source coordinates
    md_source['ulx'], md_source['uly'], z = tx.TransformPoint(md_source['ulx'], md_source['uly'])
    md_source['lrx'], md_source['lry'], z = tx.TransformPoint(md_source['lrx'], md_source['lry'])   
    
    out_of_tile =  md_source['ulx'] >= md_dest['lrx'] or \
                   md_source['lrx'] <= md_dest['ulx'] or \
                   md_source['uly'] <= md_dest['lry'] or \
                   md_source['lry'] >= md_dest['uly']
    
    return out_of_tile


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
    
    # Add good data to data_out array   
    if action == 'sum':
        data_out += data_resampled
    elif action == 'max':
        data_out[data_resampled < data_out] = data_resampled[data_resampled < data_out]
    elif action == 'min':
        data_out[np.logical_or(data_resampled < data_out, data_out == 0)] = data_resampled[np.logical_or(data_resampled < data_out, data_out == 0)]

    return data_out


def getSourceMetadata(S1_file, pol = 'VV'):
    '''
    Function to extract georefence info from 
    
    Args:
        S1_file: A Sentinel-1 .dim file
        pol: Polarisation of an input file to extract metadata from. Defaults to 'VV'.
    Returns:
        A list describing the extent of the file, in the format [xmin, ymin, xmax, ymax].
        EPSG code of the coordinate reference system of the file.
    '''
    
    from osgeo import gdal, osr
    
    # Remove trailing / from safe files if present 
    S1_file = S1_file.rstrip('/')
    
    assert S1_file[-4:] == '.dim', "S1_file must be a .dim file."
    
    # Find the xml file that contains file metadata
    input_file = glob.glob(S1_file[:-4] + '.data/*_%s*.img'%pol)[0]
        
    ds = gdal.Open(input_file,0)
    geo_t = ds.GetGeoTransform()
    
    # Get array size
    nrows = ds.RasterYSize
    ncols = ds.RasterXSize
    
    # Get extent data
    ulx = geo_t[0]
    uly = geo_t[3]
    xres = geo_t[1]
    yres = geo_t[5]
    lrx = ulx + (xres * ncols)
    lry = uly + (yres * nrows)
    extent = [ulx, lry, lrx, uly]
    
    res = abs(xres)
    
    # Find EPSG code to define projection
    srs = osr.SpatialReference(wkt=ds.GetProjection())
    srs.AutoIdentifyEPSG()
    EPSG = int(srs.GetAttrValue("AUTHORITY", 1))
    
    # Extract date string from filename
    date = S1_file.split('/')[-1].split('_')[-5]
    
    return extent, res, EPSG, date


def buildMetadataDictionary(extent_dest, res, EPSG):
    '''
    Build a metadata dictionary to describe the destination georeference info.
    
    Args:
        extent_dest: List desciribing corner coordinate points in destination CRS [xmin, ymin, xmax, ymax]
        res: Integer describing pixel size in m
        EPSG: EPSG code of destination coordinate reference system. Must be a UTM projection. See: https://www.epsg-registry.org/ for codes.
    
    Returns:
        A dictionary containg projection info.
    '''
    
    from osgeo import osr
    
    # Set up an empty dictionary
    md = {}
    
    # Define projection from EPSG code
    md['EPSG_code'] = EPSG

    # Get GDAL projection string
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(EPSG)
    md['proj'] = proj
    
    # Get image extent data
    md['ulx'] = float(extent_dest[0])
    md['lry'] = float(extent_dest[1])
    md['lrx'] = float(extent_dest[2])
    md['uly'] = float(extent_dest[3])
    md['xres'] = float(res)
    md['yres'] = float(-res)

    # Save current resolution for future reference
    md['res'] = res
    
    # Calculate array size
    md['nrows'] = int(round((md['lry'] - md['uly']) / md['yres']))
    md['ncols'] = int(round((md['lrx'] - md['ulx']) / md['xres']))
    
    # Define gdal geotransform (Affine)
    md['geo_t'] = (md['ulx'], md['xres'], 0, md['uly'], 0, md['yres'])
    
    return md


def getFilesInTile(source_files, pol, md_dest, verbose = False):
    """
    Takes a list of source files as input, and determines where each falls within extent of output tile and contains data for input polarisation.
    
    Args:
        source_files: A list of input files.
        pol: Polarisation. Set to either 'VV' or 'VH'.
        md_dest: Dictionary from buildMetaDataDictionary() containing output projection details.
        verbose: Set True to print progress.

    Returns:
        A reduced list of source_files containing only files that will contribute to each tile.
    """
    
    # Sort source files alphabetically by tile reference.    
    source_files = _sortSourceFiles(source_files)
          
    # Determine which L3A images are within specified tile bounds
    if verbose: print 'Searching for source files within specified tile...'
    
    do_tile = []
    
    for S1_file in source_files:
                                           
        # Get source file metadata
        extent_source, res, EPSG_source, date = getSourceMetadata(S1_file)
        
        # Define source file metadata dictionary
        md_source = buildMetadataDictionary(extent_source, res, EPSG_source)
                
        # Skip processing the file if image falls outside of tile area
        if _testOutsideTile(md_source, md_dest):
            do_tile.append(False)
            continue
        
        if len(glob.glob(S1_file.rstrip('/')[:-4] + '.data/*_%s*.img'%pol)) == 0:
            do_tile.append(False)
            continue
        
        if verbose: print '    Found one: %s'%S1_file
        do_tile.append(True)
    
    # Get subset of source_files in specified tile
    source_files_tile = list(np.array(source_files)[np.array(do_tile)])
    
    return source_files_tile


def generateDataArray(source_files, pol, md_dest, output_dir = os.getcwd(), output_name = 'S1_output', verbose = False):
    """generateDataArray(source_files, pol, md_dest, output_dir = os.getcwd(), output_name = 'S1_output', verbose = False)
    
    Function which generates an output GeoTiff file from list of pre-processed S1 source files for a specified output polarisation and extent.

    Args:
        source_files: A list of pre-processed S1 input files.
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
    
    # Reduce the pool of source_files to only those that overlap with output tile.
    source_files_tile = getFilesInTile(source_files, pol, md_dest, verbose = verbose)
        
    # It's only worth processing anything if at least one input image is inside tile
    if len(source_files_tile) == 0:
        print "No data inside specified tile for polarisation %s. Not processing this tile."%pol
        return 'NODATA'
                
    # For each source file
    for n, source_file in enumerate(source_files_tile):
        
        if verbose: print '    Adding pixels from %s'%source_file.split('/')[-1]
               
        # Get source file metadata
        extent_source, res, EPSG_source, date = getSourceMetadata(source_file)
        
        # Define source file metadata dictionary
        md_source = buildMetadataDictionary(extent_source, res, EPSG_source)
        
        # Update output arrays if we're finished with the previous date. Skip on first iteration as there's no data yet.
        if n != 0:
            if date != last_date:
                
                data_num = _updateDataArray(data_num, (data_date != 0) * 1, action = 'sum')
                data_sum = _updateDataArray(data_sum, data_date, action = 'sum')
                data_var = _updateDataArray(data_var, data_date ** 2, action = 'sum')
                data_min = _updateDataArray(data_min, data_date, action = 'min')
                data_max = _updateDataArray(data_max, data_date, action = 'max')

        # Update date for next loop
        last_date = date
        
        # Load source data for the band
        data = _loadSourceFile(source_file, pol)
        
        # Write array to a gdal dataset
        ds_source = _createGdalDataset(md_source, dtype = 6, data_out = data)                

        # Create an empty gdal dataset for destination
        ds_dest = _createGdalDataset(md_dest, dtype = 6)
                
        # Reproject source to destination projection and extent
        data_resampled = _reprojectImage(ds_source, ds_dest, md_source, md_dest)
        
        # Update array for this date (allowing only 1 measurement per date to be included in sum)
        data_date = _updateDataArray(data_date, data_resampled, action = 'max')
        
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
    data_std[data_num > 1] = ((data_num * data_var - data_sum * data_sum)[data_num > 1]) / ((data_sum * (data_num - 1))[data_num > 1])
    data_std[data_std < 0] = 0.
    data_std = np.sqrt(data_std)
    
    if verbose: print 'Outputting polarisation %s'%pol
    
    # Generate default output filename
    filename = '%s/%s_%s_%s_R%sm.tif'%(output_dir, output_name, '%s', pol, str(md_dest['res']))
    
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


    

def main(source_files, extent_dest, EPSG_dest, output_res = 20, pol = 'both', output_dir = os.getcwd(), output_name = 'S1_output', verbose = False):
    """main(source_files, extent_dest, EPSG_dest, output_res = 20, pol = 'both', output_dir = os.getcwd(), output_name = 'S1_output', verbose = False)
    
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
    
    # Remove trailing / from output directory if present 
    output_dir = output_dir.rstrip('/')
    
    # Convert band and res list to numpy arrays for indexing
    if pol == 'both':
        pol_list = np.array(['VV', 'VH'])
    else:
        pol_list = np.array([pol])
    
    # Build a dictionary with output projection metadata
    md_dest = buildMetadataDictionary(extent_dest, output_res, EPSG_dest)    
        
    # Keep track of output filenames
    filenames = []
    
    # Process images for each polarisation
    for pol in pol_list:
        
        if verbose: print 'Doing polarisation %s'%pol
                
        # Combine pixels into output images for each band
        filename = generateDataArray(source_files, pol, md_dest, output_dir = output_dir, output_name = output_name, verbose = verbose)
        
        # Keep track of output filenames
        filenames.append(filename)
        
    # Build VRT output files for straightforward visualisation
    if verbose: print 'Building .VRT images for visualisation.'
    
    if pol_list.tolist() == ['VV','VH'] or pol_list.tolist() == ['VH', 'VV'] and 'NODATA' not in pol_list:
                
        # Build a VV/VH image
        filename_VVVH = buildVVVH(filenames[pol_list.tolist().index('VV')]%'mean', filenames[pol_list.tolist().index('VH')]%'mean', output_dir = output_dir, output_name = output_name)
        
        # Build a false colour composite image
        # False colour image (VV, VH, VV/VH)
        buildVRT(filenames[pol_list.tolist().index('VV')]%'mean', filenames[pol_list.tolist().index('VH')]%'mean', filename_VVVH,'%s/%s_FCC.vrt'%(output_dir, output_name))
        
        # Build an alternative false colour composite image
        buildVRT(filenames[pol_list.tolist().index('VV')]%'min', filenames[pol_list.tolist().index('VH')]%'min', filenames[pol_list.tolist().index('VV')]%'stdev','%s/%s_FCC2.vrt'%(output_dir, output_name))
    
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
    optional.add_argument('-r', '--resolution', type = int, metavar = 'RES', default = 20, help=  "Output resolution in metres. Defaults to 20 m.")
    optional.add_argument('-o', '--output_dir', type = str, metavar = 'PATH', default = os.getcwd(), help = "Output directory. If nothing specified, downloads will output to the present working directory, given a standard filename.")
    optional.add_argument('-n', '--output_name', type=str, metavar = 'NAME', default = 'S1_output', help="Optionally specify a string to precede output filename.")
    optional.add_argument('-p', '--pol', type=str, metavar = 'POL', default = 'both', help="Specify a single polarisation ('VV' or 'VH') or 'both'. Defaults to processing both.")
    optional.add_argument('-v', '--verbose', action = 'store_true', help = "Print script progress.")

    # Get arguments
    args = parser.parse_args()
    
    # Extract all eligible input files (.dim, or directory containing .dim)
    infiles = _prepInfiles(args.infiles)
    
    # Convert arguments to absolute paths    
    infiles = sorted([os.path.abspath(i) for i in infiles])
    output_dir = os.path.abspath(args.output_dir)

    main(infiles, args.target_extent, args.epsg, output_dir = output_dir, output_name = args.output_name, output_res = args.resolution, pol = args.pol, verbose = args.verbose)
    
    
