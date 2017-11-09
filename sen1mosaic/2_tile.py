#!/usr/bin/env python

import argparse
import glob
import glymur
import numpy as np
import os
from scipy import ndimage
import subprocess

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


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
    When building a large mosaic, it's necessary for input tiles to be in a consistent order to avoid strange overlap artefacts. This function sorts a list of files in alphabetical order by their tile reference.
    
    Args:
        source_files: A list of processed Sentinel-1 files
    
    Returns:
        A list of source_files alphabetised by date.
    '''
    
    source_files.sort() #TODO: ensure this is sorted by date
    
    return source_files


def _loadSourceFile(S1_file, pol):
    '''
    Loads a Sentinel-1 pre-processed file of a given polarisation into a numpy array.
    
    Args:
        S1_file: /path/to/a/ pre-processed S1 file
        pol: 
    
    Returns:
        A numpy array.
    '''
    
    from osgeo import gdal
    
    # Remove trailing slash from input filename, if it exists
    S1_file = S1_file.rstrip('/')
    
    # Identify source file following the standardised file pattern
    image_path = glob.glob(S1_file[:-4] + '.data/*_%s.img'%pol)[0]
       
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
        image_n: A numpy array representing the image number from _updateMaskArrays().
        n: An integer describing the image number (first image = 1, second image = 2 etc.).
        scl_out: A numpy array representing the SCL mask from _updateMaskArrays().
    
    Returns:
        The data_out array with pixels from data_resampled added.
        
    '''
    
    # Add good data to data_out array   
    if action == 'sum':
        data_out += data_resampled
    elif action == 'max':
        data_out[data_resampled > data_out] = data_resampled[data_resampled > data_out]
    else:
        print 'ERROR in update band array action selection,'

    return data_out


def getSourceMetadata(S1_file):
    '''
    Function to extract georefence info from 
    
    Args:
        S1_file:
    Returns:
        A list describing the extent of the file, in the format [xmin, ymin, xmax, ymax].
        EPSG code of the coordinate reference system of the file.
    '''
    
    from osgeo import gdal, osr
             
    # Remove trailing / from safe files if present 
    S1_file = S1_file.rstrip('/')
    
    # Find the xml file that contains file metadata
    input_file = glob.glob(S1_file[:-4] + '.data/*_VV.img')[0]
        
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
    
    # Find EPSG code to define projection
    srs = osr.SpatialReference(wkt=ds.GetProjection())
    srs.AutoIdentifyEPSG()
    EPSG = int(srs.GetAttrValue("AUTHORITY", 1))
    
    # Extract date string from filename
    date = S1_file.split('/')[-1].split('_')[-5]
    
    return extent, EPSG, date


def buildMetadataDictionary(extent_dest, res, EPSG):
    '''
    Build a metadata dictionary to describe the destination georeference info
    
    Args:
        extent_dest: List desciribing corner coordinate points in destination CRS [xmin, ymin, xmax, ymax]
        res: Integer describing pixel size in m (10, 20, or 60)
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
    md['nrows'] = int((md['lry'] - md['uly']) / md['yres'])
    md['ncols'] = int((md['lrx'] - md['ulx']) / md['xres'])
    
    # Define gdal geotransform (Affine)
    md['geo_t'] = (md['ulx'], md['xres'], 0, md['uly'], 0, md['yres'])
    
    return md


def getFilesInTile(source_files, md_dest):
    '''
    Takes a list of source files as input, and determines where each falls within extent of output tile.
    
    Args:
        source_files: A list of input files.
        md_dest: Dictionary from buildMetaDataDictionary() containing output projection details.

    Returns:
        A reduced list of source_files containing only files that will contribute to each tile.
    '''
    
    # Sort source files alphabetically by tile reference.    
    source_files = _sortSourceFiles(source_files)
    
    # Extract this image's resolution from md_dest
    res = md_dest['res']
                
    # Determine which L3A images are within specified tile bounds
    print 'Searching for source files within specified tile...'
    
    do_tile = []
    
    for S1_file in source_files:
                                           
        # Get source file metadata
        extent_source, EPSG_source, date = getSourceMetadata(S1_file)
        
        # Define source file metadata dictionary
        md_source = buildMetadataDictionary(extent_source, res, EPSG_source)
                
        # Skip processing the file if image falls outside of tile area
        if _testOutsideTile(md_source, md_dest):
            do_tile.append(False)
            continue
        
        print '    Found one: %s'%S1_file
        do_tile.append(True)
    
    # Get subset of source_files in specified tile
    source_files_tile = list(np.array(source_files)[np.array(do_tile)])
    
    return source_files_tile



def generateDataArray(source_files, pol, md_dest, output_dir = os.getcwd(), output_name = 'S1_output'):
    """
    
    Function which generates an output GeoTiff file from list of pre-processed S1 source files for a specified output polarisation and extent.

    Args:
        source_files: A list of pre-processed S1 input files.
        pol: 
        md_dest: Dictionary from buildMetaDataDictionary() containing output projection details.
        output_dir: Optionally specify directory for output file. Defaults to current working directory.
        output_name: Optionally specify a string to prepend to output files. Defaults to 'L3B_output'.
        
    Returns:
        A numpy array containing mosaic data for the input band.
    """
    
    # Create array to contain output array for this band
    data_out = _createOutputArray(md_dest, dtype = np.float32) # To add output data to
    n_images = _createOutputArray(md_dest, dtype = np.int16) # To track number of images for calculating mean
    data_date = _createOutputArray(md_dest, dtype = np.float32) # To add new data for each date (taking max value forward)
    
    # Extract this image's resolution from md_dest
    res = md_dest['res']
        
    # For each source file
    for n, source_file in enumerate(source_files):
                
        print '    Adding pixels from %s'%source_file.split('/')[-1]
        
        # Get source file metadata
        extent_source, EPSG_source, date = getSourceMetadata(source_file)
                
        # Define source file metadata dictionary
        md_source = buildMetadataDictionary(extent_source, res, EPSG_source)
        
        # Update output arrays if we're finished with a previous overpass
        if n != 0:
            if date != last_date:
                data_out = _updateDataArray(data_out, data_date, action = 'sum')
                n_images = _updateDataArray(n_images, (data_date != 0) * 1, action = 'sum')

        # Update date for next loop
        last_date = date
        
        # Load source data for the band
        data = _loadSourceFile(source_file, pol)
        
        # Write array to a gdal dataset
        ds_source = _createGdalDataset(md_source, data_out = data)                

        # Create an empty gdal dataset for destination
        ds_dest = _createGdalDataset(md_dest, dtype = 1)
                
        # Reproject source to destination projection and extent
        data_resampled = _reprojectImage(ds_source, ds_dest, md_source, md_dest)
        
        # Update array for this date (allowing only 1 measurement per date to be included in mean)
        data_date = _updateDataArray(data_date, data_resampled, action = 'max')
        
        # Tidy up
        ds_source = None
        ds_dest = None
   
    # Update output arrays on final loop
    data_out = _updateDataArray(data_out, data_date, action = 'sum')
    n_images = _updateDataArray(n_images, (data_date != 0) * 1, action = 'sum')
    
    # Get rid of zeros in cases of no data
    n_images[n_images==0] = 1 
    
    # Change data_out to a mean
    data_out = data_out / n_images.astype(np.float32)
    
    print 'Outputting polarisation %s'%pol

    # Write output for this band to disk
    ds_out = _createGdalDataset(md_dest, data_out = data_out,
                        filename = '%s/%s_%s_R%sm.tif'%(output_dir, output_name, pol, str(res)),
                        driver='GTiff', dtype = 6, options = ['COMPRESS=LZW'])

    return data_out


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


    

def main(source_files, extent_dest, EPSG_dest, output_res = 10,
    pol_list = ['VV', 'VH'],
    output_dir = os.getcwd(), output_name = 'S1_output'):
    """
    
    Function to run through the entire chain for converting output of sen2Three into custom mosaics. This is the function that is initiated from the command line.
    
    Args:
        source_files: A list of level 3A input files.
        extent_dest: List desciribing corner coordinate points in destination CRS [xmin, ymin, xmax, ymax].
        EPSG_dest: EPSG code of destination coordinate reference system. Must be a UTM projection. See: https://www.epsg-registry.org/ for codes.
        res_list: Optionally specify a list of integers describing pixel size in m (10, 20, or 60). Must be accompanied by a band_list of the same size. Defaults to native resolution of each band.
        band_list: Optionally specify a list of output band names. Must be accompanied by a res_list of the same size. Defaults to processing all 10 and 20 m bands.
        output_dir: Optionally specify an output directory.
        output_name: Optionally specify a string to precede output file names.
    """
    
    # TEST SECTION #
    #source_files = sorted(glob.glob('/home/sbowers3/DATA/S1_testdata/*.dim'))
    #extent_dest = [600000,7900000,700000,8100000]
    #EPSG_dest = 32736
    #output_res = 10
    #pol_list = ['VV','VH']
    #output_dir = '/home/sbowers3/DATA/S1_testdata/'
    #output_name = 'S1_output'
    # TEST SECTION #

    
    assert len(extent_dest) == 4, "Output extent must be specified in the format [xmin, ymin, xmax, ymax]"
    assert len(source_files) >= 1, "No source files in specified location."
    
    # Convert band and res list to numpy arrays for indexing
    pol_list = np.array(pol_list)
    
    # Remove trailing / from output directory if present 
    output_dir = output_dir.rstrip('/')
    
    # Build a dictionary with output projection metadata
    md_dest = buildMetadataDictionary(extent_dest, output_res, EPSG_dest)    
        
    # Reduce the pool of source_files to only those that overlap with output tile
    source_files_tile = getFilesInTile(source_files, md_dest)
        
    # It's only worth processing a tile if at least one input image is inside tile
    assert len(source_files_tile) >= 1, "No data inside specified tile. Not processing this tile."
        
    # Process images for each polarisation
    for pol in pol_list:
            
        print 'Doing polarisation %s'%pol
        
        # Using image_n, combine pixels into outputs images for each band
        band_out = generateDataArray(source_files_tile, pol, md_dest, output_dir = output_dir, output_name = output_name)
    
    # Build VRT output files for straightforward visualisation
    print 'Building .VRT images for visualisation'
    
    # NOT IMPLEMENTED YET
    # False colour image (VV, VH, VV/VH)
    #buildVRT('%s/%s_B04_R10m.tif'%(output_dir, output_name), '%s/%s_B03_R10m.tif'%(output_dir, output_name),
    #          '%s/%s_B02_R10m.tif'%(output_dir, output_name), '%s/%s_RGB.vrt'%(output_dir, output_name))
    
    
    print 'Processing complete!'



if __name__ == "__main__":
    
    # Set up command line parser
    parser = argparse.ArgumentParser(description = "Process Sentinel-2 level 3A data to unofficial 'level 3B'. This script mosaics L3A into a customisable grid square, based on specified UTM coordinate bounds. Files are output as GeoTiffs, which are easier to work with than JPEG2000 files.")

    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # Required arguments
    required.add_argument('infiles', metavar = 'S1_FILES', type = str, nargs = '+', help = 'Sentinel-1 processed input files in .dim format. Specify a valid S1 input file or multiple files through wildcards (e.g. PATH/TO/*.dim).')
    required.add_argument('-te', '--target_extent', nargs = 4, metavar = ('XMIN', 'YMIN', 'XMAX', 'YMAX'), type = float, help = "Extent of output image tile, in format <xmin, ymin, xmax, ymax>.")
    required.add_argument('-e', '--epsg', type=int, help="EPSG code for output image tile CRS. This must be UTM. Find the EPSG code of your output CRS as https://www.epsg-registry.org/.")

    # Optional arguments
    optional.add_argument('-o', '--output_dir', type=str, metavar = 'DIR', default = os.getcwd(), help="Optionally specify an output directory. If nothing specified, downloads will output to the present working directory, given a standard filename.")
    optional.add_argument('-n', '--output_name', type=str, metavar = 'NAME', default = 'S1_output', help="Optionally specify a string to precede output filename.")

    # Get arguments
    args = parser.parse_args()

    # Get absolute path of input .safe files.
    args.infiles = [os.path.abspath(i) for i in args.infiles]

    main(args.infiles, args.target_extent, args.epsg, output_dir = args.output_dir, output_name = args.output_name)
    
    
