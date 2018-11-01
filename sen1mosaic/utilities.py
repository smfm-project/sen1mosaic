
import datetime
import glob
import lxml.etree as ET
import numpy as np
import os

import sen2mosaic.utilities

import pdb

"""
class Metadata(object):
    '''
    This is a generic metadata class for Geosptial data
    '''
    
    def __init__(self, extent, res, EPSG):
        '''
        Args:
            extent: A list in the form [xmin. ymin, xmax, ymax]
            res: Pixel resolution
            EPSG: The EPSG code of the desired resolution
        '''
           
        
        # Define projection from EPSG code
        self.EPSG_code = EPSG
        
        # Define resolution
        self.res = res
        
        self.xres = float(res)
        self.yres = float(-res)
        
        # Define image extent data
        self.extent = extent
        
        self.ulx = float(extent[0])
        self.lry = float(extent[1])
        self.lrx = float(extent[2])
        self.uly = float(extent[3])
        
        # Get projection
        self.proj = self.__getProjection()
                
        # Calculate array size
        self.nrows = self.__getNRows()
        self.ncols = self.__getNCols()
        
        # Define gdal geotransform (Affine)
        self.geo_t = self.__getGeoT()
        
        
    def __getProjection(self):
        '''
        '''
        
        from osgeo import osr
        
        # Get GDAL projection string
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(self.EPSG_code)
        
        return proj
    
    def __getNRows(self):
        '''
        '''
        
        return int(round((self.lry - self.uly) / self.yres))
    
    def __getNCols(self):
        '''
        '''
        
        return int(round((self.lrx - self.ulx) / self.xres))
    
    def __getGeoT(self):
        '''
        '''
        
        geo_t = (self.ulx, self.xres, 0, self.uly, 0, self.yres)
        
        return geo_t
"""



class LoadScene(object):
    '''
    Load a Sentinel-1 (pre-processed) scene
    '''
        
    def __init__(self, filename):
        '''
        Args:
            filename: The path to a Sentinel-1 .dim file
        '''
                
        # Format filename, and check that it exists
        self.filename = self.__checkFilename(filename)
                
        # Get file format
        self.file_format = self.__getFormat()
        
        # Save satellite name
        self.satellite = 'S1'
        
        # Save image type (S1_single, S1_dual, S2)
        self.image_type = self.__getImageType()
        
        self.tile = self.__getTileID()
        
        self.__getMetadata()
        
        # Define source metadata
        self.metadata = sen2mosaic.utilities.Metadata(self.extent, self.resolution, self.EPSG)
        
        
    def __checkFilename(self, filename):
        '''
        Test that the granule exists
        '''
        
        # Get rid of trailing '/' if present
        filename = filename.rstrip()
        
        # Test that file exists
        assert os.path.exists(filename),"Cannot find file %s "%filename
        
        return filename
    
    def __getFormat(self):
        '''
        Test that the file of of an appropriate format
        '''
        
        if self.filename.split('/')[-1].split('.')[-1] == 'dim':
            return 'BEAM-DIMAP'

        else:
            print 'File %s does not match any expected file pattern'%self.filename
            raise IOError
    
    def __getImageType(self):
        '''
        Test whethere S1 file is single ('S1single') or dual ('S1dual') polarised 
        '''
        
        if len(glob.glob(self.filename.split('.dim')[0] +'.data/*0_VH*.img')) > 0:
            image_type = 'S1dual'
        else:
            image_type = 'S1single'
        
        return image_type
    
    def __getTileID(self):
        '''
        '''
        
        return '_'.join(self.filename.split('/')[-1].split('_')[-4:])
    
    def __getMetadata(self):
        '''
        Extract metadata from the Sentinel-1 file.
        '''
        
        self.extent, self.EPSG, self.resolution, self.datetime, self.overpass = getS1Metadata(self.filename)
    
    def __getImagePath(self, pol = 'VV'):
        '''
        Get the path to a mask or polarisation
        '''

        # Identify source file following the standardised file pattern
        
        image_path = glob.glob(self.filename.split('.dim')[0]+'.data/*0_%s*.img'%pol)
        
        if len(image_path) == 0:
            raise IOError
        
        return image_path[0]

    
    def getMask(self, md = None):
        '''
        Load the mask to a numpy array.
        '''
        
        from osgeo import gdal
                
        # Load mask
        image_path = self.__getImagePath()
        
        mask = gdal.Open(image_path, 0).ReadAsArray()
        
        # Keep track of pixels where data are contained
        data = mask != 0
        
        # Reproject?
        if md is not None:
            data = sen2mosaic.utilities.reprojectBand(self, data, md, dtype = 1, resampling = gdal.GRA_Mode).astype(np.bool)
        
        # Return pixels without data
        return data == 0
    
    def getBand(self, pol, md = None):
        '''
        Load a single polarisation
        
        Args:
            pol: 'VV' or 'VH'
        '''

        from osgeo import gdal

        assert pol in ['VV', 'VH'], "Polarisation can only be 'VV' or 'VH'."
        if self.image_type == 'S1single':
            assert pol != 'VH', "The Sentinel-1 image must be dual polarised to load the VH polarisation."
        
        image_path = self.__getImagePath(pol = pol)
        
        # Load the image
        data = gdal.Open(image_path, 0).ReadAsArray()

        # Reproject?
        if md is not None:
            data = sen2mosaic.utilities.reprojectBand(self, data, md, dtype = 6, resampling = gdal.GRA_Average) 
        
        return data



def getS1Metadata(dim_file):
    '''
    Function to extract georefence info from level 2A Sentinel 2 data in .SAFE format.
    
    Args:
        dim_file: 
        
    Returns:
        A list describing the extent of the .dim file, in the format [xmin, ymin, xmax, ymax].
        EPSG code of the coordinate reference system of the image
        The image resolution
    '''
    
    from osgeo import gdal, osr
    
    assert os.path.exists(dim_file), "The location %s does not contain a Sentinel-1 .dim file."%dim_file
    
    tree = ET.ElementTree(file = dim_file)
    root = tree.getroot()
    
    # Get array size
    size = root.find("Raster_Dimensions")  
    nrows = int(size.find('NROWS').text)
    ncols = int(size.find('NCOLS').text)
    
    geopos = root.find("Geoposition/IMAGE_TO_MODEL_TRANSFORM").text.split(',')
    ulx = float(geopos[4])
    uly = float(geopos[5])
    xres = float(geopos[0])
    yres = float(geopos[3])
    lrx = ulx + (xres * ncols)
    lry = uly + (yres * nrows)
    extent = [ulx, lry, lrx, uly]
    
    res = abs(xres)
        
    wkt = root.find("Coordinate_Reference_System/WKT").text
    
    srs = osr.SpatialReference(wkt = wkt)
    srs.AutoIdentifyEPSG()
    EPSG = int(srs.GetAttrValue("AUTHORITY", 1))
    
    # Extract date string from filename
    datestring = root.find("Production/PRODUCT_SCENE_RASTER_START_TIME").text.split('.')[0]
    this_datetime = datetime.datetime.strptime(datestring, '%d-%b-%Y %H:%M:%S')
    
    # Get ascending/descending overpass
    overpass = root.find("Dataset_Sources/MDElem/MDElem/MDATTR[@name='PASS']").text
    
    return extent, EPSG, res, this_datetime, overpass


def prepInfiles(infiles):
    """
    Function to identify valid input files for processing chain
    
    Args:
        infiles: A list of input files, directories, or tiles for Sentinel-1 inputs.
    Returns:
        A list of all Sentinel-1 .dim files in infiles.
    """
    
    # Make interable if only one item
    if not isinstance(infiles, list):
        infiles = [infiles]
        
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



def getSourceFilesInTile(source_files, md_dest, pol = 'VV', start = '20140101', end = datetime.datetime.today().strftime('%Y%m%d'), verbose = False):
    """getSourceFilesInTile(source_files, md_dest, pol = 'VV', start = '20140101', end = datetime.datetime.today().strftime('%Y%m%d'), verbose = False)
    
    Takes a list of source files as input, and determines where each falls within extent of output tile and contains data for input polarisation.
    
    Args:
        source_files: A list of input files.
        md_dest: Dictionary from buildMetaDataDictionary() containing output projection details.
        pol: Polarisation to process (i.e. 'VV' or 'VH). Defaults to 'VV'.
        start: First date to consider (in form YYYYMMDD). Defaults to 20140101.
        end: Last date to consider (in form YYYYMMDD). Defaults to today's date.
        verbose: Set True to print progress.

    Returns:
        A reduced list of source_files containing only files that will contribute to each tile.
    """
    
    def testOutsideTile(md_source, md_dest):
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
        tx = osr.CoordinateTransformation(md_source.proj, md_dest.proj)
        
        # And translate the source coordinates
        md_source.ulx, md_source.uly, z = tx.TransformPoint(md_source.ulx, md_source.uly)
        md_source.lrx, md_source.lry, z = tx.TransformPoint(md_source.lrx, md_source.lry)   
        
        out_of_tile =  md_source.ulx >= md_dest.lrx or \
                    md_source.lrx <= md_dest.ulx or \
                    md_source.uly <= md_dest.lry or \
                    md_source.lry >= md_dest.uly
        
        return out_of_tile
    
    def testOutsideDate(scene, start = '20140101', end = datetime.datetime.today().strftime('%Y%m%d')):
        '''
        Function that uses LoadScene class to test whether a tile falls within the specified time range.
        
        Args:
            scene: Object from utilties.LoadScene()
            start: Start date to process, in format 'YYYYMMDD' Defaults to start of Sentinel-2 era.
            end: End date to process, in format 'YYYYMMDD' Defaults to today's date.
            
        Returns:
            A boolean (True/False) value.
        '''
                
        start = datetime.datetime.strptime(start,'%Y%m%d')
        end = datetime.datetime.strptime(end,'%Y%m%d')
        
        if scene.datetime > end:
            return True
        if scene.datetime < start:
            return True
        
        return False
                
    def testPolarisation(scene, pol):
        '''
        Function to test whether polarisation is available in the tile.
        
        Args:
            scene: Object from utilties.LoadScene()
            pol: 'VV' or 'VH'
        '''
        
        if scene.image_type == 'S1dual':
            return False
        if scene.image_type == 'S1single' and pol == 'VV':
            return False
        
        return True
    
    # Determine which images are within specified tile bounds
    if verbose: print 'Searching for source files within specified tile...'
    
    do_tile = []
    
    for source_file in source_files:
                        
        # Skip processing the file if image falls outside of tile area
        if testOutsideTile(source_file.metadata, md_dest):
            do_tile.append(False)
            continue
        
        if testOutsideDate(source_file, start, end):
            do_tile.append(False)
            continue
        
        if testPolarisation(source_file, pol):
            do_tile.append(False)
            continue
        
        if verbose: print '    Found one: %s'%source_file.filename.split('/')[-1]
        do_tile.append(True)
    
    # Get subset of source_files in specified tile
    source_files_tile = list(np.array(source_files)[np.array(do_tile)])
    
    return source_files_tile


def sortScenes(scenes):
    '''
    Function to sort a list of scenes by tile, then by date. This reduces some artefacts in mosaics.
    
    Args:
        scenes: A list of utilitites.LoadScene() Sentinel-2 objects
    Returns:
        A sorted list of scenes
    '''
    
    scenes_out = []
    
    scenes = np.array(scenes)
    
    dates = np.array([scene.datetime for scene in scenes])
    
    for date in np.unique(dates):
        scenes_out.extend(scenes[dates == date].tolist())
    
    return scenes_out
  


if __name__ == '__main__':
    '''
    '''
    
    import argparse
    
    # Set up command line parser
    parser = argparse.ArgumentParser(description = "This file contains functions to assist in the mosaicking and masking of Sentinel-1 data. A command line interface for image mosaicking is provided in mosaic.py.")
    
    args = parser.parse_args()
