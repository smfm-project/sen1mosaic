import glob
from setuptools import setup

setup(name='sen1mosaic',
      packages = ['sen1mosaic'],
      version='0.2',
      data_files=[('./cfg/',glob.glob('./cfg/*'))],
      description='Tools to generate mosaics of Sentinel-1 data.',
      url='https://bitbucket.org/sambowers/sen1mosaic',
      author='Samuel Bowers',
      author_email='sam.bowers@ed.ac.uk',
      license='GNU General Public License',
      zip_safe=False)

