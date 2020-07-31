#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created by Jeremy. J. Harwood (2019)
# Contact: Jeremy.Harwood@physics.org or J.Harwood3@herts.ac.uk

# Please cite Harwood, Vernstrom & Stroe 2019, MNRAS 491 803 if you have made use of this class.

# Automatically fits Gaussians and aligns images in pixel space. Particularly useful for spectral index and spectral age fitting on resolved scales.
# Change the var_ variables at the top of the script to suit your setup
# The script should be run in CASA using either exec() or the -C command line argument

import os, glob, scipy, sys, subprocess
import numpy as np
from astropy.io import fits
#from astropy.wcs import utils
#import regions


class _Constants:
    _VERSION = "v1.0.0"
    _MAX_MODES = 2
    _REFERENCE_ROUNDING_ACCURACY = 2
    _CASA_COMMANDS_FILE_NAME = "brats_casa_commands.py"
    _DEFAULT_OUTPUT_PATH = "./"
    _DEFAULT_CASA_PATH = "/soft/casapy/bin/casa"
    _CASA_HEADER = '''\
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This temporary script was created by the BRATS alignment tool and can be safely deleted once a run has completed.

# Created by Jeremy. J. Harwood (2019)
# Contact: Jeremy.Harwood@physics.org or J.Harwood3@herts.ac.uk

# Please cite Harwood, Vernstrom & Stroe 2019, MNRAS 491 803 if you have made use of this script.'''


class Setup(_Constants):

    def __init__(self, input_files:list, region_files:list, output_dir=_Constants._DEFAULT_OUTPUT_PATH, mode=2, reference_image=0, reference_location=[1024.0, 1024.0], casa_path=_Constants._DEFAULT_CASA_PATH, overwrite_files=False):
        _Setters._set_input(self, input_files)
        _Setters._set_regions(self, region_files)
        _Setters._set_output(self, output_dir)
        _Setters._set_mode(self,mode)
        _Setters._set_reference_image(self,reference_image)
        _Setters._set_reference_location(self,reference_location)
        _Setters._set_casa_path(self,casa_path)
        self._overwrite = overwrite_files

    def set_mode(self, mode:int):
        _Setters._set_mode(self, mode)

    def set_reference_image(self, mode:int):
        _Setters._set_reference_image(self, reference_image)

    def set_reference_location(self, reference_location:int):
        _Setters._set_reference_location(self, reference_location)

    def set_casa_path(self, casa_path:str):
        self._casa_path = _set_casa_path(self, casa_path)

    def overwrite(self, overwrite_files:bool):
        self._overwrite = overwrite_files
    
    def align(self):
        _Core._align(self)
  

# Setter class for setting and updating parameters
class _Setters(Setup):

    def _set_input(self, __input_files):
        if len(__input_files) == 0:
            raise IndexError('Input files list cannot be empty!')
        else:
            self._input = []
            for path in __input_files:
                if not path.strip():
                   raise ValueError('The input files list contains a blank element. All elements must contain a valid value.')
                elif len(glob.glob(path)) == 0:
                    raise FileNotFoundError('Unable to locate any valid files matching %s.' % (file))
                else:
                    self._input.append(sorted(glob.glob(path)))

    def _set_output(self, __output_directory):
        if not __output_directory:
            raise EOFError('Output directory string cannot be empty!')
        elif not os.path.isdir(__output_directory):
            print('The output directory %s does not exist. Attempting to create it...' % (__output_directory))
            os.mkdir(__output_directory)
            #raise NotADirectoryError('The output directory %s does not exist. Please ensure the path is correct.' % (__output_directory))

        if __output_directory[len(__output_directory)-1] != '/':   # Standardise the path format to make maniulation safer and easier to manage
            __output_directory += '/'

        self._output = __output_directory

    def _set_regions(self, __region_files):
        if len(__region_files) == 1:
            if not __region_files[0].strip():
                   raise ValueError('The region files list contains a blank element. All elements must contain a valid value.')
            elif not os.path.isfile(__region_files[0]):
                   raise FileNotFoundError('The region file %s does not exist. Please ensure all paths are correct.' % (__region_files[0]))
            self._regions = [__region_files[0] * len(self._input)]
        elif len(__region_files) == len(self._input):
            for file in __region_files:
               if not file.strip():
                   raise ValueError('The region files list contains a blank element. All elements must contain a valid value.')
               elif not os.path.isfile(file):
                   raise FileNotFoundError('The region file %s does not exist. Please ensure all paths are correct.' % (file))
            self._regions = [__region_files]
        elif len(__region_files) == 0:
            raise IndexError('Region files list cannot be empty!')
        else:
            raise IndexError('Region files list must be equal to either the number of images or to one [Regions: %i Images: %i]' % (len(__region_files), len(self._input)) )

    def _set_mode(self, __mode):
        if not 0 <= __mode <= self._MAX_MODES:
            raise ValueError('Modes must be in the range 0 and %i' % (self._MAX_MODES))
        self._mode = __mode

    def _set_reference_image(self,__reference_image):
        if self._mode == 2 and not -1 <= __reference_image < len(self._input):
            raise IndexError('Reference image must be an index between -1 and %d' % (len(self._input)-1))
        elif not 0 <= __reference_image < len(self._input):
            raise IndexError('Reference image must be an index between 0 and %d' % (len(self._input)-1))
        else:
            self._reference_image = __reference_image

    def _set_reference_location(self,__reference_location):
        if not len(__reference_location) == 2:
            raise IndexError('Reference location must be a list containing 2 element (currently %d)' % (len(__reference_location)))
        self._reference_location = __reference_location

    def  _set_casa_path(self, __casa_path): # This needs error checking!
        self._casa_path = __casa_path


class _Core(Setup):

    def _align (self):
        self.__commands = []
        _Core._check_mode(self)
        _Core._define_parameters(self)
        _Core._define_reference(self)
        #_Core._define_align(self)

        _Core._casa_write(self)

        __casa_log_file = " > " + self._output + "/brats_align_casa.log"
        #subprocess.run([self._casa_path, "--nologger", "--nogui", "-c", self._CASA_COMMANDS_FILE_NAME, ">> " + self._output + "brats_align_casa.log"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        #subprocess.run([self._casa_path, "--nologger ", "--nogui ", "-c ", self._CASA_COMMANDS_FILE_NAME])
        print(self._casa_path + " --nologger --nogui -c " + self._CASA_COMMANDS_FILE_NAME, __casa_log_file)
        subprocess.run([self._casa_path, " --nologger ", "--nogui ", "-c ", self._CASA_COMMANDS_FILE_NAME, __casa_log_file])

    def _check_mode(self):
        if self._mode == 0 and self._reference_image == -1:
            raise IndexError('A reference image of -1 cannot be used for mode 1 (predefined image)')

    def _define_parameters(self):
        self.__commands.append("__input = " + str(self._input))

    def _define_reference(self):
        
        if self._mode == 0:
            print('Using a predefined image for the reference coordinates')
            __reference_fitting_file = self._input[self._reference_image][0]
            self.__commands.append("__gauss_fit = imfit(imagename=\'" + __reference_fitting_file + "\', region=\'" + str(self._regions[self._reference_image]) + "\', logfile=\'" + self._output + "tmp_reference_coords.txt\', dooff=True, append=False)")
            self.__commands.append("__reference_x, __reference_y = __gauss_fit['results']['component0']['pixelcoords'].round(_REFERENCE_ROUNDING_ACCURACY")           
        elif self._mode == 1:
            print('Using a fixed reference of {:.2f}'.format(_reference_location[0]) + ', {:.2f}'.format(_reference_location[1]) + " for alignment")
            self.__commands.append("__reference_x, __reference_y = " + str(_reference_location))
        elif self._mode == 2:
            print('Using the mean value of all images for the reference coordinates')
            self.__commands.extend(["__sum_x = 0", "__sum_y = 0"])
            self.__commands.append("__sum_images = 0")

            if (self.set_reference_image == -1):
                self.__commands.append("for fits_location in __input:")
                self.__commands.append("\t__sum_images+=len(fits_location)")
                self.__commands.append("\tfor file in fits_location:")
                self.__commands.append("\t\t__gauss_fit = imfit(imagename=file, region=\'" + str(self._regions[self._reference_image]) + "\', dooff=True)")
                self.__commands.append("\t\t__sum_x += __gauss_fit[\'results\'][\'component0'][\'pixelcoords\'][0]")
                self.__commands.append("\t\t__sum_y += __gauss_fit[\'results\'][\'component0'][\'pixelcoords\'][1]")
            else:
                self.__commands.append("__sum_images+=len(__input[" + str(self._reference_image) + "])")
                self.__commands.append("for file in __input[" + str(self._reference_image) + "]:")
                self.__commands.append("\t__gauss_fit = imfit(imagename=file, region=\'" + str(self._regions[self._reference_image]) + "\', dooff=True)")
                self.__commands.append("\t__sum_x += __gauss_fit[\'results\'][\'component0\'][\'pixelcoords\'][0]")
                self.__commands.append("\t__sum_y += __gauss_fit[\'results\'][\'component0\'][\'pixelcoords\'][1]")
                
            self.__commands.extend(["__reference_x = (__sum_x/__sum_images)", "__reference_y = (__sum_y/__sum_images)"])
    
        self.__commands.append("print(\'Aligning to reference coordinates: {:.2f}\'.format(__reference_x) + ', {:.2f}\'.format(__reference_y))")

    def _define_get_position(self) -> list:
            self.__commands.append("def _get_position(self, file_name, region_name):")
            self.__commands.append("\treturn imfit(imagename=file_name, region=region_name, dooff=True)")



    def _define_align(self):
        self.__commands.append("print(\'Performing Gaussian fitting and shifting images to the reference pixel...\')")
        self.__commands.append("for fits_location in __input:")
        self.__commands.append("\tfor file in fits_location:")
        self.__commands.append("\t\tif file.upper().endswith(\'.FITS\'):")
        self.__commands.append("\t\t\tprint(\'The file %s does not contain a valid FITS extension. Skipping...\' % (file))")
        self.__commands.append("\t\t\tcontinue")
        self.__commands.append("")
        self.__commands.append("")
        self.__commands.append("")
        self.__commands.append("")
        self.__commands.append("")
        self.__commands.append("")
        self.__commands.append("")
        self.__commands.append("")
        self.__commands.append("")
        self.__commands.append("")
        self.__commands.append("")
        self.__commands.append("")
        self.__commands.append("")
        self.__commands.append("")

#for j in range(len(var_input)):

    
#    imgnum = 0
#    res_mean_offsetx = 0
#    res_mean_offsety = 0
#    output_list= []

#    print('Performing Gaussian fitting and shifting images to the reference pixel...')

#    for file in glob.glob(var_input[j]):
        
#        #Fit a gaussian and determine the offset
#        dict_gaussfit=imfit(imagename=file, region= var_region[j],dooff=True)
#        par_xoffset = par_ref[0] - dict_gaussfit['results']['component0']['pixelcoords'][0]
#        par_yoffset = par_ref[1] - dict_gaussfit['results']['component0']['pixelcoords'][1]

#        #print(file, str(par_xoffset), str(par_yoffset))
    
#        res_mean_offsetx += abs(par_xoffset)
#        res_mean_offsety += abs(par_yoffset)

#        # Get the image data
#        image_data, image_hdr = fits.getdata(file, 0, header=True)

#        # See how many degenerate axis we have so we can add them back later
#        par_degaxis = image_data.ndim - 2

#        # We need to squeeze any degenerate axis. This should be a problem for general use but may in specific cases.
#        image_squeezed = image_data.squeeze()

#        # Do the shift
#        image_aligned = np.zeros(image_squeezed.shape)

#        # Note the x,y inversion here
#        scipy.ndimage.interpolation.shift(image_squeezed, np.array([par_yoffset, par_xoffset]), image_aligned)

#        #Add back any degenerate axis
#        for i in range(par_degaxis):
#            image_aligned = np.expand_dims(image_aligned, 0)
#        #print(image_aligned.shape)

#        # Export the new FITS file
#        #var_output_app = var_output + '.aligned'
#        extensionloc = file.rfind('.')
#        directoryloc = file.rfind('/')
#        var_output_app = var_output + file[directoryloc+1:extensionloc] + '_aligned' + file[extensionloc:]
#        output_list.append(var_output_app)

#        #image_out = fits.PrimaryHDU(image_aligned)
#        fits.writeto(filename=var_output_app, data=image_aligned, header=image_hdr, overwrite=True)

#        imgnum+=1

#    res_mean_offsetx /= imgnum
#    res_mean_offsety /= imgnum
    
#    print('Mean offset (before shift) {:.3g}:'.format(res_mean_offsetx) + ', {:.3g}'.format(res_mean_offsety))


    def _casa_write(self):
        
        if os.path.isfile(self._CASA_COMMANDS_FILE_NAME) and self._overwrite != True:
            raise FileExistsError('The temporary CASA command file %s already exists (overwrite_files=False).' % (self._CASA_COMMANDS_FILE_NAME))

        self.__commands_file = open(self._CASA_COMMANDS_FILE_NAME, "w")
        
        self.__commands_file.write(self._CASA_HEADER + "\n")

        for command in self.__commands:
            #print(command)
            self.__commands_file.write("\n" + command)

        self.__commands_file.close()


