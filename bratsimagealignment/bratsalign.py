#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created by Jeremy. J. Harwood (2019)
# Contact: Jeremy.Harwood@physics.org or J.Harwood3@herts.ac.uk

# Please cite Harwood, Vernstrom & Stroe 2019, MNRAS 491 803 if you have made use of this script.

# Automatically fits Gaussians and aligns images in pixel space. Particularly useful for spectral index and spectral age fitting on resolved scales.
# Change the var_ variables at the top of the script to suit your setup
# The script should be run in CASA using either exec() or the -C command line argument

import os, glob, scipy, sys
import numpy as np
#from astropy.io import fits


class _Constants:
    _MAX_MODES = 2
    _REFERENCE_ROUNDING_ACCURACY = 2

class Setup(_Constants):

    def __init__(self, input_files:list, output_dir:str, region_files:list, mode=2, reference_image=0, reference_location=[1024.0, 1024.0]):
        _Setters._set_input(self, input_files)
        _Setters._set_output(self, output_dir)
        _Setters._set_regions(self, region_files)
        _Setters._set_mode(mode)
        _Setters._set_reference_image(reference_image)
        _Setters._set_reference_location(reference_location)

    def set_mode(self, mode:int):
        _Setters._set_mode(mode)

    def set_reference_image(self, mode:int):
        _Setters._set_reference_image(reference_image)

    def set_reference_location(self, reference_location:int):
        _Setters._set_reference_location(reference_location)

# Setter class for setting and updating parameters
class _Setters(Setup):

    #def __init__(self):
    #    super(__Setters, self).__init__()

    def _set_input(self, __input_files):
        if len(__input_files) == 0:
            raise IndexError('Input files list cannot be empty!')
        else:
            for file in __input_files:
                if not file.strip():
                   raise ValueError('The input files list contains a blank element. All elements must contain a valid value.')
                elif len(glob.glob(file)) == 0:
                    raise FileNotFoundError('Unable to locate any valid files matching %s.' % (file))
                elif file.rfind('/') == -1:
                    var_input[i] = './' + var_input[i] # Standardise the path format to make maniulation safer and easier to manage

            self._input = [__input_files]

    def _set_output(self, __output_directory):
        if not __output_directory:
            raise EOFError('Output directory string cannot be empty!')
        elif not os.path.isdir(__output_directory):
            raise NotADirectoryError('The output directory %s does not exist. Please ensure the path is correct.' % (file)
        elif __output_directory[len(__output_directory)-1] != '/':   # Standardise the path format to make maniulation safer and easier to manage
            __output_directory += '/'

        self._output = [__output_directory]

    def _set_regions(self, __region_files):
        if len(__region_files) == 1:
            self._regions = [__region_files[0]] * len(self._input)
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

    def _set_mode(__mode):
        if not 0 <= __mode < self._MAX_MODES:
            raise ValueError('Modes must be in the range 0 and %i' % (self._MAX_MODES))
        self._mode = __mode

    def _set_reference_image(__reference_image):
        if not 0 <= __reference_image < len(self._input):
            raise IndexError('Reference image must be an index between 0 and %d' % (len(self._input)-1))
        self._reference_image = __reference_image

    def _set_reference_location(__reference_location):
        if not 0 <= len(__reference_location) < 2:
            raise IndexError('Reference location must be a list containing 2 element (currently %d)' % (len(__reference_location)))
        self._reference_location = __reference_location

class _core(Setup):

    def _align(self):

        if self._mode == 0:

            print('Determining the reference coordinates (predefined image)...')
            __reference_fitting_file =  glob.glob(self._input[_reference_image])[0]
            __gaussfit=imfit(imagename=__reference_fitting_file, region=self._regions[_reference_image], dooff=True)
            __reference_x, __reference_y = dict_gaussfit['results']['component0']['pixelcoords'].round(_REFERENCE_ROUNDING_ACCURACY)
            print('Aligning to reference coordinates: {:.2f}'.format(__reference_x) + ', {:.2f}'.format(__reference_y))


################################################################################################################
########                                                                                                 #######
########     Do not change anything past this point unless you are sure you know what you are doing!     #######
########                                                                                                 #######
################################################################################################################


#imgnum = 0

#if var_reftype == 0:

#    print('Determining the reference coordinates (predefined image)...')
    
#    for file in glob.glob(var_input[0]):
#        if imgnum == var_refimage:
#            dict_gaussfit=imfit(imagename=file, region= var_region[0],dooff=True)
#            par_ref = dict_gaussfit['results']['component0']['pixelcoords'].round(2)
#            print('Aligning to reference coordinates: {:.2f}'.format(par_ref[0]) + ', {:.2f}'.format(par_ref[1]))

#        imgnum+=1

#elif var_reftype == 1:
#    print('Aligning to a fixed reference of: ' + str(var_ref[0]) + ', ' + str(var_ref[0]))
#    par_ref = var_ref

#elif var_reftype == 2:
#    print('Determining the reference coordinates (mean)...')

#    tmp_totalx = 0
#    tmp_totaly = 0
    
#    for file in glob.glob(var_input[0]):
#        dict_gaussfit=imfit(imagename=file, region= var_region[0], dooff=True)
#        tmp_totalx += dict_gaussfit['results']['component0']['pixelcoords'][0]
#        tmp_totaly += dict_gaussfit['results']['component0']['pixelcoords'][1]

#        imgnum+=1

#    par_ref = [(tmp_totalx/imgnum), (tmp_totaly/imgnum)]
        
#    print('Aligning to reference coordinates: {:.2f}'.format(par_ref[0]) + ', {:.2f}'.format(par_ref[1]))


#else:
#    sys.exit('Error: Unknown reference type of ' + str(var_reftype) + '. Please check var_reftype and try again. Exiting...')


#for j in range(len(var_input)):
    
#    imgnum = 0
#    res_mean_offsetx = 0
#    res_mean_offsety = 0
#    output_list= []

#    print('Performing Gaussian fitting and shifting images to the reference pixel...')

#    for file in glob.glob(var_input[j]):

#        #Check the file has a valid extension associated with it
#        if var_input[j].rfind('.') == -1:
#            sys.exit('Error: The image ' + file + 'appears to have no file extension associated with it (e.g. .fits). This can cause bad thing to happen. Please check the files names and try again. Exiting...')
        
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


#    print('Performing second Gaussian fitting to check alignment...')

#    # Check if the alignment worked
#    res_mean_offsetx = 0
#    res_mean_offsety = 0

#    #print('*** Aligned values ***')
#    for file in output_list:

#        #Fit a gaussian and determine the offset
#        dict_gaussfit=imfit(imagename=file, region= var_region[j],dooff=True)
#        par_xoffset = par_ref[0] - dict_gaussfit['results']['component0']['pixelcoords'][0]
#        par_yoffset = par_ref[1] - dict_gaussfit['results']['component0']['pixelcoords'][1]

#        res_mean_offsetx += abs(par_xoffset)
#        res_mean_offsety += abs(par_yoffset)
    
#        #print(file, str(par_xoffset), str(par_yoffset))

#    res_mean_offsetx /= imgnum
#    res_mean_offsety /= imgnum
    
#    print('Mean offset (after shift) {:.3g}:'.format(res_mean_offsetx) + ', {:.3g}'.format(res_mean_offsety))