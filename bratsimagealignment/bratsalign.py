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


class Setup:

    def __init__(self, input_files: list, output_dir: str, region_files:list):
        try :
            self.__set_input(input_files)
            self.__set_output(output_dir)
            self.__set_regions(region_files)
            self.mode = 2
            self.ref_image = 0
            self.ref_location = [1024.0, 1024.0]
        except ValueError as e:
            print('\n*** ValueError: ' + str(e) + ' ***\n')
        except FileNotFoundError as e:
            print('\n*** FileNotFoundError: ' + str(e) + ' ***\n')
        except EOFError as e:
            print('\n*** EOFError: ' + str(e) + ' ***\n')
        except NotADirectoryError as e:
            print('\n*** NotADirectoryError: ' + str(e) + ' ***\n')
        except (KeyboardInterrupt, SystemExit):
            print('\n*** KeyboardInterrupt: User interupt detected. Exiting the program.... ***\n')
            raise
        except Exception as e:
            print('\n*** Exception: ' + str(e) + ' ***\n')

    # Error check and set the input file list
    def __set_input(self, __input_files):
        if len(__input_files) == 0:
            raise ValueError('Input files list cannot be empty!')
        else:
            for file in __input_files:
               if not file.strip():
                   raise ValueError('The input files list contains a blank element. All elements must contain a valid value.')
               elif not os.path.isfile(file):
                   raise FileNotFoundError('The input file %s does not exist. Please ensure all paths are correct.' % (file))
            self._input = [__input_files]

    def __set_output(self, __output_directory):
        if not __output_directory:
            raise EOFError('Output directory string cannot be empty!')
        else:
            if not os.path.isdir(__output_directory):
                   raise NotADirectoryError('The output directory %s does not exist. Please ensure the path is correct.' % (file))
            self._output = [__output_directory]

    # Error check and set the region files list
    def __set_regions(self, __region_files):
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
            raise ValueError('Region files list cannot be empty!')
        else:
            raise ValueError('Region files list must be equal to either the number of images or to one [Regions: %i Images: %i]' % (len(__region_files), len(self._input)) )

    
    ### TODO ###

    #Set mode

    #Set reference image

    #Set reference location



    #def execfile(self, file_loc: str) -> dict:
    #    return _Script._execfile(self, file_loc)

    #def multiexec(self, file_list: list, num_proc: int) -> dict:
    #    return _Multiexec._multiexec(self, file_list, num_proc)



#class align(Setup):

#    if len(var_input) == len(var_region):

################################################################################################################
########                                                                                                 #######
########     Do not change anything past this point unless you are sure you know what you are doing!     #######
########                                                                                                 #######
################################################################################################################


## Do some error checking and correction before we begin
#if len(var_input) != len(var_region):
#    sys.exit('Error: The list length of var_input and var_region do not match. Please check and try again. Exiting...')

#if os.path.isdir(var_output):
#    print('Sucessfully found output directory...')
    
#else:
#    try:
#        print('Output directory not found. Attempting to create it...')
#        os.makedirs(var_output)
#        print('Success!')
#    except:
#        traceback.print_exc()
#        sys.exit('Error: Unable to find or create the output directory. Please check permission and try again. Exiting...')

#for i in range(len(var_input)):
#    if var_input[i].rfind('/') == -1:
#        var_input[i] = './' + var_input[i]
#        print('Adjusting input')

#if var_output[len(var_output)-1] != '/':   
#    var_output = var_output + '/'
#    print('Adjusting output')

## Check all the input directories contain valid files
#for file in var_input:
#    if len(glob.glob(file)) == 0:
#        sys.exit('Error: Unable to locate any valid file in ' + file + '. Please check var_input and try again. Exiting...')

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

