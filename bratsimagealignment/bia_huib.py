#!/usr/bin/python3
# -*- coding: utf-8 -*-

# List of raised exceptions:
# ValueError
# FileNotFoundError
# IndexError
# NotADirectoryError
# ValueError


# Modes
# 0 = Predefined image
# 1 = Fixed reference
# 2 = Mean 

# Reference image
# Mode 0 = Not used
# Mode 1 = The index of the input list to use a reference. If the index points to more than one file e.g. ./*.fits the first file will be used, ordered alphabetically. Default 0.
# Mode 2 = The index of the input list for which the mean should be calculated. To use all index use -1. Default 0.

import traceback
import bratsalign as bia

# List of paths and file names of the images to align. Use '*' for all files. First element of list used for reference pixel determination.
var_input = ['./huib_unaligned_1/*.FITS', './huib_unaligned_2/*.FITS']
var_output = './huib_aligned' # Path to directory for the output images
var_region = ['huib_1.crtf', 'huib_2.crtf']

var_reftype = 0 # Choose if we should a reference point set by a predefined image (0), a fixed value (1), or the mean (2) 
var_refimage = 0 # Choose which image to align to. Only used if var_reftype = 0
var_ref = [1024.0, 1024.0] # List in [x, y] coodinates to align to. Only used if var_fixedref = 1

try:
    #example = bia.Setup(var_input, var_output, var_region, overwrite_files=True, mode=0)
    #example = bia.Setup(var_input, var_output, var_region, overwrite_files=True, mode=1, reference_image=0)
    example = bia.Setup(var_input, var_region, var_output, overwrite_files=True, mode=var_reftype, reference_image=0)
    example.align()
    #bia.Setup().align()
except Exception as e:
    traceback.print_exc()
    print(str(e))
