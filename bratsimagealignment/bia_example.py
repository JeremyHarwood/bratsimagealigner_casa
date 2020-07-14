#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import bratsalign as bia

# List of paths and file names of the images to align. Use '*' for all files. First element of list used for reference pixel determination.
var_input = ['./path/to/some/fits/files/*.fits',  './path/to/some/more/fits/files/*.fits', './path/to/some/even/more/fits/files/*.fits'] 
#var_input = ['./path/to/some/fits/files/*.fits',  './path/to/some/more/fits/files/*.fits'] 
var_output = './path/to/out/folder/' # Path to directory for the output images
#var_region = ['aregionfile.crtf', './pathtoanother/regionfile.crtf', 'aregionfile.crtf'] # List of regions files. Casa format, one per input (can be the same, or different if significantly shifted)
var_region = ['aregionfile.crtf', './pathtoanother/regionfile.crtf'] # List of regions files. Casa format, one per input (can be the same, or different if significantly shifted)

var_reftype = 2 # Choose if we should a reference point set by a predefined image (0), a fixed value (1), or the mean (2) 
var_refimage = 4 # Choose which image to align to. Only used if var_fixedref = 0
var_ref = [1024.0, 1024.0] # List in [x, y] coodinates to align to. Only used if var_fixedref = 1


bia.Setup(var_input, var_output, var_region)