#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This temporary script was created by the BRATS alignment tool and can be safely deleted once a run has completed.

# Created by Jeremy. J. Harwood (2019)
# Contact: Jeremy.Harwood@physics.org or J.Harwood3@herts.ac.uk

# Please cite Harwood, Vernstrom & Stroe 2019, MNRAS 491 803 if you have made use of this script.

__input = [['./testimages\\J1206_503_L_AB_final_deep_SPW16_2048px.fits', './testimages\\J1206_503_L_AB_final_deep_SPW20_2048px.fits', './testimages\\J1206_503_L_AB_final_deep_SPW21_2048px.fits', './testimages\\J1206_503_L_AB_final_deep_SPW22_2048px.fits', './testimages\\J1206_503_L_AB_final_deep_SPW23_2048px.fits', './testimages\\J1206_503_L_AB_final_deep_SPW26_2048px.fits', './testimages\\J1206_503_L_AB_final_deep_SPW27_2048px.fits', './testimages\\J1206_503_L_AB_final_deep_SPW28_2048px.fits', './testimages\\J1206_503_L_AB_final_deep_SPW29_2048px.fits', './testimages\\J1206_503_L_AB_final_deep_SPW30_2048px.fits', './testimages\\J1206_503_L_AB_final_deep_SPW31_2048px.fits']]
__sum_x = 0
__sum_y = 0
__sum_images = 0
__sum_images+=sizeof(__input[0])
for file in __input[0]:
	__gauss_fit = imfit(imagename=file, region=var_region[0], dooff=True)
	__sum_x += __gauss_fit['results']['component0']['pixelcoords'][0])
	__sum_y += __gauss_fit['results']['component0']['pixelcoords'][1])
__reference_x = (__sum_x/__sum_images)
__reference_y = (__sum_y/__sum_images)
print('Aligning to reference coordinates: {:.2f}'.format(__reference_x) + ', {:.2f}'.format(__reference_y))
