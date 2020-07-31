#subprocess.run([self._casa_path, "--nologger", "--nogui", "-c", self._CASA_COMMANDS_FILE_NAME, ">> " + self._output + "brats_align_casa.log"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)

            
################################################################################################################
########                                                                                                 #######
########     Do not change anything past this point unless you are sure you know what you are doing!     #######
########                                                                                                 #######
################################################################################################################


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