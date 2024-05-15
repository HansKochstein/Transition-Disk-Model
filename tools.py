import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.colors import LogNorm
import os
import re
import numpy as np


''' Redefining inner Grid '''


def grid_refine_inner_edge(x_orig,nlev,nspan):
    x     = x_orig.copy()
    rev   = x[0]>x[1]
    for ilev in range(nlev):
        x_new = 0.5 * ( x[1:nspan+1] + x[:nspan] )
        x_ref = np.hstack((x,x_new))
        x_ref.sort()
        x     = x_ref
        if rev:
            x = x[::-1]
    return x

def manage_image_files(total_frames, delete_files=False):

    directory = "/Users/pbo/Desktop/Bosse_Internship_Project/radmc_for_ppdisk_English_ver1.0/inc_scattering/movie"
    
    # Delete the specified range of image files if delete_files is True
    if delete_files:
        frame_numbers = [f"image_{i:04d}.out" for i in range(1, total_frames + 1)]
        for filename in frame_numbers:
            full_path = os.path.join(directory, filename)
            if os.path.exists(full_path):  # Check if the file exists
                os.remove(full_path)       # Delete the file
                print(f"Deleted {directory}")  # Optional: print a confirmation message

    # Check for any remaining image files outside the specified range
    all_files = os.listdir(directory)
    pattern = re.compile(r"image_\d{4}\.out$")
    remaining_files = [file for file in all_files if pattern.match(file)]

    if remaining_files:
        print('Image files outside the specified range have been found:')
        for file in remaining_files:
            print(file)
    else:
        print('All specified image-files successfully deleted, no additional files found.')


import numpy as np

def circular_mask(data_array, radius):

    new_data_array = data_array.copy()

    array_size = data_array.shape
    center = (array_size[0] // 2, array_size[1] // 2)

    # Create the mask efficiently
    x_grid, y_grid = np.ogrid[:array_size[0], :array_size[1]]
    distance_from_center = np.sqrt((x_grid - center[0])**2 + (y_grid - center[1])**2)
    circmask = distance_from_center > radius

    # Apply the mask
    new_data_array[~circmask] = data_array.min()  # Set values inside the circle to 0

    return new_data_array


def trace_outer_ring(intensity_data, degrees_step=10, ax=None, grid_x=None, grid_y=None, int_threshold=0.1, color = 'white'):
    center_y, center_x = intensity_data.shape[0] // 2, intensity_data.shape[1] // 2
    angles = np.arange(0, 360, degrees_step)

    max_intensity_points = []
    for angle in angles:
        radians = np.deg2rad(angle)
        y_values = np.round(center_y + np.arange(center_y) * np.sin(radians)).astype(int)
        x_values = np.round(center_x + np.arange(center_x) * np.cos(radians)).astype(int)

        line_intensities = intensity_data[y_values, x_values]

        # Find indices where intensity is above the threshold
        valid_indices = np.where(line_intensities > int_threshold)[0]

        if valid_indices.size > 0:
            # Get the index of the maximum intensity within the valid range
            max_intensity_index = valid_indices[np.argmax(line_intensities[valid_indices])]
            max_y, max_x = y_values[max_intensity_index], x_values[max_intensity_index]
            max_intensity_points.append((max_y, max_x))

    max_y, max_x = zip(*max_intensity_points)
    if ax==None:
        plt.scatter([(x-100)*1.23 for x in max_x], [(y-100)*1.23 for y in max_y], 
                                 color=color, s=30, marker='+')
    else:
        ax[grid_x, grid_y].scatter([(x-100)*1.23 for x in max_x], [(y-100)*1.23 for y in max_y], 
                                 color=color, s=30, marker='+')
    return max_intensity_points

# au  = 1.49598e13     # Astronomical Unit       [cm]

# dpc = 138

# frame_numbers = [f"movie/image_{i:04d}.out" for i in range(1, total_frames + 1)] 
# # Set up the figure and axis for the animation
# fig, ax = plt.subplots()  # Simplified subplot creation

# def update_frame(i):
#     ax.clear()  # Clear the previous frame
#     fname = frame_numbers[i]
#     imrangeau      =  [-100,100,-100,100]
#     with open(fname,'r') as f:
#         s        = f.readline()
#         im_n      = f.readline()
#         im_nx     = int(im_n.split()[0])
#         im_ny     = int(im_n.split()[1])
#         nlam      = int(f.readline())
#         pixsize   = f.readline()
#         pixsize_x = float(pixsize.split()[0])
#         pixsize_y = float(pixsize.split()[1])
#         wavelength_mic = np.zeros(nlam)
#         im    = np.zeros((im_nx,im_ny,4,nlam)) #[MODIFIED] array dimension modified for Stokes I,Q,U,V
#         rr_au = np.zeros((im_nx,im_ny,4,nlam)) #[MODIFIED] array dimension modified for Stokes I,Q,U,V
#         phi   = np.zeros((im_nx,im_ny,4,nlam)) #[MODIFIED] array dimension modified for Stokes I,Q,U,V
#         for ilam in range(nlam): wavelength_mic[ilam] = f.readline()
#         for ilam in range(nlam):
#             dummy = f.readline() #blank line
#             for iy in range(im_ny):
#                 for ix in range(im_nx):
#                     image_temp = f.readline()
#                     im[ix,iy,:,ilam] = [image_temp.split()[istokes] for istokes in [0,1,2,3]] #[MODIFIED] array dimension modified for Stokes I,Q,U,V
#                     rr_au[ix,iy,ilam] = np.sqrt(((ix-im_nx/2)*pixsize_x)**2+((iy-im_ny/2)*pixsize_y)**2)/au
#                     phi[ix,iy,ilam]   = np.arctan2((iy-im_ny/2)*pixsize_y,(ix-im_nx/2)*pixsize_x) 
#     sizex_au = im_nx*pixsize_x/au
#     sizey_au = im_ny*pixsize_y/au

#     #polarized intensity is stored in the array im_PI
#     im_PI = np.zeros((im_nx,im_ny,nlam,4))
#     im_PI = ((im[:,:,1,0])**2.0 + (im[:,:,2,0])**2.0)**0.5

#     d_pc = 138.0

#     conv          = pixsize_x * pixsize_y / (d_pc*pc)**2. * 1e23 
#     im_PI_forplot = im_PI.copy()*conv
#     pltunit       = '[Jy/pix]'

#     # fontsize       = 12
#     # fontfamily     = 'serif'
#     interpolation  = 'bicubic'
#     cmap           = 'magma'
#     intmax         =  np.max(im_PI_forplot)
#     intmin         =  np.max(im_PI_forplot)*0.001
#     tickcolor      = 'white'
#     xlabel         = 'x [AU]'
#     ylabel         = 'y [AU]'
#     norm           =  None
#     # drawpolvect    =  False
#     title          = f'Polarized Intensity of {float(wavelength_mic[0])}$\mu m$' 

#     # fig2=plt.figure(num=None, figsize=(8,6), dpi=200, facecolor='w', edgecolor='k')
#     # ax2=fig2.add_subplot(111)
#     norm=LogNorm(vmin=intmin, vmax=intmax)
#     plt.imshow(im_PI_forplot[:,:].T,interpolation=interpolation,cmap=cmap,norm=norm,origin='lower',extent= [-sizex_au/2.0,sizex_au/2.0,-sizey_au/2.0,sizey_au/2.0])
#     plt.xlabel(xlabel)
#     plt.ylabel(ylabel)
#     plt.xlim(imrangeau[0],imrangeau[1])
#     plt.ylim(imrangeau[2],imrangeau[3])
#     plt.title(title)
#     # plt.tick_params(which='major', color=tickcolor,width=2.0,length=6.0)
#     # plt.text(x=30, y=-90, s=f"Inclination: {incl_value}°, Position Angle: {posang_value}°", fontsize=6, bbox=dict(facecolor='white', alpha=1))
#     #cbar=plt.colorbar()
#     #cbar.set_label(pltunit)

#     return fig,


'''PARAMETER FOR TWO OVERLAYING OUTER CRESENDENT SHAPES'''

# nr       = 120                       # Nr of radial grid points
# nphi     = 100                       # Nr of azimuthual grid points 
# ntheta   = 100                       # Nr of polar grid points

# rin      = 1*au                     # Inner radius
# rout     = 120*au                   # Outer radius

# #INNER DISK
# sigmad0   = 0.01                     # Sigma dust (in the middle of the Ring)  [*0.05 assuming a higher surface density]
# ringcau1  = 6                       # R_{ring, center}: the radial location of the ring center in [au]
# ringwau1  = 1                     # R_{ring, width}: the radial width of the ring in [au]
# elevD     = 50                      # inner DISK elevation in DEGREES
# rotaD1    = 70                     # inner DISK rotation in DEGREES (azimuthal rotation)

# #TESTING
# INNERASSYMETRIC = False
# azipos1   = pi                      # inner DISK rotation of cresendent structur in RAD
# azimwid1  = 0.45                    # inner DISK azimuthual width of cresendent structure in RAD

# #OUTER DISK
# sigmadc2  = 0.6                      # Sigma_{dust, center}: dust surface density at the center of the ring [g/cm^2]
# ringcau2  = 88                      # R_{ring, center}: the radial location of the ring center in [au]
# ringwau2  = 8.5                       # R_{ring, width}: the radial width of the ring in [au]
# rotaD2    = 100                      # outer DISK rotation in (azimuthal rotation) & outer DISK rotation of cresendent structur in DEG

# OUTERASSYMETRIC = True
# azimwid2  = 0.39                    # outer DISK azimuthual width of cresendent structure in RAD

# ELLIPSE = True
# a_semi  = ringcau2 * au             # [cm] semi-major axis
# e       = 0.27                      # eccentricity


'''OFFSET'''
# OFFSET = False
# offset_pos = [18*au, np.deg2rad(90-rotaD2)]     # offset position (r0, gamma) | [au, deg] polar coordinates, [x0,y0]-> r0=(x0^2+y0^2)**1/2 -> 
# offset_posxy= [18*au,7*au]                # offset position (x,y) -> (y,x-) | [au, au] cartesian  
# calculating new postion of ring with respect to offsetted origin 
# ringc2_offset = offset_pos[0]*np.cos(pp-offset_pos[1])+np.sqrt(ringc2**2 - offset_pos[0]**2 * np.sin(pp-offset_pos[1])**2)      #polar
# ringc2_offset = offset_posxy[0]*np.cos(pp)+offset_posxy[1]*np.sin(pp)+np.sqrt(ringc2**2-(offset_posxy[0]*np.sin(pp)-offset_posxy[1]*np.cos(pp))**2) #cartesian

# if OFFSET == True:
#     ringc2_mod = offset_pos[0]*np.cos(pp-offset_pos[1])+np.sqrt(ringc2**2 - offset_pos[0]**2 * np.sin(pp-offset_pos[1])**2)      #polar
# else:
#     ringc2_mod = ringc2

''' SECOND RING '''

# if ELLIPSE == True:
#     ringc3_mod = (85*au)*(1-e**2)/(1-e*np.cos(pp))
# else:
#     ringc3_mod = ringc2


# sigmad3    = 0.005  * np.exp(-((newr-ringc3_mod)**2/(11*au)**2)/2.0) * np.exp(-((pp-pi)**2/(2)**2)/2.0)          # overlay of 2nd outer ring  thing ring 

# rhod3     = ( sigmad3 / (np.sqrt(2.e0*np.pi)*hp2) ) * np.exp(-(zr**2/hpr2**2)/2.e0)                                # volume density

# rhod3     = np.roll(rhod3, rotation2, axis=2)
# rhod     += rhod3  
