""" 
**************************************NOTES****************************************

Script last updated 08/07/18

**************************************INSTRUCTIONS**********************************

This script is for running contrast curve automatically by loading in data we have
Planet subtraction is excluded as this will need a student to tell if the candiate/background star is 
fully subtracted or not. This can prevent over subtraction. 

star_info.txt = cent_x, cent_y, ref_x, ref_y, llsg_rank
name.txt = Folder name in Danger Zone 
name_file = file names for angle and list of frames
date.txt = epochs from spreadsheet

**************************************UPDATES**************************************

08/07/18
	1. 

01/06/17
	1. Option to break from loop when looping over different centre postions without
	   having to wait till all positions tried
	2. Stops users enter more six planets at planet injection stage, recomends use of 6
	
05/06/17
	1. Includes allowing for coordinate re-entry after looping over centre positions
	2. allow user to enter what multiple of background flux is used for planet brightness

06/06/17
	1. Clearer instructions and error checking when entering planet flux
	2. Background flux check now takes four 100x100 pixel grids from each image in the 
	   data cube. If the average flux in any of each of this grids is greater than 1 stdev
	   away from the population mean that grid is rejected and not used to build a mean
	   background flux
	3. Minor adjustment to Checkfuntion() allowing for more input options
	4. Changed Checkfunction to allow for exit to be called.	

07/06/17
	1. Addition of loop to allow planet injection phase to run mutliple times.
	2. Loop to redo entirly the re-centring if ds9 inspection shows problem
	3. Made if more obvious when you need to enter reference coordinates
	4. fixed bug in coordinate re-entry
	5. Made Image removal function more efficient and removed the necessity to know how many
	   images were going to be removed before entering numbers.
	6. Background flux now takes grids within 0.5 stdev of mean
08/06/17
	1. Removed Planet injection loop due to bug
	2. The .png files of each reduction were outputing blank files. No way was found to 
	   fix this therefore the code was removed. still saved as .fits but no .png.
	   This appears to be beacuse it is saving wrong figure but this is due to a VIP
	   hence no solution was found.
14/06/17
	1. Fixed the ininate loop when asking to re enter coodinates
	2. Fixed the re entered coordinates not being correctly over writtten
	   
**************************************************************************************   
	
"""

"""
Script designed to build a data cube from fits files with the aid of the VIP module
in python, reduce the images using the PCA and LLSG method of imaging reduction. 
Script also then calculates the flux from the exposure times of the fits files and 
then generates contrast curves based on the data cube produced and this value of flux. 
""" 
""" Input variables here """
psf_xy =[95,111]			#Centre of psf.
averageflux=0				#Initialises the background flux.				#Just a range of flux's for the first guess algorithm .
flvlup = 10000
flvldown = 1000				#to loop over.
#Initialises the loop counter, loop counter is used to prevent error warning on first trial
#when looping through integer check functions.
pxscale_keck=0.00953 		#Pixel scale for keck 2.
seperation=60 				#Radial seperation between synthetic planets in terms of pixels.
starphot=0					#Initialises the Starphot parametre.
sigma = 5					#Sets the value of sigma for the contrast curves.
""" End of variables """ 

"""
All packages used in this script 
"""
from matplotlib.pyplot import  *
import numpy as np
import vip_hci as vip
from vip_hci.preproc import cube_recenter_2dfit, cube_recenter_dft_upsampling, cosmetics, badpixremoval
from vip_hci.negfc import firstguess
from vip_hci.negfc import cube_planet_free
from vip_hci.negfc import show_corner_plot, show_walk_plot
from vip_hci.negfc import mcmc_negfc_sampling
from vip_hci.phot import cube_inject_companions
from vip_hci.pca import pca_adi_annular
plots=vip.var.pp_subplots
import math
import pandas as pd
from math import cos
from math import sin
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import figure
import os
"""
End of package imports 
"""

"""
Defining functions start 
"""

"""
Starphot calculation calculates the flux from the psf, reads the number of
coadds and the integration time from the headers. The function then uses these
values to calculate the parametre starphot the flux from the images in the cube.
"""

def StarphotCalculation(Number_images, psf_xy, starphot):
			
	for a in range(0, Number_images):	
		hdulist=fits.open(file_names[a],ignore_missing_end=True,verbose=False) #Creates list of open files
		Coadd_cube=hdulist[0].header['COADDS']				#Reads the number of coadds from the first fits file
		int_time_cube=hdulist[0].header['ITIME']			#Reads the integration time for each coadd from the first fits file
	
	hdulist2=fits.open(initialpsf,ignore_missing_end=True,verbose=False)
	int_time_psf=hdulist2[0].header['ITIME']			#Reads the integration time for each coadd from the psf fits file
	Coadd_psf=hdulist2[0].header['COADDS']				#Reads the number of coadds from the psf fits file
	
	ycord=np.full((Number_images),psf_xy[1], dtype=int)
	xcord=np.full((Number_images),psf_xy[0], dtype=int) 
	
	flux=vip.phot.contrcurve.aperture_flux(psf, ycord, xcord, fwhm, ap_factor=1,
	mean=False, verbose=False)
	
	psfflux=flux[1]
	starphot=psfflux * ((int_time_cube * Coadd_cube)/(int_time_psf * Coadd_psf))
	return starphot

				
"""
Backgroundflux uses a 100 by 100 square away from the centre 
of the images to get a rough guess of the background noise. 
It takes four 100x100 squares for each image in the data cube 
at oposing points of the image. It then checks that each cube has an
average flux within 0.5 stdev of the population mean. If this criteria
is failed the square is thown out. The average flux is calcuated from 
the remaining squares
"""

def backgroundflux(cube,averageflux):
	number_images=len(cube)										#calculates number of image in cube
	totalflux=np.zeros(shape=[4, number_images])				#Initialises the total flux array variable.
	averageflux_array=np.zeros(shape=[4, number_images])		#Initialises an array to contain average of each 100x100 grid
													
	for k in range (0, (number_images)):
		for a in range (100, 201): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(100, 201):
				totalflux[0,k]=totalflux[0,k] + cube[k,a,j]
				
		for a in range (100, 201): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(800, 901):
				totalflux[1,k]=totalflux[1,k] + cube[k,a,j]
			
		for a in range (800, 901): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(800, 901):
				totalflux[2,k]=totalflux[2,k] + cube[k,a,j]
			
		for a in range (900, 801): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(100, 201):
				totalflux[3,k]=totalflux[3,k] + cube[k,a,j]
		
		for a in range(0, 4):
			averageflux_array[a,k]=totalflux[a,k]/10000 #Calculates the average flux in each grid by dividing by number of pixels (10000)
			
	upper_bound = np.mean(averageflux_array) + (np.std(averageflux_array)/2) #Creates upper and lower bounds of 0.5 stdev above and below population mean
	lower_bound = np.mean(averageflux_array) - (np.std(averageflux_array)/2)
	
	
			
	check_array=np.ones(shape=[4, number_images]) #Initialises array of ones so Boolean logic will accept as 'true'. Same size and shape as averageflux_array
	for k in range (0, (number_images)):
		for a in range(0, 4):
			if averageflux_array[a,k] > upper_bound or averageflux_array[a,k] < lower_bound: 
				check_array[a,k]=0			#If flux out of bounds changes check_array data point to zero so boolean logic would give false
				
	Loop_count=1 #counts how many images are used
	totalflux=0 #Initialises total flux variable that counts all flux
	for k in range (0, (number_images)): 
		for a in range(0, 4):
			if check_array[a,k]:		#If check_array point [a,k] is true, it adds the corresponding averageflux_array point and adds to totalflux
				totalflux=totalflux + averageflux_array[a,k]
				Loop_count=Loop_count + 1
			
	averageflux=totalflux/Loop_count    #Divides the totalflux by the number of images used.
	print "Number of grids used= \n", Loop_count
	print "Average flux= \n", averageflux
	print "Standard Deviation= \n", np.std(averageflux_array)
	Loop_count=0
	
	return averageflux
	
"""
Readangles reads in the filenames and their parralactic 
angles from their headers and writes them to a text file.
"""

def Readangles(Number_images,name_input):
	
	#Loads the filenames into python from the textfile containing them.
	file_names = np.loadtxt('{star_name}_filenames.txt'.format(star_name=name[a]),dtype = 'S') 	#Loads filenames into a text file.
	Number_images=len(file_names)			#Counts the length of the list.
	angles=np.zeros((1, Number_images))		#Initialises an array for the angles.							 
	
	#Loop creates a list of open files and reads the 'PARANG' parametre from the headers.
	for a in range(0, Number_images):
		hdulist=fits.open(file_names[a], ignore_missing_end=True, verbose=False)
		angles[0,a]=hdulist[0].header['PARANG']
	
	
	#Alters the angles from a 2D array into a 1D list for ease of use.
	np.savetxt('{star_name}_angles.txt'.format(star_name=name[a]), angles)		#Saves the angles into a text file.
	angles=np.loadtxt('{star_name}_angles.txt'.format(star_name=name[a])) 		#Loads the angles back into python.


"""	
Buildcube initialises an empty array and loads the images into 
this array, removing the badpixels and also subtracting the flatfield.
"""	

def Buildcube(Number_images,name_input):
	
	file_names= np.loadtxt('{star_name}_filenames.txt'.format(star_name=name[a]),dtype='S') #Load the file names
	Number_images=len(file_names)						#Counts the length of the list.
	
	#Initialises two cubes to loop through.
	cube0=np.zeros((Number_images, 1024, 1024))				
	cube=np.zeros((Number_images, 1024, 1024))				
	
	flatfield0='./flat_Kp.fits'		#Loads in the flat-field.	
		
	flatfield=vip.fits.open_fits(flatfield0, n=0, header=False, 
	ignore_missing_end=True, verbose=True)		#Opens the flat_field.
	
	Images_loaded =0	#Initialises a counter.
	
	#Loop opens the images and loads them into the data cube format required.
	#Loop also divides by the flatfield. 
	for a in range (0, Number_images):
		Images_loaded = Images_loaded + 1	#Counter.
		print "Processing %d/%d." %(Images_loaded, Number_images)
		
		cube0=vip.fits.open_fits(file_names[a], n=0, header=False,
		ignore_missing_end=True, verbose=False)	#Opens the fits file using VIP.
		
		for j in range(0, 1024):
			for k in range(0, 1024):
				cube[a][j][k]= cube0[j][k]/flatfield[j][k]
		if Images_loaded == Number_images:
			print "\nImages loaded in.\n"
	
	#Removes the 'bad pixels' in the image using a VIP function to do so.
	print "Now removing the badpixels in the Cube:\n"
	cube=vip.preproc.badpixremoval.cube_fix_badpix_isolated(cube, bpm_mask=None, sigma_clip=3,
	num_neig=5, size=5, protect_mask=False, radius=30, verbose=True, debug=False)		
		
	#Writes the created cube to a fits file.	
	vip.fits.write_fits('Newcube.fits', cube, dtype32=True, verbose=True)	
	
	#Loads the angles and appends them to the cube, then overwrites the previous cube
	#with the appended one.       
	angles=np.loadtxt('{star_name}_angles.txt'.format(star_name=name[a]))
	cube=vip.fits.append_extension('Newcube.fits', angles)             
	vip.fits.info_fits('Newcube.fits')
	
				
"""
Backgroundflux uses a 100 by 100 square away from the centre 
of the images to get a rough guess of the background noise. 
It takes four 100x100 squares for each image in the data cube 
at oposing points of the image. It then checks that each cube has an
average flux within 0.5 stdev of the population mean. If this criteria
is failed the square is thown out. The average flux is calcuated from 
the remaining squares
"""

def backgroundflux(cube,averageflux):
	number_images=len(cube)										#calculates number of image in cube
	totalflux=np.zeros(shape=[4, number_images])				#Initialises the total flux array variable.
	averageflux_array=np.zeros(shape=[4, number_images])		#Initialises an array to contain average of each 100x100 grid
													
	for k in range (0, (number_images)):
		for a in range (100, 201): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(100, 201):
				totalflux[0,k]=totalflux[0,k] + cube[k,a,j]
				
		for a in range (100, 201): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(800, 901):
				totalflux[1,k]=totalflux[1,k] + cube[k,a,j]
			
		for a in range (800, 901): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(800, 901):
				totalflux[2,k]=totalflux[2,k] + cube[k,a,j]
			
		for a in range (900, 801): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(100, 201):
				totalflux[3,k]=totalflux[3,k] + cube[k,a,j]
		
		for a in range(0, 4):
			averageflux_array[a,k]=totalflux[a,k]/10000 #Calculates the average flux in each grid by dividing by number of pixels (10000)
			
	upper_bound = np.mean(averageflux_array) + (np.std(averageflux_array)/2) #Creates upper and lower bounds of 0.5 stdev above and below population mean
	lower_bound = np.mean(averageflux_array) - (np.std(averageflux_array)/2)
	
	
			
	check_array=np.ones(shape=[4, number_images]) #Initialises array of ones so Boolean logic will accept as 'true'. Same size and shape as averageflux_array
	for k in range (0, (number_images)):
		for a in range(0, 4):
			if averageflux_array[a,k] > upper_bound or averageflux_array[a,k] < lower_bound: 
				check_array[a,k]=0			#If flux out of bounds changes check_array data point to zero so boolean logic would give false
				
	Loop_count=1 #counts how many images are used
	totalflux=0 #Initialises total flux variable that counts all flux
	for k in range (0, (number_images)): 
		for a in range(0, 4):
			if check_array[a,k]:		#If check_array point [a,k] is true, it adds the corresponding averageflux_array point and adds to totalflux
				totalflux=totalflux + averageflux_array[a,k]
				Loop_count=Loop_count + 1
			
	averageflux=totalflux/Loop_count    #Divides the totalflux by the number of images used.
	print "Number of grids used= \n", Loop_count
	print "Average flux= \n", averageflux
	print "Standard Deviation= \n", np.std(averageflux_array)
	Loop_count=0
	
	return averageflux
	
"""
Removeimage removes the bad images by erasing 
the name of the file from the list of filenames 
"""

def Removeimage(Number_images,name_input):

	#Loads in the filenames.
	file_names=np.loadtxt('{star_name}_filenames.txt'.format(star_name=name[a]),dtype='S') 
	user_input=raw_input("Type the file numbers (eg. 1,4,5.. ) of the images you would like to remove\n Seperated by single commas\n")
	Bad_images = [int(a) for a in user_input.split(',') if a.isdigit()]
	Bad_images=sorted(Bad_images, reverse =True)	#Sorts the numbers in descending order.
	
	#Loop removes the filenames corresponding to the numbers entered from the list of filenames.
	for a in range(0,len(Bad_images)):
		file_names=np.delete(file_names, (Bad_images[a]-1))
	
	return file_names

"""		
Checkfunction allows for yes no answers and prompts the 
user to re-enter if initially enetered incorrectly.

Yes = 0
No = 1

"""

def Checkfunction():

	check=2		
	while check == 2:	#Ensures that this will be looped over unless a yes or no is entered.
	
		yes=set(['yes', 'y', 'Y','Yes','YES'])	#Sets any of these entries to equal yes.
		no=set(['no', 'n', 'N','No','NO'])	#Sets any of these entries to equal no.
		exit=set(['EXIT','exit','exit()','Exit()','Exit','EXIT()'])
		question=raw_input('Yes (y) or no (n)?\n')
		if question in yes:
			check=0
			
		if question in no:
			check=1	
		
		if question in exit:
			raise SystemExit()
		
		#An error if loop to prompt the user to change their response.	
		if question not in yes and question not in no:
			print '\nUnrecognised answer\n'
		
	return check

"""
Contrastcurvedata writes the contrast curve results into 
a text file, saves and then loads it back, plots the LLSG
and the PCA data on the same graph and saves it to a file
"""
"""
Note: This function is terribly inefficient, could be much 
better looped over but haven't had the time to do so.
"""
	
def Contrastcurvedata(PCA_contrast_curve, LLSG_contrast_curve,ADI_contrast_curve, pxscale_keck, sigma, pt_sub,name):
	
	#Saves the PCA,ADI and LLSG curve outputs and reads them back in.
	from astropy.io import fits
	import matplotlib.pyplot as plt
	import matplotlib.patches as mpatches
	from pylab import figure
	
	np.savetxt('new_PCA_{star_name}_curve_outputs'.format(star_name=name), PCA_contrast_curve)
	PCA_data=np.loadtxt('new_PCA_{star_name}_curve_outputs'.format(star_name=name))
	np.savetxt('new_LLSG_{star_name}_curve_outputs'.format(star_name=name), LLSG_contrast_curve)
	LLSG_data=np.loadtxt('new_LLSG_{star_name}_curve_outputs'.format(star_name=name))
	np.savetxt('new_ADI_{star_name}_curve_outputs'.format(star_name=name), ADI_contrast_curve)
	ADI_data=np.loadtxt('new_ADI_{star_name}_curve_outputs'.format(star_name=name))
	
	#Initialises the PCA variables (500 being the size).
	PCA_Distance=np.zeros(500)
	PCA_Sensitivity_Gauss=np.zeros(500)
	PCA_Sensitivity_Student=np.zeros(500)
	LLSG_Distance=np.zeros(500)
	LLSG_Sensitivity_Gauss=np.zeros(500)
	LLSG_Sensitivity_Student=np.zeros(500)
	ADI_Distance=np.zeros(500)
	ADI_Sensitivity_Gauss=np.zeros(500)
	ADI_Sensitivity_Student=np.zeros(500)

	#Loop loads the PCA outputs into the PCA variables.
	for a in range (0,500):
		PCA_Distance[a] =PCA_data[a,0] * pxscale_keck
		PCA_Sensitivity_Gauss[a]=PCA_data[a,2]
		PCA_Sensitivity_Student[a]= PCA_data[a,3]	

	#Loop loads the LLSG outputs into the LLSG variables.
		LLSG_Distance[a] =LLSG_data[a,0] * pxscale_keck
		LLSG_Sensitivity_Gauss[a]=LLSG_data[a,2]
		LLSG_Sensitivity_Student[a]= LLSG_data[a,3]
	
	#Loop loads the ADI outputs into the ADI variables.
		ADI_Distance[a] =ADI_data[a,0] * pxscale_keck
		ADI_Sensitivity_Gauss[a]=ADI_data[a,2]
		ADI_Sensitivity_Student[a]= ADI_data[a,3]
	
	#Plotting all 3 on one curve: 
		
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	PCA_Legend = mpatches.Patch(color='#fcb141', label='PCA Contrast Curve')
	LLSG_Legend = mpatches.Patch(color='#2f6fac', label='LLSG Contrast Curve')
	ADI_Legend = mpatches.Patch(color='black', label='ADI Contrast Curve')
	plt.legend(handles=[PCA_Legend,LLSG_Legend,ADI_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.
	
	#Plots a grid (The logarithmic lines present in the background of the plot).
	plt.grid('on', which='both', alpha=0.2, linestyle='solid')
	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	plt.ylim((0.000001,0.01))
	
	#Creates the variables that will be plotted.
	Curve = [0,0,0]	#Initialises the Curve variable.
	Curve[2] = plt.plot(ADI_Distance, ADI_Sensitivity_Gauss, linewidth =2, color='black')
	Curve[1] = plt.plot(PCA_Distance, PCA_Sensitivity_Gauss, linewidth =2, color='#fcb141')
	Curve[0] = plt.plot(LLSG_Distance, LLSG_Sensitivity_Gauss, linewidth =2, color='#2f6fac')
	
	savefig('new_Contrast_curves_{star_name}.png'.format(star_name=name))
	
	#LLSG contrast curve 
	
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	LLSG_Legend = mpatches.Patch(color='#2f6fac', label='LLSG Contrast Curve')
	plt.legend(handles=[LLSG_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.
	
	#Plots a grid (The logarithmic lines present in the background of the plot).
	plt.grid('on', which='both', alpha=0.2, linestyle='solid')
	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	plt.ylim((0.000001,0.01))
	
	#Creates the variables that will be plotted.
	Curve = 0	#Initialises the Curve variable.
	Curve = plt.plot(LLSG_Distance, LLSG_Sensitivity_Gauss, linewidth =2, color='#2f6fac')
	
	#Saves the figure and then shows the plot. 
	savefig('new_LLSG_curves_{star_name}.png'.format(star_name=name))
	
	#The ADI curve:
	
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	ADI_Legend = mpatches.Patch(color='black', label='ADI Contrast Curve')
	plt.legend(handles=[ADI_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.
	
	#Plots a grid (The logarithmic lines present in the background of the plot).
	plt.grid('on', which='both', alpha=0.2, linestyle='solid')
	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	plt.ylim((0.000001,0.01))
	
	#Creates the variables that will be plotted.
	Curve = 0	#Initialises the Curve variable.
	Curve = plt.plot(ADI_Distance, ADI_Sensitivity_Gauss, linewidth =2, color='black')
	
	#Saves the figure and then shows the plot. 
	savefig('new_ADI_curve_{star_name}.png'.format(star_name=name))
	
	#PCA contrast curve. 
	
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	PCA_Legend = mpatches.Patch(color='#fcb141', label='PCA Contrast Curve')
	plt.legend(handles=[PCA_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.
	
	#Plots a grid (The logarithmic lines present in the background of the plot).
	plt.grid('on', which='both', alpha=0.2, linestyle='solid')
	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	plt.ylim((0.000001,0.01))
	
	#Creates the variables that will be plotted.
	Curve = 0	#Initialises the Curve variable.
	Curve = plt.plot(PCA_Distance, PCA_Sensitivity_Gauss, linewidth =2, color='#fcb141')
	
	#Saves the figure and then shows the plot. 
	savefig('new_PCA_curve_{star_name}.png'.format(star_name=name))
""" 
End of defining functions 
"""

name = [line.rstrip('\n') for line in open('name.txt')]
print name

date = [line.rstrip('\n') for line in open('date.txt')]
print date

name_file = [line.rstrip('\n') for line in open('name_file.txt')]

info = np.loadtxt('star_info.txt')
print info

coor_x = np.zeros(len(info))
coor_y = np.zeros(len(info))
ref_x = np.zeros(len(info))
ref_y = np.zeros(len(info))
llsg_rank = np.zeros(len(info))

for a in range (0,len(name)):
	os.chdir('{star_name}/{epoch}'.format(star_name=name[a],epoch=date[a]))
	# run the script in here
	file_names=np.loadtxt('{star}_filenames.txt'.format(star=name_file[a]), dtype='S')	
	Number_images=len(file_names)
	coor_x[a] = info[a,0]
	coor_y[a]= info[a,1]
	ref_x[a]= info[a,2]
	ref_y[a]= info[a,3]
	llsg_rank[a]= info[a,4]
	#Loads the psf and the cube fits into python.
	initialpsf='./psf.fits'			
	cube='./Newcube.fits'
				
	#Opens the cube (using VIP) with cube_orig as HDU:0 and calls it cube_orig, 
	#and the parallactic angles as HDU:1, HDU:1 being the previously appended part. 
	#Uses VIP to also open the point spread function previously loaded in.
	cube_orig, angs=vip.fits.open_adicube(cube)
	psf=vip.fits.open_fits(initialpsf, n=0, header=False, ignore_missing_end=True, verbose=True)
	
	print "Fitting a 2D gaussian to centre the images\n"

	#Uses VIP's 2D gaussian fitting algorithm to centre the cube.
	gauss=vip.var.fit_2dgaussian(psf, crop=True, cropsize=30, cent=(psf_xy[0], psf_xy[1]),
	full_output=True, debug=False)		
	#Calculates the fwhm from the VIP gaussian function.
	print gauss[0:1]
	fwhm_x=gauss.iloc[0,3]
	fwhm_y=gauss.iloc[0,4]
	fwhm=np.mean([fwhm_x, fwhm_y])
	print "fwhm=",fwhm

	fwhm_int = math.ceil(fwhm)
	fwhm_int = int(fwhm_int)
	print "round up =", fwhm_int
	
	#Shy1 = Shift in y, shx1 = shift in x. Calculates the shifts needed to centre the image.
	#Negative = True (Bright spot at the centre) Negative = False (No bright spot at centre).
	cube1, shy1, shx1=cube_recenter_2dfit(cube_orig, (int(coor_x[a]), int(coor_y[a])), fwhm,
	model='gauss',nproc=1, subi_size=fwhm_int, negative=False,full_output=True, debug=False)

	#Plots the shifts in x and y used to centre the images.
	#figure(figsize=(8,3))								
	#plot(shy1, 'o-', label='shifts in y', alpha=0.5)	
	#plot(shx1, 'o-', label='shifts in x', alpha=0.5)
	#legend(loc='best')
	#grid('on')
	#suptitle('Positive 2d Gaussian fit')
	#plots()
	
	#Writes the values of the centered cube into a fits file.
	vip.fits.write_fits('centeredcube_{star_name}.fits'.format(star_name=name[a]), cube1, dtype32=True, verbose=True)	

	cube=cube1	#Loads up the centered cube.
	#Plots the original cube vs the new centered cube.

	im1=vip.preproc.cosmetics.frame_crop(cube_orig[0], 1000, verbose=False)
	im2=vip.preproc.cosmetics.frame_crop(cube[0], 1000, verbose=False)

	#plots(im1, im2, grid=True, title='Original first frame / First frame after re-centering')

	print "\n\n\n======================REFERENCE COORDINATE============================="

	print 'Will now reduce the data using ADI median subtraction\n'

	#Reduces the image using ADI median subtraction then plots the image.
	fr_adi = vip.madi.adi(cube, angs, mode='fullfr')
	#plots(fr_adi, title='ADI median reduction of {star_name}'.format(star_name=name[a]), dpi=100, vmin=-10, vmax=100, colorb=True)
	vip.fits.write_fits('new_ADI_{star_name}.fits'.format(star_name=name[a]), fr_adi, dtype32=True, verbose=True) 

	print 'Will now reduce the data using full frame PCA\n'

	#Calculates the optimum number of principle components to use in PCA reduction.

	opt_pcs=vip.pca.pca_optimize_snr(cube, angs, fwhm=fwhm, source_xy=(int(ref_x[a]), int(ref_y[a])),
	mask_center_px=None, fmerit='px', range_pcs=(1,20))

	#Reduces the image using full frame PCA then plots the image.
	fr_pca1=vip.pca.pca(cube, angs, ncomp=opt_pcs, mask_center_px=None)			
	#plots(fr_pca1, title='FFPCA reduction of {star_name}'.format(star_name=name[a]), dpi=100, vmin=-10, vmax=100, colorb=True)
	#Loop asks user if they would like to save the image.
	vip.fits.write_fits('new_FFPCA_{star_name}.fits'.format(star_name=name[a]), fr_pca1, dtype32=True, verbose=True) 
	
	#Reduces the image using the LLSG algorithm then plots the image.
	fr_llsg=vip.llsg.llsg(cube,angs,fwhm,rank=int(llsg_rank[a]),thresh=2.5,max_iter=20,random_seed=10)
	#plots(fr_llsg,title ='LLSG reduced image of {star_name}'.format(star_name=name[a]), dpi= 100, vmin= -10,vmax=100, colorb=True)

	vip.fits.write_fits('new_LLSG_{star_name}.fits'.format(star_name=name[a]), fr_llsg, dtype32=True, verbose=True) 

	#Starphot is the  total aperture flux from the star during
	#the integration time (if no coronograph).
	starphot=StarphotCalculation(Number_images, psf_xy, starphot)
	starphot_array=np.zeros(500)
	for i in range (0,500):
		starphot_array[i] = starphot

	np.savetxt('{name}_starphot'.format(name=name_file[a]), starphot_array)

	number_planets=5 #Initialises the variable as an odd number to use the integer check loop.

	#Loop for if artificial planets are used, loop calculates where to place them,
	#injects them and then removes them again.
	if number_planets > 0:
		rad_dists=np.empty(number_planets)

	#Calculates the orbital radius to place the injected planets at.	
	for i in range(number_planets):
		rad_dists[i]=seperation * (i+1)
	
	n_branches=3

	#Calculates the locations of the branches with their respective injected planets.
	if n_branches > 0:
		theta=360 / n_branches	#The angle to separate each branch by.

	#Checks if the user would like to input their own 
	#brightness value for the injected planets.
	brightness_check = 2
	flux_multiple=3
						
	averageflux=backgroundflux(cube, averageflux)
	flvl=float(averageflux)*float(flux_multiple)

	#Creates a normalised PSF and a data cube with the injected planets.
	psf=vip.phot.fakecomp.psf_norm(psf, fwhm, size=100, threshold=None, mask_core=None)
	cubefc=cube_inject_companions(cube, psf, angs, flevel=flvl, plsc=pxscale_keck,rad_dists=rad_dists,theta=theta, n_branches=n_branches,verbose=False) 

	#Initialises a bunch of variables to do with calculating the location
	#of the injected planets.

	average=np.zeros(2)
	sum=np.zeros(2)
	newcent=np.zeros(2)
	synplanetlocation=np.zeros(2)

	#From the shift of frames in pixels, find the center of the image
	Nx=np.prod(shx1.shape)

	for i in range(0,Nx):
		sum[0]=sum[0] +shx1[i]
	
	average[0]=sum[0]/Nx	

	Ny=np.prod(shy1.shape)
	for i in range(0,Ny):
		sum[1]=sum[1] +shy1[i]
	
	average[1]=sum[1]/Ny	

	newcent[0]=coor_x[a] + average[0] 
	newcent[1]=coor_y[a] + average[1]


	synplanetlocation[0]=newcent[0] + (seperation * sin(-theta))
	synplanetlocation[1]=newcent[1] + (seperation * cos(-theta))
	synplanetlocation[0]=round(synplanetlocation[0])			#There will be a noticeable error if the shift between images is large
	synplanetlocation[1]=round(synplanetlocation[1])			#These are only rough guesses, tends to be out by a couple of pixels

	print "synplanetlocation = ", synplanetlocation

	#optimises the number of principle components with the new cubefc.
	opt_pcs=vip.pca.pca_optimize_snr(cubefc, angs, fwhm=fwhm,source_xy=(synplanetlocation[0], synplanetlocation[1]),
	mask_center_px=None, fmerit='mean', range_pcs=(1,20))

	#Plots the data cube with the injected planets
	fr_pca3=vip.pca.pca(cubefc, angs, ncomp=opt_pcs)

	vip.fits.write_fits('new_Planet_injectionPCA.fits', fr_pca3, dtype32=True, verbose=True) 

	#_=vip.phot.frame_quick_report(fr_pca3, fwhm=fwhm, source_xy=(synplanetlocation[0],synplanetlocation[1])) #give nan in here

	source_xy=[(ref_x[a], ref_y[a])]

	print "source_xy = ", source_xy

	r_0, theta_0, f_0=firstguess(cubefc, angs, psf, ncomp=10, plsc=pxscale_keck,planets_xy_coord=source_xy,
	fwhm=fwhm,annulus_width=3, aperture_radius=3,f_range=np.linspace(100,2000,20),  simplex=True, display=False, verbose=True)

	plpar_bpicb=[(r_0, theta_0, f_0)]
	cube_emp=cube_planet_free(plpar_bpicb, cubefc, angs, psf, pxscale_keck)
	fr_pca_emp=vip.pca.pca(cube_emp, angs, ncomp=opt_pcs, verbose=False)
	pt_sub = 0.1

	PCA_contrast_curve=vip.phot.contrast_curve(cube_emp, angs, psf, fwhm, pxscale_keck, starphot,
	sigma=sigma, nbranch=n_branches, algo=vip.pca.pca, ncomp=opt_pcs, debug=True,save_plot='PCA', plot=False)

	LLSG_contrast_curve=vip.phot.contrast_curve(cube_emp, angs, psf, fwhm, pxscale_keck, starphot,
	sigma=sigma, nbranch=n_branches, algo=vip.llsg.llsg, debug=True,save_plot='LLSG', plot=False)

	ADI_contrast_curve=vip.phot.contrast_curve(cube_emp, angs, psf, fwhm, pxscale_keck, starphot,
	sigma=sigma, nbranch=n_branches, algo=vip.madi.adi, debug=True,save_plot='ADI', plot=False) 

	Contrastcurvedata(PCA_contrast_curve, LLSG_contrast_curve, ADI_contrast_curve, pxscale_keck, sigma,pt_sub, name=name[a])

	# change dir 
	os.chdir('../..')
	a = a + 1



					   
