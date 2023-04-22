import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import ndimage
from scipy.io import loadmat
import tkinter as tk
from tkinter import filedialog

#########################################
# Initialize tkinter (Tk)
# used for browsing files in explorer
#
root = tk.Tk()
root.withdraw()
#
#########################################

file_path = filedialog.askopenfilename()	# abre explorer
file_name = os.path.split(file_path)[1] # file path without '.txt'
main_folder = os.path.split(file_path)[0]+'/'
image_folder = main_folder+'IMAGES/raw_render/'
isExist = os.path.exists(image_folder)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(image_folder)

def get_img(voxels):
	r, c, frames = np.shape(voxels)
	B = 0

	# this loop makes B equal to the higher pixel value in all frames from 'Imagens'
	for it in np.arange(frames):
		A = voxels[:,:,it].max()
		if A>B:
			B = A

	for fr in np.arange(frames):
		plt.imsave(f"{image_folder}{file_name}_frame{1+fr}.png", voxels[:,:,fr], cmap='gray',vmin=0,vmax=B)

image_file = loadmat(file_path)
Imagem = image_file['Imagens']
get_img(Imagem)