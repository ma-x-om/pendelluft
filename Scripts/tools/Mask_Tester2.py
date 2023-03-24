### This script aims to automatically segment the video into the Lung Mask and also return the horizontal center of the lung
###		it will receive a .mat file containing the video frames. This skips the process of reprocessing the video into frames.
### obs: always with an axial cut in mind

from scipy.io import loadmat
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
import cv2

def createContour(Image):
	nx, ny =  np.shape(Image)
	Image_Contour = np.zeros((nx,ny))
	Image_Contour_ixs, _ = cv2.findContours(Image, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
	Image_Contour_squeeze = np.squeeze(Image_Contour_ixs)

	for pixel in contour:
		Image_Contour[pixel[1],pixel[0]] = 1

	return Image_Contour

root = tk.Tk()
root.withdraw()

file_path = filedialog.askopenfilename(title = "Abra o arquivo gerado em Code_1.py (Images_U_V.mat)", filetypes=(("Arquivos MAT","*.mat"),("Todos os arquivos","*.*")))	# open explorer
main_folder = os.path.split(file_path)[0]+'/'
Mat_File = loadmat(file_path)

mask_path = filedialog.askopenfilename(title = "Abra o arquivo de LungMask", filetypes=(("Arquivos MAT","*.mat"),("Todos os arquivos","*.*")))	# open explorer
mask_folder = os.path.split(mask_path)[0]+'/'
mask_file = loadmat(mask_path)

loaded_mask = mask_file['BW']
video = Mat_File['Imagens']

sum_of_all_frames = video.sum(2)
normalized_sum = sum_of_all_frames/sum_of_all_frames.max()
sum_uint8 = (normalized_sum*255).astype(np.uint8)

# The following lines will compare the calculated mask with the loaded mask using a variable threshold starting at 1 (max) and comparing how much "True" pixels there are
# We must also set a tolerance in which we can assume the calculated_mask is the same as the loaded_mask

# first by number of pixels:

#tolerance = 0.02 # in %
#threshold = 1 # in %
#step = 0.001

tolerance = 0.02 # in %
threshold_p = 1 # in %
step = 0.001

n_of_true_pixels = np.count_nonzero(loaded_mask)
similarity_by_n_of_pixels = 0
while abs(1-similarity_by_n_of_pixels) > tolerance:
	threshold_p -= step
	calculated_mask_p = normalized_sum >= threshold_p
	calc_pixels = np.count_nonzero(calculated_mask_p)
	similarity_by_n_of_pixels = calc_pixels/n_of_true_pixels

# now by pixel location:

threshold_l = 1 # in %

mask_indexes = np.where(loaded_mask == 1)
similarity_by_location = 0
while abs(1-similarity_by_location) > tolerance:
	threshold_l -= step
	calculated_mask_l = normalized_sum >= threshold_l
	pixels_correlation_array = calculated_mask_l[mask_indexes]
	pixels_correlation = pixels_correlation_array.sum()
	similarity_by_location = pixels_correlation/n_of_true_pixels

plt.figure()
plt.subplot(2,2,(1,2))
plt.imshow(loaded_mask, cmap='gray')
plt.axis('off')
plt.title('Loaded Mask')
plt.subplot(223)
plt.imshow(calculated_mask_p, cmap='gray')
plt.axis('off')
plt.title(f'Similarity of True pixels number\nSimilarity of {similarity_by_n_of_pixels*100:.2f}%\nThreshold of {threshold_p*100:.2f}%')
plt.subplot(224)
plt.imshow(calculated_mask_l, cmap='gray')
plt.axis('off')
plt.title(f'Similarity of True pixels location\nSimilarity of {similarity_by_location*100:.2f}%\nThreshold of {threshold_l*100:.2f}%')
plt.suptitle(f'CALCULATED MASKS\nTolerance of {tolerance*100}%\n step of {step*100}%')
plt.tight_layout()
plt.show()

