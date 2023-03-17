### This script aims to automatically segment the video into the Lung Mask and also return the horizontal center of the lung
###		it will receive a .mat file containing the video frames. This skips the process of reprocessing the video into frames.
### obs: always with an axial cut in mind

from scipy.io import loadmat
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
import skimage.morphology as skim

root = tk.Tk()
root.withdraw()

file_path = filedialog.askopenfilename(title = "Abra o arquivo gerado em Code_1.py (Images_U_V.mat)", filetypes=(("Arquivos MAT","*.mat"),("Todos os arquivos","*.*")))	# open explorer
main_folder = os.path.split(file_path)[0]+'/' #	pasta que contém o arquivo escolhido
Mat_File = loadmat(file_path)	# open the file in python

video = Mat_File['Imagens']

# The idea for the automatic lung mask is we will sum all the frames and then apply a threshold, which will give us a binary matrix.

sum_of_all_frames = video.sum(2)
normalized_sum = sum_of_all_frames/sum_of_all_frames.max()

threshold = 0.135 # in %
lung_mask = normalized_sum >= threshold

# For the lung center we will sum up all the columns and then find the local minimum and hope for the best

lung_center_sum = lung_mask.sum(0)
center_candidates = skim.local_minima(lung_center_sum)
center_indexes = np.where(center_candidates)[0]
image_center = np.size(lung_center_sum)//2

minimum = 200 # random great number
for i in center_indexes:
	if abs(i-image_center)<minimum:
		minimum = abs(i-image_center) # distance between index and image center
		center_index = i


plt.imshow(lung_mask, cmap='gray')
plt.axvline(center_index, color='red')
plt.axis('off')
plt.title(f'Threshold of {threshold*100}%\nCenter at pixel {center_index+1} of {np.size(lung_center_sum)}')
plt.show()