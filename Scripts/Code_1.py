## Ties together the following .m scripts:
#	1. Gera_Imagens_E_Campos_Vetoriais.m
#	x. HSoptical.m (not used, instead Farneback is used through a cv2 function)
#	x. quiversDHHD2D.m (not needed, instead there is a pyplot function for that)

## Receives:	txt video containing EIT frames
## Outputs:	images with drawn Motion Field Vectors for each frame
##			.mat file with the images, u and v MFV	

# load libraries
#
import os
from scipy.io import loadmat
from scipy.io import savemat
import matplotlib.pyplot as plt
import numpy as np
import skimage.transform as skt
import cv2
import tkinter as tk
from tkinter import filedialog
import SVD_applier as svd

# load external scripts
#	(as of 21/02/2023 they aren't used)
#import quivers as qv
#import quiver2 as qv2

#########################################
# Initialize tkinter (Tk)
# used for browsing files in explorer
#
root = tk.Tk()
root.withdraw()
#
#########################################

def matlab_style_gauss2D(shape=(3,3),sigma=0.5):
    """
    2D gaussian mask - should give the same result as MATLAB's
    fspecial('gaussian',[shape],[sigma])
    """
    m,n = [(ss-1.)/2. for ss in shape]
    y,x = np.ogrid[-m:m+1,-n:n+1]
    h = np.exp( -(x*x + y*y) / (2.*sigma*sigma) )
    h[ h < np.finfo(h.dtype).eps*h.max() ] = 0
    sumh = h.sum()
    if sumh != 0:
        h /= sumh
    return h

# carregar o video (arquivo com frames)
file_path = filedialog.askopenfilename()	# abre explorer
file_name = os.path.split(file_path)[1] # file path without '.txt'
main_folder = os.path.splitext(file_path)[0]+'/'
isExist = os.path.exists(main_folder)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(main_folder)
image_folder = main_folder+'IMAGES/'
isExist = os.path.exists(image_folder)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(image_folder)
dados = np.loadtxt(file_path)				# atribui o texto do arquivo à variável "dados"
(N_frames, N_pixels) = np.shape(dados)		# atribui as dimensões nº de frames e de pixels por frame

minimo = np.amin(dados[9:-10,:],axis=0)		# vetor com o frame "minimo"

r = c = int(np.sqrt(N_pixels)) # r = c = 32		# r e c tem as dimensoes do frame

Imagens = np.empty([r,c,N_frames])			# great array to hold the future images
tempoIm = []								# initialize empty lists
tempocv = []

# o loop abaixo itera, frame by frame, o video original aplicando a filtro I e salvando
# 	o resultado no array "Imagens"
for fr in np.arange(N_frames):
	imagem2d = np.reshape(dados[fr,:]-minimo,(r,c),order='F')
	imagem2d = skt.rotate(imagem2d, -90)
	Imagens[:,:,fr] = imagem2d

Imagens = svd.apply_video_svd(Imagens)

# here we initialize the empty arrays u and v, which will hold the motion field vectors
u = np.empty([r,c,N_frames-1])	# x axis motion field vector
v = np.empty([r,c,N_frames-1])	# y axis motion field vector

# here we calculate the optical flow
# let's use Farneback instead of H&S
for fr in np.arange(N_frames-1):
	# Farneback Optical Flow Calc with default parameters
	flow = cv2.calcOpticalFlowFarneback(Imagens[:,:,fr], Imagens[:,:,fr+1], None, 0.5, 3, 15, 3, 5, 1.2, 0) #1.2 is used as 0 in john's code
	# we then split "flow" array components into u's and v's corresponding "frame"
	u[:,:,fr], v[:,:,fr] = flow[:,:,0], flow[:,:,1]

savemat(main_folder+'Images_U_V.mat',{'Imagens': Imagens, 'u': u, 'v': v})	# saves what we calculated so far
															#	into a tidy ".mat" file
															#	where the video was loaded
# we define B as 0 so it can be used as our maximum comparison variable
B = 0

# this loop makes B equal to the higher pixel value in all frames from 'Imagens'
for it in np.arange(N_frames):
	A = Imagens[:,:,it].max()
	if A>B:
		B = A

X, Y = np.meshgrid(np.arange(r),np.arange(c))	# creates a meshgrid with dimensions (r, c)
												# 	and atributes it to variables X and Y

# loop a seguir itera pelos frames de 'Imagens'
print(f"Rendering and saving pictures to {image_folder}. If you wish to break/stop it, do it the hard way (CTRL C).")
for it in np.arange(N_frames-1): #np.arange(N_frames-1)
	plt.figure(dpi=500)
	plt.imshow(Imagens[:,:,it],vmin=0,vmax=B)	# pega o frame e equaliza a partir de B=vmax
	plt.axis('off')
	plt.colorbar()
	#qv2.plot_quivers(X, Y, u[:,:,it],v[:,:,it], 1, 3, '')
	plt.quiver(X, Y, u[:,:,it],-v[:,:,it])		# pyplot's function for plotting quivers
												# receives u and v for the respective frame
												#
	plt.savefig(f'{image_folder}Image{it:04.0f}.png',bbox_inches='tight')
	#plt.show()
	plt.close()
