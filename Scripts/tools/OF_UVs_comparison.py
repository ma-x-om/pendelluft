# this script is intended to compare HS and Farneback methods by plotting U and V vectors from two .mat files on the same frame image
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import os
import tkinter as tk
from tkinter import filedialog
from matplotlib.patches import Rectangle

#########################################
# Initialize tkinter (Tk)
# used for browsing files in explorer
#
root = tk.Tk()
root.withdraw()
#
#########################################

# here we load .mat files with Images, u and v for the normal pig
#   make sure that those files are in the same folder as this script
hs_file_path = filedialog.askopenfilename(title = "Abra o arquivo Images_U_V.mat para H&S", filetypes=(("Arquivos MAT","*.mat"),("Todos os arquivos","*.*")))	# open explorer
fb_file_path = filedialog.askopenfilename(title = "Abra o arquivo Images_U_V.mat para Farneback", filetypes=(("Arquivos MAT","*.mat"),("Todos os arquivos","*.*")))	# open explorer
HS = loadmat(hs_file_path)
FB = loadmat(fb_file_path)

frame_ML = int(input("Fala o frame que você quer, comandante (o primeiro frame é o 1)\n")) #ML stands for MatLab
frame_py = frame_ML-1 #py stands for python

hs_file_name = os.path.split(hs_file_path)[1] # file path without '.txt'
main_folder = os.path.splitext(hs_file_path)[0]+'/'
isExist = os.path.exists(main_folder)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(main_folder)
image_folder = main_folder+f'Comparison_renders_f{frame_ML}/'
isExist = os.path.exists(image_folder)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(image_folder)


img_hs = HS['Imagens'][:,:,frame_py]
u_hs = HS['u'][:,:,frame_py]
v_hs = HS['v'][:,:,frame_py]

u_fb = FB['u'][:,:,frame_py]
v_fb = FB['v'][:,:,frame_py]

r, c = np.shape(img_hs)

x, y = np.meshgrid(np.arange(r),np.arange(c))

#plot both u and v for HS and Farneback on the same frame

def draw_vectors(vector, Vec_type, OF_type):
	plt.figure(dpi=300)
	plt.imshow(vector)
	plt.axis('off')
	plt.colorbar()
	#plt.quiver(X, Y, u,-v)		# pyplot's function for plotting quivers
												# receives u and v for the respective frame
												#
	#plt.title(f'{OF_type} - {Vec_type}')
	plt.savefig(f'{image_folder}{OF_type}_flow_{Vec_type}.png',bbox_inches='tight')
	plt.close()

def draw_quivers(image, u, v,OF_type):
	r, c = np.shape(u)

	reta_x = c*0.0625
	reta_y = 0
	reta_L = c*0.3
	reta_H = r*0.06375

	arr_au = abs(u.max())
	arr_iu = abs(u.min())
	arr_av = abs(v.max())
	arr_iv = abs(v.min())

	arr_u = max(arr_au,arr_iu)
	arr_v = max(arr_av,arr_iv)

	arr_m = (arr_u+arr_v)/2

	fig1, ax1 = plt.subplots()
	ax1.imshow(image, cmap='gray')
	ax1.set_axis_off()
	Q = ax1.quiver(x,y,u,-v,color='red')
	ax1.add_patch(Rectangle((0, 0), reta_L, reta_H, facecolor='white',edgecolor='black',fill=True,lw=2))
	qk = ax1.quiverkey(Q, 0.1, 0.95, arr_m, f'{arr_m:.1E} px', labelpos='E', coordinates='axes') 
	plt.tight_layout()
	plt.savefig(f'{image_folder}{OF_type}_flow_quivers.png',bbox_inches='tight')
	plt.show()

def draw_quivers_old(image, u, v,OF_type):
	plt.figure()
	#plt.title(f'Frame {frame_ML}')
	plt.imshow(image, cmap='gray')
	#plt.colorbar(shrink=0.6)
	plt.axis('off')
	plt.quiver(x,y,u,-v,color='red')
	plt.tight_layout()
	plt.savefig(f'{image_folder}{OF_type}_flow_quivers.png',bbox_inches='tight')
	plt.close()

draw_quivers(img_hs, u_fb, v_fb, "Farneback")
draw_quivers(img_hs, u_hs, v_hs, "HS")

draw_vectors(u_fb, 'U', 'Farneback')
draw_vectors(v_fb, 'V', 'Farneback')
draw_vectors(u_hs, 'U', 'HornSchunck')
draw_vectors(v_hs, 'V', 'HornSchunck')