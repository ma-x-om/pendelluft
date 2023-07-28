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

def get_arrow_value(u,v):

	arr_au = abs(u.max())
	arr_iu = abs(u.min())
	arr_av = abs(v.max())
	arr_iv = abs(v.min())

	arr_u = max(arr_au,arr_iu)
	arr_v = max(arr_av,arr_iv)

	arr_m = (arr_u+arr_v)/2

	return arr_m

# here we load .mat files with Images, u and v for the normal pig
#   make sure that those files are in the same folder as this script
hs_file_path = filedialog.askopenfilename(title = "Abra o arquivo Images_U_V.mat para H&S", filetypes=(("Arquivos MAT","*.mat"),("Todos os arquivos","*.*")))	# open explorer
fb_file_path = filedialog.askopenfilename(title = "Abra o arquivo Images_U_V.mat para Farneback", filetypes=(("Arquivos MAT","*.mat"),("Todos os arquivos","*.*")))	# open explorer
HS = loadmat(hs_file_path)
FB = loadmat(fb_file_path)

frame_ML = int(input("Fala o frame que você quer, comandante (o primeiro frame é o 1)\n")) #ML stands for MatLab
frame_py = frame_ML-1 #py stands for python

img_hs = HS['Imagens'][:,:,frame_py]
u_hs = HS['u'][:,:,frame_py]
v_hs = HS['v'][:,:,frame_py]

u_fb = FB['u'][:,:,frame_py]
v_fb = FB['v'][:,:,frame_py]

r, c = np.shape(img_hs)

x, y = np.meshgrid(np.arange(r),np.arange(c))

reta_x = c*0.0625
reta_y = 0
reta_L = c*0.3
reta_H = r*0.05375*2

#plot both u and v for HS and Farneback on the same frame

arr_m_hs = get_arrow_value(u_hs,v_hs)
arr_m_fb = get_arrow_value(u_fb,v_fb)

fig1, ax1 = plt.subplots(dpi=300)
ax1.set_title(f'Frame {frame_ML}')
ax1.imshow(img_hs, cmap='summer')
ax1.set_axis_off()
qv_hs = ax1.quiver(x,y,u_hs,-v_hs,color='blue')
qv_fb = ax1.quiver(x,y,u_fb,-v_fb,color='red')
ax1.add_patch(Rectangle((0, 0), reta_L, reta_H, facecolor='white',edgecolor='black',fill=True,lw=0.5))
qk_hs = ax1.quiverkey(qv_hs, 0.1, 0.95, arr_m_hs, f'{arr_m_hs:.1E} px', labelpos='E', coordinates='axes', fontproperties=dict(size=6))
qk_fb = ax1.quiverkey(qv_fb, 0.1, 0.9, arr_m_fb, f'{arr_m_fb:.1E} px', labelpos='E', coordinates='axes', fontproperties=dict(size=6))
ax1.plot([],[],color='blue', label='Horn & Schunck OF')
ax1.plot([],[],color='red', label='Farneback OF')
ax1.legend(loc="upper center", bbox_to_anchor=(0.5, 0.15), ncol=1, fontsize = 'x-small', fancybox= True)
plt.tight_layout()
plt.show()