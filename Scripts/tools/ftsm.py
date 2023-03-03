from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt

HS = loadmat('Normal_HS.mat')
FB = loadmat('Normal_FB.mat')

frame_ML = 36 #ML stands for MatLab
frame_py = frame_ML-1 #py stands for phuck you

img_hs = HS['Imagens'][:,:,frame_py]
u_hs = HS['u'][:,:,frame_py]
v_hs = HS['v'][:,:,frame_py]

u_fb = FB['u'][:,:,frame_py]
v_fb = FB['v'][:,:,frame_py]

r, c = np.shape(img_hs)

x, y = np.meshgrid(np.arange(r),np.arange(c))

#plot the shit

plt.figure(dpi=300)
plt.title(f'Frame {frame_ML}')
plt.imshow(img_hs, cmap='summer')
plt.colorbar(shrink=0.6)
plt.axis('off')
plt.quiver(x,y,u_hs,-v_hs,color='blue')
plt.quiver(x,y,u_fb,-v_fb,color='red')
plt.plot([],[],color='blue', label='Horn & Schunck OF')
plt.plot([],[],color='red', label='Farneback OF')
plt.legend(loc="upper center", bbox_to_anchor=(0.5, 0.15), ncol=1, fontsize = 'x-small', fancybox= True)
plt.show()