### little code to generate a video file in greyscale from images contained in a .mat file

import cv2
import os
from scipy.io import loadmat
import numpy as np
import tkinter as tk
from tkinter import filedialog
 

root = tk.Tk()
root.withdraw()

file_path = filedialog.askopenfilename(title = "Abra o video dentro do .MAT (Images_U_V.mat)", filetype=(("Arquivos .MAT","*.mat"),("Todos os arquivos","*.*")))    # abre explorer
main_folder = os.path.split(file_path)[0]+'/'
dados = loadmat(file_path)
Imagem = dados['Imagens']
r, c, N_frames = np.shape(Imagem)
print(f"{r}x{c}, {N_frames} frames")

frameSize = (r, c)
FPS = 50

# cv2 can't work with float64 data type
#   and for some reason I can't get it to work with float32
#   so I convert it to uint8, which works perfectly for our case

np_Imagens = np.array(Imagem)
np_Imagens = np_Imagens.astype(np.float64) / np_Imagens.max()
np_Imagens = 255 * np_Imagens
np_Imagens = np_Imagens.astype(np.uint8)

out = cv2.VideoWriter((main_folder+'video_render.avi'),cv2.VideoWriter_fourcc(*'DIVX'), FPS, frameSize)
print(f"Video saved at {main_folder}")

for i in range(N_frames):
    img = np_Imagens[:,:,i]
    bgr = cv2.cvtColor(img, cv2.COLOR_GRAY2RGB) 
    out.write(bgr)

out.release()
