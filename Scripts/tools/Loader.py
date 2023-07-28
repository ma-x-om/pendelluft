import tkinter as tk
from tkinter import filedialog
import os

root = tk.Tk()
root.withdraw()

file_path_freqmap = filedialog.askopenfilename(title = "Abra o arquivo do mapa de frequências (calculos_freqmap.mat)", filetype=(("Arquivos MAT","*.mat"),("Todos os arquivos","*.*")))	# open explorer
freqmap_folder = os.path.split(file_path_freqmap)[0]+'/' #	pasta que contém o arquivo escolhido
from scipy.io import loadmat
freqmap_file = loadmat(file_path_freqmap)	# open the file in python
phi_file = loadmat(freqmap_folder+'calculos_phi.mat')

arr_LE = freqmap_file['arr_LE'][0]
arr_LW = freqmap_file['arr_LW'][0]
arr_E = freqmap_file['arr_E']
arr_W = freqmap_file['arr_W']
arr_LEm = freqmap_file['arr_LEm'][0]
arr_LWm = freqmap_file['arr_LWm'][0]

from numpy import shape
Quad, r, c = shape(arr_E) # Quad = quadrants

# Here lies image showing time again, kids!
import matplotlib.pyplot as plt

plt.figure()
plt.subplot(221)
PLot = arr_E[0]+arr_E[1]
plt.imshow(PLot)
plt.plot(arr_LE[1][:,0],arr_LE[1][:,1],'rx')
plt.plot(arr_LE[0][:,0],arr_LE[0][:,1],'rx')
plt.axis('off')
plt.colorbar()
plt.plot([c/2-1, c/2-1],[0, r-1],'k')
plt.title('E max')

plt.subplot(222)
PLot = arr_W[0]+arr_W[1]
plt.imshow(PLot)
plt.plot(arr_LW[1][:,0],arr_LW[1][:,1],'rx')
plt.plot(arr_LW[0][:,0],arr_LW[0][:,1],'rx')
plt.axis('off')
plt.colorbar()
plt.plot([c/2-1, c/2-1],[0, r-1],'k')
plt.title('W max')

plt.subplot(223)
PLot = arr_E[2]+arr_E[3]
plt.imshow(PLot)
plt.plot(arr_LE[2][:,0],arr_LE[2][:,1],'rx')
plt.plot(arr_LE[3][:,0],arr_LE[3][:,1],'rx')
plt.axis('off')
plt.colorbar()
plt.plot([c/2-1, c/2-1],[0, r-1],'k')
plt.title('E min')

plt.subplot(224)
PLot = arr_W[2]+arr_W[3]
plt.imshow(PLot)
plt.plot(arr_LW[2][:,0],arr_LW[2][:,1],'rx')
plt.plot(arr_LW[3][:,0],arr_LW[3][:,1],'rx')
plt.axis('off')
plt.colorbar()
plt.plot([c/2-1, c/2-1],[0, r-1],'k')
plt.title('W min')
plt.show()

### Cáculo de FI

E_phi = (arr_E[0]+arr_E[1])*(arr_E[2]+arr_E[3])
W_phi = (arr_W[0]+arr_W[1])*(arr_W[2]+arr_W[3])

plt.figure()

plt.subplot(221)
plt.imshow(E_phi)
plt.axis('off')
plt.colorbar()
plt.plot([c/2-1, c/2-1],[0, r-1],'k')

plt.subplot(222)
plt.imshow(W_phi)
plt.axis('off')
plt.colorbar()
plt.plot([c/2-1, c/2-1],[0, r-1],'k')
plt.show()

##
LE_phi = phi_file['LE_phi'][0][0]
RE_phi = phi_file['RE_phi'][0][0]
LW_phi = phi_file['LW_phi'][0][0]
RW_phi = phi_file['RW_phi'][0][0]


print('----------------------------------------------------')
print(f'\nDivergence-Free Potentials: {LE_phi} | {RE_phi}\n')
print(f'\nCurl-Free Potentials: {LW_phi} | {RW_phi}\n')
