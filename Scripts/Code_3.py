## Ties together the following .m script:
#	1. Principal.m
#		1.1 Posicoes.m
#		1.2 freqmap.m

##	Receives:	'.mat' files from 'Code_1.py' and 'Code_2.py'
#					'Images_U_V.mat'
#					'potenciaisE_W.mat'
#				'LungMask.mat' -> máscara para seccionar a área do pulmão nas imagens e centro do pulmão
#					Calculados automaticamente em "Code_2.py" através do script "Automatic_Segmentation.py"
##	Calls:	'Posicoes.py'
##			'freqmap.py'
##	Outputs:	calculated data into a file

## CHANGES FROM ITS MATLAB COUNTERPART
#	1. Instead of simultaneously processing pigs N and P, this code will process only the 
#		chosen file (from whatever pig or being it may come from)
#	2. It seems that Principal.m does not save its highly precious calculated values in a
#		specific file. That was changed in here. 

# load libraries
#
from scipy.io import loadmat
from scipy.io import savemat
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os
import matplotlib.pyplot as plt
import skimage.morphology as skim
# load external scripts
#
import posicoes as CalcPos
import freqmap as CalcFreq

#########################################
# Initialize tkinter (Tk)
# used for browsing files in explorer
#
root = tk.Tk()
root.withdraw()
#
#########################################

# let's load that precious '.mat' file
file_path_pot = filedialog.askopenfilename(title = "Abra o arquivo dos potenciais (potenciaisE_W.mat)", filetypes=(("Arquivos MAT","*.mat"),("Todos os arquivos","*.*")))	# open explorer
pot_folder = os.path.split(file_path_pot)[0]+'/' #	pasta que contém o arquivo escolhido
pot_file = loadmat(file_path_pot)	# open the file in python

potE = pot_file['potE']
potW = pot_file['potW']

r, c, N_length = np.shape(potE) # r e c = dimensoes do array de potenciais e n_length = nº de frames - 1

''' Since right now we have a lonely file for the lung mask, I deemed it pratical to let such file in the same folder
print('Abra o arquivo da máscara do pulmão\n\n')
file_path_mask = filedialog.askopenfilename()	# open explorer
mask_folder = os.path.splitext(file_path_mask)[0]+'/' #	pasta que contém o arquivo escolhido
mask_file = loadmat(file_path_mask)	# open the file in python
'''

mask_file = loadmat(pot_folder+'LungMask.mat')
mask = mask_file['BW']
mask_center = int(mask_file['Center'])
Right = mask.copy()
Right[:,0:mask_center] = 0
Left = mask.copy()
Left[:,((mask_center+1)):r] = 0

######### here we begin to retrieve maximi and minimi from PotE and PotW
#########	and separate it in Left and Right
arr_pot = [np.zeros((r,c,N_length)), np.zeros((r,c,N_length)), np.zeros((r,c,N_length)), np.zeros((r,c,N_length))]

arr_potE = [np.zeros((r,c,N_length)), np.zeros((r,c,N_length)), np.zeros((r,c,N_length)), np.zeros((r,c,N_length))]
arr_potW = [np.zeros((r,c,N_length)), np.zeros((r,c,N_length)), np.zeros((r,c,N_length)), np.zeros((r,c,N_length))]

for i in np.arange(N_length):
	arr_pot[0][:,:,i] = skim.local_maxima(potE[:,:,i],connectivity=8)
	arr_pot[1][:,:,i] = skim.local_minima(potE[:,:,i],connectivity=8)
	arr_pot[2][:,:,i] = skim.local_maxima(potW[:,:,i],connectivity=8)
	arr_pot[3][:,:,i] = skim.local_minima(potW[:,:,i],connectivity=8)

	arr_potE[0][:,:,i] = arr_pot[0][:,:,i]*Left
	arr_potE[1][:,:,i] = arr_pot[0][:,:,i]*Right
	arr_potE[2][:,:,i] = arr_pot[1][:,:,i]*Left
	arr_potE[3][:,:,i] = arr_pot[1][:,:,i]*Right

	arr_potW[0][:,:,i] = arr_pot[2][:,:,i]*Left
	arr_potW[1][:,:,i] = arr_pot[2][:,:,i]*Right
	arr_potW[2][:,:,i] = arr_pot[3][:,:,i]*Left
	arr_potW[3][:,:,i] = arr_pot[3][:,:,i]*Right


# posE and posW lists, which are indexes for maxima and minima of each potential
list_posE = []
list_posW = []

for i in range(4):
	tempE = CalcPos.calc_PosEmin(potE,arr_potE[i])
	tempW = CalcPos.calc_PosEmin(potW,arr_potW[i])
	list_posE.append(tempE.copy())
	list_posW.append(tempW.copy())
	print('*', end='')


#####
#####		Great tool for checking if the plot is the same as Matlab's
'''
print('Q\n 0 = L E MAX | 1 = R E MAX | 2 = L E MIN | 3 = R E MIN\n')
print(f'frame: 0 to {N_length-1}\n')
def check_mm(q, f):
	print(list_posE[q][:,:,f])
	print(f'Q{q} frame {f+1}\n\n')

loop = 'y'
while loop=='y' or loop=='':
	qq, ff = input('QQ FF: ').split()
	qq=int(qq)
	ff=int(ff)
	check_mm(qq, ff)
	loop = input("Proceed? (y/n) ")
'''
#####
#####

### Frequency Maps

def go_min(posmin, pot, D=6, rad=5):
	minimo = 10
	L, Cmin = CalcFreq.calc_freq(posmin, D, rad, r, c)
	Lmin = np.zeros(np.shape(L)[0])
	length = np.shape(pot)[2]
	for i in np.arange(np.shape(L)[0]):
		for j in np.arange(length):
			if pot[L[i,0],L[i,1],j] < minimo:
				minimo = pot[L[i,0],L[i,1],j]
		Lmin[i] = minimo
		minimo = 10

	print('*', end='')
	return L, Cmin, Lmin

def go_max(posmax, pot, D=6, rad=5):
	maximo = -5
	L, Cmax = CalcFreq.calc_freq(posmax, D, rad, r, c)
	Lmax = np.zeros(np.shape(L)[0])
	length = np.shape(pot)[2]
	for i in np.arange(np.shape(L)[0]):
		for j in np.arange(length):
			if pot[L[i,0],L[i,1],j] > maximo:
				maximo = pot[L[i,0],L[i,1],j]
		Lmax[i] = maximo
		maximo = -5

	print('*', end='')
	return L, Cmax, Lmax


# another family of arrays stored in lists (WARNING: they will get very big)
arr_LE = []
arr_LW = []
arr_E = []
arr_W = []
arr_LEm = [] #  this and the next declared variable are used only to print the highest
arr_LWm = [] #		and the lowest value in the plot (Used only in Victor's matlab code "Principal.m")

#Q1 <-> i=0 <-> L MAX
#Q2 <-> i=1 <-> R MAX
#Q3 <-> i=2 <-> L MIN
#Q4 <-> i=3 <-> R MIN

for i in range(2):
	tempE = go_max(list_posE[i],potE)
	tempW = go_max(list_posW[i],potW)

	arr_LE.append(tempE[0].copy())
	arr_E.append(tempE[1].copy())
	arr_LEm.append(tempE[2].copy())
	arr_LW.append(tempW[0].copy())
	arr_W.append(tempW[1].copy())
	arr_LWm.append(tempW[2].copy())

	print(']')

for i in range(2,4):
	tempE = go_min(list_posE[i],potE)
	tempW = go_min(list_posW[i],potW)

	arr_LE.append(tempE[0].copy())
	arr_E.append(tempE[1].copy())
	arr_LEm.append(tempE[2].copy())
	arr_LW.append(tempW[0].copy())
	arr_W.append(tempW[1].copy())
	arr_LWm.append(tempW[2].copy())

	print(']')
'''
good_day = np.array(arr_LW, dtype=object)
print(f"Good day! Of a shape of {np.shape(good_day)}")
print(good_day)
'''
# Here lies image showing time again, kids!

plt.figure()
plt.subplot(221)
PLot = arr_E[0]+arr_E[1]
plt.imshow(PLot)
plt.plot(arr_LE[1][:,0],arr_LE[1][:,1],'kx')
plt.plot(arr_LE[0][:,0],arr_LE[0][:,1],'kx')
plt.axis('off')
plt.colorbar()
plt.axvline(mask_center,color='k')
plt.title('E max')

plt.subplot(222)
PLot = arr_W[0]+arr_W[1]
plt.imshow(PLot)
plt.plot(arr_LW[1][:,0],arr_LW[1][:,1],'kx')
plt.plot(arr_LW[0][:,0],arr_LW[0][:,1],'kx')
plt.axis('off')
plt.colorbar()
plt.axvline(mask_center,color='k')
plt.title('W max')

plt.subplot(223)
PLot = arr_E[2]+arr_E[3]
plt.imshow(PLot)
plt.plot(arr_LE[2][:,0],arr_LE[2][:,1],'kx')
plt.plot(arr_LE[3][:,0],arr_LE[3][:,1],'kx')
plt.axis('off')
plt.colorbar()
plt.axvline(mask_center,color='k')
plt.title('E min')

plt.subplot(224)
PLot = arr_W[2]+arr_W[3]
plt.imshow(PLot)
plt.plot(arr_LW[2][:,0],arr_LW[2][:,1],'kx')
plt.plot(arr_LW[3][:,0],arr_LW[3][:,1],'kx')
plt.axis('off')
plt.colorbar()
plt.axvline(mask_center,color='k')
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
plt.axvline(mask_center,color='k')

plt.subplot(222)
plt.imshow(W_phi)
plt.axis('off')
plt.colorbar()
plt.axvline(mask_center,color='k')
plt.show()

##
LE_phi = sum(sum(E_phi[:,0:15]))
RE_phi = sum(sum(E_phi[:,16:31]))
LW_phi = sum(sum(W_phi[:,0:15]))
RW_phi = sum(sum(W_phi[:,16:31]))


#isExist = os.path.exists(pot_folder)
#if not isExist:
   # Create a new directory because it does not exist
#   os.makedirs(pot_folder)
savemat((pot_folder+'calculos_phi.mat'),{'LE_phi': LE_phi, 'RE_phi':RE_phi,'LW_phi': LW_phi, 'RW_phi': RW_phi})
savemat((pot_folder+'calculos_pos.mat'),{'list_posE': np.array(list_posE,dtype=object), 'list_posW':np.array(list_posW,dtype=object)})
savemat((pot_folder+'calculos_freqmap.mat'),{'arr_E': arr_E, 'arr_W':arr_W, 'arr_LWm': np.array(arr_LWm,dtype=object), 'arr_LE':np.array(arr_LE,dtype=object), 'arr_LW':np.array(arr_LW,dtype=object), 'arr_LEm':np.array(arr_LEm,dtype=object)})


print('\n\n----------------------------------------------------')
print(f'\nDivergence-Free Potentials: {LE_phi} | {RE_phi}\n')
print(f'\nCurl-Free Potentials: {LW_phi} | {RW_phi}\n')
