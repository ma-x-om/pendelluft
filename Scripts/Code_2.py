## Ties together the following .m scripts:
#	1. Gera_Potenciais.m
#	2. potential_calculation2.m <-> potential_calculation2.py (John's)
#		2.1 runDHHD.py (John's | Calls 'potential_calculation2.py')

##	Receives:	'.mat' file from 'Code_1.py'
##	Calls:	'runDHHD.py'
##			'potential_calculation2.py'
##	Outputs:	updated '.mat' file with Divergence-Free (PotE) and Curl-Free (PotW) components

# load libraries
#
from scipy.io import loadmat
from scipy.io import savemat
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os
import time

# load external scripts
#
import runDHHD as rdhhd

#########################################
# Initialize tkinter (Tk)
# used for browsing files in explorer
#
root = tk.Tk()
root.withdraw()
#
#########################################

# let's load that precious '.mat' file
file_path = filedialog.askopenfilename(title = "Abra o arquivo gerado em Code_1.py (Images_U_V.mat)", filetypes=(("Arquivos MAT","*.mat"),("Todos os arquivos","*.*")))	# open explorer
main_folder = os.path.split(file_path)[0]+'/' #	pasta que contém o arquivo escolhido
Mat_File = loadmat(file_path)	# open the file in python

# decompose that Mat_File into MFVs u and v
u = Mat_File['u']
v = Mat_File['v']

r, c, N_length = np.shape(u) # r e c = dimensoes do MFV e n_length = nº de frames - 1
#X, Y = np.meshgrid(np.arange(r),np.arange(c))	# meshgrid like in 'Code_1.py'

#N_length = 20	# we can use it as a temporary measure so that we don't have to use all of the frames

# define empty arrays for components
potE = np.empty([r,c,N_length])	# Divergence-Free potential
potW = np.empty([r,c,N_length])	# Curl-Free potentials

initial_time_seconds = time.time()

for i in np.arange(N_length):
	potE[:,:,i], potW[:,:,i] = rdhhd.runDHHD(u[:,:,i],v[:,:,i]) #[2:4] with original runDHHD; 'show_time' = True to show time to exit potential2
	print(f'Análise {(i+1)}/{N_length+1}', end='\r')

delta_time_seconds = time.time() - initial_time_seconds

print(f"\nFim das análises\n Tempo levado = {delta_time_seconds:.0f} s ou {delta_time_seconds/60:.0f} min")

savemat((main_folder+'potenciaisE_W.mat'),{'potE': potE, 'potW':potW})

print(f'Análises salvas em {main_folder}\n')