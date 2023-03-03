from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()

def load_img():
	path_normal = filedialog.askopenfilename()
	dados_normal = loadmat(path_normal)

	return dados_normal

def segment(dados_normal):
	fatia = None
	#fatia = 30

	if not fatia:
		fatia = int(input('Qual será o frame?\n'))

	img = dados_normal['Imagens'][:,:,fatia]
	u = dados_normal['u'][:,:,fatia]
	v = dados_normal['v'][:,:,fatia]

	return img, u, v

def plotar(img, u, v):
	r, c = np.shape(u)
	x, y = np.meshgrid(np.arange(r),np.arange(c))
	plt.figure(dpi=300, figsize=(2,2))
	plt.imshow(img)
	plt.axis('off')
	plt.colorbar()
	plt.quiver(x, y, u,-v)
	plt.show()

def plotar_hidden(img, x, y, u, v, frame):
	plt.figure(dpi=300, figsize=(2,2))
	plt.imshow(img)
	plt.axis('off')
	plt.colorbar()
	plt.quiver(x, y, u,-v)
	plt.savefig(f'Renders/Comparisons/Frame_{frame}.png',bbox_inches='tight')
	plt.close()

def plotar_intervalo(dados_normal):
	f0 = int(input('Insira o intervalo de frames que deseja:\n[Frame Inicial]: '))
	ff = int(input('[Frame Final]: '))
	img = dados_normal['Imagens']
	u = dados_normal['u']
	v = dados_normal['v']
	valo = np.arange(f0,ff+1)
	r, c = np.shape(u[:,:,f0])
	x, y = np.meshgrid(np.arange(r),np.arange(c))
	print(f'Ploting {ff-f0+1} frames: {f0} to {ff}\n')
	for frame in valo:
		plotar_hidden(img[:,:,frame], x, y, u[:,:,frame], v[:,:,frame], frame)
	print('Saved.\n')


def main():
	repetir = True
	new_norm=True
	while repetir == True:
		if new_norm == True:
			dados = load_img()
		img, u, v = segment(dados)
		plotar(img, u, v)
		#plotar(segment(load_img()))
		q = input('Repetir? (s/ss/n)\n')
		if q == 'n':
			repetir = False
		if q == 'f':
			new_norm = False
		if q == 's':
			new_norm = True
		if q == 'p':
			plotar_intervalo(dados)

main()