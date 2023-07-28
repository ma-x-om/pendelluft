import numpy as np
import skimage.morphology as skim

def calc_freq(Pos, D, rad, rows, columns):

	FreqMap = np.zeros((rows, columns))

	a = np.shape(Pos)[2] # n° de frames

	for i in np.arange(a):
		x = Pos[:,0,i]
		x=x[np.nonzero(~np.isnan(x))] # pega os valores de x que são diferentes de NaN
		x = np.array(x, dtype=int)	
		y = Pos[:,1,i]
		y=y[np.nonzero(~np.isnan(y))]
		y = np.array(y, dtype=int)
		if x.size != 0:
			for j in np.arange(x.size):
				FreqMap[x[j],y[j]] += 1

	for i in np.arange(rows):
		for j in np.arange(columns):
			if FreqMap[i,j] != 0 and FreqMap[i,j] < D:
				FreqMap[i,j] = 0

	B = skim.local_maxima(FreqMap,connectivity=8)
	r, c = np.nonzero(B==1)
	CP = np.zeros((rows,columns,np.shape(r)[0]))

	X, Y = np.meshgrid(np.arange(rows),np.arange(columns))

	for i in np.arange(np.shape(r)[0]):
		CP[:,:,i] = (Y - r[i])**2 + (X - c[i])**2 <= rad**2

	RP = np.zeros((rows,columns))

	for i in np.arange(rows):
		for j in np.arange(columns):
			if sum(CP[i,j,:]) != 0:
				RP[i,j] = 1

	FreqMap = FreqMap*RP
	L = np.array((c, r), dtype=np.int64).T

	print('*', end='')
	return L, FreqMap