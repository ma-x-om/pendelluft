import numpy as np

def calc_PosEmin(NpotE, NpotEm):
	A = 0
	L = np.shape(NpotE)[2] # gets the frames number

	for i in np.arange(L):
		r, c = np.nonzero(NpotEm[:,:,i]==1)
		B = np.shape(r)[0] # qtd de pontos extremos para o par de frames i
		if B > A:
			A = B # guarda em A o maior valor de pontos extremos para um frame entre todos os frames de NpotEm

	PosEm = np.zeros((A,2,L)) # matriz 3D para guardar os indices dos extremos de cada frame
	PosEm[:] = np.nan

	for i in np.arange(L):
		r, c = np.nonzero(NpotEm[:,:,i]==1)
		B = np.shape(r)[0]
		PosEm[0:B,:,i] = np.array((r,c)).T

	return PosEm