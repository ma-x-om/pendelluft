import numpy as np

def calc_PosEmin(NpotE, NpotEm):
	A = 0
	L = np.shape(NpotE)[2]

	for i in np.arange(L):
		r, c = np.nonzero(NpotEm[:,:,i]==1)
		B = np.shape(r)[0]
		if B > A:
			A = B

	PosEm = np.zeros((A,2,L))
	PosEm[:] = np.nan

	for i in np.arange(L):
		r, c = np.nonzero(NpotEm[:,:,i]==1)
		B = np.shape(r)[0]
		PosEm[0:B,:,i] = np.array((r,c)).T

	return PosEm