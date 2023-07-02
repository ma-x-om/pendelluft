## This code will calculate the Pendelluft degree in an specified 15 seconds analysis window

import cv2
import numpy as np
from numpy import shape
import timeit
from scipy.io import loadmat

def load_files():
	import tkinter as tk
	from tkinter import filedialog
	import os

	root = tk.Tk()
	root.withdraw()

	file_path_freqmap = filedialog.askopenfilename(title = "Abra o arquivo do mapa de frequências (calculos_freqmap.mat)", filetype=(("Arquivos MAT","*.mat"),("Todos os arquivos","*.*")))	# open explorer
	global start_time
	start_time = timeit.default_timer()
	
	freqmap_file = loadmat(file_path_freqmap)	# open the file in python
	freqmap_folder = os.path.split(file_path_freqmap)[0]+'/'

	mask_file = loadmat(freqmap_folder+'LungMask.mat')
	mask = mask_file['BW']


	arr_LW = freqmap_file['arr_LW'][0]
	arr_W = freqmap_file['arr_W']

	return arr_W, arr_LW, mask

def proccess_arrays(arr_W, arr_LW, mask):
	arr_Wmax = arr_W[0]+arr_W[1]
	arr_Wmin = arr_W[2]+arr_W[3]
	arr_Wmaxpoints = np.concatenate((arr_LW[0], arr_LW[1]), axis =0)
	arr_Wminpoints = np.concatenate((arr_LW[2], arr_LW[3]), axis =0)

	return arr_Wmax, arr_Wmin, arr_Wmaxpoints, arr_Wminpoints, mask


### DEFININDO FUNÇÕES

def haus_distance_table(array_A, array_B):		#	Responsável por gerar uma tabela que relaciona os pontos do array A com os pontos do array B,
	r_A, c_A = shape(array_A)					#		calculando para nós o quadrado da distância entre cada um deles e a razão do valor B/A
	r_B, c_B = shape(array_B)
	A_max = array_A[0][2]
	B_max = array_B[0][2]

	distances_array = np.zeros((r_B, r_A, 2))
	text_parameters = np.empty((r_B, r_A), dtype=object)
	for i in range(r_A):
		for j in range(r_B):
			distancia = (array_B[j][0]-array_A[i][0])**2+(array_B[j][1]-array_A[i][1])**2
			distances_array[j][i][0] = distancia
			ratio_from_max_A = 100*array_A[i][2]/A_max
			ratio_from_max_B = 100*array_B[j][2]/B_max
			ratio_BA = ratio_from_max_B/ratio_from_max_A
			distances_array[j][i][1] = ratio_BA
			desv_BA = 100*abs(ratio_BA-1)
			list_input = f"{distancia:.0f}: {desv_BA:.0f}% ({100*ratio_BA:.0f}%)"
			text_parameters[j][i] = list_input

	print(f'Tabela de parâmetros: (distância² | Desvio de valor relativo | MinRelativ/MaxRelativ)\n{text_parameters}')

	return distances_array

def Cm_calc(m_REF, M_GOAL, dist):				#	Calcula Cmax ou Cmin
	Cm_here = np.float64(M_GOAL/m_REF)/dist

	return Cm_here

def find_pairs(array_REF, array_GOAL):			#	Encontra os pares entre array_REF e array_GOAL. Descrito em mais detalhes em seu Readme.
	distance_table = haus_distance_table(array_REF, array_GOAL)
	r, c, d = shape(distance_table) # r = nº de pontos de destino (linhas da tabela de distâncias)
									# c = nº de pontos de referência (colunas da tabela)
									# d = dimensões (0 = distância e 1 = relação)

	print(f"\n Identificando os {c} pares...")
	pares_array = np.empty((c,3))
	lista_refs = np.arange(c) 	# lista com os indices dos objetos de referência
	lista_goals = np.arange(r)
	for i in lista_refs: # itera por coluna (de referencia)
		minimo_one = np.amin(distance_table[:,i,0])	#	'one' pois é a primeira tentativa de identificar os menores valores de distâncias entre o ref e o GOAL
		indexes = np.where(distance_table[:,i,0] == minimo_one) # aqui e na seguinte linha conferimos se a distância mínima encontrada se repete para outro GOAL
		_, n_of_mins = shape(indexes)
		
		minimos_one_index = np.empty((n_of_mins, 2))	# array para guardar os indices das distâncias minimas encontradas
		for k in range(n_of_mins):
			minimos_one_index[k] = [indexes[0][k],i]

		# agora vamos verificar qual GOAL é "melhor", pois REF só pode fazer par com 1 GOAL
		# baseando-se numa relação entre o valor que esse pixel tem em array_REF e array_GOAL
		melhor_minimo_one = np.array([0, 0, minimo_one, 100]) # (y (min), x (max), distancia, desvio)
		for k in minimos_one_index:
			y = int(k[0])
			x = int(k[1])
			desvio = abs(1 - distance_table[y][x][1])
			if desvio < melhor_minimo_one[3]:
				melhor_minimo_one = [y, x, minimo_one, desvio]
		y = melhor_minimo_one[0]	# este é o índice da linha que representa o melhor GOAL

		pares_array[i] = melhor_minimo_one[:3]

	return pares_array

def plot_Pd_I(pares, PD_idx, relation, differentiate = False):			#	Faz a plotagem. 'pares' = endereços dos pares; 'PD_idx' = Índice de Pendelluft;
	import matplotlib.pyplot as plt 									#		'relation' = se é mín/máx (i/a) ou máx/mín (a/i)
	Quad, r, c = shape(arr_W) # Quad = quadrants						#		'differentiate' = True or False: se ligado mostrará os pontos de Wmax e Wmin
	#vamos normalizar:													#		com cores diferentes ao realizar o plt.imshow()
	arr_Wmax_normalized = arr_Wmax/np.max(arr_Wmax)
	arr_Wmin_normalized = arr_Wmin/np.max(arr_Wmin)

	plt.figure("Correlação entre máximos e mínimos dos potenciais Curl-Free")

	if differentiate:
		link_color = 'm'
		mask_contour, _ = cv2.findContours(mask, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)

		EvalContour = np.squeeze(mask_contour) #reduce dimensionality of array
		EvalpointsList=[]
		for row in EvalContour:
		    EvalpointsList.append(tuple(row))
		arr_Wfull = arr_Wmax_normalized - arr_Wmin_normalized
		plt.imshow(arr_Wfull, cmap='BrBG', vmin=-1, vmax=1)
		for p in EvalpointsList:
		  plt.plot(p[0],p[1], 'k.', label="Limite pulmonar")
		plt.colorbar()
	else:
		link_color = 'w'
		arr_Wfull = arr_Wmax_normalized + arr_Wmin_normalized + mask*0.20
		plt.imshow(arr_Wfull, vmin=0, vmax=1)		

	plt.plot(arr_LW[1][:,0],arr_LW[1][:,1], 'rX', markerfacecolor="none", ms = 10, markeredgewidth=2, label="Máximo")
	plt.plot(arr_LW[0][:,0],arr_LW[0][:,1], 'rX', markerfacecolor="none", ms = 10, markeredgewidth=2)
	plt.plot(arr_LW[2][:,0],arr_LW[2][:,1], 'ob', markerfacecolor="none", ms = 10, markeredgewidth=2, label="Mínimo")
	plt.plot(arr_LW[3][:,0],arr_LW[3][:,1], 'ob', markerfacecolor='none', ms = 10, markeredgewidth=2)

	for i in pares:
		if np.sum((i[1]-i[0])) == 0:
			plt.plot([i[0][1], i[1][1]],[i[0][0], i[1][0]],f'o{link_color}', label="Relação máx e mín")
		else:
			plt.plot([i[0][1], i[1][1]],[i[0][0], i[1][0]],link_color, label="Relação máx e mín")

	plt.axis('off')
	plt.plot([c/2-1, c/2-1],[0, r-1],'k')
	handles, labels = plt.gca().get_legend_handles_labels()
	by_label = dict(zip(labels, handles))
	plt.legend(by_label.values(), by_label.keys())

	if relation == "a/i":
		plt.title(f'Índice de Pendelluft: {PD_idx:.2f}\nTipo de correlação: mín -> máx')
	else:
		plt.title(f'Índice de Pendelluft: {PD_idx:.2f}\nTipo de correlação: máx -> mín')

	plt.tight_layout()
	plt.show()

def Pd_I_Calc(array_max, array_min, relation="i/a"):
	np.seterr(divide='ignore')
	if relation == "i/a":
		a_REF = array_max
		a_GOAL = array_min
		pares_header_string = f"\nPares identificados:\n|   MAX    |   MIN   | Distância | Valor Relativo | Cmax |"
	elif relation == "a/i":
		a_REF = array_min
		a_GOAL = array_max
		pares_header_string = f"\nPares identificados:\n|   MIN    |   MAX   | Distância | Valor Relativo | Cmin |"

	pairs_start_time = timeit.default_timer()
	lista_pares = find_pairs(a_REF, a_GOAL)
	n_de_pares = shape(lista_pares)[0]
	Cm = np.empty((n_de_pares))
	Cm_sum = 0
	pares_address = np.empty((n_de_pares, 2, 2))

	print(pares_header_string)
	for i in range(n_de_pares):
		pair_REF = a_REF[int(lista_pares[i][1])]
		pair_GOAL = a_GOAL[int(lista_pares[i][0])]
		Cm[i] = Cm_calc(pair_REF[2], pair_GOAL[2], lista_pares[i][2])
		print(f"  {pair_REF[:2]} {pair_GOAL[:2]}      {np.sqrt(lista_pares[i][2]):.0f}           {pair_GOAL[2]/pair_REF[2]:.2f}        {Cm[i]:.3f}")
		Cm_sum += 1/Cm[i]
		pares_address[i] = pair_REF[:2], pair_GOAL[:2]

	Pd_index = Cm_sum/n_de_pares
	print(f"\nÍndice de Pendelluft: {Pd_index:.2f}")
	end_time = timeit.default_timer()
	delta_time_total = end_time - start_time
	delta_time_pares = end_time - pairs_start_time
	time_string = f"\nTempo total levado = {delta_time_total:.6f} s ou {1000*delta_time_total:.3f} ms\n\tTempo levado para processamento de {n_de_pares} pares = {1000*delta_time_pares:.3f} ms\n\t\tMédia de {1000*delta_time_pares/n_de_pares:.3f} ms/Par"
	print(time_string)

	plot_Pd_I(pares_address, Pd_index, relation)


def organize_values(arr_points, arr_W, tipo):	#	Organiza os endereços dados em 'arr_points' com suas respectivas intensidades de 'arr_W' em um array
	r,c = shape(arr_points)						#		único 'arr_points_values'. 'tipo' serve somente como String para os prints da função.
	arr_points_values = np.zeros((r,c+1))
	for i in range(r):
		arr_points_values[i] = [arr_points[i][1], arr_points[i][0], arr_W[arr_points[i][1]][arr_points[i][0]]]
	arr_points_values=arr_points_values[arr_points_values[:,2].argsort()[::-1]]

	#pretty print
	print(f"\nPontos de {tipo} e seus valores:\n {r} pontos\n| n | Y  | X  | Value |")
	for i in range(r):
		print(f"  {i}   {arr_points_values[i][0]:.0f}   {arr_points_values[i][1]:.0f}    {arr_points_values[i][2]:.0f}")
	return arr_points_values


###### AQUI COMEÇA ######
arr_W, arr_LW, mask = load_files()
#start_time = timeit.default_timer()
arr_Wmax, arr_Wmin, arr_Wmaxpoints, arr_Wminpoints, _ = proccess_arrays(arr_W, arr_LW, mask)

### Maximos de W (Curl-Free Potential)
arr_Wmaxpoints_values = organize_values(arr_Wmaxpoints, arr_Wmax, "máximo")

### Minimos de W
arr_Wminpoints_values = organize_values(arr_Wminpoints, arr_Wmin, "mínimo")

print('\n')
Pd_I_Calc(arr_Wmaxpoints_values, arr_Wminpoints_values, 'i/a')