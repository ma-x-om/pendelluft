# load libraries
#
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import convolve as filter2
import cv2
import tkinter as tk
from tkinter import filedialog

#########################################
# Initialize tkinter (Tk)
# used for browsing files in explorer
#
root = tk.Tk()
root.withdraw()
#
#########################################

def get_magnitude(u, v):
    scale = 3
    sum = 0.0
    counter = 0.0

    for i in range(0, u.shape[0], 8):
        for j in range(0, u.shape[1],8):
            counter += 1
            dy = v[i,j] * scale
            dx = u[i,j] * scale
            magnitude = (dx**2 + dy**2)**0.5
            sum += magnitude

    mag_avg = sum / counter

    return mag_avg

def draw_quiver(u,v,beforeImg,OF_type):
    scale = 3
    ax = plt.figure().gca()
    ax.imshow(beforeImg, cmap = 'gray')
    #plt.axis('off')

    magnitudeAvg = get_magnitude(u, v)
    
    for i in range(0, u.shape[0], 8):
        for j in range(0, u.shape[1],8):
            dy = v[i,j] * scale
            dx = u[i,j] * scale
            magnitude = (dx**2 + dy**2)**0.5
            #draw only significant changes
            if magnitude > magnitudeAvg:
                ax.arrow(j,i, dx, dy, length_includes_head= True, head_width=2, color = 'red')

    plt.draw()
    plt.xlim([0,c])
    plt.ylim([r,0])
    plt.title(f"Quivers from {OF_type}")
    plt.savefig(f'{image_folder}{OF_type}_flow_quivers.png',bbox_inches='tight')
    plt.close()
    #plt.show()

def draw_comparison(UF, VF, UH, VH, img):
	scale = 3
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5), dpi=288)
	fig.suptitle("Farneback OF x Horn & Schunck OF")
	#ax1 = plt.figure().gca()
	ax1.imshow(img, cmap ='gray')
	#ax2 = plt.figure().gca()
	ax2.imshow(img, cmap ='gray')

	magnitudeAvg = get_magnitude(UF, VF)
    
	for i in range(0, UF.shape[0], 8):
		for j in range(0, UF.shape[1],8):
			dy = VF[i,j] * scale
			dx = UF[i,j] * scale
			magnitude = (dx**2 + dy**2)**0.5
			#draw only significant changes
			if magnitude > magnitudeAvg:
			   ax1.arrow(j,i, dx, dy, length_includes_head= True, head_width=2, color = 'red')

	magnitudeAvg = get_magnitude(UH, VH)
    
	for i in range(0, UH.shape[0], 8):
		for j in range(0, UH.shape[1],8):
			dy = VH[i,j] * scale
			dx = UH[i,j] * scale
			magnitude = (dx**2 + dy**2)**0.5
			#draw only significant changes
			if magnitude > magnitudeAvg:
			   ax2.arrow(j,i, dx, dy, length_includes_head= True, head_width=2, color = 'red')

	plt.draw()
	ax1.set_xlim([0,c])
	ax1.set_ylim([r,0])
	ax2.set_xlim([0,c])
	ax2.set_ylim([r,0])
	ax1.set_axis_off()
	ax2.set_axis_off()
	ax1.set_title('Farneback')
	ax2.set_title('Horn & Schunck')
	fig.savefig(f'{image_folder}Comparison_flow_quivers.png',bbox_inches='tight')

	plt.show()

def get_derivatives(img1, img2):
    #derivative masks
    x_kernel = np.array([[-1, 1], [-1, 1]]) * 0.25
    y_kernel = np.array([[-1, -1], [1, 1]]) * 0.25
    t_kernel = np.ones((2, 2)) * 0.25

    fx = filter2(img1,x_kernel) + filter2(img2,x_kernel)
    fy = filter2(img1, y_kernel) + filter2(img2, y_kernel)
    ft = filter2(img1, -t_kernel) + filter2(img2, t_kernel)

    return [fx,fy, ft]



#input: images name, smoothing parameter, tolerance
#output: images variations (flow vectors u, v)
#calculates u,v vectors and draw quiver
def computeHS(beforeImg, afterImg, alpha, delta):
    #removing noise
    beforeImg  = cv2.GaussianBlur(beforeImg, (5, 5), 0)
    afterImg = cv2.GaussianBlur(afterImg, (5, 5), 0)

    # set up initial values
    u = np.zeros((beforeImg.shape[0], beforeImg.shape[1]))
    v = np.zeros((beforeImg.shape[0], beforeImg.shape[1]))
    fx, fy, ft = get_derivatives(beforeImg, afterImg)
    avg_kernel = np.array([[1 / 12, 1 / 6, 1 / 12],
                            [1 / 6, 0, 1 / 6],
                            [1 / 12, 1 / 6, 1 / 12]], float)
    iter_counter = 0
    while True:
        iter_counter += 1
        u_avg = filter2(u, avg_kernel)
        v_avg = filter2(v, avg_kernel)
        p = fx * u_avg + fy * v_avg + ft
        d = 4 * alpha**2 + fx**2 + fy**2
        prev = u

        u = u_avg - fx * (p / d)
        v = v_avg - fy * (p / d)

        diff = np.linalg.norm(u - prev, 2)
        #converges check (at most 300 iterations) 
        if  diff < delta or iter_counter > 300:
            # print("iteration number: ", iter_counter)
            break

    #draw_quiver(u, v, beforeImg)

    return u, v

def draw_vectors(vector, Vec_type, OF_type):
	plt.figure(dpi=200)
	plt.imshow(vector)
	plt.axis('off')
	plt.colorbar()
	#plt.quiver(X, Y, u,-v)		# pyplot's function for plotting quivers
												# receives u and v for the respective frame
												#
	plt.title(f'{OF_type} - {Vec_type}')
	plt.savefig(f'{image_folder}{OF_type}_flow_{Vec_type}.png',bbox_inches='tight')
	plt.close()


# carregar o video (arquivo com frames)
file_path1 = filedialog.askopenfilename()	# abre explorer
file_name1 = os.path.split(file_path1)[1] # file path without '.txt'
file_path2 = filedialog.askopenfilename()	# abre explorer
file_name2 = os.path.split(file_path2)[1] # file path without '.txt'
main_folder = os.path.splitext(file_path1)[0]+'/'
isExist = os.path.exists(main_folder)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(main_folder)
image_folder = main_folder+'IMAGES/'
isExist = os.path.exists(image_folder)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(image_folder)

print('Wait some time...')
frame1 = cv2.imread(file_path1,0).astype(float)				# atribui o texto do arquivo à variável "dados"
frame2 = cv2.imread(file_path2,0).astype(float)			# atribui o texto do arquivo à variável "dados"
(r, c) = np.shape(frame1)		# atribui as dimensões nº de frames e de pixels por frame

Imagens = np.empty([r,c,2])			# great array to hold the future images
tempoIm = []								# initialize empty lists
tempocv = []


# here we calculate the optical flow
# let's use Farneback instead of H&S
# Farneback Optical Flow Calc with default parameters
flow = cv2.calcOpticalFlowFarneback(frame1, frame2, None, 0.5, 3, 15, 3, 5, 0, 0) #1.2 is used as 0 in john's code
# we then split "flow" array components into u's and v's corresponding "frame"
u_Farneback, v_Farneback = flow[:,:,0], flow[:,:,1]
print("Farneback calculated.")

u_HS, v_HS = computeHS(frame2, frame1, alpha = 15, delta = 10**-1)
print("HornSchunck calculated.")

draw_vectors(u_Farneback, 'U', 'Farneback')
draw_vectors(v_Farneback, 'V', 'Farneback')
draw_vectors(u_HS, 'U', 'HornSchunck')
draw_vectors(v_HS, 'V', 'HornSchunck')

draw_quiver(u_Farneback, v_Farneback, frame1, 'Farneback')
draw_quiver(u_HS, v_HS, frame1, 'HornSchunck')

draw_comparison(u_Farneback, v_Farneback, u_HS, v_HS, frame1)


print("The end?")