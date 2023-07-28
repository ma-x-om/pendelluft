import matplotlib.pyplot as plt
import numpy as np
import cv2
import tkinter as tk
from tkinter import filedialog
import imageio

root = tk.Tk()
root.withdraw()

def load_files():
    file_path = filedialog.askopenfilename(multiple=True)
    
    return file_path

def get_img(file, visible = True):
    img = cv2.imread(file,0).astype(float)
    r, c = np.shape(img)
    if visible:
        plt.imshow(img, cmap='gray')
        plt.axis('off')
        plt.show()

    return img, r, c

def write_svd(video_file, max_brightness):
    im_dir = filedialog.askdirectory()+'/'
    _,_,frames = np.shape(video_file)
    for fr in np.arange(frames):
        plt.imsave(f"{im_dir}SVD_frame{1+fr:03.0f}.png", video_file[:,:,fr], cmap='gray', vmin=0, vmax= max_brightness)
    print(f'images written to {im_dir}')
    write_gif(im_dir)

def write_gif(im_dir):
    gif_file = im_dir+'SVD.gif'
    with imageio.get_writer(gif_file, mode='I') as writer:
        for filename in os.listdir(im_dir):
            if filename.endswith(".png"):
                image = imageio.v2.imread(im_dir+filename)
                writer.append_data(image)
    print(f'gif written to {gif_file}')

def find_nearest(array, value):
    array = np.asarray(array)
    values = [value]
    if np.shape(np.shape(values))[0] > 1:
        values = value
    idx = []
    for i in values:
        idx.append((np.abs(array - i)).argmin())
    return idx

def do_image_svd(percent=None, level=None):
    arquivo = load_files()
    img, r, c = get_img(arquivo[0])

    u, s, vh = np.linalg.svd(img, full_matrices=False)
    S = np.diag(s)
    cum_sum = np.cumsum(np.diag(S))/np.sum(np.diag(S))
    n_levels = np.shape(cum_sum)[0]
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,3))
    ax1.semilogy(np.diag(S))
    ax1.set_title("Singular Values")
    ax2.plot(cum_sum)
    ax2.set_title("Cummulative Sum")
    plt.tight_layout()
    plt.show()
    
    if percent:
        percent = np.array(percent)/100
        Scomps = np.array(find_nearest(cum_sum, percent))+1
    elif level:
        Scomps = level
    else:
        Scomps = np.array(find_nearest(cum_sum, [0.35, 0.60, 0.80]))+1
    
    j = 0
    for i in Scomps:
        Aq3 = u[:,:i]@S[:i,:i]@vh[:i,:]
        plt.figure(j+1)
        j += 1
        plt.imshow(Aq3, cmap='gray')
        plt.title(f'{cum_sum[i-1]*100:.2f}%, {i}/{n_levels} levels')
        plt.axis('off')
        plt.show()

def do_video_svd(percent=None, level=None, frame = 10, render = False):
    arquivos = load_files()
    n_frames = np.shape(arquivos)[0]
    if frame > n_frames:
        frame = n_frames-1
    _, r, c = get_img(arquivos[0], False)
    images_array = np.empty([r,c,n_frames])
    for i in np.arange(n_frames):
        images_array[:,:,i] , _, _ = get_img(arquivos[i], False)
        
    B = 0
    for it in np.arange(n_frames):
        A = images_array[:,:,it].max()
        if A>B:
            B = A
            
    images_reshaped = images_array.reshape([r*c,n_frames])

    u, s, vh = np.linalg.svd(images_reshaped, full_matrices=False)
    S = np.diag(s)
    cum_sum = np.cumsum(np.diag(S))/np.sum(np.diag(S))
    n_levels = np.shape(cum_sum)[0]
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,3))
    ax1.semilogy(np.diag(S))
    ax1.set_title("Singular Values")
    ax2.plot(cum_sum)
    ax2.set_title("Cummulative Sum")
    plt.tight_layout()
    plt.show()
    
    if percent:
        percent = np.array(percent)/100
        Scomps = np.array(find_nearest(cum_sum, percent))+1
    elif level:
        Scomps = level
    else:
        Scomps = np.array(find_nearest(cum_sum, [0.35, 0.60, 0.80]))+1
        
    for i in Scomps:
        Aq3 = u[:,:i]@S[:i,:i]@vh[:i,:]
        Aq32 = Aq3.reshape([r,c,n_frames])
        
        Bl = 0
        for it in np.arange(n_frames):
            Al = Aq32[:,:,it].max()
            if Al>Bl:
                Bl = Al
                
        fig2, (ax3, ax4) = plt.subplots(1,2)
        #plt.figure()
        ax3.imshow(Aq32[:,:,frame], cmap='gray',vmin=0,vmax=Bl)
        ax3.set_axis_off()
        ax3.set_title(f'{cum_sum[i-1]*100:.2f}%, {i}/{n_levels} levels\n Frame {frame}')
        ax4.imshow(images_array[:,:,frame], cmap='gray',vmin=0,vmax=B)
        ax4.set_axis_off()
        ax4.set_title(f'Input video\n Frame {frame}')
        plt.tight_layout()
        plt.show()
    
    if render:
        write_svd(Aq32, Bl)

def apply_video_svd(images_array, percent=None, level=None, Brightness = False):
    r, c, n_frames =np.shape(images_array)
            
    images_reshaped = images_array.reshape([r*c,n_frames])

    u, s, vh = np.linalg.svd(images_reshaped, full_matrices=False)
    S = np.diag(s)
    cum_sum = np.cumsum(np.diag(S))/np.sum(np.diag(S))
    
    if percent:
        percent = np.array(percent)/100
        Scomps = np.array(find_nearest(cum_sum, percent))+1
    elif level:
        Scomps = level
    else:
        Scomps = np.array(find_nearest(cum_sum, [0.90]))+1
        
    for i in Scomps:
        Aq3 = u[:,:i]@S[:i,:i]@vh[:i,:]
        Aq32 = Aq3.reshape([r,c,n_frames])
        
        Bl = 0
        for it in np.arange(n_frames):
            Al = Aq32[:,:,it].max()
            if Al>Bl:
                Bl = Al
    
    if Brightness:
        return Aq32, Bl
    return Aq32