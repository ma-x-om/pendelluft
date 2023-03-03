# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 12:31:44 2019

@author: vgmarques and John Sims
@botched by: Mateus Ramirez
"""
import numpy as np
import potential_calculation2 as pt2
#from scipy.ndimage.filters import  maximum_filter,minimum_filter
import time

def runDHHD(mv_u,mv_v):

    nr,nc = mv_u.shape
    
#    im_xy = np.zeros(mv_u.shape);
#    
    mv_x,mv_y = np.meshgrid(np.arange(0,nc),np.arange(0,nr)) #% meshgrid for velocity vectors
#
#    t_i = time.time()   
    mesh = pt2.mesh_points(mv_x,mv_y);# Achtung!!!
#    print('mesh_points: %0.3f s'%(time.time()-t_i))
    
#    t_i = time.time()   
    triangle = pt2.triangle_definition(mv_x,mv_y);# achtung
#    print('triangle_definition: %0.3f s'%(time.time()-t_i))
    
    #    generate the Basis Gradient table
#    t_i = time.time()   
    gradphi = pt2.phi_gradient(triangle,mesh)# achtung
#    print('phi_gradient: %0.3f s'%(time.time()-t_i))
    
    # Rewriten to avoid loops
#    t_i = time.time()   
    
    tri_and_points = [pt2.get_neighbor(i,triangle) for i in range(0,mesh.shape[0])]
    
    N_tri = [tri_and_points[i][0] for i in range(0,mesh.shape[0])]
    N_point = [ tri_and_points[i][0] for i in range(0,mesh.shape[0])]
    
#    print('get_neighbors: %0.3f s'%(time.time()-t_i))

    node_info = {'mesh':mesh,'triangle':triangle,
                 'gradphi':gradphi,'N_tri':N_tri,'N_point':N_point}
#    %  perform DHHD, creating curl-free and divergence-free potential functions.
    
    
    
    t_i = time.time()   

    potE,potW = pt2.potential_calculation2(mv_x,mv_y,mv_u,mv_v,node_info);
    print('Exited potential_calculation2: %0.3f s'%(time.time()-t_i))
    '''
    veldcalc = {}
    velrcalc  ={}
    
    
    lm = maximum_filter(potE,8)
    potEmax = (potE == lm) 
    
    lm = minimum_filter(potE,8)
    potEmin = (potE == lm) 
    
    lm = maximum_filter(potW,8)
    potWmax = (potW == lm) 
    
    lm = minimum_filter(potW,8)
    potWmin= (potW == lm) 
    
    
    Emax_r,Emax_c= np.where(potEmax==1);
    Emin_r,Emin_c= np.where(potEmin==1);
    Wmax_r,Wmax_c= np.where(potWmax==1);
    Wmin_r,Wmin_c= np.where(potWmin==1);

    veldcalc['u'],veldcalc['v'] = np.gradient(potE)
    velrcalc['u'],velrcalc['v'] = np.gradient(potW)
    
    
#    velH.u = mv_u - veldcalc.u - velrcalc.u;
#    velH.v = mv_v - veldcalc.v - velrcalc.v;
#    
    '''
    return potE, potW
