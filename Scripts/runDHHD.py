# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 12:31:44 2019

@author: vgmarques and John Sims
@botched by: Mateus Ramirez
"""
import numpy as np
import potential_calculation2 as pt2
import time

def runDHHD(mv_u,mv_v, show_time = False):   

    nr,nc = mv_u.shape
    
    mv_x,mv_y = np.meshgrid(np.arange(0,nc),np.arange(0,nr)) #% meshgrid for velocity vectors
 
    mesh = pt2.mesh_points(mv_x,mv_y);# Achtung!!!

    triangle = pt2.triangle_definition(mv_x,mv_y);# achtung

    gradphi = pt2.phi_gradient(triangle,mesh)# achtung

    # Rewriten to avoid loops
    tri_and_points = [pt2.get_neighbor(i,triangle) for i in range(0,mesh.shape[0])]
    
    N_tri = [tri_and_points[i][0] for i in range(0,mesh.shape[0])]
    N_point = [ tri_and_points[i][0] for i in range(0,mesh.shape[0])]

    node_info = {'mesh':mesh,'triangle':triangle,
                 'gradphi':gradphi,'N_tri':N_tri,'N_point':N_point}
#    %  perform DHHD, creating curl-free and divergence-free potential functions.
    
    if show_time: t_i = time.time()   

    potE,potW = pt2.potential_calculation2(mv_x,mv_y,mv_u,mv_v,node_info);
    if show_time: print('Exited potential_calculation2: %0.3f s'%(time.time()-t_i))

    return potE, potW
