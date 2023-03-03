# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 11:30:58 2019

@author: vgmarques and John Sims
Translating John's codes to pyzhon
"""

import numpy as np
import time
import os 
#import personal_functions as pf
#import scipy.signal as sig


def potential_calculation2(mv_x, mv_y, mv_u, mv_v,node_info):
    ''' 
    potE: the curlfree potential function E
    potW: the divergence-free potential function W
    mv_x: horizontal coords of the grids
    mv_y: vertical coords of the grids
    mv_u: horizontal components of the input motion field
    mv_v: vertical components of the input motion field.    
    '''
    M,N = mv_x.shape; 
    S1= construct_S1(mv_x, mv_y,node_info); #construct element matrix S1 # achtung
    
    [sa,sb] = np.shape(S1);

    potE = np.zeros(shape = (M,N),dtype = np.double)
    potW = np.zeros(shape = (M,N),dtype = np.double)
        
    '''
    vectors B and C are described on p 38 39 of thesis.  
    B is an L x 1  vector determined by the given velocity field and gradiens
    of the basis functions grad(phi_i)
    
    C is an an Lx1 vector, which can be computed using equ 3-23
    c_i = sum (over triangles) of a cross product where k component is
    non-zero.
    
    
    %dbstop construct_BC
    '''
    
    bcmat = np.arange(0,N*M)

    bcmat = np.reshape(bcmat,(M,N),order = 'F'); #I have matrix M by N 
    
    bclist = np.append(np.append(bcmat[:,0],bcmat[:,-1],axis = 0),
                       np.append(bcmat[0,:].T,bcmat[-1,:].T,axis = 0),
                       axis = 0)#list of border elements !
    
    bclist  = np.unique(bclist); #take unique values and sort.
    Ident= np.identity(M*N)
    
    
    # replace bclist columns in S1 with columns from in Ident
    S1[:,bclist]=Ident[:,bclist];
    
    IS1 = np.linalg.inv(S1);
    
    
    B,C = construct_BC(mv_x, mv_y, mv_u, mv_v,node_info);  # construct the two vector B and C # achtung
    
    # apply boundary conds to B and C; set values to zero.
    Br=B
    Br[bclist]=0;
    Cr=C
    Cr[bclist]=0;
    
    
    Er=np.matmul(IS1,Br);    #solve for the (L-1)x1 vector Er
    Wr=np.matmul(IS1,Cr);       #solve for the (L-1)x1 vector Wr
    
    Er[bclist]=0
    Wr[bclist]=0
    E_pot=Er; #[0;Er];       #reconstruct Lx1 vector E
    W_pot=Wr; #[0;Wr];       #reconstruct Lx1 vector W
    potE=np.reshape(E_pot,(M,N),order = 'F')    #re-organize MxN potential surface E
    potW =np.reshape(W_pot,(M,N),order = 'F')   #re-organize MxN potential surface W

    return potE, potW

'''

'''
 
def construct_S1(mv_x,mv_y,node_info):
    mesh = node_info['mesh']
    triangle = node_info['triangle']
    gradphi = node_info['gradphi']
    N_tri = node_info['N_tri']
#    N_point = node_info['N_point']

    #S1 is the element matrix

    M = mesh.shape[0] #size(mesh,1)
    S1= np.zeros((M,M),dtype = np.double)
    
    #calculate the element matrix S1
    
    #  There are M rows and N columns
#    t_i = time.time()   
    for i in range(0,M):
#        search for the neighbouring triangles and neighboring nodes of a
#        reference node
        #neighbor_triangle, neighbor_point = get_neighbor(i,triangle)# achtung
        for t in N_tri[i]:
            tri=triangle[t,:]
            order_in_triangle= np.where(tri==i)[0][0]
            del_phi01 = gradphi[t,order_in_triangle,0] # achtung
            del_phi02 = gradphi[t,order_in_triangle,1]# achtung
            for j in range(0,3):
                mm=tri[j];
                del_phi11 = gradphi[t,j,0]# achtung
                del_phi12 = gradphi[t,j,1]# achtung
                S1[i,mm] += del_phi01*del_phi11+del_phi02*del_phi12
#    print('S1 loop: %0.3f s'%(time.time()-t_i))

    return S1
'''

'''

def construct_BC(mv_x,mv_y,mv_u,mv_v,node_info):

#     mv_x: horizontal coords of grids
#     mv_y vertical coords of grids
#     mv_u: horizontal component of the input motion field
#     mv_v: vertical component of the input motion field
#     B: vector B
#     C: vector C
    

    
    mesh = node_info['mesh']
    triangle = node_info['triangle']
    gradphi = node_info['gradphi']
    N_tri = node_info['N_tri']
#    N_point = node_info['N_point']
#     generate the Average Vector table
#    %UV is the average vector table
    
    delta_x=mv_x[0,1] - mv_x[0,0]
    delta_y=mv_y[1,0] - mv_y[0,0]
    
    UV = triangle_uv(triangle,mesh,mv_x[0,0], mv_y[0,0], delta_x, delta_y, mv_u, mv_v) # achtung
    

    
    M = mesh.shape[0]; #size(mesh,1) #achtung
    B = np.zeros(shape = (M,1),dtype = np.double)
    C = np.zeros(shape = (M,1),dtype = np.double)
    
#    calculate vectors B and C
#    t_i = time.time()      
    for i in range(0,M):
        for t in N_tri[i]:
            tri = triangle[t,:];
            order_in_triangle = np.where(tri==i)[0]
            del_phi01= gradphi[t, order_in_triangle,0] # achtung
            del_phi02= gradphi[t, order_in_triangle,1] # achtung
            B[i] = B[i] + del_phi01*UV[t,0] + del_phi02*UV[t,1];
            C[i] = C[i] - del_phi01*UV[t,1] + del_phi02*UV[t,0];
#    print('BC loop: %0.3f s'%(time.time()-t_i))

    return B,C

'''

'''
def mesh_points(mv_x,mv_y):
#    generate Node Table
#    mv_x: horizontal coords of the grids
#     mv_y: vertical coords of the grids
#    xy: the coords of all nodes
    
    M,N = mv_x.shape;
#    MN = np.prod(sz);
#    xx = np.reshape(mv_x,(MN,1));
#    
#    yy = np.reshape(mv_y,(MN,1));
    xy = np.zeros((M*N,2));
    k = 0
    for i in range(0,M):
        for j in range(0,N):
            xy[k,0] = mv_y[i,j]
            xy[k,1] = mv_x[i,j]
            k += 1
            
    xy = np.asarray(sorted(xy,key=lambda x: (x[1],x[0]))) # Make sorting like matlab TODO
    
    return xy
'''

'''
def triangle_definition(mv_x,mv_y): 
#    %this function originally caused an error.  I adapted it to create value xx
#    
#    generate Grid Table, which defines the grid topology
#     mv_x: horizontal coords of the grids
#     mv_y: vertical coords of the grids
#     triangles: Nodes of all triangular meshes, i.e. the grid topology
    
    [M,N] = mv_x.shape;
    
    period = (M-1);
    total_tri= 2*period*(N-1);
    triangles = np.zeros((total_tri,3))
    
    for i in range(0,N-1):
        for j in range(0,period):
            k= (i)*period + j
            start = np.floor((k)/period)*M + np.mod(k,period)
            triangles[2*k,:] = [start,start+M+1,start+M]
            triangles[2*k+1,:] = [start,start+1,start+M+1]
    
    return np.uint(triangles)
'''

'''
def triangle_uv(triangle,mesh,start_x, start_y,delta_x,delta_y,uu,vv):
    
#     generate the Average Vector Table
#     triangle: the Grid Table
#     mesh: the Node Table
#    start_x: start coord in horizontal dir
#    start_y: start coord in vertical dir
#    delta_x: step in horiz dir
#    delta_y: step in vertical dir
#    uu: horizontal component of the input motion field
#    vv: vertical component of the input motion field
#    UV: the resulted Average Vector table
    
    M = triangle.shape[0] #size(triangle,1);
    UV = np.zeros((M,2))
     
    # Remove loop to speed up computation
#    t_i = time.time()
    # Locs
#    loc0 = mesh[triangle[:,0],:]
#    loc1 = mesh[triangle[:,1],:]
#    loc2 = mesh[triangle[:,2],:]
#    
#    #Inds
#    ind_y0 = np.asarray((loc0[:,0]-start_y)/delta_y,dtype = np.int)
#    ind_x0 = np.asarray((loc0[:,1]-start_x)/delta_x,dtype = np.int)
#    ind_y1 = np.asarray((loc1[:,0]-start_y)/delta_y,dtype = np.int)
#    ind_x1 = np.asarray((loc1[:,1]-start_x)/delta_x,dtype = np.int)    
#    ind_y2 = np.asarray((loc2[:,0]-start_y)/delta_y,dtype = np.int)
#    ind_x2 = np.asarray((loc2[:,1]-start_x)/delta_x,dtype = np.int)
#    
#    
#    UV[:,0] = (uu[ind_y0,ind_x0]+uu[ind_y1,ind_x1]+uu[ind_y2,ind_x2])/3
#    UV[:,1] = (vv[ind_y0,ind_x0]+vv[ind_y0,ind_x0]+vv[ind_y0,ind_x0])/3
    

    
    for i in range(0,M):
        loc1 = mesh[triangle[i,0],:]
        loc2 = mesh[triangle[i,1],:]
        loc3 = mesh[triangle[i,2],:]
        ind_y1 = int(np.floor((loc1[0]-start_y)/delta_y))  # Tirei os +1 desses índices
        ind_x1 = int(np.floor((loc1[1]-start_x)/delta_x))
        ind_y2 = int(np.floor((loc2[0]-start_y)/delta_y))
        ind_x2 = int(np.floor((loc2[1]-start_x)/delta_x))
        ind_y3 = int(np.floor((loc3[0]-start_y)/delta_y))
        ind_x3 = int(np.floor((loc3[1]-start_x)/delta_x))
        UV[i,0] = (uu[ind_y1,ind_x1]+uu[ind_y2,ind_x2]+uu[ind_y3,ind_x3])/3
        UV[i,1] = (vv[ind_y1,ind_x1]+vv[ind_y2,ind_x2]+vv[ind_y3,ind_x3])/3
#    print('triangle_uv: %0.3f s'%(time.time()-t_i))

    return UV

'''

'''
def phi_gradient(triangle,mesh):
#    generate basis gradient table
#    triangle: the grid table
#    mesh: the node table
#    gradphi: the resulted basis gradient table
    
    M = triangle.shape[0]#size(triangle,1);
    gradphi = np.zeros((M,3,2));
    area2= get_area2(triangle[0,:],mesh); # calc the 2*triangle_area
    AB = np.zeros((3,2));
    for i in range(0,M):
        y0=mesh[triangle[i,0],0]
        x0=mesh[triangle[i,0],1]
        y1=mesh[triangle[i,1],0]
        x1=mesh[triangle[i,1],1]
        y2=mesh[triangle[i,2],0]
        x2=mesh[triangle[i,2],1]
        
        AB[0,0]=y1-y2
        AB[0,1]=x2-x1
        AB[1,0]=y2-y0
        AB[1,1]=x0-x2
        AB[2,0]=y0-y1
        AB[2,1]=x1-x0
        
        for j in range(0,3):
            for k in range(0,2):
                gradphi[i,j,k]=-AB[j,k]
    
        gradphi=gradphi/area2;
    return gradphi 
'''

'''

def  get_neighbor(node,triangle):
#     get the neighboring triangles and neighboring nodes of a node
#    node: the reference node
#     triangle: the grid table
#     neighbor_triangle: the neighbor triangle of the reference node
#    neighbor_point: the neighbor points of the reference node
    
#    M = triangle.shape[0]#size(triangle,1);
    
    neighbor_triangle=[]
    
#    t_i = time.time()
    # Loop changed to reduce comp. time
    check = np.sum(triangle==node,axis = 1)
    for i,C in enumerate(check):
        if C:
            neighbor_triangle=np.append(neighbor_triangle,i)
#    print('neighbor_triangle: %0.3f'%(time.time()-t_i))  

    neighbor_point=[]
#    N = len(neighbor_triangle)
#    t_i = time.time()
    for n_tri in neighbor_triangle:
        temp = triangle[int(n_tri),:]
        for j in range(0,3):
            if np.sum(neighbor_point==temp[j])>0:
                continue
            else:
                neighbor_point = np.append(neighbor_point, temp[j])

    neighbor_point= np.sort(neighbor_point) #   sort in ascending order
#    print('neighbor_point: %0.3f'%(time.time()-t_i)) 
    
    return np.uint(neighbor_triangle),np.uint(neighbor_point)

'''

'''
def get_area2(triangle,mesh):
    
    
    # calc the 2*traingle_area of a triangle
    # triangle: the grid table
    # mesh: the node table
    
    
    # area2: the resulted 2*triangle_area
    x0 = mesh[triangle[0],0];
    y0 = mesh[triangle[0],1];
    x1 = mesh[triangle[1],0];
    y1 = mesh[triangle[1],1];
    x2 = mesh[triangle[2],0];
    y2 = mesh[triangle[2],1];
    
    area2= np.abs((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0));
    return area2
