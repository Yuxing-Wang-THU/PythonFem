from matplotlib import pyplot as plt
from matplotlib.patches import Wedge
from matplotlib import cm 
import numpy as np
import imageio
import shutil
import os 

model_list = {0:'Linear', 1:'STVK', 2:'Co-rotated', 3: 'Neohookean'}

def rect_collision(point, rect_pos):
    if point[0] < rect_pos[0][0]:
        return False
    if point[1] < rect_pos[0][1]:
        return False
    if point[0] > rect_pos[0][0] + rect_pos[1][0]:
        return False
    if point[1] > rect_pos[0][1] +rect_pos[1][1]:
        return False
    return True

# Mesh
def get_mesh(args):
    # Element的索引，每一个element包含三个点的索引
    element_idx = np.zeros((args.element_num,3),dtype = int)

    # 为每个节点构建对应的index
    cnt = 0
    for j in range(args.N_y-1):
        for i in range(args.N_x-1):
            # 第1个三角形，下三角
            idx = j * args.N_x + i
            element_idx[cnt,0] = idx
            element_idx[cnt,1] = idx + 1
            element_idx[cnt,2] = idx + args.N_x
            # 第2个三角形，上三角
            cnt += 1
            element_idx[cnt,0] = idx + 1
            element_idx[cnt,1] = idx + args.N_x + 1
            element_idx[cnt,2] = idx + args.N_x
            cnt += 1
    return element_idx

# Plot
def plot_voxel(node_pos, element_num, element_idx, empty_voxel_index):
    for ie in range(element_num):
        if ie in empty_voxel_index:
            continue
        element_x = node_pos[element_idx[ie,:],0]
        element_y = node_pos[element_idx[ie,:],1]
        dx = np.array([[element_x[1] - element_x[0],element_x[2] - element_x[0]],
                       [element_y[1] - element_y[0],element_y[2] - element_y[0]]])
        current_area = np.linalg.det(dx)*0.5
        plt.fill(element_x,element_y, c=cm.Blues(0.5/current_area))

# Deformation
def Afine_2D(points, theta, tx, ty, sx, sy):
    rad = np.radians(theta)
    # Rotation
    rotation_matrix = np.array([[np.cos(rad), -np.sin(rad), 0], [np.sin(rad), np.cos(rad), 0], [0, 0, 1]])
    # Translation
    translation_matrix = np.array([[1, 0, tx], [0, 1, ty], [0, 0, 1]])
    # Stretching
    scale_matrix = np.array([[sx, 0, 0], [0, sy, 0], [0, 0, 1]])
    transform_matrix = np.dot(np.dot(translation_matrix, rotation_matrix), scale_matrix)
    rts = np.hstack((points, np.ones((points.shape[0], 1)))) 
    rts = np.dot(rts, transform_matrix.T)[:,:-1]
    return rts

# Calculate Area
def get_area(points):
    return 0.5  *  (points[0,0] * (points[1,1] - points[2,1])
                  + points[1,0] * (points[2,1] - points[0,1]) 
                  + points[2,0] * (points[0,1] - points[1,1]))
    
def vec_to_array(vec, array):
    for i in range(vec.shape[0]):
        array[:,i] = vec[i].flatten("F")
    return array 

# Linear elasity
def linear_elasity(F,mu,lam):
    # Calculate strain tensor of Linear Elasticticy
    # Small strain tensor  e = 1/2(F+F^T) -I
    strain = (F + F.T) * 0.5 - np.identity(2)
    # Calculate Energy density 
    doubleInner = strain[0,0]*strain[0,0] + strain[1,0]*strain[1,0] + strain[0,1]*strain[0,1] + strain[1,1]*strain[1,1]
    # Ψ(F) = µe:e + λ/2 x trace^2(strain)
    energy = doubleInner * mu + lam * 0.5 * np.trace(strain) ** 2
    # Get PK1 stress
    # piola = µ(F + F^T - 2I) + λtrace(F-I)I
    piola = mu * (F + F.T - 2 * np.identity(2)) + lam * (F[0,0] - 1 + F[1,1] - 1) * np.identity(2)
    return energy, piola, strain

# St. Venant-Kirchhoff elasity
def stvk(F,mu,lam):
    # Green strain tensor  e = 1/2(F^TF-I) 
    strain = (np.dot(F.T,F) - np.identity(2)) * 0.5
    # Calculate Energy density function for updating position
    doubleInner = strain[0,0]*strain[0,0] + strain[1,0]*strain[1,0] + strain[0,1]*strain[0,1] + strain[1,1]*strain[1,1]
    # Ψ(F) = µe:e + λ/2 x trace^2(strain)
    energy = doubleInner * mu + lam * 0.5 * (strain[0,0] + strain[1,1]) ** 2
    # Get PK1 stress
    # piola = F[2µE + λtrace(E)I]
    piola = np.dot(F, 2 * mu * strain + lam * (strain[0,0] + strain[1,1]) * np.identity(2))
    return energy, piola, strain

# Corotated linear elasticity
def corotated(F,mu,lam):
    # SVD decomposition of F
    U, sigma, Vt = np.linalg.svd(F)

    # Calculate Energy density function
    # Get polar decomposition
    R = np.dot(U,Vt)
    S = np.dot(np.linalg.inv(R),F)
    e_c = S - np.identity(2)

    # Ψ(F) = µe_c:e_c + λ/2 x trace^2(e_c)
    doubleInner = e_c[0,0]*e_c[0,0] + e_c[1,0]*e_c[1,0] + e_c[0,1]*e_c[0,1] + e_c[1,1]*e_c[1,1]
    energy = mu * doubleInner + 0.5 * lam * np.trace(e_c)**2
    
    # Get PK1 stress
    # piola = = 2µ(F − R) + λtr(RTF − I)R
    RTF_I = np.dot(R.T,F) - np.identity(2)
    piola = 2 * mu * (F - R) + lam * np.trace(RTF_I) * R
    return energy, piola, e_c

# Neohookean elasticity
def neohookean(F,mu,lam):
    # Calculate strain tensor of Linear Elasticticy
    # Green strain tensor  e = 1/2(F^TF-I)
    strain = (np.dot(F.T,F)- np.identity(2)) * 0.5
   
    # Get Log J
    logJ = np.log(max(np.linalg.det(F),0.01))
    # logJ = np.log(np.linalg.det(F))
    # log_I3 = np.log(max(np.linalg.det(F)**2,0.01))

    I_1 = np.trace(np.dot(F.T,F))
    # Calculate Energy density function for updating position

    # Ψ(I1, J) = µ/2(I1 − 2) − µlog(J) + λ/2log^2(J)
    energy = 0.5 * mu * (I_1 - 2) - mu * logJ + 0.5 * lam * logJ**2
    
    # Ψ(I1, J) = µ/2(I1 − 2) − µlog(J) + λ/2log^2(J)
    # energy = 0.5 * mu * (I_1 - log_I3 - 2) + lam/8 * log_I3**2
    
    # Get PK1 stress
    # piola = µ(F − µF−T) + λ log(J)F−T
    F_inv = np.linalg.inv(F)
    F_inv_T = F_inv.T
    
    piola = mu * (F - F_inv_T) + lam * logJ * F_inv_T
    # piola = mu * (F - F_inv_T) + lam * 0.5 * log_I3 * F_inv_T
    return energy, piola, strain

# Implict Euler method for STVK
def implicit_euler_stvk(args, B_m, F, strain, nodal_force,nodal_vel, multiple=False):
    # Implicit Euler
    # ∂D_s/∂x，其中D_s是一个2x2的矩阵，x是一个3x2的矩阵（仨点俩维度）, 因此计算出来是一个2x2x3x2的四阶张量
    # 使用向量化进行计算的话可以得到每一列,向量化后D_s是一个长度为4的向量，x是一个长度为6的向量，求Jacobin，分子主序，最后矩阵是4x6

    dDsdx = np.zeros((6,2,2))
    dDsdx[0,:,:] = np.array([[-1,-1],[0,0]])
    dDsdx[1,:,:] = np.array([[0,0],[-1,-1]])

    dDsdx[2,:,:] = np.array([[1,0],[0,0]])
    dDsdx[3,:,:] = np.array([[0,0],[1,0]])

    dDsdx[4,:,:] = np.array([[0,1],[0,0]])
    dDsdx[5,:,:] = np.array([[0,0],[0,1]])
   
    # delta deformation gradient
    dF, dE, dP, dH  = np.zeros_like(dDsdx), np.zeros_like(dDsdx),np.zeros_like(dDsdx),np.zeros_like(dDsdx)

    # 𝛿F, 𝛿P, 𝛿H
    for i in range(6):
        dF[i,:,:] = np.dot(dDsdx[i,:,:],B_m) 
        d_F = dF[i,:,:]
        dE[i,:,:] = (np.dot(d_F.T,F) + np.dot(F.T,d_F))*0.5
        d_E = dE[i,:,:]
        dP[i,:,:] = np.dot(d_F,2 * args.mu * strain + args.lam * (strain[0,0] + strain[1,1]) * np.identity(2))
        dP[i,:,:] += np.dot(F,2 * args.mu * d_E + args.lam * (d_E[0,0] + d_E[1,1]) * np.identity(2))
        dH[i,:,:] = -args.init_area * np.dot(dP[i,:,:],B_m.T)
    
    if multiple:
        return dH
    # K
    array_H = vec_to_array(vec=dH,array=np.zeros((4,6)))

    K = np.zeros((6,6))
    # 3 个顶点
    for n in range(3):
        # 2 个维度
        for d in range(2):
            # 第 idx 列，每列3 x 2 个元素
            idx = n * 2 + d
            # 先填写第一第二个顶点，第零个顶点之后填
            K[2,idx] += dH[idx,0,0]
            K[3,idx] += dH[idx,1,0]
            K[4,idx] += dH[idx,0,1]
            K[5,idx] += dH[idx,1,1]
            
            K[0,idx] += - dH[idx,0,0] - dH[idx,0,1]
            K[1,idx] += - dH[idx,1,0] - dH[idx,1,1]
    # A
    A = args.mass * np.identity(6) -  K  * args.dt * args.dt
    # b
    b = np.zeros((6))
    for n in range(3):
        for d in range(2):
            b[n*2+d] = args.mass * nodal_vel[n,d] + args.dt * nodal_force[n,d]
    # X       
    X = np.dot(np.linalg.inv(A), b)
    for n in range(3):
        if args.fix and n == 0:
            continue
        for d in range(2):
            nodal_vel[n,d] = X[n*2+d]

    return nodal_vel      

# Implict Euler method for Neohookean
def implicit_euler_neohookean(args, B_m, F, nodal_force,nodal_vel,multiple=False):
    # Implicit Euler
    # ∂D_s/∂x，其中D_s是一个2x2的矩阵，x是一个3x2的矩阵（仨点俩维度）, 因此计算出来是一个2x2x3x2的四阶张量
    # 使用向量化进行计算的话可以得到每一列,向量化后D_s是一个长度为4的向量，x是一个长度为6的向量，求Jacobin，分子主序，最后矩阵是4x6

    dDsdx = np.zeros((6,2,2))
    dDsdx[0,:,:] = np.array([[-1,-1],[0,0]])
    dDsdx[1,:,:] = np.array([[0,0],[-1,-1]])

    dDsdx[2,:,:] = np.array([[1,0],[0,0]])
    dDsdx[3,:,:] = np.array([[0,0],[1,0]])

    dDsdx[4,:,:] = np.array([[0,1],[0,0]])
    dDsdx[5,:,:] = np.array([[0,0],[0,1]])
   
    # delta deformation gradient
    dF, dP, dH  = np.zeros_like(dDsdx), np.zeros_like(dDsdx),np.zeros_like(dDsdx)

    # 𝛿F, 𝛿P, 𝛿H
    F_inv = np.linalg.inv(F)
    logJ = np.log(max(np.linalg.det(F),0.01))
    for i in range(6):
        dF[i,:,:] = np.dot(dDsdx[i,:,:],B_m) 
        d_F = dF[i,:,:]
        FTDFTF_T = np.dot(np.dot(F_inv.T,d_F.T),F_inv.T)
        dP[i,:,:] = args.mu * d_F + (args.mu-args.lam*logJ) * FTDFTF_T + args.lam * np.trace(np.dot(F_inv,d_F)) * F_inv.T
        dH[i,:,:] = -args.init_area * np.dot(dP[i,:,:],B_m.T)
    
    if multiple:
        return dH
    # K
    array_H = vec_to_array(vec=dH,array=np.zeros((4,6)))

    K = np.zeros((6,6))
    # 3 个顶点
    for n in range(3):
        # 2 个维度
        for d in range(2):
            # 第 idx 列，每列3 x 2 个元素
            idx = n * 2 + d
            # 先填写第一第二个顶点，第零个顶点之后填
            K[2,idx] += dH[idx,0,0]
            K[3,idx] += dH[idx,1,0]
            K[4,idx] += dH[idx,0,1]
            K[5,idx] += dH[idx,1,1]
            
            K[0,idx] += - dH[idx,0,0] - dH[idx,0,1]
            K[1,idx] += - dH[idx,1,0] - dH[idx,1,1]
    # A
    A = args.mass * np.identity(6) -  K  * args.dt * args.dt
    # b
    b = np.zeros((6))
    for n in range(3):
        for d in range(2):
            b[n*2+d] = args.mass * nodal_vel[n,d] + args.dt * nodal_force[n,d]
    # X       
    X = np.dot(np.linalg.inv(A), b)
    for n in range(3):
        if args.fix and n == 0:
            continue
        for d in range(2):
            nodal_vel[n,d] = X[n*2+d]

    return nodal_vel      

# Make Ax=b for Implict Euler method
def solver(args,A,K,b,node_vel,node_force):
    for i in range(args.krow):
        for j in range(args.krow):
            if i == j:
                A[i,j] = args.mass - K[i,j]*args.dt*args.dt
            else:
                A[i,j] = - K[i,j]*args.dt*args.dt
    for i in range(args.node_num):
        b[i*2 + 0] = args.mass * node_vel[i,0] + args.dt * node_force[i,0]
        b[i*2 + 1] = args.mass * node_vel[i,1] + args.dt * node_force[i,1]
        
    X = np.dot(np.linalg.inv(A),b)
    return X

# Simulation for a single tetrahedron
def run_fem_simulation_single(args):
    # Dir for saving imgs 
    try:
        shutil.rmtree("./images")
    except:
        pass
    if not os.path.exists('images'):
        os.makedirs('images')

    # Figure
    plt.ion()
    plt.figure(figsize=(9, 8))
    plt.xlim(0, 7)
    plt.ylim(0, 7)

    # Two Lame
    mu = args.Y / ( 2 * (1 + args.V))
    lam = args.Y * args.V / (1 + args.V) / (1 - 2 * args.V)
    args.mu = mu
    args.lam = lam

    # Gravity
    gravity = np.array([0, -9.8]) * args.mass

    # Velocity
    vel = np.zeros_like(args.ref_node_pos)

    # Init Area
    init_area = get_area(args.ref_node_pos)
    args.init_area = init_area

    # Position after deformation 
    x = Afine_2D(args.ref_node_pos, args.theta, args.tx, args.ty, args.sx, args.sy)
   
    # Plot the triangular after deformation
    plt.fill(x[:,0], x[:,1], 'b')
    plt.text(0.01, 6.8, f'Init Area: {init_area}', fontsize=15, color='red')
    plt.text(1.8, 6.8, f'mu: {round(mu,2)}, lamda: {round(lam,2)}', fontsize=15, color='blue')
    plt.text(5.0, 6.8, f'Iteration: {0}', fontsize=15, color='green')
    plt.text(4.0, 6.5, f'Constitutive model: {model_list[args.model]}', fontsize=15, color='black')
    plt.text(0.01, 6.2, f'Implict time integration: {args.implict}', fontsize=15, color='black')
    plt.text(x[0][0], x[0][1], f'{round(x[0][0],2), round(x[0][1],2)}', fontsize=15, color='black')
    plt.text(x[1][0], x[1][1], f'{round(x[1][0],2), round(x[1][1],2)}', fontsize=15, color='black')
    plt.text(x[2][0], x[2][1], f'{round(x[2][0],2), round(x[2][1],2)}', fontsize=15, color='black')
    plt.draw()
    plt.savefig('images/plot_{}.png'.format(0))
    plt.pause(1)

    '''Step 1. Pre-calculate D_m (2 x 2)'''
    # [[X_1-X_0, X_2-X_0]
    #  [Y_1-Y_0, Y_2-Y_0]]
    D_m = np.array([[args.ref_node_pos[1,0] - args.ref_node_pos[0,0],
                     args.ref_node_pos[2,0] - args.ref_node_pos[0,0]],
                    [args.ref_node_pos[1,1] - args.ref_node_pos[0,1],
                     args.ref_node_pos[2,1] - args.ref_node_pos[0,1]]])
    assert init_area == 0.5 * np.linalg.det(D_m)

    # B_m = D_m^-1, (2 x 2)
    B_m = np.linalg.inv(D_m)
    
    # Iteration
    iter = 0
    while(iter < args.max_iter):
        plt.clf()
        iter += 1
        for substep in range(10):
            nodal_force = np.zeros((3,2))

            '''Step 2. Calculate D_s (2 x 2)'''
            D_s = np.array([[x[1,0] - x[0,0],
                             x[2,0] - x[0,0]],
                            [x[1,1] - x[0,1],
                             x[2,1] - x[0,1]]])
            
            '''Step 3. Calculate deformation gradient'''
            # F = DsDm^-1
            F = np.dot(D_s,B_m)
            
            '''Step 4. Calculate PK1 stress and Energy''' 
            if args.model == 0:
                energy, piola, strain = linear_elasity(F,mu,lam)
            elif args.model == 1:
                energy, piola, strain = stvk(F,mu,lam)
            elif args.model == 2:
                energy, piola, strain = corotated(F,mu,lam)
            elif args.model == 3:
                energy, piola, strain = neohookean(F,mu,lam)
            else:
                print("Model is not implemented")
                raise    

            '''Step 5. Calculate Force''' 
            # H = - WP(B_m)^T
            H = - init_area * np.dot(piola, B_m.T)
            if args.gravity:
                # Force 1, the first column of H
                nodal_force[1,:] = np.array([H[0,0],H[1,0]]) + gravity
                # Force 2, the second column of H
                nodal_force[2,:] = np.array([H[0,1],H[1,1]]) + gravity 
                # Force 0, balanced force
                nodal_force[0,:] = - nodal_force[1,:] - nodal_force[2,:] + gravity
            else:
                # Force 1, the first column of H
                nodal_force[1,:] = np.array([H[0,0],H[1,0]])
                # Force 2, the second column of H
                nodal_force[2,:] = np.array([H[0,1],H[1,1]])
                # Force 0, balanced force
                nodal_force[0,:] = - nodal_force[1,:] - nodal_force[2,:]

            '''Step 6. Update positions of three points''' 
            if not args.implict:
                # Force-based update
                if not args.fix:
                    vel[0,:] += args.dt * nodal_force[0,:] / args.mass 
                vel[1,:] += args.dt * nodal_force[1,:] / args.mass 
                vel[2,:] += args.dt * nodal_force[2,:] / args.mass 

                x[0,:] += args.dt * vel[0,:]
                x[1,:] += args.dt * vel[1,:]
                x[2,:] += args.dt * vel[2,:]
                
                # Damping
                vel[0,:] *= np.exp(-args.dt * 5.0)
                vel[1,:] *= np.exp(-args.dt * 5.0)
                vel[2,:] *= np.exp(-args.dt * 5.0)
            else:
                if args.model == 1:
                    vel = implicit_euler_stvk(args, B_m, F, strain, nodal_force,vel)
                    for n in range(3):
                        for d in range(2):
                            x[n,d] += vel[n,d]*args.dt
                            vel[n,d] *= np.exp(-args.dt * 5.0)
                elif args.model == 3:
                    vel = implicit_euler_neohookean(args, B_m, F, nodal_force,vel)
                    for n in range(3):
                        for d in range(2):
                            x[n,d] += vel[n,d]*args.dt
                            vel[n,d] *= np.exp(-args.dt * 5.0)
                else:
                    raise
            
            current_area = get_area(x)
        
            # Boundary Condition
            for i in range(3):
                for j in range(2):
                    if x[i][j] < 0:
                        x[i][j] = 0
                        vel[i][j] = -0.99*vel[i][j]  
                    if x[i][j] > 7:
                        x[i][j] = 7
                        vel[i][j] = -0.99*vel[i][j]
                    
        plt.text(0.01, 6.5, f'Energy: {round(energy,2)}', fontsize=15, color='purple')
        plt.fill(x[:,0], x[:,1], 'b')
        plt.xlim(0, 7)
        plt.ylim(0, 7)
        plt.text(0.01, 6.8, f'Init Area: {init_area}', fontsize=15, color='red')
        plt.text(1.8, 6.8, f'mu: {round(mu,2)}, lamda: {round(lam,2)}', fontsize=15, color='blue')
        plt.text(5.0, 6.8, f'Iteration: {iter}', fontsize=15, color='green')
        plt.text(4.0, 6.5, f'Constitutive model: {model_list[args.model]}', fontsize=15, color='black')
        plt.text(0.01, 6.2, f'Implict time integration: {args.implict}', fontsize=15, color='black')
        plt.text(x[0][0], x[0][1], f'{round(x[0][0],2), round(x[0][1],2)}', fontsize=15, color='black')
        plt.text(x[1][0], x[1][1], f'{round(x[1][0],2), round(x[1][1],2)}', fontsize=15, color='black')
        plt.text(x[2][0], x[2][1], f'{round(x[2][0],2), round(x[2][1],2)}', fontsize=15, color='black')
        plt.text(1.8, 6.5, f'Current Area: {round(current_area,2)}', fontsize=15, color='red')
        plt.draw()
        plt.savefig('images/plot_{}.png'.format(iter))
        plt.pause(0.0001)
    plt.ioff()

    # Make Gifs
    images = []
    for i in range(iter):
        filename = 'images/plot_{}.png'.format(i)
        images.append(imageio.imread(filename))
    imageio.mimsave(f'./gifs/single_2D_{model_list[args.model]}_implicit_{args.implict}.gif', images, duration=0.01)
    shutil.rmtree("./images")


    # Simulation

# Simulation for multiple tetrahedrons
def run_fem_simulation_multiple(args):
    # Make dirs
    try:
        shutil.rmtree("./images")
    except:
        pass
    if not os.path.exists('images'):
        os.makedirs('images')
    
    # Figure
    plt.ion()
    fig, ax = plt.subplots(figsize=(15, 9))
    plt.xlim(0, 15)
    plt.ylim(0, 9)
    
    mu = args.Y / ( 2 * (1 + args.V))
    lam = args.Y * args.V / (1 + args.V) / (1 - 2 * args.V)
    args.mu = mu
    args.lam = lam

    # Gravity
    Gravity = np.array([0, -9.8]) * args.mass

    # Nodes
    # Init position，left down coroner
    init_x, init_y = 0.0, 4.0
    # dx
    cube_len = 1.0
    # Number of nodes (x-axis)
    N_x = args.voxelx+1
    # Number of nodes (y-axis)
    N_y = args.voxely+1
    args.N_x = N_x
    args.N_y = N_y
    node_num = N_x * N_y
    args.node_num = node_num

    # Elements
    element_num = 2 * args.voxelx * args.voxely
    args.element_num = element_num

    # Mesh
    element_idx = get_mesh(args)

    # Position,velosity and force
    node_pos = np.zeros((node_num,2),dtype = float)
    node_vel = np.zeros((node_num,2),dtype = float)
    nodal_force = np.zeros((node_num,2),dtype = float)
    
    # Set position
    for j in range(N_y):
        for i in range(N_x):
            node_pos[j*N_x+i] = np.array([init_x + i * cube_len, init_y + j * cube_len])
    
    init_area = get_area(node_pos[element_idx[0]])
    args.init_area = init_area
    assert init_area == 0.5

    '''Step 1. Pre-calculate D_m: n x (2 x 2)'''
    # [[X_i-X_l, X_j-X_l]
    #  [Y_i-Y_l, Y_j-Y_l]]
    D_m = np.zeros((args.element_num,2,2))
    for ie in range(element_num):
        p0 = node_pos[element_idx[ie,0]]
        p1 = node_pos[element_idx[ie,1]]
        p2 = node_pos[element_idx[ie,2]]
        dX = np.array([[p1[0] - p0[0],p2[0] - p0[0]],
                       [p1[1] - p0[1],p2[1] - p0[1]]])
        D_m[ie] = np.linalg.inv(dX)
    
    node_pos = Afine_2D(node_pos, theta=0.0, tx=0.0, ty=0.0, sx=0.9, sy=1.0)
    
    # Add obstacles
    # ax.add_patch(args.rectangle)
    ax.add_patch(args.semicircle)
    plot_voxel(node_pos,element_num, element_idx,args.em)
    plt.draw()
    plt.tight_layout()
    plt.savefig('images/plot_{}.png'.format(0))

    # simulate
    time = 0
    if args.implict:
        krow = node_num * 2
        args.krow = krow
        K = np.zeros((krow,krow))
        A = np.zeros((krow,krow))
        b = np.zeros((krow))

    while(time < args.max_iter):
        time += 1
        plt.cla()
        # SubStep
        for i in range(50):
            nodal_force[:,:] = 0
            # Fem for each element
            '''Step 2. Calculate D_s (2 x 2)'''
            for ie in range(element_num):
                # empty voxel
                if ie in args.em:
                    continue
                p0 = node_pos[element_idx[ie,0]]
                p1 = node_pos[element_idx[ie,1]]
                p2 = node_pos[element_idx[ie,2]]
                D_s = np.array([[p1[0] - p0[0],p2[0] - p0[0]],
                                [p1[1] - p0[1],p2[1] - p0[1]]])
                
                '''Step 3. Calculate deformation gradient'''
                F = np.dot(D_s, D_m[ie])
                
                '''Step 4. Calculate PK1 stress and Energy''' 
                if args.model == 0:
                    energy, piola, strain = linear_elasity(F,mu,lam)
                elif args.model == 1:
                    energy, piola, strain = stvk(F,mu,lam)
                elif args.model == 2:
                    energy, piola, strain = corotated(F,mu,lam)
                elif args.model == 3:
                    energy, piola, strain = neohookean(F,mu,lam)
                else:
                    print("Model is not implemented")
                    raise    

                H = -init_area * np.dot(piola,D_m[ie].transpose())
                # Force 1, the first column of H
                f1 = np.array([H[0,0],H[1,0]]) 
                # Force 2, the second column of H
                f2 = np.array([H[0,1],H[1,1]]) 
                # Force 0, balanced force
                f0 = - f1 - f2 
               
                nodal_force[element_idx[ie,0],:] += f0
                nodal_force[element_idx[ie,1],:] += f1
                nodal_force[element_idx[ie,2],:] += f2
                
                if args.implict:
                    if args.model == 1:
                        dH = implicit_euler_stvk(args, D_m[ie], F, strain, nodal_force, node_vel, multiple=True)
                    elif args.model == 3:
                        dH = implicit_euler_neohookean(args, D_m[ie], F, nodal_force, node_vel,multiple=True)
                    for n in range(3):
                        nidx = element_idx[ie,n]
                        for d in range(2):
                            kidx = nidx * 2 + d
                            didx = n * 2 + d
                            idx = element_idx[ie,1] * 2
                            K[idx,kidx] += dH[didx,0,0]
                            idx = element_idx[ie,1] * 2 + 1
                            K[idx,kidx] += dH[didx,1,0]
                            idx = element_idx[ie,2] * 2
                            K[idx,kidx] += dH[didx,0,1]
                            idx = element_idx[ie,2] * 2 + 1
                            K[idx,kidx] += dH[didx,1,1]
                            idx = element_idx[ie,0] * 2
                            K[idx,kidx] += - dH[didx,0,0] - dH[didx,0,1]
                            idx = element_idx[ie,0] * 2 + 1
                            K[idx,kidx] += - dH[didx,1,0] - dH[didx,1,1]

            '''Step 6. Update positions of three points''' 
            # Update
            if not args.implict:
                for i in range(node_num):
                    acc = nodal_force[i]/args.mass + Gravity
                    node_vel[i] += args.dt*acc
                    node_vel[i] *= np.exp(-args.dt*1)
            else:
                # Add Gravity
                for i in range(node_num):
                    nodal_force[i] += Gravity
                # Implict Euler
                X = solver(args,A,K,b,node_vel,nodal_force)
                for i in range(node_num):
                    node_vel[i,0] = X[i * 2 + 0]
                    node_vel[i,1] = X[i * 2 + 1]
                    node_pos[i,:] += node_vel[i,:] * args.dt
                    node_vel[i] *= np.exp(-args.dt * 5.0)

            # bend
            if args.mode == "bend":
                for j in range(N_y):
                    ind = 0 + j * N_x
                    node_vel[ind,:] = np.array([0, 0])
                    node_pos[ind,:] = np.array([init_x, init_y + j * cube_len])  # rest pose attached to the wall

            # Boundary Condition
            for i in range(node_num):
                # Ball boundary condition
                # 圆心指向点的向量
                disp = node_pos[i] - args.ball_pos
                # 点距圆心的距离
                disp2 = np.sum(np.square(disp))
                # 如果小于半径，说明有碰撞
                if disp2 <= args.ball_radius**2:
                    # 内积，获得速度在该向量上的投影
                    NoV = node_vel[i].dot(disp)
                    # <0，说明在挤压这个圆
                    if NoV < 0:
                        # 把速度掰成相切的
                        node_vel[i] -= NoV * disp / disp2

                # Rectangle boundary condition
                # if rect_collision(node_pos[i], args.rectangle_pos):
                #     if args.rectangle_pos[0][0] < node_pos[i][0] < (args.rectangle_pos[0][0]+args.rectangle_pos[1][0]):
                #         node_vel[i][1] = -0.99*node_vel[i][1] 
                #     else:
                #         node_vel[i][0] = -0.99*node_vel[i][0] 
                 
                for j in range(2):
                    if node_pos[i][j] < 0:
                        node_pos[i][j] = 0
                        node_vel[i][j] = -0.99*node_vel[i][j]  
                    if node_pos[i][j] > 15:
                        node_pos[i][j] = 15
                        node_vel[i][j] = -0.99*node_vel[i][j]

                node_pos[i] += args.dt*node_vel[i]
    
        # ax.add_patch(args.rectangle)
        ax.add_patch(args.semicircle)
        plot_voxel(node_pos,element_num, element_idx, args.em)
        plt.xlim(0, 15)
        plt.ylim(0, 9)
        plt.draw()
        plt.tight_layout()
        plt.savefig('images/plot_{}.png'.format(time))
        plt.pause(0.01)
    plt.ioff()

    images = []
    for i in range(args.max_iter):
        filename = 'images/plot_{}.png'.format(i)
        images.append(imageio.imread(filename))
    imageio.mimsave(f'./gifs/Multiple_2D_{model_list[args.model]}_implicit_{args.implict}_{args.mode}_{len(args.em)}.gif', images, duration=0.01)
    shutil.rmtree("./images")