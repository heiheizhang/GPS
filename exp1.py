import os
import numpy as np


data_dir = os.path.join(os.path.dirname(os.path.abspath('__file__')), 'data')
txt_path = os.path.join(data_dir, 'exp1_data.txt')


def coords_read_from_txt(txt_path):
    coords = []
    delta_I_T = []
    rho = []
    with open(txt_path) as f:
        for line in f:
            if line.find('SATXYZA') == 1: #find the data line we need
                line_list = line.strip().split(';')#fist, split the str into the head and the content
                content = line_list[1].split(',')
                n_satellites = int(content[1])
                for i in range(n_satellites):
                    #print(content[2+9*i])
                    coords.append([float(content[3+9*i]), float(content[4+9*i]), float(content[5+9*i])])
                    delta_I_T.append([float(content[6 + 9 * i]), float(content[7 + 9 * i]), float(content[8 + 9 * i])])
                #print("=========================================")
            if line.find('RANGEA') == 1:
                line_list = line.strip().split(';')  # fist, split the str into the head and the content
                content = line_list[1].split(',')
                n_satellites = int(content[0]) - 1 #rangea里的数据多了一颗153号卫星
                for i in range(n_satellites):
                    #print(content[3+10*i])
                    rho.append([float(content[3 + 10*i])])
                #print("=========================================")
    #print(np.array(coords))
    X = np.array(coords) #finally, X.shape = (11, 3)
    DIT = np.array(delta_I_T)
    rho = np.array(rho)
    #print(X.shape, DIT.shape, Rho.shape)
    return X, DIT, rho


# def compute_G(X, xk):
#     G = np.ones((X.shape[0], X.shape[1] + 1)) # plus one because of equation (2)
#     for i in range(X.shape[0]):
#         d = np.linalg.norm(X[i, :] - xk)
#         for j in range(X.shape[1]):
#             G[i, j] = (X[i, j] - xk[-1, j]) / d
#     print(G)

def compute_r(X, xk):
    r = np.linalg.norm(X - xk, axis=1).reshape(-1, 1) #计算距离
    return r

def compute_G(X, xk):
    r = compute_r(X, xk)
    G = - (X - xk) / r
    G = np.column_stack((G, np.ones(G.shape[0])))
    return G

def compute_rhoc(rho, DIT):
    return rho + DIT[:, 0].reshape(-1, 1) - np.sum(DIT[:, 1:3], axis=1).reshape(-1, 1) # rhoc = rho + delta - I - T

def compute_dtu(rhoc, r):
    return rhoc - r

x0 = np.array([[1.0, 0.0, 0.0]])
dtu = 0.0
X, DIT, rho = coords_read_from_txt(txt_path)
rhoc = compute_rhoc(rho, DIT)

delta = np.ones((4, 1))
i = 0
while (np.sum(abs(delta)) > 0.00000001):
    print("第%d轮迭代"%i)
    i += 1
    r = compute_r(X, x0)
    #print(r)
    #dtu = rhoc - r
    #print(dtu)
    G = compute_G(X, x0)
    #print(G)
    b = rhoc - r - dtu
    #print(r)
    delta = np.dot(np.dot(np.linalg.inv(np.dot(G.T, G)), G.T), b)
    x0 += delta.T[0, 0:3]
    dtu += delta.T[0, 3]
    print(x0)
