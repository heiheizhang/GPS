import os
import numpy as np
import re
from math import *

data_dir = os.path.join(os.path.dirname(os.path.abspath('__file__')), 'data')
txt_path = os.path.join(data_dir, 'exp2_data.txt')


def rawephemb_read_from_txt(txt_path):
    rawephemb = []
    n = 0
    sync = "AA44121C" # more precisely, it should be "AA4412", but one data part contains the "AA4412"
    with open(txt_path) as f:
        for line in f:
            #matchobj = re.findall(sync, line)
            line_list = line.strip().split(sync)
            for s in line_list:
                if (s != ''):
                    rawephemb.append(sync + s)
                    n += 1
            break
    return rawephemb, n

def data_spilt(rawephemb):
    temp = []
    start = 80 # rawephemb[start] = '8' then "B 0E ..."
    real_start = 80 + 30*2 + 7*2 # jump the data we don't use
    for i in range(len(rawephemb)):
        temp.append(rawephemb[i]) # first split each letter
    C_rs      = ''.join(temp[real_start+0 : real_start+4])
    #print(C_rs)
    delt_n    = ''.join(temp[real_start+4 : real_start+8])
    #print(delt_n)
    M_0       = ''.join(temp[real_start+8 : real_start+16])
    #print(M_0)
    C_uc      = ''.join(temp[real_start+16 : real_start+20])
    #print(C_uc)
    e_s       = ''.join(temp[real_start+20 : real_start+28])
    #print(e_s)
    C_us      = ''.join(temp[real_start+28 : real_start+32])
    #print(C_us)
    A_sqrt    = ''.join(temp[real_start+32 : real_start+40])
    #print(A_sqrt)
    t_oe      = ''.join(temp[real_start+40 : real_start+44])
    #print(t_oe)
    t         = ''.join(temp[real_start+52 : real_start+57])
    #print(t)
    C_ic      = ''.join(temp[real_start+58 : real_start+62])
    #print(C_ic)
    OMEGA_0   = ''.join(temp[real_start+62 : real_start+70])
    #print(OMEGA_0)
    C_is      = ''.join(temp[real_start+70 : real_start+74])
    #print(C_is)
    i_0       = ''.join(temp[real_start+74 : real_start+82])
    #print(i_0)
    C_rc      = ''.join(temp[real_start+82 : real_start+86])
    #print(C_rc)
    w         = ''.join(temp[real_start+86 : real_start+94])
    #print(w)
    OMEGA_dot = ''.join(temp[real_start+94 : real_start+100])
    #print(OMEGA_dot)
    i_dot     = ''.join(temp[real_start+102 : real_start+106])
    #print(i_dot)

    #print('\n')
    return {
        'C_rs':two_complement(C_rs, -5),
        'delt_n':two_complement(delt_n, -43) * pi,
        'M_0':two_complement(M_0, -31) * pi,
        'C_uc':two_complement(C_uc, -29),
        'e_s':HexToDeci(e_s, -33),
        'C_us':two_complement(C_us, -29),
        'A_sqrt':HexToDeci(A_sqrt, -19),
        't_oe':HexToDeci(t_oe, 4),
        't':t_GPS(t, 0),
        'C_ic':two_complement(C_ic, -29),
        'OMEGA_0':two_complement(OMEGA_0, -31) * pi,
        'C_is':two_complement(C_is, -29),
        'i_0':two_complement(i_0, -31) * pi,
        'C_rc':two_complement(C_rc, -5),
        'w':two_complement(w, -31) * pi,
        'OMEGA_dot':two_complement(OMEGA_dot, -43) * pi,
        'i_dot':displacement(i_dot, -43) * pi
    }

def hex2dec(hex_str, sf): #sf means scale factor
    temp = []
    l = len(hex_str)
    fan = pow(2, 4*l)
    if (l%2 == 0): #如果字节数为2的倍数
        for i in range(l//2 - 1 , -1, -1): #自己试验一下
           temp.append(hex_str[2*i:2*i+2])
    else:
        temp.append(hex_str[l-1]) #先加入最高位，剩下的和上面一样处理
        l = l - 1
        for i in range(l // 2 - 1, -1, -1):
            temp.append(hex_str[2 * i:2 * i + 2])
    num_str = ''.join(temp)
    num = int(num_str, 16)
    if (int(num_str[0], 16) > 7):#最高位数字大于7要变成负数
        num = num - fan
    return num * pow(2, sf)

def two_complement(hex_str, sf):
    l = len(hex_str)
    fan = pow(2, 4 * l)
    sum = 0
    for i in range(l):
        sum = sum * 16 + int(hex_str[i], 16)

    if (int(hex_str[0], 16) > 7):#最高位数字大于7要变成负数
        sum = sum - fan
    return sum * pow(2, sf)

def displacement(hex_str, sf):
    l = len(hex_str)
    fan = pow(2, 4 * l)
    sum = 0
    for i in range(l):
        sum = sum * 16 + int(hex_str[i], 16)

    if (int(hex_str[0], 16) > 7):#最高位数字大于7要变成负数
        sum = sum / 4
        sum = sum - fan
    else:
        sum = sum / 4
    return sum * pow(2, sf)

def HexToDeci(hex_str, sf):
    l = len(hex_str)
    #fan = pow(2, 4 * l)
    sum = 0.0
    for i in range(l):
        sum = sum * 16 + int(hex_str[i], 16)

    return sum * pow(2, sf)

def t_GPS(hex_str, sf):
    l = len(hex_str)
    #fan = pow(2, 4 * l)
    sum = 0.0
    for i in range(l):
        sum = sum * 16 + int(hex_str[i], 16)

    sum = sum / 8
    sum = sum*6 + 20
    return sum * pow(2, sf)

# first we prepare the data
rawephemb, satnum = rawephemb_read_from_txt(txt_path)

C_rs = []#
delt_n = []#
M_0 = []#
C_uc = []#
e_s = []#
C_us = []
A_sqrt = []#
t_oe =[]#
t = []
C_ic = []
OMEGA_0 = []
C_is = []
i_0 = []
C_rc = []
w = []
OMEGA_dot = []
i_dot = []

for i in range(len(rawephemb)):
    #print(len(rawephemb[i]))
    para_dict = data_spilt(rawephemb[i])
    #print(parm_dict['C_rs'])
    C_rs.append(para_dict['C_rs'])
    delt_n.append(para_dict['delt_n'])
    M_0.append(para_dict['M_0'])
    C_uc.append(para_dict['C_uc'])
    e_s.append(para_dict['e_s'])
    C_us.append(para_dict['C_us'])
    A_sqrt.append(para_dict['A_sqrt'])
    t_oe.append(para_dict['t_oe'])
    t.append(para_dict['t'])
    C_ic.append(para_dict['C_ic'])
    OMEGA_0.append(para_dict['OMEGA_0'])
    C_is.append(para_dict['C_is'])
    i_0.append(para_dict['i_0'])
    C_rc.append(para_dict['C_rc'])
    w.append(para_dict['w'])
    OMEGA_dot.append(para_dict['OMEGA_dot'])
    i_dot.append(para_dict['i_dot'])
#rint(C_rs)
#print(delt_n)
#print(M_0)
# C_uc = []#
# e_s = []#
# C_us = []
# A_sqrt = []#
# t_oe =[]#
# t = []
# C_ic = []
# OMEGA_0 = []
# C_is = []
# i_0 = []
# C_rc = []
# w = []
# OMEGA_dot = []
# i_dot = []
#second we calculate the formulations
GM = 3.986005e14
w_e = 7.2921151467e-5
threshhold = 10e-6
pi = 3.14159265358979
tk = []
X = []
Y = []
Z = []
Vx = []
Vy = []
Vz = []

for i in range(satnum):
    tk.append(t[i] - t_oe[i])
    if (tk[i] > 302400):
        tk[i] = tk[i] - 604800
    if (tk[i] < -302400):
        tk[i] = tk[i] + 604800

    n_0 = sqrt(GM) / pow(A_sqrt[i], 3)
    n = n_0 + delt_n[i]
    M = M_0[i] + n*tk[i]

    x1 = M - e_s[i]
    x2 = M + e_s[i]
    #while (fabs(x2 - x1) > threshhold):
    for j in range(20):
        Q = (x1 + x2) / 2
        v = M + e_s[i]*sin(Q) - Q
        if (v == 0):
            E = Q
            break
        else:
            v_x1 = M + e_s[i]*sin(x1) - x1
            if (v*v_x1 > 0):
                x1 = Q
            else:
                x2 = Q

    E = (x1 + x2) / 2
    d_E = n / ( 1 - e_s[i]*cos(E) )

    temp1 = sqrt( 1 - pow(e_s[i], 2) ) * sin(E)
    temp2 = cos(E) - e_s[i]
    f = atan2(temp1, temp2)
    u_1 = w[i] + f
    d_u_1 = (sqrt( (1 - e_s[i]**2) * d_E )) / ( 1 - e_s[i]*cos(E) )

    delt_u = C_uc[i] * cos(2*u_1) + C_us[i] * sin(2*u_1)
    delt_r = C_rc[i] * cos(2*u_1) + C_rs[i] * sin(2*u_1)
    delt_i = C_ic[i] * cos(2*u_1) + C_is[i] * sin(2*u_1)
    d_delt_u = 2 * d_u_1 * ( C_us[i]*cos(2*u_1) ) - C_uc[i]*sin(2*u_1)
    d_delt_r = 2 * d_u_1 * ( C_rs[i]*cos(2*u_1) ) - C_rc[i]*sin(2*u_1)
    d_delt_i = 2 * d_u_1 * ( C_is[i]*cos(2*u_1) ) - C_ic[i]*sin(2*u_1)

    r = pow(A_sqrt[i], 2) * ( 1 - e_s[i]*cos(E) )
    uu = u_1 + delt_u
    r += delt_r
    I = i_0[i] + delt_i + i_dot[i] * tk[i]
    d_uu = d_u_1 + delt_u
    d_r = pow(A_sqrt[i], 2) * e_s[i] * d_E * sin(E) + d_delt_r
    d_I = i_dot[i] + d_delt_i

    x = r * cos(uu)
    y = r * sin(uu)
    d_x = d_r*cos(uu) - r*d_uu*sin(uu)
    d_y = d_r*sin(uu) - r*d_uu*cos(uu)

    L = OMEGA_0[i] + (OMEGA_dot[i] - w_e) * tk[i] - w_e*t_oe[i]
    d_L = OMEGA_dot[i] - w_e

    X.append( x * cos(L) - y * cos(I) * sin(L) )
    Y.append( x * sin(L) + y * cos(I) * cos(L) )
    Z.append( y * sin(I) )

    Vx.append( -Y[i]*d_L - ( d_y*cos(I) - y*sin(I)*d_I ) * sin(L) + d_x*cos(L) )
    Vy.append( X[i]*d_L + ( d_y*cos(I) - y*sin(I)*d_I ) * cos(L) + d_x*sin(L) )
    Vz.append( d_y*sin(L) + y*d_I*cos(I) )

    t1 = 0.078
    X[i] = X[i]*cos(w_e*t1) + Y[i]*sin(w_e*t1)
    Y[i] = -X[i]*sin(w_e*t1) + Y[i]*cos(w_e*t1)
    #Z[i]不变

for i in range(satnum):
    print(X[i],',', Y[i],',', Z[i])
    pass



