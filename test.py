from math import *
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
    sum = 0.0
    for i in range(l):
        sum = sum * 16 + int(hex_str[i], 16)

    if (int(hex_str[0], 16) > 7):#最高位数字大于7要变成负数
        sum = sum - fan
    return sum * pow(2, sf)

print(two_complement('123456', 0))