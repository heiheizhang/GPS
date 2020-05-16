import numpy as np
from math import *
a = 6378137.0
alpha = 1.0 / 298.257223563 #椭圆曲率
e2 = alpha * (2 - alpha)
PI = np.pi

def WGStoSpace(lon, lat, alt):
    N = a / sqrt(1 - e2*sin(lat * pi/180)**2)
    x = (N + alt) * cos(lat * pi/180) * cos(lon * pi/180)
    y = (N + alt) * cos(lat * pi/180) * sin(lon * pi/180)
    z = (N * (1-e2) + alt) * sin(lat * pi/180)
    return x, y, z

def SpacetoWGS(x, y, z):
    lat = 0
    lon = 0
    H = 0
    lon = atan(y/x) * 180/pi
    if (lon < 0 and lon > -140):
        lon = lon + 180
    p = sqrt(x**2 + y**2)
    for j in range(6):
        N = a / sqrt(1 - e2*sin(lat)**2)
        lat = atan(z / p / (1 - e2 * N/(N+H)))
        H = p / cos(lat) - N
    lat = lat * 180/pi
    return lon, lat, H

print(SpacetoWGS(-2486346.951481741, 4819685.447984258, 3346005.2513101604))

print(WGStoSpace(117.288, 31.8459, 145.85))