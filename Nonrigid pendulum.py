# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 16:57:33 2018

@author: Rishi Acharya
"""


import math
import matplotlib

G = 9.80665
DURATION = 5
STEP = 0.00001

L = 0.5
R_IN = 0.5
R_DOT_IN = 0
M = 0.2
THETA_IN = 2
THETA_DOT_IN = 0
K = 10000


def r_2dot(r:float, r_dot:float,
               theta:float, theta_dot:float,
               l: float,
               m:float, k:float
               )->float:
    out = ((G*m*math.cos(theta))+(m*r*(theta_dot**2)))/m
    if l < r:
        out += k*(L-r)/m
    return out


def theta_2dot(r:float, r_dot:float,
               theta:float, theta_dot:float,
               l: float,
               m:float, k:float
               )->float:
    if r != 0:
        out = (-G*math.sin(theta)-2*r_dot*theta_dot)/r
    else:
        out = False
    return out

def get_energy(r:float, r_dot:float,
               theta:float, theta_dot:float,
               l: float,
               m:float, k:float
               )->float:
    e_k = 0.5 * m * (
            ((math.sin(theta)*r_dot)+(math.cos(theta)*r*theta)**2)\
            +((math.cos(theta)*r_dot)-(r*math.sin(theta)*theta_dot)**2)
    )
    e_p = m*G*(l-(r*math.cos(theta)))
    if l < r:
        e_p += 0.5*k*((L-r)**2)
    e_t = e_k + e_p
    return e_t

t = 0

thetas = []
theta_dots = []
rs = []
r_dots = []
xs = []
ys = []
times = []

thetas.append(THETA_IN)
rs.append(R_IN)
theta_dots.append(THETA_DOT_IN)
r_dots.append(R_DOT_IN)
xs.append(R_IN*math.sin(THETA_IN))
ys.append(R_IN*math.cos(THETA_IN))
times.append(t)

print("Start loop")

while t <= DURATION:
    theta_dot = theta_dots[-1] + (theta_2dot(rs[-1], r_dots[-1],
                          thetas[-1], theta_dots[-1],
                          L, M, K)* STEP)
    theta = (thetas[-1] + theta_dot * STEP)
    while theta > math.pi:
        theta -= (2*math.pi)
    while theta < -1*math.pi:
        theta += 2*math.pi
    r_dot = r_dots[-1] + (r_2dot(rs[-1], r_dots[-1],
                          thetas[-1], theta_dots[-1],
                          L, M, K)* STEP)
    r = rs[-1] + (r_dot * STEP)
    theta_dots.append(theta_dot)
    thetas.append(theta)
    r_dots.append(r_dot)
    rs.append(r)
    xs.append(r*math.sin(theta))
    ys.append(r*math.sin(theta))
    times.append(t)
    t += STEP

energy_change = get_energy(rs[-1], r_dots[-1],
                          thetas[-1], theta_dots[-1],
                          L, M, K)\
               - get_energy(rs[0], r_dots[0],
                          thetas[0], theta_dots[0],
                          L, M, K)

print("Energy change = ")
print(str(100*energy_change/get_energy(rs[0], r_dots[0],
                          thetas[0], theta_dots[0],
                          L, M, K))+"%")


fig_theta_time = matplotlib.pyplot.figure()
axes_theta_time = fig_theta_time.add_subplot(111)
axes_theta_time.plot(thetas, rs, c='tab:orange')
axes_theta_time.grid()
axes_theta_time.set(
        xlabel='', ylabel='',
        title='Time theta'
        )