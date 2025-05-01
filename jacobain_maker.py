import sympy as sp 
from sympy import Symbol
import numpy as np



def Rz(theta):
    return sp.Matrix([
        [sp.cos(theta), -sp.sin(theta), 0, 0],
        [sp.sin(theta), sp.cos(theta), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])

def Tz(d):
    return sp.Matrix([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, d],
        [0, 0, 0, 1]
    ])

def Tx(a):
    return sp.Matrix([
        [1, 0, 0, a],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])

def Rx(alpha):
    return sp.Matrix([
        [1, 0, 0, 0],
        [0, sp.cos(alpha), -sp.sin(alpha), 0],
        [0, sp.sin(alpha), sp.cos(alpha), 0],
        [0, 0, 0, 1]
    ])


def T_i(theta, n, dh):
        return Rz(theta) @ Tz(dh[n]['d']) @ Tx(dh[n]['a']) @ Rx(dh[n]['alpha'])




dh = [{'d':.672,       'a': 0,         'alpha': -sp.pi/2, },
      {'d':.139,       'a': .431,      'alpha': 0},
      {'d': -.139,     'a': 0,         'alpha': sp.pi/2},
      {'d': .431,      'a': 0,         'alpha': 0},
      {'d': 0,         'a': 0,         'alpha': -sp.pi/2},
      {'d': 0,         'a': 0,         'alpha': sp.pi/2},
      {'d': .056,      'a': 0,         'alpha': 0}]



q = [Symbol('q' + str(i+1)) for i in range(6)]

T = sp.Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])

T = T @ T_i(q[0], 0, dh)
T = T @ T_i(q[1], 1, dh)
T = T @ T_i(q[2], 2, dh)
T = T @ T_i(   0, 3, dh)
T = T @ T_i(q[3], 4, dh)
T = T @ T_i(q[4], 5, dh)
T = T @ T_i(q[5], 6, dh)

T.simplify()

x = T[3]
y = T[7]
z = T[11]

x_dot = [x.diff(i).simplify() for i in q]
y_dot = [y.diff(i).simplify() for i in q]
z_dot = [z.diff(i).simplify() for i in q]



# R = [[T[0], T[1], T[2]],
#      [T[4], T[5], T[6]],
#      [T[8], T[9], T[10]]


fi = sp.atan2(T[6], T[2]).simplify()
theta = sp.atan2(sp.sqrt(T[6]**2 + T[2]**2), T[10]).simplify()
psi = sp.atan2(T[9], -T[8]).simplify()

fi_dot = [fi.diff(i).simplify() for i in q]
theta_dot = [theta.diff(i).simplify() for i in q]
psi_dot = [psi.diff(i).simplify() for i in q]

J = [x_dot, y_dot, z_dot, fi_dot, theta_dot, psi_dot]