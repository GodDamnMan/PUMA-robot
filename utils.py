import numpy as np



# Матрицы преобразований dh
def Rz(theta):
    return np.array([
        [np.cos(theta), -np.sin(theta), 0, 0],
        [np.sin(theta), np.cos(theta), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])

def Tz(d):
    return np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, d],
        [0, 0, 0, 1]
    ])

def Tx(a):
    return np.array([
        [1, 0, 0, a],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])

def Rx(alpha):
    return np.array([
        [1, 0, 0, 0],
        [0, np.cos(alpha), -np.sin(alpha), 0],
        [0, np.sin(alpha), np.cos(alpha), 0],
        [0, 0, 0, 1]
    ])

# Матрицы поворота через углы эйлера zyz
def Rzyz(fi,theta,psi):
    cf, ct, cp = np.cos([fi, theta, psi])
    sf, st, sp = np.sin([fi, theta, psi])
    # TODO если нужны другие углы
    R = np.array([[cf*ct*cp - sf*sp,    -cf*ct*sp - sf*cp,    cf*st],
                  [sf*ct*cp + cf*sp,    -sf*ct*sp + cf*cp,    sf*st],
                  [-st*cp          ,    st*sp            ,    ct   ]])
    return R


def T06(traget):
    R = Rzyz(traget[3],traget[4],traget[5])
    return np.array([[R[0,0],   R[0,1],   R[0,2],   traget[0]],
                     [R[1,0],   R[1,1],   R[1,2],   traget[1]],
                     [R[2,0],   R[2,1],   R[2,2],   traget[2]],
                     [     0,        0,        0,          1]])
