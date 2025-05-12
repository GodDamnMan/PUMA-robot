from sympy import symbols, cos, sin, Matrix, pi, diff, lambdify

q1, q2, q3, q4, q5, q6 = symbols('q1 q2 q3 q4 q5 q6')
dq1, dq2, dq3, dq4, dq5, dq6 = symbols('dq1 dq2 dq3 dq4 dq5 dq6')
ddq1, ddq2, ddq3, ddq4, ddq5, ddq6 = symbols('ddq1 ddq2 ddq3 ddq4 ddq5 ddq6')

theta = [q1, q2, q3, q4, q5, q6]
dtheta = [dq1, dq2, dq3, dq4, dq5, dq6]
ddtheta = [ddq1, ddq2, ddq3, ddq4, ddq5, ddq6]

m1 = 10.0
m2 = 17.4
m3 = 4.8
m4 = 0.82
m5 = 0.34
m6 = 0.09
m7 = 0.09

L1 = 0.672
L2 = 0.431
L3 = 0.431
L4 = 0.056
L5 = 0.056
L6 = 0.056
L7 = 0.056


r1 = 0.05
r2 = 0.05
r3 = 0.05
r4 = 0.05 
r5 = 0.05
r6 = 0.05
r7 = 0.05


g = 9.81

com_offsets = [
    {'x': 0, 'y': 0, 'z': 0.336},
    {'x': 0.2155, 'y': 0, 'z': 0},
    {'x': 0, 'y': 0, 'z': 0}, 
    {'x': 0, 'y': 0, 'z': 0.2155},
    {'x': 0, 'y': 0, 'z': 0},
    {'x': 0, 'y': 0, 'z': 0.028},
    {'x': 0, 'y': 0, 'z': 0.028}
]

dh = [{'d':.672,       'a': 0,         'alpha': -pi/2, },
                   {'d':.139,       'a': .431,      'alpha': 0},
                   {'d': -.139,     'a': 0,         'alpha': pi/2},
                   {'d': .431,      'a': 0,         'alpha': 0},
                   {'d': 0,         'a': 0,         'alpha': -pi/2},
                   {'d': 0,         'a': 0,         'alpha': pi/2},
                   {'d': .056,      'a': 0,         'alpha': 0}]

def Rz(theta):
    return Matrix([
        [cos(theta), -sin(theta), 0, 0],
        [sin(theta), cos(theta), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])

def Tz(d):
    return Matrix([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, d],
        [0, 0, 0, 1]
    ])

def Tx(a):
    return Matrix([
        [1, 0, 0, a],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])

def Rx(alpha):
    return Matrix([
        [1, 0, 0, 0],
        [0, cos(alpha), -sin(alpha), 0],
        [0, sin(alpha), cos(alpha), 0],
        [0, 0, 0, 1]
    ])

def Ti(theta, n):
    return Rz(theta) @ Tz(dh[n]['d']) @ Tx(dh[n]['a']) @ Rx(dh[n]['alpha'])


def compute_ode_matrices(q_vals, dq_vals):
    T1 =      Ti(q1, 0)
    T2 = T1 @ Ti(q2, 1)
    T3 = T2 @ Ti(q3, 2)
    Ts = T3 @ Ti(0, 3)
    T4 = Ts @ Ti(q4, 4)
    T5 = T4 @ Ti(q5, 5)
    T6 = T5 @ Ti(q6, 6)
    T7 = T6 @ Ti(0, 6)

    p0 = Matrix([0, 0, 0, 1])
    p1 = T1 @ p0
    p2 = T2 @ p0
    p3 = T3 @ p0
    ps = Ts @ p0
    p4 = T4 @ p0
    p5 = T5 @ p0
    p6 = T6 @ p0
    p7 = T7 @ p0

    com1 = (p0 + p1) / 2
    com2 = (p1 + p2) / 2
    com3 = (p2 + p3) / 2
    com4 = (p3 + p4) / 2
    com5 = (p4 + p5) / 2
    com6 = (p5 + p6) / 2
    com7 = (p6 + p7) / 2

    # Calculate velocities
    vx1 = (diff(p1[0], q1) * dq1)
    vy1 = (diff(p1[1], q1) * dq1)
    vz1 = (diff(p1[2], q1) * dq1)

    vx2 = (diff(p2[0], q1) * dq1 + diff(p2[0], q2) * dq2)
    vy2 = (diff(p2[1], q1) * dq1 + diff(p2[1], q2) * dq2)
    vz2 = (diff(p2[2], q1) * dq1 + diff(p2[2], q2) * dq2)

    vx3 = (diff(p3[0], q1) * dq1 + diff(p3[0], q2) * dq2 + diff(p3[0], q3) * dq3)
    vy3 = (diff(p3[1], q1) * dq1 + diff(p3[1], q2) * dq2 + diff(p3[1], q3) * dq3)
    vz3 = (diff(p3[2], q1) * dq1 + diff(p3[2], q2) * dq2 + diff(p3[2], q3) * dq3)

    vx4 = (diff(p4[0], q1) * dq1 + diff(p4[0], q2) * dq2 + diff(p4[0], q3) * dq3 + diff(p4[0], q4) * dq4)
    vy4 = (diff(p4[1], q1) * dq1 + diff(p4[1], q2) * dq2 + diff(p4[1], q3) * dq3 + diff(p4[1], q4) * dq4)
    vz4 = (diff(p4[2], q1) * dq1 + diff(p4[2], q2) * dq2 + diff(p4[2], q3) * dq3 + diff(p4[2], q4) * dq4)

    vx5 = (diff(p5[0], q1) * dq1 + diff(p5[0], q2) * dq2 + diff(p5[0], q3) * dq3 + diff(p5[0], q4) * dq4 + diff(p5[0], q5) * dq5)
    vy5 = (diff(p5[1], q1) * dq1 + diff(p5[1], q2) * dq2 + diff(p5[1], q3) * dq3 + diff(p5[1], q4) * dq4 + diff(p5[1], q5) * dq5)
    vz5 = (diff(p5[2], q1) * dq1 + diff(p5[2], q2) * dq2 + diff(p5[2], q3) * dq3 + diff(p5[2], q4) * dq4 + diff(p5[2], q5) * dq5)

    vx6 = (diff(p6[0], q1) * dq1 + diff(p6[0], q2) * dq2 + diff(p6[0], q3) * dq3 + diff(p6[0], q4) * dq4 + diff(p6[0], q5) * dq5 + diff(p6[0], q6) * dq6)
    vy6 = (diff(p6[1], q1) * dq1 + diff(p6[1], q2) * dq2 + diff(p6[1], q3) * dq3 + diff(p6[1], q4) * dq4 + diff(p6[1], q5) * dq5 + diff(p6[1], q6) * dq6)
    vz6 = (diff(p6[2], q1) * dq1 + diff(p6[2], q2) * dq2 + diff(p6[2], q3) * dq3 + diff(p6[2], q4) * dq4 + diff(p6[2], q5) * dq5 + diff(p6[2], q6) * dq6)

    vx7 = (diff(com7[0], q1) * dq1 + diff(com7[0], q2) * dq2 + diff(com7[0], q3) * dq3 + diff(com7[0], q4) * dq4 + diff(com7[0], q5) * dq5 + diff(com7[0], q6) * dq6)
    vy7 = (diff(com7[1], q1) * dq1 + diff(com7[1], q2) * dq2 + diff(com7[1], q3) * dq3 + diff(com7[1], q4) * dq4 + diff(com7[1], q5) * dq5 + diff(com7[1], q6) * dq6)
    vz7 = (diff(com7[2], q1) * dq1 + diff(com7[2], q2) * dq2 + diff(com7[2], q3) * dq3 + diff(com7[2], q4) * dq4 + diff(com7[2], q5) * dq5 + diff(com7[2], q6) * dq6)


    U1 = m1 * g * com1[2]   
    U2 = m2 * g * com2[2]
    U3 = m3 * g * com3[2]
    U4 = m4 * g * com4[2]
    U5 = m5 * g * com5[2]
    U6 = m6 * g * com6[2]
    U7 = m7 * g * com7[2]

    U = U1 + U2 + U3 + U4 + U5 + U6 + U7

    G = Matrix([
    [(diff(U, q1))],
    [(diff(U, q2))],
    [(diff(U, q3))],
    [(diff(U, q4))],
    [(diff(U, q5))],
    [(diff(U, q6))],
    # [(diff(U, q7))],
    ])

    z0 = Matrix([0, 0, 1])
    z1 = T1[:3, 2]
    z2 = T2[:3, 2]
    z3 = T3[:3, 2]
    z4 = T4[:3, 2]
    z5 = T5[:3, 2]

    omega1 = dq1 * z0
    omega2 = dq1 * z0 + dq2 * z1
    omega3 = dq1 * z0 + dq2 * z1 + dq3 * z2
    omega4 = dq1 * z0 + dq2 * z1 + dq3 * z2 + dq4 * z3
    omega5 = dq1 * z0 + dq2 * z1 + dq3 * z2 + dq4 * z3 + dq5 * z4
    omega6 = dq1 * z0 + dq2 * z1 + dq3 * z2 + dq4 * z3 + dq5 * z4 + dq6 * z5
    omega7 = dq1 * z0 + dq2 * z1 + dq3 * z2 + dq4 * z3 + dq5 * z4 + dq6 * z5



    I1 = Matrix([
    [(1/12)*m1*(3*r1**2 + L1**2), 0, 0],
    [0, (1/12)*m1*(3*r1**2 + L1**2), 0],
    [0, 0, (1/2)*m1*r1**2]
    ])
    I2 = Matrix([
    [(1/12)*m2*(3*r2**2 + L2**2), 0, 0],
    [0, (1/12)*m2*(3*r2**2 + L2**2), 0],
    [0, 0, (1/2)*m2*r2**2]
    ])
    I3 = Matrix([
    [(1/12)*m3*(3*r3**2 + L3**2), 0, 0],
    [0, (1/12)*m3*(3*r3**2 + L3**2), 0],
    [0, 0, (1/2)*m3*r3**2]
    ])
    I4 = Matrix([
    [(1/12)*m4*(3*r4**2 + L4**2), 0, 0],
    [0, (1/12)*m4*(3*r4**2 + L4**2), 0],
    [0, 0, (1/2)*m4*r4**2]
    ])
    I5 = Matrix([
    [(1/12)*m5*(3*r5**2 + L5**2), 0, 0],
    [0, (1/12)*m5*(3*r5**2 + L5**2), 0],
    [0, 0, (1/2)*m5*r5**2]
    ])
    I6 = Matrix([
    [(1/12)*m6*(3*r6**2 + L6**2), 0, 0],
    [0, (1/12)*m6*(3*r6**2 + L6**2), 0],
    [0, 0, (1/2)*m6*r6**2]
    ])
    I7 = Matrix([
    [(1/12)*m7*(3*r7**2 + L7**2), 0, 0],
    [0, (1/12)*m7*(3*r7**2 + L7**2), 0],
    [0, 0, (1/2)*m7*r7**2]
    ])


    R1 = T1[:3, :3]
    R2 = T2[:3, :3]
    R3 = T3[:3, :3]
    R4 = T4[:3, :3]
    R5 = T5[:3, :3]
    R6 = T6[:3, :3]
    R7 = T7[:3, :3]

    I1_world = R1 * I1 * R1.T
    I2_world = R2 * I2 * R2.T
    I3_world = R3 * I3 * R3.T
    I4_world = R4 * I4 * R4.T
    I5_world = R5 * I5 * R5.T
    I6_world = R6 * I6 * R6.T
    I7_world = R7 * I7 * R7.T

    # Kinetic energy
    v_com1 = Matrix([vx1, vy1, vz1])
    v_com2 = Matrix([vx2, vy2, vz2])
    v_com3 = Matrix([vx3, vy3, vz3])
    v_com4 = Matrix([vx4, vy4, vz4])
    v_com5 = Matrix([vx5, vy5, vz5])
    v_com6 = Matrix([vx6, vy6, vz6])
    v_com7 = Matrix([vx7, vy7, vz7])

    v_com1_squared = (v_com1.T * v_com1)[0]
    v_com2_squared = (v_com2.T * v_com2)[0]
    v_com3_squared = (v_com3.T * v_com3)[0]
    v_com4_squared = (v_com4.T * v_com4)[0]
    v_com5_squared = (v_com5.T * v_com5)[0]
    v_com6_squared = (v_com6.T * v_com6)[0]
    v_com7_squared = (v_com7.T * v_com7)[0]

    omega1_I_omega1 = (omega1.T * I1_world * omega1)[0]
    omega2_I_omega2 = (omega2.T * I2_world * omega2)[0]
    omega3_I_omega3 = (omega3.T * I3_world * omega3)[0]
    omega4_I_omega4 = (omega4.T * I4_world * omega4)[0]
    omega5_I_omega5 = (omega5.T * I5_world * omega5)[0]
    omega6_I_omega6 = (omega6.T * I6_world * omega6)[0]
    omega7_I_omega7 = (omega7.T * I7_world * omega7)[0]

    E1 = (1/2) * m1 * v_com1_squared + (1/2) * omega1_I_omega1
    E2 = (1/2) * m2 * v_com2_squared + (1/2) * omega2_I_omega2
    E3 = (1/2) * m3 * v_com3_squared + (1/2) * omega3_I_omega3
    E4 = (1/2) * m4 * v_com4_squared + (1/2) * omega4_I_omega4
    E5 = (1/2) * m5 * v_com5_squared + (1/2) * omega5_I_omega5
    E6 = (1/2) * m6 * v_com6_squared + (1/2) * omega6_I_omega6
    E7 = (1/2) * m7 * v_com7_squared + (1/2) * omega7_I_omega7

    E_total = (E1 + E2 + E3 + E4 + E5 + E6 + E7)

    # M
    M11 = (diff(diff(E_total, dq1), dq1))
    M12 = (diff(diff(E_total, dq1), dq2))
    M13 = (diff(diff(E_total, dq1), dq3))
    M14 = (diff(diff(E_total, dq1), dq4))
    M15 = (diff(diff(E_total, dq1), dq5))
    M16 = (diff(diff(E_total, dq1), dq6))

    M21 = (diff(diff(E_total, dq2), dq1))
    M22 = (diff(diff(E_total, dq2), dq2))
    M23 = (diff(diff(E_total, dq2), dq3))
    M24 = (diff(diff(E_total, dq2), dq4))
    M25 = (diff(diff(E_total, dq2), dq5))
    M26 = (diff(diff(E_total, dq2), dq6))

    M31 = (diff(diff(E_total, dq3), dq1))
    M32 = (diff(diff(E_total, dq3), dq2))
    M33 = (diff(diff(E_total, dq3), dq3))
    M34 = (diff(diff(E_total, dq3), dq4))
    M35 = (diff(diff(E_total, dq3), dq5))
    M36 = (diff(diff(E_total, dq3), dq6))

    M41 = (diff(diff(E_total, dq4), dq1))
    M42 = (diff(diff(E_total, dq4), dq2))
    M43 = (diff(diff(E_total, dq4), dq3))
    M44 = (diff(diff(E_total, dq4), dq4))
    M45 = (diff(diff(E_total, dq4), dq5))
    M46 = (diff(diff(E_total, dq4), dq6))

    M51 = (diff(diff(E_total, dq5), dq1))
    M52 = (diff(diff(E_total, dq5), dq2))
    M53 = (diff(diff(E_total, dq5), dq3))
    M54 = (diff(diff(E_total, dq5), dq4))
    M55 = (diff(diff(E_total, dq5), dq5))
    M56 = (diff(diff(E_total, dq5), dq6))

    M61 = (diff(diff(E_total, dq6), dq1))
    M62 = (diff(diff(E_total, dq6), dq2))
    M63 = (diff(diff(E_total, dq6), dq3))
    M64 = (diff(diff(E_total, dq6), dq4))
    M65 = (diff(diff(E_total, dq6), dq5))
    M66 = (diff(diff(E_total, dq6), dq6))

    M = Matrix([
    [M11, M12, M13, M14, M15, M16],
    [M21, M22, M23, M24, M25, M26],
    [M31, M32, M33, M34, M35, M36],
    [M41, M42, M43, M44, M45, M46],
    [M51, M52, M53, M54, M55, M56],
    [M61, M62, M63, M64, M65, M66],
    ])

    

    def compute_coriolis_matrix(M, theta, dtheta) -> Matrix:
        n = len(theta)
        C = Matrix.zeros(n, n)
        for k in range(n):
            for j in range(n):
                C_kj = 0
                for i in range(n):
                    c_ijk = (1/2) * (
                        diff(M[k, j], theta[i]) +
                        diff(M[k, i], theta[j]) -
                        diff(M[i, j], theta[k])
                    )
                    C_kj += c_ijk * dtheta[i]
                C[k, j] = (C_kj)
        return C

    C = compute_coriolis_matrix(M, theta, dtheta)


    with open("ode.txt", 'w') as f:
        f.write("M = ")
        f.write(str(M.tolist()))
        f.write("\n")

        f.write("C = ")
        f.write(str(C.tolist()))
        f.write("\n")

        f.write("G = ")
        f.write(str(G.tolist()))
        f.write("\n")


