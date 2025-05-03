
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

dh = [
    {'d': 0.672, 'a': 0,     'alpha': -np.pi/2},
    {'d': 0.139, 'a': 0.431, 'alpha': 0},
    {'d': -0.139,'a': 0,     'alpha': np.pi/2},
    {'d': 0.431, 'a': 0,     'alpha': 0},
    {'d': 0,     'a': 0,     'alpha': -np.pi/2},
    {'d': 0,     'a': 0,     'alpha': np.pi/2},
    {'d': 0.056, 'a': 0,     'alpha': 0}
]

def dh_transform(a, alpha, d, theta):
    ct, st = np.cos(theta), np.sin(theta)
    ca, sa = np.cos(alpha), np.sin(alpha)
    return np.array([
        [ct, -st * ca,  st * sa, a * ct],
        [st,  ct * ca, -ct * sa, a * st],
        [0,       sa,      ca,     d],
        [0,        0,       0,     1]
    ])

def forward_kinematics(joint_angles):
    T = np.eye(4)
    positions = [T[:3, 3]]
    for i, angle in enumerate(joint_angles):
        a = dh[i]['a']
        alpha = dh[i]['alpha']
        d = dh[i]['d']
        theta = angle
        A = dh_transform(a, alpha, d, theta)
        T = T @ A
        positions.append(T[:3, 3])
    return np.array(positions)

def compute_septic_coeffs(q0, qf, v0, vf, a0, af, j0, jf, t0, tf):
    t = t0
    T = tf
    M = np.array([
        [1, t,   t**2,   t**3,    t**4,     t**5,      t**6,       t**7],
        [0, 1,    2*t,    3*t**2,  4*t**3,   5*t**4,    6*t**5,     7*t**6],
        [0, 0,    2,       6*t,     12*t**2,  20*t**3,   30*t**4,    42*t**5],
        [0, 0,    0,       6,        24*t,     60*t**2,   120*t**3,   210*t**4],
        [1, T,   T**2,   T**3,    T**4,     T**5,      T**6,       T**7],
        [0, 1,    2*T,    3*T**2,  4*T**3,   5*T**4,    6*T**5,     7*T**6],
        [0, 0,    2,       6*T,     12*T**2,  20*T**3,   30*T**4,    42*T**5],
        [0, 0,    0,       6,        24*T,     60*T**2,   120*T**3,   210*T**4]
    ])
    b = np.array([q0, v0, a0, j0, qf, vf, af, jf], dtype=np.float64)
    coeffs = np.linalg.solve(M, b)
    return coeffs

def evaluate_poly_derivatives(coeffs, t):
    q = np.polyval(coeffs[::-1], t)
    v = np.polyval(np.polyder(coeffs[::-1], 1), t)
    a = np.polyval(np.polyder(coeffs[::-1], 2), t)
    j = np.polyval(np.polyder(coeffs[::-1], 3), t)
    return q, v, a, j

q_start = [0, 0, 0, 0, 0, 0]
q_middle1 = [0.5, -0.3, 0.2, 0.1, 0.2, -0.1]
q_middle2 = [1.0, 0.0, 0.5, 0.2, 0.4, 0.0]
q_final = [0.6, 0.2, 0.1, 0.4, 0.3, 0.2]
waypoints = [q_start, q_middle1, q_middle2, q_final]
t_points = [0, 2, 4, 6]

v0 = a0 = j0 = vf = af = jf = 0

trajectories = []
for j in range(6):
    q0 = q_start[j]
    q1 = q_middle1[j]
    q2 = q_middle2[j]
    q3 = q_final[j]

    coeffs1 = compute_septic_coeffs(q0, q1, v0, 0, a0, 0, j0, 0, t_points[0], t_points[1])
    _, v1, a1, j1 = evaluate_poly_derivatives(coeffs1, t_points[1])
    coeffs2 = compute_septic_coeffs(q1, q2, v1, 0, a1, 0, j1, 0, t_points[1], t_points[2])
    _, v2, a2, j2 = evaluate_poly_derivatives(coeffs2, t_points[2])

    coeffs3 = compute_septic_coeffs(q2, q3, v2, vf, a2, af, j2, jf, t_points[2], t_points[3])

    trajectories.append([coeffs1, coeffs2, coeffs3])

frames = 150
times = np.linspace(0, 6, frames)
all_positions = []



for t in times:
    q = []
    for j in range(6):
        if t <= t_points[1]:
            coeffs = trajectories[j][0]
        elif t <= t_points[2]:
            coeffs = trajectories[j][1]
        else:
            coeffs = trajectories[j][2]


        qj, _, _, _ = evaluate_poly_derivatives(coeffs, t)
        q.append(qj)
    pos = forward_kinematics(q)
    all_positions.append(pos)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(0, 1.5)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.view_init(elev=30, azim=130)

line, = ax.plot([], [], [], '-o', lw=2)

for q_point, color, label in zip(waypoints, ['red', 'green', 'blue', 'orange'], ['Start', 'Mid1', 'Mid2', 'End']):
    pos = forward_kinematics(q_point)[-1]
    ax.scatter(*pos, color=color, s=50, label=label)
ax.legend()

def update(frame):
    pos = all_positions[frame]
    line.set_data(pos[:, 0], pos[:, 1])
    line.set_3d_properties(pos[:, 2])
    return line,

ani = FuncAnimation(fig, update, frames=len(all_positions), blit=False)
HTML(ani.to_jshtml())
