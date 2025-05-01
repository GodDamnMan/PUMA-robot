import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, FFMpegWriter

# === DH and FK Functions ===
dh = [{'d': .672,      'a': 0,     'alpha': -np.pi/2},
      {'d': .139,      'a': .431,  'alpha': 0},
      {'d': -.139,     'a': 0,     'alpha': np.pi/2},
      {'d': .431,      'a': 0,     'alpha': 0},
      {'d': 0,         'a': 0,     'alpha': -np.pi/2},
      {'d': 0,         'a': 0,     'alpha': np.pi/2},
      {'d': .056,      'a': 0,     'alpha': 0}]

def dh_transform(a, alpha, d, theta):
    ct, st = np.cos(theta), np.sin(theta)
    ca, sa = np.cos(alpha), np.sin(alpha)
    return np.array([
        [ct, -st * ca,  st * sa, a * ct],
        [st,  ct * ca, -ct * sa, a * st],
        [0,       sa,      ca,     d],
        [0,        0,       0,     1]
    ])

def forward_kinematics(joint_angles):
    T = np.eye(4)
    positions = [T[:3, 3]]
    for i, angle in enumerate(joint_angles):
        a = dh[i]['a']
        alpha = dh[i]['alpha']
        d = dh[i]['d']
        theta = angle
        A = dh_transform(a, alpha, d, theta)
        T = T @ A
        positions.append(T[:3, 3])
    return np.array(positions)

# === Trajectory Planning ===
def compute_quintic_coeffs(q0, qf, v0, vf, a0, af, t0, tf):
    M = np.array([
        [1, t0, t0**2, t0**3, t0**4, t0**5],
        [0, 1, 2*t0, 3*t0**2, 4*t0**3, 5*t0**4],
        [0, 0, 2, 6*t0, 12*t0**2, 20*t0**3],
        [1, tf, tf**2, tf**3, tf**4, tf**5],
        [0, 1, 2*tf, 3*tf**2, 4*tf**3, 5*tf**4],
        [0, 0, 2, 6*tf, 12*tf**2, 20*tf**3]
    ])
    b = np.array([q0, v0, a0, qf, vf, af])
    return np.linalg.solve(M, b)

def evaluate_polynomial(coeffs, t):
    return np.polyval(coeffs[::-1], t)

# Key waypoints
q_start = [0, 0, 0, 0, 0, 0]
q_middle = [0.5, -0.3, 0.2, 0.1, 0.2, -0.1]
q_end = [1.0, 0.0, 0.5, 0.2, 0.4, 0.0]
t_points = [0, 2, 4]
v0 = vf = a0 = af = 0

# Compute trajectory for each joint
trajectories = []
for j in range(6):
    seg1 = compute_quintic_coeffs(q_start[j], q_middle[j], v0, vf, a0, af, 0, 2)
    seg2 = compute_quintic_coeffs(q_middle[j], q_end[j], v0, vf, a0, af, 2, 4)
    trajectories.append((seg1, seg2))

# Time settings
frames = 100
times = np.linspace(0, 4, frames)
all_positions = []

for t in times:
    q = []
    for j in range(6):
        coeffs = trajectories[j][0] if t <= 2 else trajectories[j][1]
        q.append(evaluate_polynomial(coeffs, t))
    positions = forward_kinematics(q)
    all_positions.append(positions)

# === Animation ===
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(0, 1.5)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.view_init(elev=30, azim=130)

line, = ax.plot([], [], [], '-o', lw=2)

# Mark key end-effector points
for q, c, lbl in zip([q_start, q_middle, q_end], ['red', 'green', 'blue'], ['Start', 'Middle', 'End']):
    pos = forward_kinematics(q)[-1]
    ax.scatter(*pos, color=c, s=50, label=lbl)
ax.legend()

def update(frame):
    pos = all_positions[frame]
    line.set_data(pos[:, 0], pos[:, 1])
    line.set_3d_properties(pos[:, 2])
    return line,

ani = FuncAnimation(fig, update, frames=len(all_positions), blit=False)

from IPython.display import HTML

ani = FuncAnimation(fig, update, frames=len(all_positions), blit=False)
HTML(ani.to_jshtml())
