

import numpy as np
import matplotlib.pyplot as plt

def compute_septic_coeffs(q0, qf, v0, vf, a0, af, j0, jf, T):
    
    M = np.array([
        [1, 0, 0, 0,     0,      0,       0,        0],
        [0, 1, 0, 0,     0,      0,       0,        0],
        [0, 0, 2, 0,     0,      0,       0,        0],
        [0, 0, 0, 6,     0,      0,       0,        0],
        [1, T, T**2, T**3, T**4, T**5, T**6, T**7],
        [0, 1, 2*T, 3*T**2, 4*T**3, 5*T**4, 6*T**5, 7*T**6],
        [0, 0, 2, 6*T, 12*T**2, 20*T**3, 30*T**4, 42*T**5],
        [0, 0, 0, 6, 24*T, 60*T**2, 120*T**3, 210*T**4]
    ])
    b = np.array([q0, v0, a0, j0, qf, vf, af, jf])
    return np.linalg.solve(M, b)

def evaluate_poly_derivatives(coeffs, t):
   
    q = np.polyval(coeffs[::-1], t)
    v = np.polyval(np.polyder(coeffs[::-1], 1), t)
    a = np.polyval(np.polyder(coeffs[::-1], 2), t)
    return q, v, a

q_start   = [0, 0, 0, 0, 0, 0]
q_middle1 = [0.5, -0.3, 0.2, 0.1, 0.2, -0.1]
q_middle2 = [1.0, 0.0, 0.5, 0.2, 0.4, 0.0]
q_final   = [0.6, 0.2, 0.1, 0.4, 0.3, 0.2]
waypoints = [q_start, q_middle1, q_middle2, q_final]
t_points = [0, 2, 4, 6]


v0 = a0 = j0 = vf = af = jf = 0


trajectories = []

for j in range(6):
    T1 = t_points[1] - t_points[0]
    coeffs1 = compute_septic_coeffs(q_start[j], q_middle1[j], v0, 0, a0, 0, j0, 0, T1)
    v1 = np.polyval(np.polyder(coeffs1[::-1], 1), T1)
    a1 = np.polyval(np.polyder(coeffs1[::-1], 2), T1)
    j1 = np.polyval(np.polyder(coeffs1[::-1], 3), T1)


    T2 = t_points[2] - t_points[1]
    coeffs2 = compute_septic_coeffs(q_middle1[j], q_middle2[j], v1, 0, a1, 0, j1, 0, T2)
    v2 = np.polyval(np.polyder(coeffs2[::-1], 1), T2)
    a2 = np.polyval(np.polyder(coeffs2[::-1], 2), T2)
    j2 = np.polyval(np.polyder(coeffs2[::-1], 3), T2)

    T3 = t_points[3] - t_points[2]
    coeffs3 = compute_septic_coeffs(q_middle2[j], q_final[j], v2, vf, a2, af, j2, jf, T3)

    trajectories.append([coeffs1, coeffs2, coeffs3])

print("Waypoint States (t, q, qdot, qddot):")
for idx, t in enumerate(t_points):
    q = []
    qdot = []
    qddot = []
    for j in range(6):
        if t <= t_points[1]:
            coeffs = trajectories[j][0]
            local_t = t - t_points[0]
        elif t <= t_points[2]:
            coeffs = trajectories[j][1]
            local_t = t - t_points[1]
        else:
            coeffs = trajectories[j][2]
            local_t = t - t_points[2]
        qj, vj, aj = evaluate_poly_derivatives(coeffs, local_t)
        q.append(round(qj, 5))
        qdot.append(round(vj, 5))
        qddot.append(round(aj, 5))
    print(f"t = {t:.1f}")
    print(f"q     = {q}")
    print(f"qdot  = {qdot}")
    print(f"qddot = {qddot}\n")

time_vec = np.linspace(0, 6, 300)
q_mat = np.zeros((6, len(time_vec)))
qdot_mat = np.zeros((6, len(time_vec)))
qddot_mat = np.zeros((6, len(time_vec)))

for i, t in enumerate(time_vec):
    for j in range(6):
        if t <= t_points[1]:
            coeffs = trajectories[j][0]
            local_t = t - t_points[0]
        elif t <= t_points[2]:
            coeffs = trajectories[j][1]
            local_t = t - t_points[1]
        else:
            coeffs = trajectories[j][2]
            local_t = t - t_points[2]
        q, v, a = evaluate_poly_derivatives(coeffs, local_t)
        q_mat[j, i] = q
        qdot_mat[j, i] = v
        qddot_mat[j, i] = a

titles = ['Joint 1', 'Joint 2', 'Joint 3', 'Joint 4', 'Joint 5', 'Joint 6']
fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

for j in range(6):
    axes[0].plot(time_vec, q_mat[j], label=titles[j])
    axes[1].plot(time_vec, qdot_mat[j])
    axes[2].plot(time_vec, qddot_mat[j])

axes[0].set_ylabel("Position (rad)")
axes[1].set_ylabel("Velocity (rad/s)")
axes[2].set_ylabel("Acceleration (rad/s²)")
axes[2].set_xlabel("Time (s)")
axes[0].set_title("Joint Positions")
axes[1].set_title("Joint Velocities")
axes[2].set_title("Joint Accelerations")
axes[0].legend()
axes[0].grid(True)
axes[1].grid(True)
axes[2].grid(True)
plt.tight_layout()
plt.show()
