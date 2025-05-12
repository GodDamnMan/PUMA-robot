from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

from PUMA import PUMA, StatefulPUMA
import Poly as pl
from symbols import compute_tau_numeric


import warnings
warnings.filterwarnings("ignore")
# warnings.resetwarnings()




#first h/w
robot_line = None
frame_artists = []
def first_task_render(*args, basis_visible:bool = True, Robot:PUMA = None):
    # Настройка 3D-графика
    fig, ax = None, None
    args = list(args)
    if len(args) != 0:
        fig = args[0]
        ax = args[1]
    else:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([0, 1.5])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('PUMA 560 Robot Animation')

    def init():
        global robot_line
        robot_line, = ax.plot([], [], [], 'ko-', lw=3, markersize=8)
        return [robot_line]

    x_start = np.array([0.4, 0.0, 0.5, 0, np.pi/2, 0])
    x_end = np.array([0.2, 0.3, 0.7, 0, np.pi/2, 0])
    
    ax.plot([x_start[0], x_end[0]],
            [x_start[1], x_end[1]],
            [x_start[2], x_end[2]],
            'b-', linewidth=1.5, label='Goal Line')
    
    theta_start = Robot.inverse_kinematics(x_start)
    theta_end = Robot.inverse_kinematics(x_end)
    
    theta_dot_max = [1.0] * 6
    theta_ddot_max = [2.0] * 6
    
    v_max = 0.5
    a_max = 1.0
    dt = 1/120
    
    distance = np.linalg.norm(x_end[:3] - x_start[:3])
    t_acc = v_max / a_max
    d_acc = 0.5 * a_max * t_acc**2
    
    if 2 * d_acc < distance:
        t_const = (distance - 2 * d_acc) / v_max
        t_total = 2 * t_acc + t_const
    else:
        t_acc = np.sqrt(distance / a_max)
        t_const = 0
        t_total = 2 * t_acc
    
    t = np.arange(0, t_total, dt)
    
    cartesian_trajectory = []
    joint_trajectory = []
    
    for ti in t:
        if ti < t_acc:
            s = 0.5 * a_max * ti**2
        elif ti < (t_acc + t_const):
            s = d_acc + v_max * (ti - t_acc)
        else:
            t_dec = ti - t_acc - t_const
            s = d_acc + (t_const * v_max) + v_max * t_dec - 0.5 * a_max * t_dec**2
        
        ratio = s / distance
        pos = x_start[:3] + (x_end[:3] - x_start[:3]) * ratio
        target = np.concatenate([pos, x_start[3:]])

        theta = Robot.inverse_kinematics(target)
        cartesian_trajectory.append(target)
        joint_trajectory.append(theta)
    
    ee_traj = np.array([Robot.forward_kinematics(theta)[:3] for theta in joint_trajectory])
    ax.plot(ee_traj[:, 0], ee_traj[:, 1], ee_traj[:, 2], 'r--', linewidth=1.5, label='Actual')
    ax.legend()

    def update(frame):
        global robot_line, frame_artists
        theta = joint_trajectory[frame]
        points = Robot.forward_kinematics_points(theta)

        x = [p[0] for p in points]
        y = [p[1] for p in points]
        z = [p[2] for p in points]
        robot_line.set_data(x, y)
        robot_line.set_3d_properties(z)

        if basis_visible:
            for artist in frame_artists:
                artist.remove()
            frame_artists = []
            basises = Robot._get_basises(theta)
            scale = 0.2
            for i, basis in enumerate(basises):
                origin = basis[:3, 3]
                x_line, = ax.plot([origin[0], origin[0] + basis[0,0]*scale],
                                 [origin[1], origin[1] + basis[1,0]*scale],
                                 [origin[2], origin[2] + basis[2,0]*scale], 'r-', lw=2)
                y_line, = ax.plot([origin[0], origin[0] + basis[0,1]*scale],
                                 [origin[1], origin[1] + basis[1,1]*scale],
                                 [origin[2], origin[2] + basis[2,1]*scale], 'g-', lw=2)
                z_line, = ax.plot([origin[0], origin[0] + basis[0,2]*scale],
                                 [origin[1], origin[1] + basis[1,2]*scale],
                                 [origin[2], origin[2] + basis[2,2]*scale], 'b-', lw=2)
                label = ax.text(origin[0], origin[1], origin[2], f'J{i}', fontsize=10)
                frame_artists.extend([x_line, y_line, z_line, label])

        fig.canvas.draw()
        return [robot_line] + frame_artists

    ani = FuncAnimation(fig, update, frames=len(joint_trajectory),
                        init_func=init, blit=True, interval=50, repeat=True)

    plt.tight_layout()
    plt.show()



#second h/w
def second_task_render(*args, basis_visible:bool = True, Robot:StatefulPUMA = None):
    # Настройка 3D-графика
    fig, ax = None, None
    args = list(args)
    if len(args) != 0:
        fig = args[0]
        ax = args[1]
    else:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([0, 1.5])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('PUMA 560 Robot Animation')


    def init():
        global robot_line
        robot_line, = ax.plot([], [], [], 'ko-', lw=3, markersize=8)
        return [robot_line]
    
    Robot.set_joints([0,0,0,0,1,0])
    # Функция анимации
    def update(frame):
        global robot_line, frame_artists

        Robot.move_end_effector([-.35 * np.sin(frame / 100 * np.pi), 0, 0, 0, 0, 0], 0.05)
        theta = Robot.theta
        points = Robot.forward_kinematics_points(theta)


        # Соединяем точки
        x = [p[0] for p in points]
        y = [p[1] for p in points]
        z = [p[2] for p in points]

        robot_line.set_data(x, y)
        robot_line.set_3d_properties(z)



        if basis_visible:
            # Очистка предыдущих кадров
            for artist in frame_artists:
                artist.remove()
            frame_artists = []
            basises = Robot._get_basises(theta)
            # Отрисовка систем координат
            scale = 0.2
            for i, basis in enumerate(basises):
                origin = basis[:3, 3]

                # Создаем новые линии для осей
                x_line, = ax.plot([origin[0], origin[0] + basis[0,0]*scale],
                                 [origin[1], origin[1] + basis[1,0]*scale],
                                 [origin[2], origin[2] + basis[2,0]*scale], 'r-', lw=2)
                y_line, = ax.plot([origin[0], origin[0] + basis[0,1]*scale],
                                 [origin[1], origin[1] + basis[1,1]*scale],
                                 [origin[2], origin[2] + basis[2,1]*scale], 'g-', lw=2)
                z_line, = ax.plot([origin[0], origin[0] + basis[0,2]*scale],
                                 [origin[1], origin[1] + basis[1,2]*scale],
                                 [origin[2], origin[2] + basis[2,2]*scale], 'b-', lw=2)

                # Добавляем подписи
                label = ax.text(origin[0], origin[1], origin[2], f'J{i}', fontsize=10)

                frame_artists.extend([x_line, y_line, z_line, label])


        fig.canvas.draw()
        return [robot_line] + frame_artists

    # Запуск анимации
    ani = FuncAnimation(fig, update, frames=200, init_func=init,
                        blit=True, interval=50, repeat=True)

    plt.tight_layout()
    plt.show()

def plot_by_joint_movement(Robot:StatefulPUMA = None, theta_init:list = None, theta_dot = None):
    if theta_init is None:
        theta_init = np.deg2rad([30, -45, 60, 45, -30, 10])   # default theta0

    if theta_dot is None:
        theta_dot = np.array([-3, 1, -.1, 0, 0, 0])    # default joints_velocity
    else:
        theta_dot = np.array(theta_dot)


    total = 20                                # sec in simulation
    dts = [1, 0.5, 0.01]                      # sec per iter
    timesteps = [int(total/i) for i in dts]   # iters in simulation

    fig1, axs1 = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
    axs1 = axs1.flatten() 
    fig2, axs2 = plt.subplots(6, 2, figsize=(12, 9), sharex=True)
    axs2 = axs2.flatten()

    for dt, timestep in zip(dts, timesteps):
        Robot.set_joints(theta_init)
        positions = [Robot.ee]       # XYZ FTP
        velocities = []      # vel of XYZ FTP
        joints = [theta_init]          # thetas

        for _ in range(timestep):
            ee_vel = Robot.move_joints(theta_dot, dt)
            theta = Robot.theta.copy()
            pos = Robot.ee.copy()

            joints.append(theta)
            positions.append(pos)
            velocities.append(ee_vel)

        joints     = np.array(joints[:-1])     # for same initial condition and length as t
        positions  = np.array(positions[:-1])    # same
        velocities = np.array(velocities)
        t = np.arange(timestep) * dt
        ticks = range(total + 1)

        # plot joint angles
        for i, ax in enumerate(axs1):
            angles_deg = np.rad2deg(joints[:, i])
            ax.plot(t, angles_deg, label=f"dt={dt}")
            ax.set_ylabel(f"q{i+1} deg")
            ax.grid(True)
            ax.legend()
            ax.set_xticks(ticks)

        axs1[-2].set_xlabel("Time [s]")
        axs1[-1].set_xlabel("Time [s]")

        # plot ee position
        for i in range(6):
            axs2[2*i].plot(t, positions[:, i], label=f"Position dt={dt}")
            axs2[2*i].set_ylabel(["X", "Y", "Z"][i] + " [m]" if i < 3 else ["Fi", "Theta", "Psi"][i%3] + " [rad]")
            axs2[2*i].grid(True)
            axs2[2*i].legend()
            axs2[2*i].set_xticks(ticks)

            axs2[2*i+1].plot(t, velocities[:, i], label=f"Velocity dt={dt}")
            axs2[2*i+1].set_ylabel("v_" + ["X", "Y", "Z"][i] + " [m/s]" if i < 3 else "w_" + ["Fi", "Theta", "Psi"][i%3] + " [rad/s]")
            axs2[2*i+1].grid(True)
            axs2[2*i+1].legend()
            axs2[2*i+1].set_xticks(ticks)

    axs2[-2].set_xlabel("Time [s]")
    axs2[-1].set_xlabel("Time [s]")


    for axs in (axs1, axs2):
        axs[-1].set_xlabel("Time [s]")

    fig1.suptitle("Joint angles over time")
    fig2.suptitle("End-effector position over time")
    plt.tight_layout()
    plt.show()

def plot_by_ee_movement(Robot:StatefulPUMA = None, theta_init:list = None, ee_dot = None):
    if theta_init is None:
        theta_init = np.deg2rad([30, -45, 60, 45, -30, 10])   # default theta0

    if ee_dot is None:
        ee_dot = np.array([-.005, .01, -.008, .0, .1, .0])    # default ee_velocity
    else:
        ee_dot = np.array(ee_dot)

    total = 20                                # sec in simulation
    dts = [1, 0.5, 0.01]                      # sec per iter
    timesteps = [int(total/i) for i in dts]   # iters in simulation

    fig1, axs1 = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
    axs1=axs1.flatten()
    fig2, axs2 = plt.subplots(6, 2, figsize=(12, 9), sharex=True)
    axs2 = axs2.flatten()

    for dt, timestep in zip(dts, timesteps):
        Robot.set_joints(theta_init)
        positions = [Robot.ee]       # XYZ
        velocities = []      # vel of thetas
        joints = [theta_init]          # thetas

        for _ in range(timestep):
            q_vel = Robot.move_end_effector(ee_dot, dt)
            theta = Robot.theta.copy()
            pos = Robot.ee.copy()

            joints.append(theta)
            positions.append(pos)
            velocities.append(q_vel)

        joints     = np.array(joints[:-1])     # for same initial condition and length as t
        positions  = np.array(positions[:-1])    # same
        velocities = np.array(velocities)
        t = np.arange(timestep) * dt
        ticks = range(total + 1)

        # plot ee pos
        for i, ax in enumerate(axs1):
            # angles_deg = np.rad2deg(joints[:, i])
            ax.plot(t, positions[:, i], label=f"dt={dt}")
            ax.set_ylabel(["X", "Y", "Z"][i] + " [m]" if i < 3 else ["Fi", "Theta", "Psi"][i%3] + " [rad]")
            ax.set_xticks(ticks)
            ax.grid(True)
            ax.legend()

        axs1[-2].set_xlabel("Time [s]")
        axs1[-1].set_xlabel("Time [s]")

        # plot joint ang, vel
        for i in range(6):
            angles_deg = np.rad2deg(joints[:, i])
            axs2[2*i].plot(t, angles_deg, label=f"Position dt={dt}")
            axs2[2*i].set_ylabel(f"q{i+1} [deg]")
            axs2[2*i].grid(True)
            axs2[2*i].legend()
            axs2[2*i].set_xticks(ticks)

            speed_deg = np.rad2deg(velocities[:, i])
            axs2[2*i+1].plot(t, speed_deg, label=f"Angular speed dt={dt}")
            axs2[2*i+1].set_ylabel(f"q{i+1} [deg/s]")
            axs2[2*i+1].grid(True)
            axs2[2*i+1].legend()
            axs2[2*i+1].set_xticks(ticks)
    
    axs2[-2].set_xlabel("Time [s]")
    axs2[-1].set_xlabel("Time [s]")

    for axs in (axs1, axs2):
        axs[-1].set_xlabel("Time [s]")

    fig2.suptitle("Joint angles and angular velocities over time")
    fig1.suptitle("End-effector position over time")
    plt.tight_layout()
    plt.show()





#third h/w
def trapezoidal_in_joint_space(Robot:PUMA, theta_init:list = None, theta_final:list = None, theta_dot_max:list = None, theta_ddot_max:list = None, compute_error:bool = True, all_joints_0_to_pi:bool = False):
    if all_joints_0_to_pi:
        theta_init = [0] * 6
        theta_final = [np.pi] * 6
    else:
        if theta_init is None:
            theta_init = [0, -np.pi/4, np.pi/4, 0, np.pi/6, 0]
        if theta_final is None:
            theta_final = [np.pi/6, -np.pi/6, np.pi, -np.pi/4, np.pi/3, np.pi/4]

    if theta_dot_max is None:
        theta_dot_max = [1] * 6

    if theta_ddot_max is None:
        theta_ddot_max = [2] * 6

    t, q_disc, Q_dot, Q_ddot = Robot.pos_vel_acc(theta_init, theta_final, theta_dot_max, theta_ddot_max, dt = 1/120)
    Robot.plot_trajectories(t, q_disc, Q_dot, Q_ddot)

    if compute_error:
        _, q_cont, _, _ = Robot.pos_vel_acc(theta_init, theta_final, theta_dot_max, theta_ddot_max, dt = 1/120/100)
        x_cont = [Robot.forward_kinematics(qi) for qi in q_cont]
        x_disc = [Robot.forward_kinematics(qi) for qi in q_disc]

        x_cont = np.array(list(filter((0).__ne__,[x if i%100==0 else 0 for i,x in enumerate(x_cont)])))
        x_disc = np.array(x_disc)

        min_len = min(len(x_cont), len(x_disc))

        x_cont_trimmed = np.array(x_cont[:min_len])[:, :3]
        x_disc_trimmed = np.array(x_disc[:min_len])[:, :3]
        error = np.linalg.norm(x_cont_trimmed - x_disc_trimmed, axis=1)

        plt.plot(t[:min_len], error)
        plt.title("Error b/w cont and disc")
        plt.xlabel("t(s)")
        plt.ylabel("error(m)")
        plt.grid(True)
        plt.show()

def polynomial_in_joint_space(Robot:PUMA, theta_init:list = None, theta_final:list = None):
    if theta_init is None:
        theta_init = [0, -np.pi/4, np.pi/4, 0, np.pi/6, 0]
    
    if theta_final is None:
        theta_final = [np.pi/6, -np.pi/6, np.pi, -np.pi/4, np.pi/3, np.pi/4]

    a = pl.compute_trajectory_4_points(theta_init=theta_init, theta_final=theta_final)
    Robot.plot_trajectories(a['time'], a['joint_trajectory'], a['joint_velocity'], a['joint_acceleration'])


if __name__ == '__main__':
    Robot = StatefulPUMA()

    # fig, ax = Robot.plot_workspace(samples=1000, show_plot=False)
    # first_task_render(fig, ax, Robot = Robot)
    

    # Robot.inv_kinematics_tester()
    # Robot.joint_motion_sanity_check()

    # plot_by_joint_movement(Robot)
    # plot_by_ee_movement(Robot)
    

    # trapezoidal_in_joint_space(Robot, all_joints_0_to_pi=True)
    # polynomial_in_joint_space(Robot)

    theta_dot_max = [1.0] * 6
    theta_ddot_max = [2.0] * 6
    t, Q, Q_dot, Q_ddot = Robot.pos_vel_acc([0]*6, [np.pi]*6, theta_dot_max, theta_ddot_max, dt=1/120)

    torques = []
    for i in range(len(t)):
        q = Q[i]
        dq = Q_dot[i]
        ddq = Q_ddot[i]

        q7 = np.zeros(7)
        q7[:6] = q
        dq7 = np.zeros(7)
        dq7[:6] = dq
        ddq7 = np.zeros(7)
        ddq7[:6] = ddq
        tau = compute_tau_numeric(q7, dq7, ddq7)
        torques.append(tau[:6])
    torques = np.array(torques)

    plt.figure(figsize=(12, 8))
    for i in range(torques.shape[1]):
        plt.plot(t, torques[:, i], label=f'Joint {i+1}')
    plt.xlabel('Time (s)')
    plt.ylabel('Torque (Nm)')
    plt.title('Joint Torques Over Time')
    plt.legend()
    plt.grid(True)
    plt.show()

    
    ''' choose one type of render '''
    ''' moving chaoticly, w/ workspace '''

    ''' moving chaoticly, w/o workspace '''
    # first_task_render(Robot = Robot)

    ''' moving by Jacobian '''
    # second_task_render(Robot=Robot)

