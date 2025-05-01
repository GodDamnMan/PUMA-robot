from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

from PUMA import PUMA, StatefulPUMA

import warnings
warnings.filterwarnings("ignore")
# warnings.resetwarnings()





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

    # Функция анимации
    def update(frame):
        global robot_line, frame_artists
        dt = frame * 0.05
        q2 = np.pi/2 * np.sin(dt)
        q1 = np.pi * np.sin(dt)
        q3 = np.pi/2 + np.sin(5*dt)
        theta = [0, 0, 0, dt, dt, 0]
        # theta = [0, 0, np.pi/2, 0, 0, 0]
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


def trapezoidal_in_joint_space(Robot:PUMA, theta_init:list = None, theta_final:list = None, theta_dot_max:list = None, theta_ddot_max:list = None):
    if theta_init is None:
        theta_init = [0, -np.pi/4, np.pi/4, 0, np.pi/6, 0]
    
    if theta_final is None:
        theta_final = [np.pi/6, -np.pi/6, np.pi/3, -np.pi/4, np.pi/3, np.pi/4]

    if theta_dot_max is None:
        theta_dot_max = [1] * 6

    if theta_ddot_max is None:
        theta_ddot_max = [2] * 6

    _, q_const, _, _ = PUMA.pos_vel_acc(theta_init, theta_final, theta_dot_max, theta_ddot_max, dt = 0.0001)
    t, q_disc, Q_dot, Q_ddot = PUMA.pos_vel_acc(theta_init, theta_final, theta_dot_max, theta_ddot_max, dt = 1/120)

    PUMA.plot_trajectories(t, q_disc, Q_dot, Q_ddot)

    x_cont = [Robot.forward_kinematics(qi) for qi in q_const]
    x_disc = [Robot.forward_kinematics(qi) for qi in q_disc]

    x_cont = np.array(x_cont)
    x_disc = np.array(x_disc)

    min_len = min(len(x_cont), len(x_disc))

    x_cont_trimmed = np.array(x_cont[:min_len])[:, :3]
    x_disc_trimmed = np.array(x_disc[:min_len])[:, :3]

    error = np.linalg.norm(x_cont_trimmed - x_disc_trimmed, axis=1)

    plt.plot(error)
    plt.title("Error b/w cont and disc")
    plt.xlabel("t(s)")
    plt.ylabel("error(m)")
    plt.grid(True)
    plt.show()






if __name__ == '__main__':
    Robot = StatefulPUMA()
    

    # Robot.inv_kinematics_tester()
    # Robot.joint_motion_sanity_check()

    # plot_by_joint_movement(Robot)
    # plot_by_ee_movement(Robot)
    

    trapezoidal_in_joint_space(Robot)

    
    ''' choose one type of render '''
    ''' moving chaoticly, w/ workspace '''
    # fig, ax = Robot.plot_workspace(samples=1000, show_plot=False)
    # first_task_render(fig, ax, Robot = Robot)

    ''' moving chaoticly, w/o workspace '''
    # first_task_render(Robot = Robot)

    ''' moving by Jacobian '''
    # second_task_render(Robot=Robot)


    
    