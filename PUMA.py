from utils import *
import matplotlib.pyplot as plt
from numpy.random import rand
from numpy import sin, cos, sqrt


class PUMA:
    def __init__(self):  
        self.dh = [{'d':.672,       'a': 0,         'alpha': -np.pi/2, },
                   {'d':.139,       'a': .431,      'alpha': 0},
                   {'d': -.139,     'a': 0,         'alpha': np.pi/2},
                   {'d': .431,      'a': 0,         'alpha': 0},
                   {'d': 0,         'a': 0,         'alpha': -np.pi/2},
                   {'d': 0,         'a': 0,         'alpha': np.pi/2},
                   {'d': .056,      'a': 0,         'alpha': 0}]




    def _Ti(self, theta, n):
        return Rz(theta) @ Tz(self.dh[n]['d']) @ Tx(self.dh[n]['a']) @ Rx(self.dh[n]['alpha'])

    

    def _get_basises(self, theta):
        T1 = self._Ti(theta[0], 0)
        T2 = T1 @ self._Ti(theta[1], 1)
        T3 = T2 @ self._Ti(theta[2], 2)
        Ts = T3 @ self._Ti(0,        3)
        T4 = Ts @ self._Ti(theta[3], 4)
        T5 = T4 @ self._Ti(theta[4], 5)
        T6 = T5 @ self._Ti(theta[5], 6)
        return [np.eye(4), T1, T2, T3, Ts, T4, T5, T6]
    

    def _get_T06(self, theta):
        T = self._get_basises(theta)
        return T[-1]


    def forward_kinematics(self, theta):
        T = self._get_T06(theta)
        origin = T[:3, 3]
        # TODO если нужны другие углы
        e_theta = np.atan2(np.sqrt(T[1,2]**2 + T[0,2]**2), T[2,2]) 
        e_fi = 0
        e_psi = 0

        if e_theta == 0:
            e_fi = np.atan2(-T[0,1], T[0,0])
        else:
            e_fi = np.atan2(T[1,2], T[0,2])
            e_psi = np.atan2(T[2,1], -T[2,0])


        end_effector = list(map(lambda x: float(x), list(origin) + ([e_fi, e_theta, e_psi])))
        return end_effector



    # Прямая кинематика с возвратом всех промежуточных точек
    def forward_kinematics_points(self, theta):
        _, T1, T2, T3, Ts, T4, T5, T6 = self._get_basises(theta)
    
        p1 = np.array([0, 0, 0, 1]).transpose()
        p2 = T1 @ p1
        p3 = T2 @ p1
        ps = T3 @ p1
        p4 = Ts @ p1
        p5 = T4 @ p1
        p6 = T5 @ p1
        p7 = T6 @ p1

        # Возвращаем координаты всех узлов
        return [
        p1,  # База
        p2,  # Плечо
        p3,  # Локоть
        ps,  # Неподвижный joint
        p4,  # Запястье 1
        p5,  # Запястье 2
        p6,  # Запястье 3
        p7   # Конец эффектора
        ]



    def inverse_kinematics(self, target):
        x, y, z, e_fi, e_theta, e_psi = target
        R = Rzyz(e_fi, e_theta, e_psi)

        # Wrist center 
        W = np.array([x, y, z]).transpose() - R @ np.array([0,0,self.dh[6]['d']]).transpose()


        theta1 = np.atan2(W[1],W[0])
        
        l1 = self.dh[1]['a']
        l2 = self.dh[3]['d']
        xe = np.sqrt(W[0]**2 + W[1]**2)
        ze = W[2] - self.dh[0]['d']
        L = np.sqrt(xe**2 + ze ** 2) 
        p = (l1+l2+L)/2
        S = np.sqrt(p * (p - l1) * (p - l2) * (p - L))


        Calpha = (L**2 + l1**2 - l2 **2) / (2*l1*L)
        Salpha = (2*S)/(l1*L)

        Cbeta = (-L**2 + l1**2 + l2 **2) / (2*l1*l2)
        Sbeta = (2*S)/(l1*l2)

        alpha = np.atan2(Salpha,Calpha)
        beta = np.atan2(Sbeta,Cbeta)
        

        theta3 = np.pi/2 + np.pi - beta
        A = np.arctan(ze/xe)
        theta2 = -(alpha + A)


        # Wrist
        _T06 = T06(traget=target)

        T1 =      self._Ti(theta1, 0)
        T2 = T1 @ self._Ti(theta2, 1)
        T3 = T2 @ self._Ti(theta3, 2)
        Ts = T3 @ self._Ti(0,      3)

        theta4 = 0
        theta5 = 0
        theta6 = 0
        # TODO если нужны другие углы
        T36 = Ts.transpose() @ _T06
        theta5 = np.atan2(np.sqrt(T36[0,2]**2 + T36[1,2]**2),T36[2,2])
        if theta5 == 0:
            theta4 = np.atan2(T36[1,0], T36[0,0])
        else:
            theta4 = np.atan2(T36[1,2], T36[0,2])
            theta6 = np.atan2(T36[2,1], -T36[2,0])

        return np.array([theta1, theta2, theta3, theta4 , theta5, theta6])


    def plot_workspace(self, samples=10000, show_plot=True):
        """
        Plot the reachable workspace of the PUMA robot.

        Parameters:
        - samples: number of random joint configurations to sample
        - show_plot: whether to display the plot immediately
        - save_path: if provided, save the plot to this path
        """
        robot = PUMA()
        points = []

        print(f"Calculating workspace with {samples} samples...")

        # Sample random joint angles and collect end effector positions
        for _ in range(samples):
            # Generate random joint angles within reasonable limits
            theta = (np.random.rand(6) * 2 - 1) * np.array([
                np.pi,       # Joint 1: full rotation
                np.pi/2,     # Joint 2: limited range
                np.pi/2,     # Joint 3: limited range
                np.pi,       # Joint 4: full rotation
                np.pi/2,     # Joint 5: limited range
                np.pi        # Joint 6: full rotation
            ])

            # Get end effector position
            pos = robot.forward_kinematics(theta)[:3]
            points.append(pos)

        points = np.array(points)

        # Create 3D plot
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot the points
        ax.scatter(points[:, 0], points[:, 1], points[:, 2], 
                   s=1, alpha=0.5)

        # Plot settings
        ax.set_title('PUMA 560')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_xlim([-1.5, 1.5])
        ax.set_ylim([-1.5, 1.5])
        ax.set_zlim([0, 2])
        ax.grid(True)



        if show_plot:
            plt.show()

        return fig, ax
    


    def tester(self, n = 100):
        cnt = 0
        for _ in range(n):
            theta = (rand(6)*2-1) * np.pi
            theta_inv = self.inverse_kinematics(self.forward_kinematics(theta))

            pos = self.forward_kinematics(theta)
            pos_inv = self.forward_kinematics(theta_inv)

            err = [abs(i-j) for i,j in zip(pos, pos_inv)]

            if max(err) > 1e-10:
                print(pos)
                print(pos_inv)
                print('-')
                print(theta)
                print(theta_inv)
                print('-----------------')
            else:
                cnt += 1

        print(cnt, 'out of', n, 'tests passed')
    

    
    def get_jacobian(self, theta:list):
        q1, q2, q3, q4, q5, _ = theta


        J = [[-0.056*sin(q1)*sin(q5)*cos(q4)*cos(q2 + q3) - 0.056*sin(q1)*sin(q2 + q3)*cos(q5) - 0.431*sin(q1)*sin(q2 + q3) - 0.431*sin(q1)*cos(q2) - 0.056*sin(q4)*sin(q5)*cos(q1), (-0.431*sin(q2) - 0.056*sin(q5)*sin(q2 + q3)*cos(q4) + 0.056*cos(q5)*cos(q2 + q3) + 0.431*cos(q2 + q3))*cos(q1), (-0.056*sin(q5)*sin(q2 + q3)*cos(q4) + 0.056*cos(q5)*cos(q2 + q3) + 0.431*cos(q2 + q3))*cos(q1), -0.056*(sin(q1)*cos(q4) + sin(q4)*cos(q1)*cos(q2 + q3))*sin(q5), -0.056*sin(q1)*sin(q4)*cos(q5) - 0.056*sin(q5)*sin(q2 + q3)*cos(q1) + 0.056*cos(q1)*cos(q4)*cos(q5)*cos(q2 + q3), 0],
            [-0.056*sin(q1)*sin(q4)*sin(q5) + 0.056*sin(q5)*cos(q1)*cos(q4)*cos(q2 + q3) + 0.056*sin(q2 + q3)*cos(q1)*cos(q5) + 0.431*sin(q2 + q3)*cos(q1) + 0.431*cos(q1)*cos(q2), (-0.431*sin(q2) - 0.056*sin(q5)*sin(q2 + q3)*cos(q4) + 0.056*cos(q5)*cos(q2 + q3) + 0.431*cos(q2 + q3))*sin(q1), (-0.056*sin(q5)*sin(q2 + q3)*cos(q4) + 0.056*cos(q5)*cos(q2 + q3) + 0.431*cos(q2 + q3))*sin(q1), 0.056*(-sin(q1)*sin(q4)*cos(q2 + q3) + cos(q1)*cos(q4))*sin(q5), -0.056*sin(q1)*sin(q5)*sin(q2 + q3) + 0.056*sin(q1)*cos(q4)*cos(q5)*cos(q2 + q3) + 0.056*sin(q4)*cos(q1)*cos(q5), 0],
            [0, -0.056*sin(q5)*cos(q4)*cos(q2 + q3) - 0.056*sin(q2 + q3)*cos(q5) - 0.431*sin(q2 + q3) - 0.431*cos(q2), -0.056*sin(q5)*cos(q4)*cos(q2 + q3) - 0.056*sin(q2 + q3)*cos(q5) - 0.431*sin(q2 + q3), 0.056*sin(q4)*sin(q5)*sin(q2 + q3), -0.056*sin(q5)*cos(q2 + q3) - 0.056*sin(q2 + q3)*cos(q4)*cos(q5), 0],
            [1, (sin(q5)*sin(q2 + q3)*cos(q4) - cos(q5)*cos(q2 + q3))*sin(q4)*sin(q5)/(sin(q4)**2*sin(q5)**2 + sin(q5)**2*cos(q4)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), (sin(q5)*sin(q2 + q3)*cos(q4) - cos(q5)*cos(q2 + q3))*sin(q4)*sin(q5)/(sin(q4)**2*sin(q5)**2 + sin(q5)**2*cos(q4)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), (sin(q5)*cos(q2 + q3) + sin(q2 + q3)*cos(q4)*cos(q5))*sin(q5)/(sin(q4)**2*sin(q5)**2 + sin(q5)**2*cos(q4)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), sin(q4)*sin(q2 + q3)/(sin(q4)**2*sin(q5)**2 + sin(q5)**2*cos(q4)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), 0],
            [0, sqrt(2)*(2*sin(q2 + q3 - q5) + 2*sin(q2 + q3 + q5) - sin(q2 + q3 - q4 - q5) + sin(q2 + q3 - q4 + q5) - sin(q2 + q3 + q4 - q5) + sin(q2 + q3 + q4 + q5))/sqrt(-4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20), sqrt(2)*(2*sin(q2 + q3 - q5) + 2*sin(q2 + q3 + q5) - sin(q2 + q3 - q4 - q5) + sin(q2 + q3 - q4 + q5) - sin(q2 + q3 + q4 - q5) + sin(q2 + q3 + q4 + q5))/sqrt(-4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20), 2*sqrt(2)*(-cos(q2 + q3 - q4) + cos(q2 + q3 + q4))*sin(q5)/sqrt(-4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20), sqrt(2)*(-2*sin(q2 + q3 - q5) + 2*sin(q2 + q3 + q5) + sin(q2 + q3 - q4 - q5) + sin(q2 + q3 - q4 + q5) + sin(q2 + q3 + q4 - q5) + sin(q2 + q3 + q4 + q5))/sqrt(-4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20), 0],
            [0, sin(q4)*sin(q5)/(sin(q4)**2*sin(q2 + q3)**2 + sin(q5)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q4)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), sin(q4)*sin(q5)/(sin(q4)**2*sin(q2 + q3)**2 + sin(q5)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q4)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), (sin(q5)*cos(q4)*cos(q2 + q3) + sin(q2 + q3)*cos(q5))*sin(q2 + q3)/(sin(q4)**2*sin(q2 + q3)**2 + sin(q5)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q4)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), (sin(q5)*sin(q2 + q3)*cos(q4) - cos(q5)*cos(q2 + q3))*sin(q4)*sin(q2 + q3)/(sin(q4)**2*sin(q2 + q3)**2 + sin(q5)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q4)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), 1]] 
        
        for i,T in enumerate(J):
            for j,a in enumerate(T):
                if np.isnan(a):
                    J[i][j] = 0

        return J













class StatefulPUMA(PUMA):
    def __init__(self):
        super().__init__()

        self.theta = [0,0,0,0,0,0]
        self.ee = self.forward_kinematics(self.theta)

        self.is_unlocked = True
        

    def get_jacobian(self) -> list:
        return super().get_jacobian(self.theta)      # get the jacobian for current state
    
    def set_joints(self, theta) -> None:
        self.theta = theta.copy()
        self.ee = self.forward_kinematics(self.theta)     # automatically set the ee postion

    def set_ee(self, ee_pos) -> None:
        self.ee = ee_pos.copy()
        self.theta = self.inverse_kinematics(ee_pos)     # automatically set all the angles


    


    def move_joints(self, vel, dt) -> list:          # simulate the movement of robot by joints angular speed 
        if vel is not np.array:
            vel = np.array(vel)
        J = self.get_jacobian()
        self.theta += vel * dt
        self.ee += J @ vel * dt
        
        return J @ vel


    def move_end_effector(self, vel_target, dt) -> list:    # simulate the movement of robot by ee velocity
        if vel_target is not np.array:
            vel_target = np.array(vel_target)

        J = self.get_jacobian()
        J_inv = np.linalg.inv(J)
        vel_q = J_inv @ vel_target

        if np.max(np.abs(vel_q)) > 2*np.pi:         # if angular velocity of any joint is too high (higher than 2pi)
            vel_q = np.array([0,0,0,0,0,0])         # then definitely smth went wrong (maybe ee is out of bound), so disable the robot from moving
                                                    
        vel_ee = J @ vel_q
        self.ee += vel_ee * dt
        self.theta += vel_q * dt

        return vel_q



    def joint_motion_sanity_check(self, n = 1000):
        theta_check = self.theta.copy()
        p, f = 0,0
        for i in range(n):
            t = i/100
            self.move_joints(np.array([1,1,1,1,1,1]), 1/100)
            pos = np.array(self.forward_kinematics(np.array(theta_check) + np.array([t,t,t,t,t,t]))) 
            p = pos[:3] - self.ee[:3] 
            f = Rzyz(*pos[3:]) -  Rzyz(*self.ee[3:])

            
            if max(abs(p)) > 1 or any(list(map(lambda x: any(list(map(lambda x: abs(x) > np.pi, x))), f))):
                print('Jacobian motion failed')
                print('iteration:', i)
                print('teoretical theta', np.array(theta_check) + np.array([t,t,t,t,t,t]))
                print('practical theta', self.theta)
                print('teoretical ee pos', self.forward_kinematics(np.array(theta_check) + np.array([t,t,t,t,t,t])))
                print('practical ee pos', self.ee)
                print('error', p, f)
   
                
                return False
            
        print('Jacobian motion passed', n, 'iterations')
        print('last error', p, f)
        return True
