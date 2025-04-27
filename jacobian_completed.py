from sympy import sin, cos, sqrt, Matrix, Symbol
import sympy as sp
from numpy import pi, linalg



q1,q2,q3,q4,q5,q6 = [Symbol('q' + str(i+1)) for i in range(6)]




J = Matrix( 
        [[-0.056*sin(q1)*sin(q5)*cos(q4)*cos(q2 + q3) - 0.056*sin(q1)*sin(q2 + q3)*cos(q5) - 0.431*sin(q1)*sin(q2 + q3) - 0.431*sin(q1)*cos(q2) - 0.056*sin(q4)*sin(q5)*cos(q1), (-0.431*sin(q2) - 0.056*sin(q5)*sin(q2 + q3)*cos(q4) + 0.056*cos(q5)*cos(q2 + q3) + 0.431*cos(q2 + q3))*cos(q1), (-0.056*sin(q5)*sin(q2 + q3)*cos(q4) + 0.056*cos(q5)*cos(q2 + q3) + 0.431*cos(q2 + q3))*cos(q1), -0.056*(sin(q1)*cos(q4) + sin(q4)*cos(q1)*cos(q2 + q3))*sin(q5), -0.056*sin(q1)*sin(q4)*cos(q5) - 0.056*sin(q5)*sin(q2 + q3)*cos(q1) + 0.056*cos(q1)*cos(q4)*cos(q5)*cos(q2 + q3), 0],
        [-0.056*sin(q1)*sin(q4)*sin(q5) + 0.056*sin(q5)*cos(q1)*cos(q4)*cos(q2 + q3) + 0.056*sin(q2 + q3)*cos(q1)*cos(q5) + 0.431*sin(q2 + q3)*cos(q1) + 0.431*cos(q1)*cos(q2), (-0.431*sin(q2) - 0.056*sin(q5)*sin(q2 + q3)*cos(q4) + 0.056*cos(q5)*cos(q2 + q3) + 0.431*cos(q2 + q3))*sin(q1), (-0.056*sin(q5)*sin(q2 + q3)*cos(q4) + 0.056*cos(q5)*cos(q2 + q3) + 0.431*cos(q2 + q3))*sin(q1), 0.056*(-sin(q1)*sin(q4)*cos(q2 + q3) + cos(q1)*cos(q4))*sin(q5), -0.056*sin(q1)*sin(q5)*sin(q2 + q3) + 0.056*sin(q1)*cos(q4)*cos(q5)*cos(q2 + q3) + 0.056*sin(q4)*cos(q1)*cos(q5), 0],
        [0, -0.056*sin(q5)*cos(q4)*cos(q2 + q3) - 0.056*sin(q2 + q3)*cos(q5) - 0.431*sin(q2 + q3) - 0.431*cos(q2), -0.056*sin(q5)*cos(q4)*cos(q2 + q3) - 0.056*sin(q2 + q3)*cos(q5) - 0.431*sin(q2 + q3), 0.056*sin(q4)*sin(q5)*sin(q2 + q3), -0.056*sin(q5)*cos(q2 + q3) - 0.056*sin(q2 + q3)*cos(q4)*cos(q5), 0],
        [1, (sin(q5)*sin(q2 + q3)*cos(q4) - cos(q5)*cos(q2 + q3))*sin(q4)*sin(q5)/(sin(q4)**2*sin(q5)**2 + sin(q5)**2*cos(q4)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), (sin(q5)*sin(q2 + q3)*cos(q4) - cos(q5)*cos(q2 + q3))*sin(q4)*sin(q5)/(sin(q4)**2*sin(q5)**2 + sin(q5)**2*cos(q4)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), (sin(q5)*cos(q2 + q3) + sin(q2 + q3)*cos(q4)*cos(q5))*sin(q5)/(sin(q4)**2*sin(q5)**2 + sin(q5)**2*cos(q4)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), sin(q4)*sin(q2 + q3)/(sin(q4)**2*sin(q5)**2 + sin(q5)**2*cos(q4)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), 0],
        [0, sqrt(2)*(2*sin(q2 + q3 - q5) + 2*sin(q2 + q3 + q5) - sin(q2 + q3 - q4 - q5) + sin(q2 + q3 - q4 + q5) - sin(q2 + q3 + q4 - q5) + sin(q2 + q3 + q4 + q5))/sqrt(-4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20), sqrt(2)*(2*sin(q2 + q3 - q5) + 2*sin(q2 + q3 + q5) - sin(q2 + q3 - q4 - q5) + sin(q2 + q3 - q4 + q5) - sin(q2 + q3 + q4 - q5) + sin(q2 + q3 + q4 + q5))/sqrt(-4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20), 2*sqrt(2)*(-cos(q2 + q3 - q4) + cos(q2 + q3 + q4))*sin(q5)/sqrt(-4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20), sqrt(2)*(-2*sin(q2 + q3 - q5) + 2*sin(q2 + q3 + q5) + sin(q2 + q3 - q4 - q5) + sin(q2 + q3 - q4 + q5) + sin(q2 + q3 + q4 - q5) + sin(q2 + q3 + q4 + q5))/sqrt(-4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20), 0],
        [0, sin(q4)*sin(q5)/(sin(q4)**2*sin(q2 + q3)**2 + sin(q5)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q4)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), sin(q4)*sin(q5)/(sin(q4)**2*sin(q2 + q3)**2 + sin(q5)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q4)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), (sin(q5)*cos(q4)*cos(q2 + q3) + sin(q2 + q3)*cos(q5))*sin(q2 + q3)/(sin(q4)**2*sin(q2 + q3)**2 + sin(q5)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q4)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), (sin(q5)*sin(q2 + q3)*cos(q4) - cos(q5)*cos(q2 + q3))*sin(q4)*sin(q2 + q3)/(sin(q4)**2*sin(q2 + q3)**2 + sin(q5)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q4)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), 1]] 
        )



J = J.subs(q1, 0)
J = J.subs(q2, 0)
J = J.subs(q3, pi/2)
J = J.subs(q4, 0)
J = J.subs(q5, 1)
J = J.subs(q6, 0)


# print(J[25].simplify())


# a = -4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20
# print(sp.solve(a))
# exit(0)


for i in range(36):
    if J[i] is sp.nan:
        J[i] = 0


for i in range(6):
    for j in range(6):
        print(round(J[i * 6 + j], 4), end=' ')
    print('0')


