from sympy import sin, cos, sqrt, Matrix, Symbol

from numpy import pi



q1,q2,q3,q4,q5,q6 = [Symbol('q' + str(i+1)) for i in range(6)]




J = Matrix( 
        [[-0.056*sin(q1)*sin(q5)*cos(q4)*cos(q2 + q3) - 0.056*sin(q1)*sin(q2 + q3)*cos(q5) - 0.431*sin(q1)*sin(q2 + q3) - 0.431*sin(q1)*cos(q2) - 0.056*sin(q4)*sin(q5)*cos(q1), (-0.431*sin(q2) - 0.056*sin(q5)*sin(q2 + q3)*cos(q4) + 0.056*cos(q5)*cos(q2 + q3) + 0.431*cos(q2 + q3))*cos(q1), (-0.056*sin(q5)*sin(q2 + q3)*cos(q4) + 0.056*cos(q5)*cos(q2 + q3) + 0.431*cos(q2 + q3))*cos(q1), -0.056*(sin(q1)*cos(q4) + sin(q4)*cos(q1)*cos(q2 + q3))*sin(q5), -0.056*sin(q1)*sin(q4)*cos(q5) - 0.056*sin(q5)*sin(q2 + q3)*cos(q1) + 0.056*cos(q1)*cos(q4)*cos(q5)*cos(q2 + q3), 0],
        [-0.056*sin(q1)*sin(q4)*sin(q5) + 0.056*sin(q5)*cos(q1)*cos(q4)*cos(q2 + q3) + 0.056*sin(q2 + q3)*cos(q1)*cos(q5) + 0.431*sin(q2 + q3)*cos(q1) + 0.431*cos(q1)*cos(q2), (-0.431*sin(q2) - 0.056*sin(q5)*sin(q2 + q3)*cos(q4) + 0.056*cos(q5)*cos(q2 + q3) + 0.431*cos(q2 + q3))*sin(q1), (-0.056*sin(q5)*sin(q2 + q3)*cos(q4) + 0.056*cos(q5)*cos(q2 + q3) + 0.431*cos(q2 + q3))*sin(q1), 0.056*(-sin(q1)*sin(q4)*cos(q2 + q3) + cos(q1)*cos(q4))*sin(q5), -0.056*sin(q1)*sin(q5)*sin(q2 + q3) + 0.056*sin(q1)*cos(q4)*cos(q5)*cos(q2 + q3) + 0.056*sin(q4)*cos(q1)*cos(q5), 0],
        [0, -0.056*sin(q5)*cos(q4)*cos(q2 + q3) - 0.056*sin(q2 + q3)*cos(q5) - 0.431*sin(q2 + q3) - 0.431*cos(q2), -0.056*sin(q5)*cos(q4)*cos(q2 + q3) - 0.056*sin(q2 + q3)*cos(q5) - 0.431*sin(q2 + q3), 0.056*sin(q4)*sin(q5)*sin(q2 + q3), -0.056*sin(q5)*cos(q2 + q3) - 0.056*sin(q2 + q3)*cos(q4)*cos(q5), 0],
        [1, (sin(q5)*sin(q2 + q3)*cos(q4) - cos(q5)*cos(q2 + q3))*sin(q4)*sin(q5)/(sin(q4)**2*sin(q5)**2 + sin(q5)**2*cos(q4)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), (sin(q5)*sin(q2 + q3)*cos(q4) - cos(q5)*cos(q2 + q3))*sin(q4)*sin(q5)/(sin(q4)**2*sin(q5)**2 + sin(q5)**2*cos(q4)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), (sin(q5)*cos(q2 + q3) + sin(q2 + q3)*cos(q4)*cos(q5))*sin(q5)/(sin(q4)**2*sin(q5)**2 + sin(q5)**2*cos(q4)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), sin(q4)*sin(q2 + q3)/(sin(q4)**2*sin(q5)**2 + sin(q5)**2*cos(q4)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), 0],
        [0, sqrt(2)*(2*sin(q2 + q3 - q5) + 2*sin(q2 + q3 + q5) - sin(q2 + q3 - q4 - q5) + sin(q2 + q3 - q4 + q5) - sin(q2 + q3 + q4 - q5) + sin(q2 + q3 + q4 + q5))/sqrt(-4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20), sqrt(2)*(2*sin(q2 + q3 - q5) + 2*sin(q2 + q3 + q5) - sin(q2 + q3 - q4 - q5) + sin(q2 + q3 - q4 + q5) - sin(q2 + q3 + q4 - q5) + sin(q2 + q3 + q4 + q5))/sqrt(-4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20), 2*sqrt(2)*(-cos(q2 + q3 - q4) + cos(q2 + q3 + q4))*sin(q5)/sqrt(-4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20), sqrt(2)*(-2*sin(q2 + q3 - q5) + 2*sin(q2 + q3 + q5) + sin(q2 + q3 - q4 - q5) + sin(q2 + q3 - q4 + q5) + sin(q2 + q3 + q4 - q5) + sin(q2 + q3 + q4 + q5))/sqrt(-4*cos(2*q4) - 4*cos(2*q5) - 4*cos(2*q2 + 2*q3) + 2*cos(2*q4 - 2*q5) + 2*cos(2*q4 + 2*q5) + 2*cos(2*q2 + 2*q3 - 2*q4) + 2*cos(2*q2 + 2*q3 + 2*q4) - 6*cos(2*q2 + 2*q3 - 2*q5) - 6*cos(2*q2 + 2*q3 + 2*q5) - cos(2*q2 + 2*q3 - 2*q4 - 2*q5) - cos(2*q2 + 2*q3 - 2*q4 + 2*q5) + 4*cos(2*q2 + 2*q3 - q4 - 2*q5) - 4*cos(2*q2 + 2*q3 - q4 + 2*q5) + 4*cos(2*q2 + 2*q3 + q4 - 2*q5) - 4*cos(2*q2 + 2*q3 + q4 + 2*q5) - cos(2*q2 + 2*q3 + 2*q4 - 2*q5) - cos(2*q2 + 2*q3 + 2*q4 + 2*q5) + 20), 0],
        [0, sin(q4)*sin(q5)/(sin(q4)**2*sin(q2 + q3)**2 + sin(q5)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q4)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), sin(q4)*sin(q5)/(sin(q4)**2*sin(q2 + q3)**2 + sin(q5)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q4)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), (sin(q5)*cos(q4)*cos(q2 + q3) + sin(q2 + q3)*cos(q5))*sin(q2 + q3)/(sin(q4)**2*sin(q2 + q3)**2 + sin(q5)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q4)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), (sin(q5)*sin(q2 + q3)*cos(q4) - cos(q5)*cos(q2 + q3))*sin(q4)*sin(q2 + q3)/(sin(q4)**2*sin(q2 + q3)**2 + sin(q5)**2*cos(q2 + q3)**2 + sin(q2 + q3)**2*cos(q4)**2*cos(q5)**2 + cos(2*q2 + 2*q3 - q4 - 2*q5)/8 - cos(2*q2 + 2*q3 - q4 + 2*q5)/8 + cos(2*q2 + 2*q3 + q4 - 2*q5)/8 - cos(2*q2 + 2*q3 + q4 + 2*q5)/8), 1]] 
        )


# print(J == J.simplify())
problem = J[27]

# problem = problem.subs(q1, 0)
# problem = problem.subs(q2, 0)
# problem = problem.subs(q3, 0)
# problem = problem.subs(q4, 0)
# problem = problem.subs(q5, 0)
# problem = problem.subs(q6, 0)

print(problem.simplify())

exit(0)

J = J.subs(q1, 0)
J = J.subs(q2, 0)
J = J.subs(q3, 0)
J = J.subs(q4, 0)
J = J.subs(q5, 0)
J = J.subs(q6, 0)


for i in range(6):
    for j in range(6):
        print(round(J[i * 6 + j], 4), ',   ', end='')
    print()


