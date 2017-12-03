import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt

u = [8, 10]


def system(x):
    return [0.6 * x[0] * x[1] + 0.4 * x[1] + 0.3 * u[0],
            0.3 * x[0] * x[1] + 0.7 * x[1] - 0.1 * u[1]]


def plt_1(x):
    x2 = (-0.3 * u[0]) / (0.6 * x + 0.4)
    plt.plot(x,  x2, color="pink")


def plt_2(x):
    x2 = (0.1 * u[1] - 0.7 * x) / (0.3 * x)
    plt.plot(x,  x2, color="green")


plt_1(np.linspace(-3, -0.3, 30))
plt_1(np.linspace(0.1, 3, 30))

plt_2(np.linspace(-3, 0, 23))
plt_2(np.linspace(-0.6, 3, 37))

x0 = optimize.fsolve(system, [-2., 0.])
x1 = optimize.fsolve(system, [1., 0.])
x2 = optimize.fsolve(system, [-1., -40.])

print("Balance position:")
print("x0 = ", x0)
print("x1 = ", x1)
print("x2 = ", x2)

print("Check:")
print(system(x0))
print(system(x1))
print(system(x2))

plt.scatter(x0[0], x0[1])
plt.scatter(x1[0], x1[1])
plt.scatter(x2[0], x2[1])

plt.xlim(-5, 5)
plt.ylim(-20, 20)
plt.grid(True)

plt.show()
