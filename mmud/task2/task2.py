import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt

u = [1, 1]


def system(x):
    return [0.6 * x[0] * x[1] + 0.4 * x[1] + 0.3 * u[0],
            0.3 * x[0] * x[1] + 0.7 * x[1] - 0.1 * u[1]]


def plt_1(x):
    x2 = (-0.3 * u[0]) / (0.6 * x + 0.4)
    plt.plot(x,  x2, color="pink")


def plt_2(x):
    x2 = (0.1 * u[1]) / (0.3 * x + 0.7)
    plt.plot(x,  x2, color="green")


plt_1(np.linspace(-4.5, -0.7, 50))
plt_1(np.linspace(-0.6, 4.5, 50))

plt_2(np.linspace(-4.5, -2.4, 50))
plt_2(np.linspace(-2.3, 4.5, 50))

p0 = optimize.fsolve(system, [-10, 0])

print("Balance position:")
print("p0 = ", p0)

print("Check:")
print(system(p0))

plt.scatter(p0[0], p0[1])

plt.xlim(-5, 5)
plt.ylim(-20, 20)
plt.grid(True)
plt.show()
