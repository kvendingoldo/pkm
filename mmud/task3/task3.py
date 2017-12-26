# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

"""
Michaylov stability criterion
"""

def ƒ(p, i):
    return [(p[0] * i ** 4).real, (p[1] * i ** 3).imag, (p[2] * i ** 2).real, (p[3] * i).imag, p[4]]


λ = [1, 0.4, 1.5, 0.8, 4]
iω = complex(0, 1)

p_iω = ƒ(λ, iω)

ω = np.arange(0, 5, 1e-3)
x_ω = p_iω[0] * ω ** 4 + p_iω[2] * ω ** 2 + p_iω[4]
y_ω = p_iω[1] * ω ** 2 + p_iω[3] * ω

plt.plot(x_ω, y_ω)
plt.xlabel('$x(\omega)$')
plt.ylabel('$y(\omega)$')
plt.grid(True)
plt.show()


