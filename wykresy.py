from matplotlib import pyplot as plt, colors, cm
import numpy as np
from math import sqrt, pi, exp

speed = 10

def f(dx, n):
    return (sqrt(speed / (2*pi*dx*n) ))

# X = np.logspace(-6, 2, 4000)
# Y = np.logspace(3, 6, 4000)

# Z = np.zeros(len(X)*len(Y))
# i = 0
# for x in X:
#     for y in Y:
#         Z[i] = f(x, y)
#         i += 1

# z = Z.reshape(len(Y), len(X))

# vmin=1e-6
# vmax=1.0
# levels = []
# N_LEVELS = 3
# for E in range(-4,2):
#     levels = np.concatenate((levels[:-1],np.linspace(10**E,10**(E+1),N_LEVELS)))

# plot = plt.contourf(X, Y, z, levels, 
#     norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max())
# )


# plt.colorbar(plot, ticks=[10**(-i) for i in range(-5, 5)])
# plt.xscale('log')
# plt.yscale('log')
# plt.title('$r_{eff}$')
# plt.xlabel('dx')
# plt.ylabel('number of particles')
# plt.grid(True, which='minor', linestyle='--', linewidth=0.2, color='black')
# plt.grid(True, which='major', linestyle='--')
# plt.show()

# plt.close()

if __name__ == '__main__':

    def maxwell(x, a):
        temp = sqrt(2.0/pi)/a**3
        return temp * x**2 * exp(-x**2/ (2.0*a**2) )

    x = np.linspace(0, 3.5*speed, 2000)
    y = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = maxwell(x[i], 6.0)
    plt.plot(x, y)
    plt.show()