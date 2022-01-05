import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import maxwell

def maxw(size = None, scale = 1.0):
    """Generates size samples of maxwell"""
    vx = np.random.normal(size=size)
    vy = np.random.normal(size=size)
    # vz = np.random.normal(size=size)
    return scale*(np.sqrt(vx*vx + vy*vy))
	
def maxwell2d(x, a=1.0):
	return 2*a*x*np.exp(-a*x**2)

mdata = maxw(100000, 4)
h, bins = np.histogram(mdata, bins = 201, range=(0.0, 30.0))

x = np.linspace(0.0, 30.0, 1000)
rv = maxwell()

fig, ax = plt.subplots(1, 1)

z = zip(rv.pdf(x), x)
# y = [1.2* aa / bb for aa,bb in z]

y1 = [maxwell2d(xx, a=0.006) for xx in x]
ax.plot(x, y1, 'k-.', lw=1, label='$a_1 = 0.006$')

y2 = [maxwell2d(xx, a=0.03) for xx in x]

plt.rcParams.update({'font.size': 14})
# ax.hist(mdata, bins = bins, density=True)
ax.plot(x, y2, 'b--', lw=1, label='$a_2 = 0.03$')

y3 = [maxwell2d(xx, a=0.06) for xx in x]
ax.plot(x, y3, 'g-.', lw=1, label='$a_3 = 0.06$')

y4 = [maxwell2d(xx, a=0.2) for xx in x]
ax.plot(x, y4, 'r', lw=1, label='$a_4 = 0.2$')

plt.title("Rozkład Maxwella-Boltzmanna")
plt.grid(True)
plt.legend()
plt.show()