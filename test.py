import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import maxwell

def maxw(size = None, scale = 1):
    """Generates size samples of maxwell"""
    vx = np.random.normal(size=size)
    vy = np.random.normal(size=size)
    # vz = np.random.normal(size=size)
    return scale*(np.sqrt(vx*vx + vy*vy))
	
def maxwell2d(x, a=0.50):
	return 2*a*x*np.exp(-a*x**2)

mdata = maxw(100000, 4)
h, bins = np.histogram(mdata, bins = 201, range=(0.0, 30.0))

x = np.linspace(0.0, 30.0, 1000)
rv = maxwell()

fig, ax = plt.subplots(1, 1)

z = zip(rv.pdf(x), x)
# y = [1.2* aa / bb for aa,bb in z]
y = [maxwell2d(xx, a=0.031) for xx in x]

ax.hist(mdata, bins = bins, density=True)
ax.plot(x, y, 'k-', lw=1, label='Maxwell pdf')
plt.title("Maxwell")
plt.show()