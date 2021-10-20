from random import random
from Particle import Particle, plot_particles
from Velocity import Velocity
from MC import MC
import time
start_time = time.time()

n = int(4*1e01)
x_max=4.0
y_max=4.0
n_c = 10
N_it = 26

particles = []
for i in range(n):
	particles.append(Particle(x_max, y_max, speed=1))


MC(n_c_x=n_c, n_c_y=n_c, N_it=N_it, x_max=x_max, y_max=y_max, particles=particles)

print("--- execution time: %s seconds ---" % (time.time() - start_time))