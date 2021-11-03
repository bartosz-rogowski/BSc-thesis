from random import random
from Particle import Particle, plot_particles
from Velocity import Velocity
from MC import MC
import time

if __name__ == '__main__':
	Start_timex = time.time()

	n = int(2*1e02)
	x_max=4.0
	y_max=4.0
	n_c = 10
	N_it = 1

	particles = []
	for i in range(n):
		particles.append(Particle(x_max, y_max, speed=10))
		
	for i in range(1):
		print(i+1)
		start_time = time.time()
		MC(n_c_x=n_c, n_c_y=n_c, N_it=N_it, x_max=x_max, y_max=y_max, particles=particles)
		print(time.time() - start_time)

	print("--- execution time: %s seconds ---" % (time.time() - Start_timex))