from MC import MC
import time
from math import sqrt, pi
from Particle import Particle

if __name__ == '__main__':
	Start_timex = time.time()

	n = int(2*1e04)
	x_max = 4.0
	y_max = 4.0
	n_c_x = 400
	n_c_y = 400
	N_it = 30
	speed = 10
	mfp_coeff = 1e1
	bins_number = 100
	d_x = x_max/n_c_x
	d_y = y_max/n_c_y

	r_eff = sqrt(speed / (2*pi*mfp_coeff*min(d_x, d_y)*n) )
	print(f"{r_eff = }")
	particles = [None for i in range(n)]
	for i in range(n):
		particles[i] = Particle(x_max, y_max, speed=speed, r_eff=r_eff)

	start_time = time.perf_counter()
	MC(
		particles=particles,
		n_c_x=n_c_x, 
		n_c_y=n_c_y, 
		N_it=N_it, 
		x_max=x_max, 
		y_max=y_max,
		speed=speed,
		bins_number=bins_number
	)

	end_time = time.perf_counter()
	print("--- execution time: %s seconds ---" % round(end_time - start_time, 4))