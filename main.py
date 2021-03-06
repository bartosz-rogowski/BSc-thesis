from MC import MC
import time
from math import sqrt, pi
from Particle import Particle

if __name__ == '__main__':
	Start_timex = time.time()

	n = int(4*1e04)
	x_max = 0.8
	y_max = 0.8
	n_c_x = 200
	n_c_y = 200
	N_it = 1
	speed = 10
	mfp_coeff = 1e-2
	bins_number = 200
	d_x = x_max/n_c_x
	d_y = y_max/n_c_y
	modes = ["BASIC", "SOD", "NOH"]
	mode = modes[2]

	r_eff = x_max*y_max / (2*pi*mfp_coeff*min(d_x, d_y)*n)
	print(f"{r_eff = }")
	particles = [None for i in range(n)]
	for i in range(n):
		if mode == "SOD":
			if i < n/10:
				particles[i] = Particle(x_max/2, 0, x_max, y_max, speed=speed, r_eff=r_eff, mode=mode)
			else:
				particles[i] = Particle(0, 0, x_max/2, y_max, speed=speed, r_eff=r_eff, mode=mode)
		else:
			particles[i] = Particle(0, 0, x_max, y_max, speed=speed, r_eff=r_eff, mode=mode)

	start_time = time.perf_counter()
	MC(
		particles=particles,
		n_c_x=n_c_x, 
		n_c_y=n_c_y, 
		N_it=N_it, 
		x_max=x_max, 
		y_max=y_max,
		speed=speed,
		bins_number=bins_number,
		mode=mode
	)

	end_time = time.perf_counter()
	print("--- execution time: %s seconds ---" % round(end_time - start_time, 4))