from MC import MC
import time

if __name__ == '__main__':
	Start_timex = time.time()

	n = int(9*1e02)
	x_max = 4.0
	y_max = 4.0
	n_c_x = 10
	n_c_y = 10
	N_it = 25
	speed = 10

	start_time = time.perf_counter()
	MC(
		n_c_x=n_c_x, 
		n_c_y=n_c_y, 
		N_it=N_it, 
		x_max=x_max, 
		y_max=y_max, 
		n=n,
		speed=speed
	)

	end_time = time.perf_counter()
	print("--- execution time: %s seconds ---" % round(end_time - start_time, 4))