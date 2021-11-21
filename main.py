from MC import MC
import time

if __name__ == '__main__':
	Start_timex = time.time()

	n = int(2*1e04)
	x_max = 4.0
	y_max = 4.0
	n_c_x = 400
	n_c_y = 400
	N_it = 25
	speed = 10
	mfp_coeff = 1e1
	bins_number = 100

	start_time = time.perf_counter()
	MC(
		n_c_x=n_c_x, 
		n_c_y=n_c_y, 
		N_it=N_it, 
		x_max=x_max, 
		y_max=y_max, 
		n=n,
		speed=speed,
		mfp_coeff=mfp_coeff,
		bins_number=bins_number
	)

	end_time = time.perf_counter()
	print("--- execution time: %s seconds ---" % round(end_time - start_time, 4))