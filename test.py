from threading import Thread, Lock
from Particle import Particle
from math import floor
import multiprocessing as mp
import numpy as np
from time import perf_counter

pi = 3.14

def fun(p, delta_t, lock):
	for p_a in particles:
		for p_b in particles:
			if p_a != p_b:
				p_a.position.x = pi
				p_b.position.y = -pi

def paraller_plz(particles, i, j, n_c_x, n_c_y, lock):
	I_k, I_l = max(0,i-1), min(n_c_x-1, i+1)
	J_k, J_l = max(0,j-1), min(n_c_y-1, j+1)
	p_c = []
	for k in range(I_k, I_l+1):
		for l in range(J_k, J_l+1):
			for p in particles[:4]:
				p_c.append(p)
	lock.acquire()
	fun(p_c, 0.01, lock)
	lock.release()
	print(i, j)

def analyse(data, i, value, lock):
	I_k, I_l = max(0,i-1), min(len(data)-1, i+1)
	p_c = []
	# for k in range(I_k, I_l+1):
	lock.acquire()
	for k in range(I_k, I_l+1):
		data[k] += value
	# for p,_ in enumerate(p_c):
	# 	p_c[p] += value
	# return data
	lock.release()

if __name__ == "__main__":
	t1 = perf_counter()
	n = 2000
	
	particles = [None for i in range(n)]
	for i in range(n):
		particles[i] = Particle(4, 4, speed=10, r_eff=9e-1)

	n_c_x = n_c_y = 10
	# results = None
	# with ProcessPoolExecutor(max_workers=mp.cpu_count()) as pool:
	# 	for i in range(n_c_x):
	# 		for j in range(n_c_y):
	# 			# print(i, j)
	# 			results = pool.map(paraller_plz, (particles, i, j, n_c_x, n_c_y))
	

	lock = Lock()
	# arggss = [ (data, i, 2, lock) for i in range(len(data)) ]

	threads = []
	for i in range(n_c_x):
		for j in range(n_c_y):
			# print(i, j)
			# results = pool.map(paraller_plz, (particles, i, j, n_c_x, n_c_y))
			thread = Thread(target=paraller_plz, args=(particles, i, j, n_c_x, n_c_y, lock))
			thread.start()
			threads.append(thread)
	for thread in threads:
		thread.join()

	bins_number = 100
	dv = 35/bins_number
	hist_v = [0 for i in range(bins_number)]
	for j, p in enumerate(particles):
		i = floor(p.get_speed()/dv)
		hist_v[i] +=1
	
	print(hist_v)

	print(perf_counter() - t1)