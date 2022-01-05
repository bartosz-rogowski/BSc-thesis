from Particle import Particle, plot_particles
from math import floor, sqrt, pi, isclose
from random import random
from matplotlib import pyplot as plt
import time
from celluloid import Camera
from tools import *
import multiprocessing as mp
from threading import Thread, Lock
import numpy as np

def solve_wall_collision(x_max, y_max, d_t, particles, literal:str):
	if literal == "R":
		for p in particles:
			if p.position.x + p.velocity.x*d_t >= x_max:
				delta_t = (x_max - p.position.x)/p.velocity.x 
				p.update_position(delta_t)
				p.velocity.x *= -1
				p.update_position(d_t - delta_t)
				p.canBeMoved = False

	if literal == "L":
		for p in particles:
			if p.position.x + p.velocity.x*d_t <= 0:
				delta_t = (0 - p.position.x)/p.velocity.x 
				p.update_position(delta_t)
				p.velocity.x *= -1
				p.update_position(d_t - delta_t)
				p.canBeMoved = False

	if literal == "B":
		for p in particles:
			if p.position.y + p.velocity.y*d_t <= 0:
				delta_t = (0 - p.position.y)/p.velocity.y 
				p.update_position(delta_t)
				p.velocity.y *= -1
				p.update_position(d_t - delta_t)
				p.canBeMoved = False

	if literal == "T":
		for p in particles:
			if p.position.y + p.velocity.y*d_t >= y_max:
				delta_t = (y_max - p.position.y)/p.velocity.y 
				p.update_position(delta_t)
				p.velocity.y *= -1
				p.update_position(d_t - delta_t)
				p.canBeMoved = False

	if literal == "TL":
		for p in particles:
			if p.position.y + p.velocity.y*d_t >= y_max and p.position.x + p.velocity.x*d_t <= 0:
				if abs(p.position.y - y_max) <= abs(p.position.x - 0):
					delta_t_1 = (y_max - p.position.y)/p.velocity.y 
					p.update_position(delta_t_1)
					p.velocity.y *= -1
					delta_t_2 = (0 - p.position.x)/p.velocity.x 
					p.update_position(delta_t_2)
					p.velocity.x *= -1
					p.update_position(d_t - delta_t_1 - delta_t_2)
					p.canBeMoved = False
				else:
					delta_t_1 = (0 - p.position.x)/p.velocity.x 
					p.update_position(delta_t_1)
					p.velocity.x *= -1
					delta_t_2 = (y_max - p.position.y)/p.velocity.y 
					p.update_position(delta_t_2)
					p.velocity.y *= -1
					p.update_position(d_t - delta_t_1 - delta_t_2)
					p.canBeMoved = False
			elif p.position.y + p.velocity.y*d_t >= y_max:
				solve_wall_collision(x_max, y_max, d_t, particles, "T")
			elif p.position.x + p.velocity.x*d_t <= 0:
				solve_wall_collision(x_max, y_max, d_t, particles, "L")

	if literal == "BR":
		for p in particles:
			if p.position.y + p.velocity.y*d_t <= 0 and p.position.x + p.velocity.x*d_t >= x_max:
				if abs(0 - p.position.y) <= abs(x_max - p.position.x):
					delta_t_1 = (0 - p.position.y)/p.velocity.y 
					p.update_position(delta_t_1)
					p.velocity.y *= -1
					delta_t_2 = (x_max - p.position.x)/p.velocity.x 
					p.update_position(delta_t_2)
					p.velocity.x *= -1
					p.update_position(d_t - delta_t_1 - delta_t_2)
					p.canBeMoved = False
				else:
					delta_t_1 = (x_max - p.position.x)/p.velocity.x 
					p.update_position(delta_t_1)
					p.velocity.x *= -1
					delta_t_2 = (0 - p.position.y)/p.velocity.y 
					p.update_position(delta_t_2)
					p.velocity.y *= -1
					p.update_position(d_t - delta_t_1 - delta_t_2)
					p.canBeMoved = False
			elif p.position.y + p.velocity.y*d_t <= 0:
				solve_wall_collision(x_max, y_max, d_t, particles, "B")
			elif p.position.x + p.velocity.x*d_t >= x_max:
				solve_wall_collision(x_max, y_max, d_t, particles, "R")


	if literal == "TR":
		for p in particles:
			if p.position.y + p.velocity.y*d_t >= y_max and p.position.x + p.velocity.x*d_t >= x_max:
				if abs(p.position.y - y_max) <= abs(x_max - p.position.x):
					delta_t_1 = (y_max - p.position.y)/p.velocity.y 
					p.update_position(delta_t_1)
					p.velocity.y *= -1
					delta_t_2 = (x_max - p.position.x)/p.velocity.x 
					p.update_position(delta_t_2)
					p.velocity.x *= -1
					p.update_position(d_t - delta_t_1 - delta_t_2)
					p.canBeMoved = False
				else:
					delta_t_1 = (x_max - p.position.x)/p.velocity.x 
					p.update_position(delta_t_1)
					p.velocity.x *= -1
					delta_t_2 = (y_max - p.position.y)/p.velocity.y 
					p.update_position(delta_t_2)
					p.velocity.y *= -1
					p.update_position(d_t - delta_t_1 - delta_t_2)
					p.canBeMoved = False
			elif p.position.y + p.velocity.y*d_t >= y_max:
				solve_wall_collision(x_max, y_max, d_t, particles, "T")
			elif p.position.x + p.velocity.x*d_t >= x_max:
				solve_wall_collision(x_max, y_max, d_t, particles, "R")
	if literal == "BL":
		for p in particles:
			if p.position.y + p.velocity.y*d_t < 0 and p.position.x + p.velocity.x*d_t < 0:
				if abs(0 - p.position.y) <= abs(p.position.x - 0):
					delta_t_1 = (0 - p.position.y)/p.velocity.y 
					p.update_position(delta_t_1)
					p.velocity.y *= -1
					delta_t_2 = (0 - p.position.x)/p.velocity.x 
					p.update_position(delta_t_2)
					p.velocity.x *= -1
					p.update_position(d_t - delta_t_1 - delta_t_2)
					p.canBeMoved = False
				else:
					delta_t_1 = (0 - p.position.x)/p.velocity.x 
					p.update_position(delta_t_1)
					p.velocity.x *= -1
					delta_t_2 = (0 - p.position.y)/p.velocity.y 
					p.update_position(delta_t_2)
					p.velocity.y *= -1
					p.update_position(d_t - delta_t_1 - delta_t_2)
					p.canBeMoved = False
			elif p.position.y + p.velocity.y*d_t < 0:
				solve_wall_collision(x_max, y_max, d_t, particles, "B")
			elif p.position.x + p.velocity.x*d_t < 0: 
				solve_wall_collision(x_max, y_max, d_t, particles, "L")


#######################################################################################


def collide_particles(p_a: Particle, p_b: Particle, delta_t, t_poca):
	if p_a.canBeMoved and p_b.canBeMoved:
		v_cm = (
			(p_a.velocity.x + p_b.velocity.x)/2, 
			(p_a.velocity.y + p_b.velocity.y)/2 
		)
		v_1c = (p_a.velocity.x - v_cm[0], p_a.velocity.y - v_cm[1])
		v_2c = (p_b.velocity.x - v_cm[0], p_b.velocity.y - v_cm[1])
		p_1c = (v_1c[0], v_1c[1])
		p_2c = (v_2c[0], v_2c[1])

		phi = 2*pi*random()

		p_11 = rotate_2D_vec(phi, p_1c)
		p_21 = (-p_11[0], -p_11[1])
		v_11 = (p_11[0], p_11[1])
		v_21 = (p_21[0], p_21[1])
		
		p_a.update_position(t_poca)
		p_b.update_position(t_poca)

		p_a.velocity.x, p_a.velocity.y = v_11[0] + v_cm[0], v_11[1] + v_cm[1]
		p_b.velocity.x, p_b.velocity.y = v_21[0] + v_cm[0], v_21[1] + v_cm[1]

		p_a.update_position(delta_t - t_poca)
		p_b.update_position(delta_t - t_poca)

		E1 = 0.5*(pythagoras(p_1c[0], p_1c[1])**2 + pythagoras(p_1c[0], p_1c[1])**2 )
		E2 = 0.5*(pythagoras(p_2c[0], p_2c[1])**2 + pythagoras(p_2c[0], p_2c[1])**2 )
		E1_prim = 0.5*(pythagoras(p_11[0], p_11[1])**2 + pythagoras(p_11[0], p_11[1])**2 )
		E2_prim = 0.5*(pythagoras(p_21[0], p_21[1])**2 + pythagoras(p_21[0], p_21[1])**2 )

		assert isclose(E1+E2, E1_prim+E2_prim, rel_tol=1e-6, abs_tol=0.0), \
			f"{E1 + E2} =/= {E1_prim + E2_prim}"

		assert isclose(p_1c[0] + p_2c[0], 0.0, \
			rel_tol=1e-6, abs_tol=1e-6), f"pęd = {p_1c[0] + p_2c[0]}"

		p_a.canBeMoved = False
		p_b.canBeMoved = False


#######################################################################################


def solve_particle_collision(particles, delta_t):
	collision_candidates = []
	for p_a in particles:
		for p_b in particles:
			if p_a != p_b:
				r_rel_t = (p_a.position.x - p_b.position.x, p_a.position.y - p_b.position.y)
				r_rel_delta_t = (
					(p_a.position.x+p_a.velocity.x*delta_t) - (p_b.position.x+p_b.velocity.x*delta_t),
					(p_a.position.y+p_a.velocity.y*delta_t) - (p_b.position.y+p_b.velocity.y*delta_t)
				)
				v_rel_t = (p_a.velocity.x - p_b.velocity.x, p_a.velocity.y - p_b.velocity.y)
				v_rel_delta_t = (v_rel_t[0], v_rel_t[1])

				chi = r_rel_t[0]*v_rel_t[0] + r_rel_t[1]*v_rel_t[1]
				chi *= r_rel_delta_t[0]*v_rel_delta_t[0] + r_rel_delta_t[1]*v_rel_delta_t[1]
				if chi < 0:
					delta = (v_rel_t[0]*r_rel_t[0]+v_rel_t[1]*r_rel_t[1])**2
					delta -= (v_rel_t[0]**2 + v_rel_t[1]**2) \
						* ((r_rel_t[0]**2 + r_rel_t[1]**2) - (p_a.r_eff + p_b.r_eff)**2)
					if delta >= 0: #only then t_o_{1,2} will be real numbers
						# t_o_1 = (-(v_rel_t[0]*r_rel_t[0]+v_rel_t[1]*r_rel_t[1])**2 - sqrt(delta)) \
						# 	/ (v_rel_t[0]**2 + v_rel_t[1]**2)
						# t_o_2 = (-(v_rel_t[0]*r_rel_t[0]+v_rel_t[1]*r_rel_t[1])**2 + sqrt(delta)) \
						# 	/ (v_rel_t[0]**2 + v_rel_t[1]**2)
						
						# t_o = -1. 
						# if (t_o_1 > 0 and t_o_1 < delta_t) or (t_o_2 > 0 and t_o_2 < delta_t):
						# 	if t_o_2 >= t_o_1:
						# 		t_o = t_o_2
						# 	else:
						# 		t_o = t_o_1

						# if t_o > 0 and t_o < delta_t:
						# 	dist = pythagoras(r_rel_t[0]+v_rel_t[0]*t_o, r_rel_t[1]+v_rel_t[1]*t_o)
						# 	collision_candidates.append({'pair': (p_a.ID, p_b.ID), 'dist': dist})

						t_poca = -( r_rel_t[0]*v_rel_t[0] + r_rel_t[1]*v_rel_t[1] ) \
							/ (v_rel_t[0]**2 + v_rel_t[1]**2)
						collision_candidates.append({'pair': (p_a, p_b), 'crit': t_poca})
	if collision_candidates:
		partners = min(collision_candidates, key = lambda x:x['crit']) 
		collide_particles(partners['pair'][0], partners['pair'][1], delta_t, partners['crit']) 
	del collision_candidates


#######################################################################################

def search_and_perform_particles_collisions(particles, i, j, n_c_x, n_c_y, indices, d_t, lock):
	I_k, I_l = max(0,i-1), min(n_c_x-1, i+1)
	J_k, J_l = max(0,j-1), min(n_c_y-1, j+1)
	# p_c = [ p for p in cell[I_k:I_l+1][J_k:J_l+1] ]
	p_c = []
	for k in range(I_k, I_l+1):
		for l in range(J_k, J_l+1):
			idx = k*n_c_y + l
			if idx == 0:
				start = 0
			else:
				start = indices[idx-1]
			end = indices[idx]
			for p in particles[start:end]:
				p_c.append(p)

	# lock.acquire()
	solve_particle_collision(p_c, d_t)
	# lock.release()


def MC(particles, n_c_x: int, n_c_y: int, N_it: float, x_max: float, y_max: float, \
	speed: float, bins_number: int, mode: str):
	'''Kinetic Monte Carlo method

	Arguments:
		particles - array of generated particles
		n_c_x - divides x-axis to n_c_x pieces
		n_c_y - divides x-axis to n_c_y pieces
		N_it - number of time iterations
		x_max - numerix box x size
		y_max - numeric box y size
		n - number of particles
		speed - starting speed of every particle
		mfp_coeff - mean free path coefficient 
		bins_number - number of histogram bins
		mode - 1 of 3 available modes: BASIC, SOD, NOH
	'''
	n = len(particles)
	d_x = x_max/n_c_x
	d_y = y_max/n_c_y

	x_min, y_min = 0, 0

	thread_dx = 3
	mp_count = 4

	# r_eff = sqrt(speed / (2*pi*mfp_coeff*min(d_x, d_y)*n) )
	# print(f"{r_eff = }")
	# particles = [None for i in range(n)]
	# for i in range(n):
	# 	particles[i] = Particle(x_max, y_max, speed=speed, r_eff=r_eff)
	
	T = 0.0
	plt.rcParams.update({'font.size': 14})
	fig = plt.figure(figsize=(8,8))
	ax = fig.add_subplot()
	camera = Camera(fig)

	ax.text(0.48, 1.01, f'T = {T:.3f}', transform=ax.transAxes)
	plot_particles(particles, x_max=x_max, y_max=y_max, n_c_x=n_c_x, n_c_y=n_c_y)
	camera.snap()

	colors = ['cyan', 'red', 'lime', 'yellow', 'blue', 'orange',  
			'purple', 'green', 'pink', 'brown', 'magenta', 'gold', 'silver']

	indices = [0 for i in range(n_c_x*n_c_y)]
	cells = [None for i in range(n_c_x*n_c_y)]
	for i in range(n_c_x):
		for j in range(n_c_y):
			cells[i*n_c_y+j] = (i, j)

	# positions = [None for i in range(N_it*n)]
	for t_it in range(N_it):
		print("time iteration:", t_it+1)
		v_max = particles[0].get_speed()
		# przyporządkowanie cząstki do komórki
		# cell = [ [ [] for y in range(n_c_y) ] for x in range(n_c_x) ]
		particles = sorted(particles, 
			key = lambda p: (
				floor(p.position.x/d_x)%n_c_x, 
				floor(p.position.y/d_y)%n_c_y 
			) 
		)
		# [i_prev, j_prev] = [floor(particles[0].position.x/d_x), floor(particles[0].position.y/d_y)]
		# print([(floor(p.position.x/d_x), floor(p.position.y/d_y) ) for p in particles])

		cell_counter = 0
		for i in range(n):
			particles[i].canBeMoved = True
			particle_cell = ( 
				floor(particles[i].position.x/d_x)%n_c_x, 
				floor(particles[i].position.y/d_y)%n_c_y 
			)
			# positions[t_it*n+i] = (particles[i].position.x, particles[i].position.y)
			try:
				while particle_cell != cells[cell_counter]:
					indices[cell_counter] = i
					cell_counter += 1
			except IndexError as e:
				print(t_it, T, i, cell_counter, particle_cell)
				# print([positions[t*n+i] for t in range(t_it+1)])
				raise e
				
			# i = floor(particles[i].position.x/d_x)
			# j = floor(particles[i].position.y/d_y)
			# try:
			# 	cell[i][j].append(particles[i])
			# except IndexError as e:
			# 	print(i,j, d_t)
			# 	raise e
				
			if i == n-1:
				if particle_cell == cells[cell_counter]:
					indices[cell_counter] = i+1
					cell_counter += 1
					while cell_counter < n_c_x*n_c_y:
						indices[cell_counter] = indices[cell_counter-1]
						cell_counter += 1

			#znalezienie maksymalnej prędkości cząsteczki
			v = particles[i].get_speed()
			v_max = v if v > v_max else v_max

		#adapting time step
		d_t = min(d_x, d_y)/v_max
		T += d_t

		if mode == "NOH":
			radii = [None for i in range(n)]
			for i, p in enumerate(particles):
				radii[i] = sqrt(
					(p.position.x - (x_min + (x_max - x_min)/2) )**2 + \
					(p.position.y - (y_min + (y_max - y_min)/2) )**2
				)
			radii = sorted(radii)

			with open("density_radial.dat", "a+") as file:
				o_bins = bins_number
				dr = sqrt(
					(x_max - (x_min + (x_max - x_min)/2) )**2 + \
					(y_max - (y_min + (y_max - y_min)/2) )**2
				)/o_bins

				i = 0
				for ir in range(o_bins):
					r = dr*(ir+1)
					count = 0
					while i < n and radii[i] <= r:
						count += 1
						i += 1
					file.write(f"{r}\t{count/(pi*dr**2*(2*ir + 1) ) / 1e7}\n")
				file.write("\n\n")

		with open("density.dat", "a+") as file1, open("bulk_velocity.dat", "a+") as file2:
			for i in range(n_c_x):
				count = 0
				v_x = 0
				for j in range(n_c_y):
					idx = i*n_c_y + j
					if idx == 0:
						start = 0
					else:
						start = indices[idx-1]
					count += indices[idx] - start

					if mode == "SOD":
						for idx in range(start, indices[idx]):
							v_x += particles[idx].velocity.x
					
				file1.write(f"{i*d_x}\t{count/n_c_y}\n")
				if mode == "SOD":
					file2.write(f"{i*d_x}\t{v_x / count}\n")
			file1.write("\n\n")
			if mode == "SOD":
				file2.write("\n\n")
		
		# cell2proc = [ [] for i in range(mp_count) ]
		
		# for A in range(0, ( (n_c_x) // (thread_dx * mp_count) ) +1):
		# 	a = A*thread_dx * mp_count
		# 	for j in range(0, n_c_y, thread_dx):
		# 		for p in range(mp_count):
		# 			i = p*thread_dx + a
		# 			I_k, I_l = max(0,i) - i, min(n_c_x-1, i+2) - i
		# 			J_k, J_l = max(0,j) - j, min(n_c_y-1, j+2) - j
		# 			# print(f'{j = }, { i = }\t {I_k = }, {I_l = }\t, {J_k = }, {J_l = }')
		# 			for I in range(I_k, I_l+1):
		# 				for J in range(J_k, J_l+1):
		# 					cell2proc[p].append((i+I+1,j+J+1))

		# lock = Lock()
		# for q in range(max(map(len, cell2proc))):
		# 	threads =[]
		# 	with mp.Pool(mp.cpu_count()) as pool:
		# 		for p in range(len(cell2proc)):
		# 			try:
		# 				if cell2proc[p][q]:
		# 					i = cell2proc[p][q][0]
		# 					j = cell2proc[p][q][1]
		# 					pool.starmap_async(
		# 						search_and_perform_particles_collisions, 
		# 						(particles, i, j, n_c_x, n_c_y, indices, d_t, lock)
		# 					)
						# thread = Thread(
						# 	target=search_and_perform_particles_collisions, 
						# 	args=(particles, i, j, n_c_x, n_c_y, indices, d_t, lock)
						# )
						# thread.start()
						# threads.append(thread)
					# except IndexError:
					# 	pass
			# for thread in threads:
			# 	thread.join()
		# del cell2proc

		for i in range(n_c_x):
			threads = []
			for j in range(n_c_y):
				search_and_perform_particles_collisions(
					particles, i, j, n_c_x, n_c_y, indices, d_t, None
				)
				# thread = Thread(
				# 	target=search_and_perform_particles_collisions, 
				# 	args=(particles, i, j, n_c_x, n_c_y, indices, d_t, lock)
				# )
				# thread.start()
				# threads.append(thread)

			# for thread in threads:
			# 	thread.join()

		for i in range(n_c_x):
			for j in range(n_c_y):
				idx = i*n_c_y + j
				if idx == 0:
					start = 0
				else:
					start = indices[idx-1]
				end = indices[idx]

				#detekcja kolizji ze sciankami

				if i > 0 and i < n_c_x - 1:
					if j == 0: 
						solve_wall_collision(x_max, y_max, d_t, particles[start:end], "B")
					if j == n_c_y - 1:
						solve_wall_collision(x_max, y_max, d_t, particles[start:end], "T")

				if j > 0 and j < n_c_y - 1:
					if i == 0:
						solve_wall_collision(x_max, y_max, d_t, particles[start:end], "L")
					if i == n_c_x - 1:
						solve_wall_collision(x_max, y_max, d_t, particles[start:end], "R")
				
				if (i == 0 and j == 0):
					solve_wall_collision(x_max, y_max, d_t, particles[start:end], "BL")
				if (i == 0 and j == n_c_y - 1):
					solve_wall_collision(x_max, y_max, d_t, particles[start:end], "TL")
				if (i == n_c_x - 1 and j == 0):
					solve_wall_collision(x_max, y_max, d_t, particles[start:end], "BR")
				if (i == n_c_x - 1 and j == n_c_y - 1):
					solve_wall_collision(x_max, y_max, d_t, particles[start:end], "TR")

		#buffer against particles located out of the numeric box
		
		for p in particles:
			p.update_position(d_t)
			if p.position.x > x_max:
				p.position.x -= 2*abs(p.position.x - x_max)
			if p.position.y > y_max:
				p.position.y -= 2*abs(p.position.y - y_max)
			if p.position.x < x_min:
				p.position.x += 2*abs(p.position.x - x_min)
			if p.position.y < y_min:
				p.position.y += 2*abs(p.position.y - y_min)
		
		ax.text(0.48, 1.01, f'T = {T:.3f}', transform=ax.transAxes)
		plot_particles(particles, x_max=x_max, y_max=y_max, n_c_x=n_c_x, n_c_y=n_c_y)
		camera.snap()
		
	animation = camera.animate(interval = 100, blit=False)
	animation.save('animacja.gif', writer = 'imagemagick')

	if mode == "BASIC":
		dv = v_max/bins_number
		# hist_v = [0 for i in range(bins_number)]
		hist_v_2 = [0 for i in range(n)]

		a = 0
		for j, p in enumerate(particles):
			# i = floor(p.get_speed()/dv) -1
			# hist_v[i] +=1
			hist_v_2[j] = p.get_speed()
			a += hist_v_2[j]**2/n
		# for i in range(bins_number):
		# 	hist_v[i] /= n
		plt.close()
		fig = plt.figure(figsize=(7,7))
		plt.hist(hist_v_2, bins=bins_number, density=True)

		x = np.linspace(0, 3.5*speed, 200)
		y_a = [maxwell2d(xx, 1.0/speed**2) for xx in x]
		y_d = [maxwell2d(xx, 1.0/a) for xx in x]
		plt.plot(x, y_a, 'r-.', markersize=0.5, label='rozw. dokładne')
		# plt.plot(x, y_d, 'g-.', markersize=0.5, label='rozw. numeryczne')
		# print(1.0/a,1.0/speed**2)
		
		plt.xlim(-0.2*speed, 3.5*speed)
		plt.title(f"Histogram prędkości dla {n:_} cząstek ($v_0 = {speed}, N_{{it}} = {N_it}$)")
		plt.xlabel("$v$")
		plt.ylabel("$f(v)$")
		plt.grid(True)
		plt.legend()
		plt.savefig('velocity_histogram.png')

