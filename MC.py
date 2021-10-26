from Particle import Particle, plot_particles
from math import floor, sqrt
from random import randint
from matplotlib import pyplot as plt
import time
from celluloid import Camera

def solve_wall_collision(x_max, y_max, d_t, particles, literal:str):
	if literal == "R":
		for p in particles:
			if p.position.x + p.velocity.x*d_t > x_max:
				delta_t = (x_max - p.position.x)/p.velocity.x 
				p.update_position(delta_t)
				p.velocity.x *= -1
				p.update_position(d_t - delta_t)
				p.canBeMoved = False

	if literal == "L":
		for p in particles:
			if p.position.x + p.velocity.x*d_t < 0:
				delta_t = (0 - p.position.x)/p.velocity.x 
				p.update_position(delta_t)
				p.velocity.x *= -1
				p.update_position(d_t - delta_t)
				p.canBeMoved = False

	if literal == "B":
		for p in particles:
			if p.position.y + p.velocity.y*d_t < 0:
				delta_t = (0 - p.position.y)/p.velocity.y 
				p.update_position(delta_t)
				p.velocity.y *= -1
				p.update_position(d_t - delta_t)
				p.canBeMoved = False

	if literal == "T":
		for p in particles:
			if p.position.y + p.velocity.y*d_t > y_max:
				delta_t = (y_max - p.position.y)/p.velocity.y 
				p.update_position(delta_t)
				p.velocity.y *= -1
				p.update_position(d_t - delta_t)
				p.canBeMoved = False

	if literal == "TL":
		for p in particles:
			if p.position.y + p.velocity.y*d_t > y_max and p.position.x + p.velocity.x*d_t < 0:
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
			elif p.position.y + p.velocity.y*d_t > y_max:
				solve_wall_collision(x_max, y_max, d_t, particles, "T")
			elif p.position.x + p.velocity.x*d_t < 0:
				solve_wall_collision(x_max, y_max, d_t, particles, "L")

	if literal == "BR":
		for p in particles:
			if p.position.y + p.velocity.y*d_t < 0 and p.position.x + p.velocity.x*d_t > x_max:
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
			elif p.position.y + p.velocity.y*d_t < 0:
				solve_wall_collision(x_max, y_max, d_t, particles, "B")
			elif p.position.x + p.velocity.x*d_t > x_max:
				solve_wall_collision(x_max, y_max, d_t, particles, "R")


	if literal == "TR":
		for p in particles:
			if p.position.y + p.velocity.y*d_t > y_max and p.position.x + p.velocity.x*d_t > x_max:
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
			elif p.position.y + p.velocity.y*d_t > y_max:
				solve_wall_collision(x_max, y_max, d_t, particles, "T")
			elif p.position.x + p.velocity.x*d_t > x_max:
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

def MC(n_c_x, n_c_y, N_it, x_max, y_max, particles):
	'''Kinetic Monte Carlo method

	Arguments:
		n_c_x - divides x-axis to n_c_x pieces
		n_c_y - divides x-axis to n_c_y pieces
		N_it - number of time iterations
		x_max - numerix box x size
		y_max - numeric box y size
		particles - array of particles
	'''
	
	d_x = x_max/n_c_x
	d_y = y_max/n_c_y
	T = 0.0
	# fig = plt.figure(figsize=(8,8))
	# ax = fig.add_subplot()
	# camera = Camera(fig)

	# ax.text(0.48, 1.01, f'T = {T:.3f}', transform=ax.transAxes)
	# plot_particles(particles, x_max=x_max, y_max=y_max, n_c=n_c_x)
	# camera.snap()

	colors = ['cyan', 'red', 'lime', 'yellow', 'blue', 'orange',  
			'purple', 'green', 'pink', 'brown', 'magenta', 'gold', 'silver']

	indices = [0 for i in range(n_c_x*n_c_y)]
	cells = [None for i in range(n_c_x*n_c_y)]
	for i in range(n_c_x):
		for j in range(n_c_y):
			cells[i*n_c_x+j] = (i, j)

	positions = [None for i in range(N_it*len(particles))]
	for t_it in range(N_it):
		v_max = sqrt(particles[0].velocity.x**2 + particles[0].velocity.y**2)
		# przyporządkowanie cząstki do komórki
		# cell = [ [ [] for y in range(n_c_y) ] for x in range(n_c_x) ]
		particles = sorted(particles, key = lambda p: (floor(p.position.x/d_x), floor(p.position.y/d_y) ) )
		# [i_prev, j_prev] = [floor(particles[0].position.x/d_x), floor(particles[0].position.y/d_y)]
		
		cell_counter = 0
		for i in range(len(particles)):
			particles[i].canBeMoved = True
			particle_cell = ( floor(particles[i].position.x/d_x), floor(particles[i].position.y/d_y) )
			positions[t_it*len(particles)+i] = (particles[i].position.x, particles[i].position.y)
			try:
				while particle_cell != cells[cell_counter]:
					indices[cell_counter] = i
					cell_counter += 1
			except IndexError as e:
				print(t_it, T, i, cell_counter, (floor(particles[i].position.x/d_x), floor(particles[i].position.y/d_y)))
				print([positions[t*len(particles)+i] for t in range(t_it+1)])
				raise e
				
			# i = floor(particles[i].position.x/d_x)
			# j = floor(particles[i].position.y/d_y)
			# try:
			# 	cell[i][j].append(particles[i])
			# except IndexError as e:
			# 	print(i,j, d_t)
			# 	raise e
				
			if i == len(particles)-1:
				if particle_cell == cells[cell_counter]:
					indices[cell_counter] = i+1
					cell_counter += 1
					while cell_counter < n_c_x*n_c_y:
						indices[cell_counter] = indices[cell_counter-1]
						cell_counter += 1

			#znalezienie maksymalnej prędkości cząsteczki
			v = sqrt(particles[i].velocity.x**2 + particles[i].velocity.y**2)
			v_max = v if v > v_max else v_max

		#adaptacja kroku czasowego
		d_t = 0.3*min(d_x, d_y)/v_max
		T += d_t
		# print(d_t)
		# sum = 0
		
		for i in range(n_c_x):
			for j in range(n_c_y):
				# I_p, I_k = max(0,i-1), min(n_c_x-1, i+1)
				# J_p, J_k = max(0,j-1), min(n_c_y-1, j+1)
				# p_c = [ p for p in cell[I_p:I_k+1][J_p:J_k+1] ]
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

				# l = len(cell[i][j])
				# print(i, j, l)
				# sum += l
				# pos_X = [p.position.x for p in cell[i][j]]
				# pos_Y = [p.position.y for p in cell[i][j]]
				# plt.plot(pos_X, pos_Y, '.', markersize=2, color=colors[(i*n_c_x+j)%13])
		
		for p in particles:
			p.update_position(d_t)
		
	# 	ax.text(0.48, 1.01, f'T = {T:.3f}', transform=ax.transAxes)
	# 	plot_particles(particles, x_max=x_max, y_max=y_max, n_c=n_c_x)
	# 	camera.snap()
		
	# animation = camera.animate(interval = 500, blit=False)
	# animation.save('animacja.gif', writer = 'imagemagick')