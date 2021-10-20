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
	fig = plt.figure(figsize=(8,8))
	ax = fig.add_subplot()
	camera = Camera(fig)

	ax.text(0.48, 1.01, f'T = {T:.3f}', transform=ax.transAxes)
	plot_particles(particles, x_max=x_max, y_max=y_max, n_c=n_c_x)
	camera.snap()

	colors = ['cyan', 'red', 'lime', 'yellow', 'blue', 'orange',  
			'purple', 'green', 'pink', 'brown', 'magenta', 'gold', 'silver']

	for t_it in range(N_it):
		v_max = sqrt(particles[0].velocity.x**2 + particles[0].velocity.y**2)
		# przyporządkowanie cząstki do komórki
		cell = [ [ [] for y in range(n_c_y) ] for x in range(n_c_x) ]
		for p in particles:
			p.canBeMoved = True
			i = floor(p.position.x/d_x)
			j = floor(p.position.y/d_y)
			try:
				cell[i][j].append(p)
			except IndexError as e:
				print(i,j, d_t)
				raise e
			#znalezienie maksymalnej prędkości cząsteczki
			v = sqrt(p.velocity.x**2 + p.velocity.y**2)
			v_max = v if v > v_max else v_max

		#adaptacja kroku czasowego
		d_t = min(d_x, d_y)/v_max
		T += d_t
		# print(d_t)
		# sum = 0
		
		for i in range(n_c_x):
			for j in range(n_c_y):
				I_p, I_k = max(0,i-1), min(n_c_x-1, i+1)
				J_p, J_k = max(0,j-1), min(n_c_y-1, j+1)
				p_c = [ p for p in cell[I_p:I_k+1][J_p:J_k+1] ]

				#detekcja kolizji ze sciankami

				if i > 0 and i < n_c_x - 1:
					if j == 0: 
						solve_wall_collision(x_max, y_max, d_t, cell[i][j], "B")
					if j == n_c_y - 1:
						solve_wall_collision(x_max, y_max, d_t, cell[i][j], "T")

				if j > 0 and j < n_c_y - 1:
					if i == 0:
						solve_wall_collision(x_max, y_max, d_t, cell[i][j], "L")
					if i == n_c_x - 1:
						solve_wall_collision(x_max, y_max, d_t, cell[i][j], "R")
				
				if (i == 0 and j == 0):
					solve_wall_collision(x_max, y_max, d_t, cell[i][j], "BL")
				if (i == 0 and j == n_c_y - 1):
					solve_wall_collision(x_max, y_max, d_t, cell[i][j], "TL")
				if (i == n_c_x - 1 and j == 0):
					solve_wall_collision(x_max, y_max, d_t, cell[i][j], "BR")
				if (i == n_c_x - 1 and j == n_c_y - 1):
					solve_wall_collision(x_max, y_max, d_t, cell[i][j], "TR")

				# l = len(cell[i][j])
				# print(i, j, l)
				# sum += l
				# pos_X = [p.position.x for p in cell[i][j]]
				# pos_Y = [p.position.y for p in cell[i][j]]
				# plt.plot(pos_X, pos_Y, '.', markersize=2, color=colors[(i*n_c_x+j)%13])
		
		for p in particles:
			p.update_position(d_t)
		# print(f"Sum: {sum}")
		ax.text(0.48, 1.01, f'T = {T:.3f}', transform=ax.transAxes)
		plot_particles(particles, x_max=x_max, y_max=y_max, n_c=n_c_x)
		camera.snap()
		# print(f'{d_t:.4f}')
	animation = camera.animate(interval = 500, blit=False)
	animation.save('animacja.gif', writer = 'imagemagick')