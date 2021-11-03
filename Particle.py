from random import random
from Position import Position
from Velocity import Velocity
from matplotlib import pyplot as plt
from math import sqrt



class Particle:
	id = 0
	def __init__(self, x_max=1, y_max=1, speed=10):
		self.position = Position(x_max, y_max)
		self.velocity = Velocity(speed)
		self.canBeMoved = True
		self.r_eff = 0.01
		self.ID = Particle.id
		Particle.id+=1

	def print(self):
		print(f"Pos({self.position.x:.3f}, {self.position.y:.3f}), \
		Vel({self.velocity.x:.4f}, {self.velocity.y:.4f}), = {sqrt(self.velocity.x**2 + self.velocity.y**2)}")

	def get_speed(self):
		return self.velocity.get_speed()

	def update_position(self, d_t):
		if self.canBeMoved:
			self.position.x += self.velocity.x*d_t
			self.position.y += self.velocity.y*d_t

def plot_particles(particles, x_max, y_max, n_c=1):
	pos_X = [p.position.x for p in particles]
	pos_Y = [p.position.y for p in particles]
	vel_X = [p.velocity.x for p in particles]
	vel_Y = [p.velocity.y for p in particles]

	# plt.quiver(pos_X, pos_Y, vel_X, vel_Y, scale=300, width=0.001)
	plt.plot(pos_X, pos_Y, 'b.', markersize=2)
	# plt.xlim(-0.2, x_max+0.2)
	# plt.ylim(-0.2, y_max+0.2)
	plt.xlim(-1, 1+x_max)
	plt.ylim(-1, 1+y_max)

	ax = plt.gca()
	# for p in particles:                                       
	# 	ax.annotate(f'{p.ID}', xy=(p.position.x, p.position.y), textcoords='data')

	cell_size = (x_max/n_c, y_max/n_c)
	for n in range(n_c+1):
		plt.plot([0, x_max], [n*cell_size[1], n*cell_size[1]], 'k--', linewidth=1)
		plt.plot([n*cell_size[0], n*cell_size[0]], [0, y_max], 'k--', linewidth=1)
	# plt.show()