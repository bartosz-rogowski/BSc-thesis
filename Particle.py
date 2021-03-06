from random import random
from Position import Position
from Velocity import Velocity
from matplotlib import pyplot as plt
from math import sqrt
from tools import pythagoras
from numpy.random import normal

class Particle:
	id = 0
	def __init__(self, x_min, y_min, x_max, y_max, speed, r_eff, mode):
		self.position = Position(x_min, y_min, x_max, y_max)
		if mode == "BASIC":
			self.velocity = Velocity(speed)

		elif mode == "SOD":
			self.velocity = Velocity(speed/sqrt(2), x=normal(), y=normal())

		elif mode == "NOH":
			self.position.x /= 2
			self.position.y /= 2
			vel = (-self.position.x, -self.position.y)
			norm = pythagoras(vel[0], vel[1])
			self.velocity = Velocity(speed, x=vel[0]/norm, y=vel[1]/norm)

		self.canBeMoved = True
		self.r_eff = r_eff
		self.ID = Particle.id
		Particle.id+=1

	def print(self):
		print(f"Pos({self.position.x:.3f}, {self.position.y:.3f}), \
		Vel({self.velocity.x:.4f}, {self.velocity.y:.4f}), = {sqrt(self.velocity.x**2 + self.velocity.y**2)}")

	def get_speed(self):
		return self.velocity.get_speed()

	def update_position(self, d_t, velocity = None):
		if self.canBeMoved:
			if not velocity:
				velocity = (self.velocity.x, self.velocity.y)
			self.position.x += velocity[0]*d_t
			self.position.y += velocity[1]*d_t

def plot_particles(particles, x_max, y_max, n_c_x, n_c_y):
	pos_X = [p.position.x for p in particles]
	pos_Y = [p.position.y for p in particles]
	vel_X = [p.velocity.x for p in particles]
	vel_Y = [p.velocity.y for p in particles]

	# plt.quiver(pos_X, pos_Y, vel_X, vel_Y, scale=300, width=0.001)
	plt.plot(pos_X, pos_Y, 'b.', markersize=2)
	plt.xlim(-0.2, x_max+0.2)
	plt.ylim(-0.2, y_max+0.2)
	# plt.xlim(0, x_max)
	# plt.ylim(0, y_max)

	ax = plt.gca()
	# for p in particles:                                       
	# 	ax.annotate(f'{p.ID}', xy=(p.position.x, p.position.y), textcoords='data')

	cell_size = (x_max/n_c_x, y_max/n_c_y)
	# for n in range(n_c_y+1):
	# 	plt.plot([0, x_max], [n*cell_size[1], n*cell_size[1]], 'k--', linewidth=1)
	# for n in range(n_c_x+1):
	# 	plt.plot([n*cell_size[0], n*cell_size[0]], [0, y_max], 'k--', linewidth=1)
	# plt.show()