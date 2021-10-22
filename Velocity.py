from random import random
from math import sin, cos, sqrt
from math import pi

class Velocity:
	def __init__(self, speed = 10):
		fi = 2*pi*random() #losowanie kata
		# self.x = speed*2*(random() - 0.5)
		# self.y = speed*2*(random() - 0.5)
		self.x = speed*cos(fi)
		self.y = speed*sin(fi)

	def get_speed(self):
		return sqrt(self.x**2 + self.y**2)
	
	def __add__(self, o):
		return self.x + o.x, self.y + o.y