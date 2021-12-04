from random import random
from math import sin, cos, pi
from tools import pythagoras

class Velocity:
	def __init__(self, speed = 10, x = 0, y = 0):
		if x*y == 0:
			fi = 2*pi*random() #losowanie kata
			# self.x = speed*2*(random() - 0.5)
			# self.y = speed*2*(random() - 0.5)
			self.x = speed*cos(fi)
			self.y = speed*sin(fi)
		else:
			self.x = speed*x
			self.y = speed*y

	def get_speed(self):
		return pythagoras(self.x, self.y)
	
	def __add__(self, o):
		return self.x + o.x, self.y + o.y