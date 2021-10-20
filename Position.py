from random import random

class Position:
	def __init__(self, x_max = 1, y_max = 1):
		self.x = x_max*random()
		self.y = y_max*random()
	
	def __add__(self, o):
		return self.x + o.x, self.y + o.y