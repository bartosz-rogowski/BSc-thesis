from random import random

class Position:
	def __init__(self, x_min = 0, y_min = 0, x_max = 1, y_max = 1):
		self.x = (x_max-x_min)*random() + x_min
		self.y = (y_max-y_min)*random() + y_min
	
	def __add__(self, o):
		return self.x + o.x, self.y + o.y