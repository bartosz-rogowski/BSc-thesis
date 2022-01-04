from math import sqrt, sin, cos, exp

def pythagoras(a, b):
	return sqrt(a**2 + b**2)

def rotate_2D_vec(phi, vec):
	return (vec[0]*cos(phi) - vec[1]*sin(phi), \
		vec[0]*sin(phi) + vec[1]*cos(phi))

def maxwell2d(x, a=0.50):
	return 2*a*x*exp(-a*x**2)