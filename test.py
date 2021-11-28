from matplotlib import pyplot as plt
from math import ceil

n_c_x = 16
n_c_y = 10

x_max = 4.0
y_max = 4.0
x_min, y_min = 0., 0.

dx = (x_max - x_min)/n_c_x
dy = (y_max - y_min)/n_c_y

for n in range(n_c_y+1):
	plt.plot([x_min, x_max], [n*dy, n*dy], 'k--', linewidth=1)
for n in range(n_c_x+1):
	plt.plot([n*dx, n*dx], [y_min, y_max], 'k--', linewidth=1)

ax = plt.gca()

# for i in range(n_c_x):
# 	for j in range(n_c_y):
# 		ax.annotate(
# 			f'{i*n_c_y+j}', 
# 			xy=( dx*(i+1) -dx/2, dy*(j+1) - dy/2), 
# 			textcoords='data',
# 			color = (0., 0., 255/255, 255/255)
# 		)

threads = 4
thread_dx = 3
color = [
	'blue',
	'red',
	'lime',
	'magenta',
	'gold',
	'black'
]

t = 0
# while(t*thread_dx*n_c_y < )
ivy = 0
for j in range(0, n_c_y, thread_dx):
	t = 0
	for i in range(0, n_c_x, thread_dx):
		I_k, I_l = max(0,i) - i, min(n_c_x-1, i+2) - i
		J_k, J_l = max(0,j) - j, min(n_c_y-1, j+2) - j
		# print(f'{j = }, { i = }\t {I_k = }, {I_l = }\t, {J_k = }, {J_l = }')
		for I in range(I_k, I_l+1):
			for J in range(J_k, J_l+1):
				ax.annotate(
					f'{( (i+I)*(J_l+1) + (J) ) % (n_c_x*thread_dx)  + ivy*n_c_x*thread_dx}', 
					xy=( dx*(i+I+1) - dx/2, dy*(j+J+1) - dy/2), 
					textcoords='data',
					color = color[t%threads],
					size = 8
				)
		t += 1
	ivy += 1



plt.show()