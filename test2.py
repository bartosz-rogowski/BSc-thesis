from matplotlib import pyplot as plt
from math import ceil

colors = ['cyan', 'red', 'lime', 'yellow', 'blue', 'orange',  
	'purple', 'green', 'pink', 'brown', 'magenta', 'gold', 'silver']
n_c_x = 19
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

thread_dx = 3
mp_count = 8


# for a in range(0, ceil(nx/mp_count)+1, dx):
# 	for b in range(0, ceil(ny/mp_count)+1, dx):
# 		for p in range(mp_count):
# 			if p <= nx:
# 				print(a+p, end=' ')
# 		print()

cell2proc = [ [] for i in range(mp_count) ]

for A in range(0, ( (n_c_x) // (thread_dx * mp_count) ) +1):
	a = A*thread_dx * mp_count
	# print(f"{a=}\t{A=}")
	for j in range(0, n_c_y, thread_dx):
		for p in range(mp_count):
			i = p*thread_dx + a
			I_k, I_l = max(0,i) - i, min(n_c_x-1, i+2) - i
			J_k, J_l = max(0,j) - j, min(n_c_y-1, j+2) - j
			# print(f'{j = }, { i = }\t {I_k = }, {I_l = }\t, {J_k = }, {J_l = }')
			for I in range(I_k, I_l+1):
				for J in range(J_k, J_l+1):
					# print(f"{p=}\t{i+I + (j+J)*n_c_x}")
					ax.annotate(
						# f'{i+I + (j+J)*n_c_x}, {p}', 
						f'{i+I+1},{j+J+1}',
						xy=( dx*(i+I+1) - dx/1.25, dy*(j+J+1) - dy/1.85), 
						textcoords='data',
						color = colors[p],
						size = 8
					)
					cell2proc[p].append((i+I+1,j+J+1))
		
for q in range(max(map(len, cell2proc))):
	for p in range(len(cell2proc)):
		print(p+1, end=':\t')
		try:
			if cell2proc[p][q]:
				print(cell2proc[p][q], end=', ')
		except IndexError:
			pass
		print()
plt.show()