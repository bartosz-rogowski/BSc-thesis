
n_c_x = 3
n_c_y = 3

# for i in range(0, nc):
#     for j in range(0, nc):
#         print(f'({i},{j}) -- I = {max(0, i-1)}:{min(nc-1, i+1)}\tJ = {max(0, j-1)}:{min(nc-1, j+1)}')

cells = []

for i in range(n_c_x):
	for j in range(n_c_y):
		cells.append((i,j))
# l = [(0,0), (0,1),(0,1),(0,1), (0, 2), (0, 2), (0, 5), (1, 2), (1, 2), (1, 5), (3,3),(5,5),(5,5)]

l = [(0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 1,), (0, 1), (0, 2), (1, 2), (1, 2)]
# l = [(0, 1), (0, 3), (0, 3), (0, 4), (0, 6), (0, 6), (1, 4), (1, 6), (1, 9), (2, 0), (2, 0), (2, 4), (2, 7), (3, 0), (3, 1), (3, 2), (3, 3), (3, 6), (4, 2), (4, 3), (4, 4), (4, 4), (4, 4), (4, 6), (5, 7), (5, 7), (6, 2), (6, 3), (6, 8), (6, 9), (7, 2), (7, 5), (8, 0), (8, 0), (8, 1), (8, 2), (8, 9), (9, 0), (9, 5), (9, 8)]
i = 0
cell_counter = 0
# idx_start = [0 for i in range(n_c_x*n_c_y)]
idx_end = [0 for i in range(n_c_x*n_c_y)]
# for i in range(len(l)):
# 	idx_start[cell_counter] = i
# 	while l[i] != cells[cell_counter]:
# 		# print(l[i], cells[cell_counter])
# 		# print(f"i={i}\tcell_counter={cell_counter}")
# 		idx_end[cell_counter] = i
# 		cell_counter += 1
# 	if i == len(l)-1:
# 		while cell_counter < 36:
# 			idx_start[cell_counter] = i
# 			idx_end[cell_counter] = i
# 			cell_counter += 1

# for cell_counter in range(36):
# 	idx_start[cell_counter] = i
# 	while l[i] == cells[cell_counter]:
# 		idx_end[cell_counter] = i
# 		if i+1 < len(l):
# 			i +=1
# 		else: 
# 			break



for i in range(len(l)):
	while l[i] != cells[cell_counter]:
		idx_end[cell_counter] = i
		cell_counter += 1	
	# if l[i] == cells[cell_counter]:
	# 	count += 1
	# else:
	# 	idx_end[cell_counter] = count
	# 	count = 0
	# 	cell_counter += 1

	if i == len(l)-1:
		if l[i] == cells[cell_counter]:
			idx_end[cell_counter] = i+1
			cell_counter += 1
			while cell_counter < n_c_x*n_c_y:
				idx_end[cell_counter] = idx_end[cell_counter-1]
				cell_counter += 1

		
print(l, '\n', idx_end)
for i in range(n_c_x):
	for j in range(n_c_y):
		idx = i*n_c_y+j
		if idx == 0:
			c = idx_end[idx]
		else:
			c = idx_end[idx] - idx_end[idx-1]
		print(f'({i}, {j}): {c}')