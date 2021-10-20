
nc = 4

for i in range(0, nc):
    for j in range(0, nc):
        print(f'({i},{j}) -- I = {max(0, i-1)}:{min(nc-1, i+1)}\tJ = {max(0, j-1)}:{min(nc-1, j+1)}')