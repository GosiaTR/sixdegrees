import sixdegrees

N = 100
k = 10
p = 0.2

N_new, edges = sixdegrees.original_small_world_network(N,k,p,use_giant_component = True)

print(N_new, edges)
