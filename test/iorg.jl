using Test

a = zeros(Bool,100,100)
a[1:2:end, 1:2:end] .= true

#@test all(compute_distances(a, 10, 10) .≈ 10)
#@test all(compute_distances(a, 1.0, 1.0) .≈ 1.0)

uniform_convec = zeros(Bool,500,500,100)
uniform_convec[1:100:end,1:100:end,:] .= true

## rand(Nx,Ny) .> 0.5 will have around 0.5*Nx*Ny positives. rand(Nx,Ny) .> 0.10 will have around 0.9*Nx*Ny so (1.0 - threshold)*Nx*Ny 
random_scarce = (rand(500,500,10) .> 0.99) # rand is in (0,1)
random_abundant =  (rand(500,500,10) .> 0.9); # rand is in (0,1)
# clustered_convec = zeros(Bool,500,500,100)
# clustered_convec[1:10,1:10, rand(1:100)] .= true
# clustered_convec[2:30,20:30,rand(1:100)] .= true
# clustered_convec[90:100,90:100,rand(1:100)] .= true

Lx = Ly = 5e5
x_coords = y_coords = range(0.0, stop = Lx, length = 500)
distances,r,x,y,iorg = compute_iorg(uniform_convec,x_coords,y_coords,Lx,Ly)
@test iorg < 0.5
#lineplot(y,x)
#distances,r,x,y,iorg = compute_iorg(clustered_convec,x_coords,y_coords,Lx,Ly)
#lineplot(y,x)
distances,r,x,y,iorg = compute_iorg(random_scarce,x_coords,y_coords,Lx,Ly)
@test isapprox(iorg, 0.5, rtol = 0.15)
# lineplot(y,x)
#distances,r,x,y,iorg = compute_iorg(random_abundant,x_coords,y_coords,Lx,Ly)
# lineplot(y,x)


#TODO it will not work without identifying individual clusters

# lineplot(y,x)


# #Find clusters

# sometrues = Bool[1 1 1 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 0 0; 0 0 0 1 0 0 0 1 0 0; 1 0 0 0 0 0 0 0 0 1; 0 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 1 0; 0 1 0 0 0 0 0 1 0 0]

# # This matrix looks like this on heatmap so if I am counting correctly, the clustering function should find 7 when not considering boundary conditions and 5 when considering periodic boundary.

# # | 1  1  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅ |
# # | ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅ |
# # | ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅ |
# # | ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  1  ⋅  ⋅ |
# # | 1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1 |
# # | ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅ |
# # | ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅ |
# # | ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅ |
# # | ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅ |
# # | ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅ |

# @test length(detect_updraft_clusters(sometrues)) == 7
# #@test length(detect_updrafts(sometrues; periodic_boundaries = true)) == 5

# clusters = detect_updraft_clusters(sometrues)

# find_cluster_centroid(clusters[1],1:10,1:10)

# centroids = [(1.0, 2.0), (5.5, 1.5), (3.5, 8.0), (4.0, 4.0), (5.0,10.0), (10.0, 2.0), (9.5, 8.5) ]

# @test Set(centroids) == Set([find_cluster_centroid(cluster,range(1.0,stop = 10, length = 10), range(1.0,stop = 10, length = 10)) for cluster in clusters]) 