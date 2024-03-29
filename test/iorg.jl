using Test


# Find clusters

sometrues = Bool[
    1 1 1 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 1 0 0
    0 0 0 1 0 0 0 1 0 0
    1 0 0 0 0 0 0 0 0 1
    0 1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 1 0
    0 1 0 0 0 0 0 1 0 0
]

# This matrix looks like this on heatmap so if I am counting correctly, the clustering function should find 7 when not considering boundary conditions and 5 when considering periodic boundary.

# | 1  1  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅ |
# | ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅ |
# | ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅ |
# | ⋅  ⋅  ⋅  1  ⋅  ⋅  ⋅  1  ⋅  ⋅ |
# | 1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1 |
# | ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅ |
# | ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅ |
# | ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅ |
# | ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅ |
# | ⋅  1  ⋅  ⋅  ⋅  ⋅  ⋅  1  ⋅  ⋅ |

# TODO
#@test length(detect_updrafts(sometrues; periodic_boundaries = true)) == 5

clusters = detect_updraft_clusters(sometrues)
centroids =
    [(1.0, 2.0), (5.5, 1.5), (3.5, 8.0), (4.0, 4.0), (5.0, 10.0), (10.0, 2.0), (9.5, 8.5)]

@test length(detect_updraft_clusters(sometrues)) == 7
@test Set(centroids) == Set([
    find_cluster_centroid(
        cluster,
        range(1.0, stop = 10, length = 10),
        range(1.0, stop = 10, length = 10),
    ) for cluster in clusters
])

a = zeros(Bool, 100, 100)
a[1:2:end, 1:2:end] .= true

#### Test iorg
# Regularly spaced convection should be < 0.5 
# Random should be around 0.5
# Clustered convection should be > 0.5

uniform_convec = zeros(Bool, 500, 500, 100)
uniform_convec[1:100:end, 1:100:end, :] .= true

## rand(Nx,Ny) .> 0.5 will have around 0.5*Nx*Ny positives. rand(Nx,Ny) .> 0.10 will have around 0.9*Nx*Ny so (1.0 - threshold)*Nx*Ny 
random_scarce = (rand(500, 500, 10) .> 0.99) # rand is in (0,1)
random_abundant = (rand(500, 500, 10) .> 0.9); # rand is in (0,1)
clustered_convec = zeros(Bool, 500, 500, 10)
clustered_convec[1:10, 1:10, :] .= true
clustered_convec[20:30, 20:30, :] .= true
clustered_convec[90:100, 90:100, :] .= true
clustered_convec[1:3:500, 400:3:500, :] .= true

Lx = Ly = 5e5
x_coords = y_coords = range(0.0, stop = Lx, length = 500)

distances, r, x, y, iorg_uniform = compute_iorg(uniform_convec, x_coords, y_coords, Lx, Ly)
distances, r, x, y, iorg_clustered =
    compute_iorg(clustered_convec, x_coords, y_coords, Lx, Ly)
distances, r, x, y, iorg_random_scarce =
    compute_iorg(random_scarce, x_coords, y_coords, Lx, Ly)
distances, r, x, y, iorg_abundant =
    compute_iorg(random_abundant, x_coords, y_coords, Lx, Ly)

distances, r, x, y, iorg_uniform_2d =
    compute_iorg(uniform_convec[:, :, end], x_coords, y_coords, Lx, Ly)
distances, r, x, y, iorg_clustered_2d =
    compute_iorg(clustered_convec[:, :, end], x_coords, y_coords, Lx, Ly)
distances, r, x, y, iorg_random_scarce_2d =
    compute_iorg(random_scarce[:, :, end], x_coords, y_coords, Lx, Ly)
distances, r, x, y, iorg_abundant_2d =
    compute_iorg(random_abundant[:, :, end], x_coords, y_coords, Lx, Ly)

@test iorg_uniform < 0.5
@test iorg_clustered > 0
@test isapprox(iorg_random_scarce, 0.5, rtol = 0.15)

@test iorg_uniform_2d < 0.5
@test iorg_clustered_2d > 0
@test isapprox(iorg_random_scarce_2d, 0.5, rtol = 0.15)
