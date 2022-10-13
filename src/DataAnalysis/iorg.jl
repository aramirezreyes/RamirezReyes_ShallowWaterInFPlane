"""
    detect_updraft_clusters(array)

Receive a 2-d array of Bool. Each `true` represents an updraft. Return an array of arrays of `CartesianIndex`. Each inner array contains the locations of adjacent updrafts, representing an updraft cluster. It uses ImageMorphology.label_components to find clusters.
"""
function detect_updraft_clusters(array)
    connected = label_components(array, ones(Bool,3,3))
    nclusters = maximum(connected)
    cluster_labels = 1:nclusters
    clustered = [findall(==(cluster), connected) for cluster in cluster_labels]
    return clustered,nclusters
end


"""
    find_cluster_centroid(cluster, x,y)

Receive an iterable of CartesianIndex{2} and vectors of coordinates. Return a tuple (xc,yc) with the coordinates of the centroid of the group of coordinates defined as the average over each direction.
"""
function find_cluster_centroid(cluster, x,y)
    npoints = length(cluster)
    x = sum(point -> x[point[1]],cluster)
    y = sum(point -> y[point[2]],cluster)
    return x/npoints,y/npoints
end

"""
    compute_distances(centroids)

Receive an array of 2-tuple with coordinatess of points. Return an array with the distance of each point to it's closes neighbohr.   
"""
function compute_distances(centroids,r)
    distances = zeros(Int,r.len)
    compute_distances!(distances,centroids,r)
end

"""
    compute_distances!(distances,centroids,r)

Receive a vector of bins (distances), an array of 2-tuples with coordinates of points and a `range` of possible nearest-neighbor distances. Counts how many occurrences of each distance appear.
"""
function compute_distances!(distances,centroids,r)
    dr = r[2] - r[1]
    for point in centroids
        distance = Inf
        for neighbor in centroids
            point == neighbor && continue
            newdistance = compute_distance(point,neighbor)
            distance = newdistance < distance ? newdistance : distance
        end
        distance == Inf && continue
        distance_index = ceil(Int,distance/dr)
        distance_index > r.len && continue
        distances[distance_index] += 1
    end
    return distances
end

"""
    compute_iorg(array::Array{Bool,3},x,y,Lx,Ly, updraft_density = count(array[:,:,end])/(Lx*Ly))

Receive a 3-d array of booleans on the form (x,y,t), vectors for the x and y coordinates, and Lx, Ly the length of the domain in each direction.
Returns the organization index, iorg, as defined by Tompkins and Semie (2017). The assumption is that each `true` value in the array represents an updraft obtained elswhere by the user.

If updraft_density is not passed, this assumes that the simulation is in steady state and that the point density of the last frame is representative of the simulation.

Tompkins, A. M., & Semie, A. G. (2017). Organization of tropical convection in low vertical wind shears: Role of updraft entrainment. Journal of Advances in Modeling Earth Systems, 9(2), 1046–1068. https://doi.org/10.1002/2016MS000802

"""
function compute_iorg(array::AbstractArray{Bool,3},x,y,Lx,Ly, updraft_density = 0.0) 
    ## Next r_max comes from trying to estimate a radius to which the NNCDF of poisson-distributed points comes close to 1.0 for a given mean density
    if iszero(updraft_density)
        updraft_density = detect_updraft_clusters(array[:,:,end])[2] / (Lx*Ly)
    end
    rmax = sqrt(-1 * log(0.001) /(updraft_density * π))
    r = range(0, stop = rmax,length = 100)
    nncdf_poisson = nncdf_poisson_func.(updraft_density,r)
    distances_count = zeros(Int,r.len)
    nclusters = 0
    for arr_2d in eachslice(array,dims=3)
        clusters,nclusters_local = detect_updraft_clusters(arr_2d)
        if nclusters_local <= 1
            continue
        end
        centroids = map(cluster -> find_cluster_centroid(cluster,x, y), clusters)
        compute_distances!(distances_count,centroids,r)
        nclusters = nclusters + nclusters_local
    end
    if nclusters <= 1
        @info "Array needs to have at least 2 convective events in at least one frame! Returning NaN"
        return distances_count, r, distances_count, nncdf_poisson, NaN
    end
    #npoints unit area
    nncdf_data_func = ecdf(r,weights = distances_count)
    nncdf_data = if iszero(distances_count)
        distances_count
    else
        nncdf_data = nncdf_data_func.(r)
    end
    iorg = 0.5*sum((nncdf_data[1:end-1] .+ nncdf_data[2:end]) .* (nncdf_poisson[2:end] .- nncdf_poisson[1:end-1]))
    return distances_count,r,nncdf_data, nncdf_poisson, iorg
end


function compute_iorg(array::AbstractArray{Bool,2},x,y,Lx,Ly, updraft_density = 0.0) 
    ## Next r_max comes from trying to estimate a radius to which the NNCDF of poisson-distributed points comes close to 1.0 for a given mean density 
    clusters,nclusters = detect_updraft_clusters(array)
    if iszero(updraft_density)
        updraft_density = nclusters / (Lx*Ly)
    end
    rmax = sqrt(-1 * log(0.001) /(updraft_density * π))
    r = range(0, stop = rmax,length = 100)
    nncdf_poisson = nncdf_poisson_func.(updraft_density,r)
    distances_count = zeros(Int,r.len)
    if nclusters <= 1
        @info "Array needs to have at least 2 convective events! Returning NaN"
        return distances_count, r, distances_count, nncdf_poisson, NaN
    end
    centroids = map(cluster -> find_cluster_centroid(cluster,x, y), clusters)
    compute_distances!(distances_count,centroids,r)
    #npoints unit area
    nncdf_data_func = ecdf(r,weights = distances_count)
    nncdf_data = if iszero(distances_count)
        distances_count
    else
        nncdf_data = nncdf_data_func.(r)
    end
    iorg = 0.5*sum((nncdf_data[1:end-1] .+ nncdf_data[2:end]) .* (nncdf_poisson[2:end] .- nncdf_poisson[1:end-1]))
    return distances_count,r,nncdf_data, nncdf_poisson, iorg
end


"""
    compute_distance(a,b)

Receives two 2-tuples of coordinates (x,y) and computes the distance between both points.
"""
@inline function compute_distance(a,b)
    hypot((b[1] - a[1]),(b[2] - a[2]))
end

@inline """
    nncdf_poisson_func(λ,r)

This function receives the point density (npoints per unit area) and a distance, and produces the Weibull distribution. This function represents the cumulative distribution function of pairwise distances between a point and its nearest neighbor if the location of each point behave as Poisson processes (Tompkins and Semie, 2017). 

Tompkins, A. M., & Semie, A. G. (2017). Organization of tropical convection in low vertical wind shears: Role of updraft entrainment. Journal of Advances in Modeling Earth Systems, 9(2), 1046–1068. https://doi.org/10.1002/2016MS000802
"""
function nncdf_poisson_func(λ,r)
    1 - exp(-λ*π*r^2)
end