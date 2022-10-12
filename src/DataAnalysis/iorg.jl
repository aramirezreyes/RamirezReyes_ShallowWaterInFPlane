"""
    detect_updraft_clusters(array)

Receive a 2-d array of Bool. Each `true` represents an updraft. Return an array of arrays of `CartesianIndex`. Each inner array contains the locations of adjacent updrafts, representing an updraft cluster. It uses ImageMorphology.label_components to find clusters.
"""
function detect_updraft_clusters(array)
    connected = label_components(array, ones(Bool,3,3))
    nclusters = 1:maximum(connected)
    clustered = [findall(==(cluster), connected) for cluster in nclusters]
    return clustered
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
function compute_distances(centroids)
    distances = Float64[]
    compute_distances!(distances,centroids)
end

"""
    compute_distances!(distances,centroids)

Receive a vector of distances and an array of 2-tuples with coordinates of points. Adds to the distances array the distance of each point to it's closes neighbohr. 
"""
function compute_distances!(distances,centroids)
    for point in centroids
        distance = Inf
        for neighbor in centroids
            point == neighbor && continue
            newdistance = compute_distance(point,neighbor)
            distance = newdistance < distance ? newdistance : distance
        end
        push!(distances,distance)
    end
    return distances
end

"""
    compute_iorg(array::Array{Bool,3},x,y,Lx,Ly, updraft_density = count(array[:,:,end])/(Lx*Ly))

Receive a 3-d array of booleans on the form (x,y,t), vectors for the x and y coordinates, and Lx, Ly the length of the domain in each direction.
Returns the organization index, iorg, as defined by Tompkins and Semie (2017). The assumption is that each `true` value in the array represents an updraft obtained elswhere by the user.

Tompkins, A. M., & Semie, A. G. (2017). Organization of tropical convection in low vertical wind shears: Role of updraft entrainment. Journal of Advances in Modeling Earth Systems, 9(2), 1046–1068. https://doi.org/10.1002/2016MS000802

"""
function compute_iorg(array::AbstractArray{Bool,3},x,y,Lx,Ly, updraft_density = count(array[:,:,end])/(Lx*Ly)) 
    ## Next r_max comes from trying to estimate a radius to which the NNCDF of poisson-distributed points comes close to 1.0 for a given mean density 
    rmax = sqrt(-1 * log(0.001) /(updraft_density * π))
    r = range(0, stop = rmax,length = 100)
    #dr = r[2] - r[1]
    nncdf_poisson = nncdf_poisson_func.(updraft_density,r)
    distances = Float64[]
    for arr_2d in eachslice(array,dims=3)
        clusters = detect_updraft_clusters(arr_2d)
        centroids = map(cluster -> find_cluster_centroid(cluster,x, y), clusters)
        compute_distances!(distances,centroids)
    end
    nncfd_data_func = ecdf(distances)
    #npoints unit area
    nncdf_data = nncfd_data_func.(r)
    iorg = 0.5*sum((nncdf_data[1:end-1] .+ nncdf_data[2:end]) .* (nncdf_poisson[2:end] .- nncdf_poisson[1:end-1]))
    return distances,r,nncdf_data, nncdf_poisson, iorg
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