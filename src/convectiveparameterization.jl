"""
    update_convective_events!(architecture,isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
Use the parameters τ_convec (field of times since last convective event started) and h_threshold to determine if one point should convect.
The version that receives architecture as first parameter is an interface for the specific implementation.

Based on:
Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
"""
function update_convective_events!(architecture :: CPU,isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
    update_convective_events_cpu!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
end

"""
    update_convective_events!(architecture,isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
Use the parameters τ_convec (field of times since last convective event started) and h_threshold to determine if one point should convect.
The version that receives architecture as first parameter is an interface for the specific implementation.

Based on:
Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
"""
function update_convective_events!(architecture :: GPU,isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
    @cuda update_convective_events_gpu!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
end

"""
    update_convective_events_cpu!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
Use the parameters τ_convec (field of times since last convective event started) and h_threshold to determine if one point should convect.
The version contains the implementation for CPU.

Based on:
Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
"""
@inline function update_convective_events_cpu!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
     @inbounds Threads.@threads for ind in eachindex(h)
        time_convecting = t - convection_triggered_time[ind]
        needs_to_convect_by_time = isconvecting[ind] && (time_convecting < τ_convec) #has been convecting less than τ_c?
        needs_to_convect_by_height = (h[ind] <= h_threshold)
        will_start_convecting = needs_to_convect_by_height && iszero(needs_to_convect_by_time) #time needs be updated?
        isconvecting[ind] = needs_to_convect_by_time || needs_to_convect_by_height 
        will_start_convecting && (convection_triggered_time[ind] = t) #Update time only if new convective event
    end
    return nothing
end

"""
    update_convective_events_gpu!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
Use the parameters τ_convec (field of times since last convective event started) and h_threshold to determine if one point should convect.
The version contains the implementation for GPU using CUDA.

Based on:
Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
"""
function update_convective_events_gpu!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
    #@show typeof(h)

    index_x = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride_x = gridDim().x * blockDim().x

    index_y = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    stride_y = gridDim().y * blockDim().y

    
    @inbounds for i in index_x:stride_x:Nx
        for j in index_y:stride_y:Ny
            time_convecting = t - convection_triggered_time[i,j]
            needs_to_convect_by_time = isconvecting[i,j] && (time_convecting < τ_convec) #has been convecting less than τ_c?
            needs_to_convect_by_height = (h[i,j] <= h_threshold)
            will_start_convecting = needs_to_convect_by_height && iszero(needs_to_convect_by_time) #time needs be updated?
            isconvecting[i,j] = needs_to_convect_by_time || needs_to_convect_by_height 
            will_start_convecting && (convection_triggered_time[i,j] = t) #Update time only if new convective event
        end
    end
    
    return nothing
end

"""
    heat(t,distance_from_conv_centersq,conv_time_triggered,q0,τ_c,R2,A0)
For a point that will be heated by convection, compute the value of the convective mass source.

Based on:
Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
"""
@inline function heat(t,distance_from_conv_centersq,conv_time_triggered,q0,τ_c,R2,A0)
    deltat   = t - conv_time_triggered
    quotient = 2.0 * (deltat - τ_c/2.0)/(τ_c)
    q        = q0*(1.0 - quotient*quotient)*(1.0 - (distance_from_conv_centersq / R2))
    return  q / (τ_c*A0)
end

"""
    nth_neighbor(i,n,N) = mod1(i + n,N)
Tells you the index of your nth neighbor considering periodic boundaries.
"""
@inline nth_neighbor(i,n,N) = mod1(i + n,N)

"""
    heat_at_point(i,j,k,clock,τ_c,convective_radius,isconvecting,convection_triggered_time,q0,Δx,Δy,Nx,Ny,numelements_to_traverse_x,numelements_to_traverse_y)
    Centered on each point it will traverese a square of numelements_to_traverse_y * numelements_to_traverse_y. If one of those neighbors are a convective center, it will heat the current point with the rules shown in:

Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
"""
function heat_at_point(i,j,k,clock,τ_c,convective_radius,isconvecting,convection_triggered_time,q0,Δx,Δy,Nx,Ny,numelements_to_traverse_x,numelements_to_traverse_y)
    forcing = 0.0
    R2 = convective_radius * convective_radius
    A0 = pi*R2
    @inbounds for neigh_y_id in (-numelements_to_traverse_y:numelements_to_traverse_y)
        y_distancesq = neigh_y_id * neigh_y_id* Δy * Δy
        #idy = nth_neighbor(j,neigh_y_id,Ny)
        idy = j + neigh_y_id
        @inbounds for neigh_x_id in (-numelements_to_traverse_x:numelements_to_traverse_x)
            x_distancesq = neigh_x_id * neigh_x_id * Δx * Δx
            distance_from_conv_centersq = x_distancesq + y_distancesq #takes time
            if distance_from_conv_centersq < R2
                #idx = nth_neighbor(i,neigh_x_id,Nx)
                idx = i + neigh_x_id
                if (isconvecting[idx , idy] == 1.0)
                    forcing += heat(clock.time,distance_from_conv_centersq,convection_triggered_time[idx,idy],q0,τ_c,R2,A0)
                end
            end
        end
    end 
    return forcing
end

"""
    u_damping(x, y, z, t, u, relaxation_parameter) = - u * relaxation_parameter
Create a linear damping function for the u field
"""
u_damping(x, y, z, t, u, relaxation_parameter) = - u * relaxation_parameter

"""
    v_damping(x, y, z, t, v, relaxation_parameter) = - v * relaxation_parameter
Create a linear damping function for the v field
"""
v_damping(x, y, z, t, v, relaxation_parameter) = - v * relaxation_parameter

"""
    model_forcing(i,j,k,grid,clock,model_fields,parameters)
This is an interface with the correct signature to register the convective parameterization to Oceananigans.jl
It also adds the radiative cooling and the relaxation in the height field with a given timescale.
"""
function model_forcing(i,j,k,grid,clock,model_fields,parameters)
    heat_at_point(i,j,k,clock,
                      parameters.τ_c,
                      parameters.R,
                      parameters.isconvecting,
                      parameters.convection_triggered_time,
                      parameters.q0,
                      grid.Δxᶜᵃᵃ,
                      grid.Δyᵃᶜᵃ,
                      grid.Nx,
                      grid.Ny,
                      parameters.nghosts_x,
                      parameters.nghosts_y,
                      ) - parameters.radiative_cooling_rate - (model_fields.h[i,j,k] - parameters.relaxation_height)*parameters.relaxation_parameter
end

