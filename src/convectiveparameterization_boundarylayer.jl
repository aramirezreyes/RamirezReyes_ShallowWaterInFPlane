"""
    update_convective_events!(architecture,isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
Use the parameters τ_convec (field of times since last convective event started) and h_threshold to determine if one point should convect.
The version that receives architecture as first parameter is an interface for the specific implementation.

Based on:
Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
"""
function update_convective_events_bl!(architecture :: CPU,isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
    update_convective_events_cpu_bl!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
end

"""
    update_convective_events!(architecture,isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
Use the parameters τ_convec (field of times since last convective event started) and h_threshold to determine if one point should convect.
The version that receives architecture as first parameter is an interface for the specific implementation.

Based on:
Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
"""
function update_convective_events_bl!(architecture :: GPU,isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
    kernel = @cuda launch=false update_convective_events_gpu_bl!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
    config = launch_configuration(kernel.fun)
    threads = min(size(isconvecting,1), config.threads)
    blocks = cld(size(isconvecting,1), threads)
    CUDA.@sync begin
        kernel(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny;threads, blocks)
    end
end

"""
    update_convective_events_cpu!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
Use the parameters τ_convec (field of times since last convective event started) and h_threshold to determine if one point should convect.
The version contains the implementation for CPU.

Based on:
Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
"""
@inline function update_convective_events_cpu_bl!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
     @inbounds Threads.@threads for ind in eachindex(h)
        time_convecting = t - convection_triggered_time[ind]
        needs_to_convect_by_time = isconvecting[ind] && (time_convecting < τ_convec) #has been convecting less than τ_c?
        needs_to_convect_by_height = (h[ind] >= h_threshold)
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
function update_convective_events_gpu_bl!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny)
    #@show typeof(h)

    index_x = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride_x = gridDim().x * blockDim().x

    index_y = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    stride_y = gridDim().y * blockDim().y

    
    @inbounds for i in index_x:stride_x:Nx
        for j in index_y:stride_y:Ny
            time_convecting = t - convection_triggered_time[i,j]
            needs_to_convect_by_time = isconvecting[i,j] && (time_convecting < τ_convec) #has been convecting less than τ_c?
            needs_to_convect_by_height = (h[i,j] >= h_threshold)
            will_start_convecting = needs_to_convect_by_height && iszero(needs_to_convect_by_time) #time needs be updated?
            isconvecting[i,j] = needs_to_convect_by_time || needs_to_convect_by_height 
            will_start_convecting && (convection_triggered_time[i,j] = t) #Update time only if new convective event
        end
    end
    
    return nothing
end

"""
    heat_at_point(i,j,k,clock,τ_c,convective_radius,isconvecting,convection_triggered_time,q0,Δx,Δy,numelements_to_traverse_x,numelements_to_traverse_y)
    Centered on each point it will traverese a square of numelements_to_traverse_y * numelements_to_traverse_y. If one of those neighbors are a convective center, it will heat the current point with the rules shown in:

Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
"""
@inline function heat_at_point_bl(i,j,k,current_time,τc,convective_radius,isconvecting,convection_triggered_time,q0,Δx2,Δy2,numelements_to_traverse, heating_stencil)
    forcing = 0.0
    @inbounds for neigh_j in (-numelements_to_traverse:numelements_to_traverse)
        @inbounds for neigh_i in (-numelements_to_traverse:numelements_to_traverse)
             if isconvecting[i + neigh_i , j + neigh_j]
                 forcing -= heat(current_time,convection_triggered_time[i + neigh_i,j + neigh_j],τc,heating_stencil[-neigh_i,-neigh_j])
             end
        end
    end 
    return forcing
end


"""
    model_forcing(i,j,k,grid,clock,model_fields,parameters)
This is an interface with the correct signature to register the convective parameterization to Oceananigans.jl
It also adds the radiative cooling and the relaxation in the height field with a given timescale.
"""
function model_forcing_bl(i,j,k,grid,clock,model_fields,parameters)
    heat_at_point_bl(i,j,k,clock.time,
                      parameters.τ_c,
                      parameters.R,
                      parameters.isconvecting,
                      parameters.convection_triggered_time,
                      parameters.q0,
                      parameters.Δx2,
                      parameters.Δy2,
                      parameters.nghosts,
                      parameters.heating_stencil) + parameters.radiative_cooling_rate - (model_fields.h[i,j,k] - parameters.relaxation_height)*parameters.relaxation_parameter
end

