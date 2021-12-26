function update_convective_events!(architecture :: CPU,isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny,nghosts_x,nghosts_y)
    update_convective_events_cpu!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny,nghosts_x,nghosts_y)
end


function update_convective_events!(architecture :: GPU,isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny,nghosts_x,nghosts_y)
    @cuda update_convective_events_gpu!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny,nghosts_x,nghosts_y)
end

## This one fails on gpu because it iterates. I need to make a kernel.
function update_convective_events_cpu!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny,nghosts_x,nghosts_y)
    @inbounds for ind in eachindex(h)
        if isconvecting[ind] == 1.0
            ((t - convection_triggered_time[ind]) < τ_convec) && continue
            isconvecting[ind] = 0.0
            convection_triggered_time[ind] = 0.0
        elseif (isconvecting[ind] == 0.0) && (h[ind] <= h_threshold)
            isconvecting[ind] = 1.0
            convection_triggered_time[ind] = t
        end
    end
    return nothing
end

function update_convective_events_gpu!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny,nghosts_x,nghosts_y)
    #@show typeof(h)

    index_x = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride_x = gridDim().x * blockDim().x

    index_y = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    stride_y = gridDim().y * blockDim().y

    
    @inbounds for i in index_x:stride_x:Nx
        for j in index_y:stride_y:Ny
           if isconvecting[i, j] == 1.0
             ((t - convection_triggered_time[i, j]) < τ_convec) && continue
             
                   isconvecting[i, j] = 0.0
                   convection_triggered_time[i, j] = 0.0
             
            elseif (isconvecting[i, j] == 0.0) && (h[i, j] <= h_threshold)
                isconvecting[i, j] = 1.0
                convection_triggered_time[i, j] = t
            end        
        end
    end
    
    return nothing
end

@inline function heat(t,distance_from_conv_centersq,conv_time_triggered,q0,τ_c,R2,A0)
    deltat   = t - conv_time_triggered
    quotient = 2.0 * (deltat - τ_c/2.0)/(τ_c)
    q        = q0*(1.0 - quotient*quotient)*(1.0 - (distance_from_conv_centersq / R2))
    return  q / (τ_c*A0)
end

@inline nth_neighbor(i,n,N) = mod1(i + n,N)

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


damping(x, y, z, t, variable, relaxation_parameter) = - variable * relaxation_parameter


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
                      ) - parameters.radiative_cooling_rate
end

# function model_forcing_gpu(i,j,k,grid,clock,model_fields,parameters)
#     forcing =  heat_at_point_gpu(i,j,k,clock,
#                        parameters.τ_c,
#                        parameters.R,
#                        parameters.isconvecting,
#                        parameters.convection_triggered_time,
#                              parameters.q0,
#                              grid.Δx,
#                              grid.Δy,
#                              grid.Nx,
#                              grid.Ny,
#                              parameters.numelements_to_traverse_x,
#                              parameters.numelements_to_traverse_y,
#                       )

# end
