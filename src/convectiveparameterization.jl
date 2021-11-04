## This one fails on gpu because it iterates. I need to make a kernel.
function update_convective_events!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny,nghosts_x,nghosts_y)
    #@show typeof(h)
    @inbounds for ind in eachindex(h)
        if isconvecting[ind]
            ((t - convection_triggered_time[ind]) < τ_convec) && continue
            isconvecting[ind] = false
            convection_triggered_time[ind] = 0.0
        elseif !isconvecting[ind] && (h[ind] <= h_threshold)
     #       @info "Inside "
            isconvecting[ind] = true
            convection_triggered_time[ind] = t
        end
    end
    isconvecting .= padarray(isconvecting[1:Nx,1:Ny], Pad(:circular, (nghosts_x,nghosts_x), (nghosts_y,nghosts_y) ))
    convection_triggered_time .= padarray(convection_triggered_time[1:Nx,1:Ny], Pad(:circular, (nghosts_x,nghosts_x), (nghosts_y,nghosts_y) ))
    #fill_ghosts!(isconvecting,Nx,Ny,nghosts_x,nghosts_y)
    #fill_ghosts!(convection_triggered_time,Nx,Ny,nghosts_x,nghosts_y)
    return nothing
end

function update_convective_events_gpu!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold,Nx,Ny,nghosts_x,nghosts_y)
    #@show typeof(h)
    
    @inbounds for ind in eachindex(h)
        if isconvecting[ind]
            ((t - convection_triggered_time[ind]) < τ_convec) && continue
            isconvecting[ind] = false
            convection_triggered_time[ind] = 0.0
        elseif !isconvecting[ind] && (h[ind] <= h_threshold)
     #       @info "Inside "
            isconvecting[ind] = true
            convection_triggered_time[ind] = t
        end
    end
    isconvecting .= padarray(isconvecting[1:Nx,1:Ny], Pad(:circular, (nghosts_x,nghosts_x), (nghosts_y,nghosts_y) ))
    convection_triggered_time .= padarray(convection_triggered_time[1:Nx,1:Ny], Pad(:circular, (nghosts_x,nghosts_x), (nghosts_y,nghosts_y) ))
    #fill_ghosts!(isconvecting,Nx,Ny,nghosts_x,nghosts_y)
    #fill_ghosts!(convection_triggered_time,Nx,Ny,nghosts_x,nghosts_y)
    return nothing
end

function update_convective_events_kernel(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold)
    #let us assume this receives one scalar and it indices
#    @inbounds for ind in eachindex(h)
    ind_x
    ind_y
    if isconvecting[ind_x,ind_y]
        ((t - convection_triggered_time[ind_x,ind_y]) < τ_convec) && continue
        isconvecting[ind_x,ind_y] = false
        convection_triggered_time[ind_x,ind_y] = 0.0
    elseif !isconvecting[ind_x,ind_y] && (h[ind_x,ind_y] <= h_threshold)
        #       @info "Inside "
        isconvecting[ind_x,ind_y] = true
        convection_triggered_time[ind_x,ind_y] = t
    end
end
end


function heat(t,distance_from_conv_centersq,conv_time_triggered,q0,τ_c,R2,A0)
    deltat   = t - conv_time_triggered
    quotient = 2.0 * (deltat - τ_c/2.0)/(τ_c)
    q        = q0*(1.0 - quotient*quotient)*(1.0 - (distance_from_conv_centersq / R2))
#    if (q/ (τ_c*A0) ) < 0
#        @info "Heating is negative. This is not right."
#    end
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
                if isconvecting[idx , idy]
                    forcing += heat(clock.time,distance_from_conv_centersq,convection_triggered_time[idx,idy],q0,τ_c,R2,A0)
                end
            end
        end
    end 
    return forcing
end


function model_forcing(i,j,k,grid,clock,model_fields,parameters)
    forcing =  heat_at_point(i,j,k,clock,
                       parameters.τ_c,
                       parameters.R,
                       parameters.isconvecting,
                       parameters.convection_triggered_time,
                             parameters.q0,
                             grid.Δx,
                             grid.Δy,
                             grid.Nx,
                             grid.Ny,
                             parameters.numelements_to_traverse_x,
                             parameters.numelements_to_traverse_y,
                      )

end
