function update_convective_events!(isconvecting,convection_triggered_time,h,t,τ_convec,h_threshold)
    for ind in eachindex(h)
        if isconvecting[ind]
            ((t - convection_triggered_time[ind]) < τ_convec) && continue
            isconvecting[ind] = false
            convection_triggered_time[ind] = 0.0
        elseif !isconvecting[ind] && (h[ind] <= h_threshold)
            @info "Inside "
            isconvecting[ind] = true
            convection_triggered_time[ind] = t
        end
    end
    return nothing
end


function heat_at_point(i,j,k,grid,h)
#    element_at_right(i) = 
    @show grid
end
