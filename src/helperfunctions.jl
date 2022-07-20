function update_convective_helper_arrays(sim, parameters)
    p = parameters
    m = sim.model
    update_convective_events!(m.architecture,p.isconvecting,p.convection_triggered_time,m.solution.h,
                              m.clock.time,p.convection_timescale,p.convection_critical_height,m.grid.Nx,m.grid.Ny, p.boundary_layer)
    Oceananigans.BoundaryConditions.fill_halo_regions!(p.isconvecting, m.architecture)
    Oceananigans.BoundaryConditions.fill_halo_regions!(p.convection_triggered_time, m.architecture)
    isnothing(p.relaxation_height) && mean!(p.mean_h,m.solution.h)
end


function build_convective_parameterization_tools(grid, parameters)
    p = parameters

    #These ones come from the external parameters
    convection_timescale       = p["convection_timescale"]
    convection_critical_height = p["convection_critical_height"]
    large_scale_forcing        = p["large_scale_forcing"]
    heating_amplitude          = p["heating_amplitude"]
    convective_radius          = p["convective_radius"]
    relaxation_parameter       = p["relaxation_parameter"]
    relaxation_height          = p["relaxation_height"]
    boundary_layer             = p["boundary_layer"]

    #We then create this ones

    Δx2                        = (grid.Δxᶜᵃᵃ)^2
    Δy2                        = (grid.Δyᵃᶜᵃ)^2
    @assert Δx2 == Δy2
    nghosts                    = Int(convective_radius ÷ grid.Δxᶜᵃᵃ) - 1

    isconvecting               = CenterField(grid,Bool)
    convection_triggered_time  = CenterField(grid)
    ## Will create heating stencil with the spatial component
    heating_stencil            = CenterField(grid,Float64; indices=(-nghosts:nghosts,-nghosts:nghosts,:))
    mean_h                     = Field{Nothing,Nothing,Center}(grid)
    fill_heating_stencil!(grid.architecture,heating_stencil,heating_amplitude,grid.Δxᶜᵃᵃ,convective_radius^2)

    parameters = (; convection_timescale, convection_critical_height, large_scale_forcing,  relaxation_parameter, relaxation_height,boundary_layer, nghosts, isconvecting, convection_triggered_time, heating_stencil, mean_h)


    return isconvecting, mean_h, parameters
end

function compute_nghosts(parameters)
    Lx = parameters["Lx"]
    Nx = parameters["Nx"]
    convective_radius = parameters["convective_radius"]
    grid_spacing_x = Lx ÷ Nx #These two need to be equal for x and y!
    nghosts = max(Int(convective_radius ÷ grid_spacing_x), 3)
    return nghosts
end

const short_parameter_names = Dict(
    "architecture" => "arch",
    "convection_timescale" => "tauc",
    "convection_critical_height" => "hc",
    "heating_amplitude" => "q0",
    "large_scale_forcing" => "r",
    "convective_radius" => "rconvec",
    "relaxation_parameter" => "taurelax",
    "relaxation_height" => "hrelax",
    "boundary_layer" => "bl"
)

function short_name(name,short_names)
    short_name = haskey(short_names,name) ? short_names[name] : name
end

function shorten_names(parameters,short_names)
    return Dict((short_name(key, short_names),value) for (key,value) in parameters)
end