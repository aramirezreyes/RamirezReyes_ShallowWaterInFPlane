function update_convective_helper_arrays(sim, parameters)
    p = parameters
    m = sim.model
    update_convective_events!(
        m.architecture,
        p.isconvecting,
        p.convection_triggered_time,
        m.solution.h,
        m.clock.time,
        p.convection_timescale,
        p.convection_critical_height,
        m.grid.Nx,
        m.grid.Ny,
        p.boundary_layer,
    )
    Oceananigans.BoundaryConditions.fill_halo_regions!(p.isconvecting, m.architecture)
    Oceananigans.BoundaryConditions.fill_halo_regions!(
        p.convection_triggered_time,
        m.architecture,
    )
    isnothing(p.relaxation_height) && mean!(p.mean_h, m.solution.h)
    return nothing
end

function update_convec_heating(sim, parameters)
    p = parameters
    m = sim.model
    set!(p.convec_heating, 0.0)
    RamirezReyes_ShallowWaterInFPlane.fill_convec_heating!(
        p.convec_heating,
        m.grid,
        m.clock,
        p,
    )
end

function convective_heating_output(model, parameters)
    p = parameters
    m = model
    set!(p.convec_heating, 0.0)
    #RamirezReyes_ShallowWaterInFPlane.fill_convec_heating!(
    #    p.convec_heating,
    #    m.grid,
    #    m.clock,
    #    p,
    #)
    event = launch!(model.architecture, model.grid, :xyz, model_forcing_only_convec!, p.convec_heating, model.clock, p)
    wait(event)
    return adapt(Array,interior(p.convec_heating))
end

function build_convective_parameterization_tools(grid, parameters)
    p = parameters

    #These ones come from the external parameters
    convection_timescale = p["convection_timescale"]
    convection_critical_height = p["convection_critical_height"]
    large_scale_forcing = p["large_scale_forcing"]
    heating_amplitude = p["heating_amplitude"]
    convective_radius = p["convective_radius"]
    relaxation_parameter = p["relaxation_parameter"]
    relaxation_height = p["relaxation_height"]
    boundary_layer = p["boundary_layer"]

    #We then create this ones

    Δx2 = (grid.Δxᶜᵃᵃ)^2
    Δy2 = (grid.Δyᵃᶜᵃ)^2
    @assert Δx2 == Δy2
    nghosts = Int(convective_radius ÷ grid.Δxᶜᵃᵃ)

    isconvecting = CenterField(grid, Bool)
    convection_triggered_time = CenterField(grid)
    ## Will create heating stencil with the spatial component
    heating_stencil =
        CenterField(grid, Float64; indices = (-nghosts:nghosts, -nghosts:nghosts, :))
    mean_h = Field{Nothing,Nothing,Center}(grid)
    convec_heating = Field{Center,Center,Center}(grid)
    fill_heating_stencil!(
        grid.architecture,
        heating_stencil,
        heating_amplitude,
        grid.Δxᶜᵃᵃ,
        convective_radius^2,
    )

    parameters = (;
        convection_timescale,
        convection_critical_height,
        large_scale_forcing,
        relaxation_parameter,
        relaxation_height,
        boundary_layer,
        nghosts,
        isconvecting,
        convection_triggered_time,
        heating_stencil,
        mean_h,
        convec_heating,
    )


    return isconvecting, mean_h, convec_heating, parameters
end

function compute_nghosts(parameters)
    Lx = parameters["Lx"]
    Nx = parameters["Nx"]
    convective_radius = parameters["convective_radius"]
    grid_spacing_x = Lx ÷ Nx #These two need to be equal for x and y!
    nghosts = max(Int(convective_radius ÷ grid_spacing_x + 1), 3)
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
    "boundary_layer" => "bl",
)

function short_name(name, short_names)
    short_name = haskey(short_names, name) ? short_names[name] : name
end

function shorten_names(parameters, short_names)
    return Dict((short_name(key, short_names), value) for (key, value) in parameters)
end


"""
    create_debug_parameters(arch :: AbstractString, ultrashort = true)
Returns a dictionary with a set of parameters for the debug simulation.

"""
function create_debug_parameters(arch::AbstractString; ultrashort = true, restart = false)

    parameters_dict = Dict{String,Union{Nothing,Real,String}}(
        "architecture" => arch,
        "output_filename" => "debug_run_" * arch,
        "f" => 5e-4, #Coriolis
        "g" => 9.8, #Gravity
        "convection_timescale" => 10800.0, #Duration of convective events (in seconds)
        "convection_critical_height" => 130.0, #Critical height that triggers convection (in meters)
        "heating_amplitude" => 3e9, #The amplitude of the convective event (q0 in paper)
        "large_scale_forcing" => (1.12 / 3) * 1.0e-8, #The amplitude of the large scale forcing
        "convective_radius" => 30000.0, #Radius of convective event (in meters)
        "relaxation_parameter" => 1.0 / (2 * 86400), # 1/τ where tau is the relaxation timescale for friction and h recovery
        "relaxation_height" => nothing, #131.0 #Target for the recovery of the h field if nothing it will relax to the mean
        "Lx" => 1.5e6, #Size of the domain (in meters)
        "Ly" => 1.5e6,
        "Lz" => 126.5, # A characteristic height
        "Nx" => 500, #Number of points
        "Ny" => 500,
        "boundary_layer" => true, #If true, convection is a mass sink, otherwise is false,
        "initialization_style" => "rand",
        "initialization_amplitude" => 4.0,
        "simulation_length_in_days" => ultrashort ? 100 / 86400 : 10000 / 86400,
        "output_interval_in_seconds" => 5,
        "timestep_in_seconds" => 5,
        "restart" => false,
        "checkpoint_interval_in_seconds" => 5,
    )

    if restart
        parameters_dict["restart"] = true
        parameters_dict["simulation_length_in_days"] = 2*parameters_dict["simulation_length_in_days"]
    end

    return parameters_dict
end


function restore_helper_fields!(isconvecting,convection_triggered_time,filename)
    jldopen(joinpath(datadir(),filename), "r") do file
        last_timestep = last(keys(file["timeseries/isconvecting"]))
        isconvecting_r = file["timeseries/isconvecting"][last_timestep]
        convection_triggered_time_r = file["timeseries/convection_triggered_time"][last_timestep]
        copyto!(isconvecting.data.parent, isconvecting_r[:,:,1])
        copyto!(convection_triggered_time.data.parent, convection_triggered_time_r[:,:,1])
    end
    nothing
end
