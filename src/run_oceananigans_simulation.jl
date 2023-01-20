"""
    This is intended to be launched from scripts/read_parameter_file_and_launch_simulation.jl
"""
function run_shallow_simulation(parameters_dict)

    architecture = if parameters_dict["architecture"] == "CPU"
        CPU()
    elseif parameters_dict["architecture"] == "GPU"
        GPU()
    end
    @info "Using architecture: $(string(architecture))..."

    simulation_length = parameters_dict["simulation_length_in_days"]
    save_every = parameters_dict["output_interval_in_seconds"]
    timestep = parameters_dict["timestep_in_seconds"]
    pickup = parameters_dict["restart"]
    checkpoint_interval = parameters_dict["checkpoint_interval_in_seconds"]

    nghosts = compute_nghosts(parameters_dict)
    grid = RectilinearGrid(
        architecture,
        size = (parameters_dict["Nx"], parameters_dict["Ny"]),
        x = (0, parameters_dict["Lx"]),
        y = (0, parameters_dict["Ly"]),
        topology = (Periodic, Periodic, Flat),
        halo = (nghosts, nghosts),
    )

    #    @info "Building grid..."

    isconvecting, mean_h, convec_heating, parameters =
        build_convective_parameterization_tools(grid, parameters_dict)

    #build forcing
    convec_forcing = Forcing(model_forcing, discrete_form = true, parameters = parameters)
    u_forcing = Forcing(
        u_damping,
        parameters = parameters.relaxation_parameter,
        field_dependencies = :u,
    )
    v_forcing = Forcing(
        v_damping,
        parameters = parameters.relaxation_parameter,
        field_dependencies = :v,
    )

    ## Build the model

    model = ShallowWaterModel(;
        grid = grid,
        timestepper = :RungeKutta3,
        momentum_advection = WENO(grid = grid),
        mass_advection = WENO(grid = grid),
        tracer_advection = WENO(grid = grid),
        gravitational_acceleration = parameters_dict["g"],
        coriolis = FPlane(f = parameters_dict["f"]),
        forcing = (h = convec_forcing, u = u_forcing, v = v_forcing),
    )



    uh, vh, h = model.solution

    ## Build velocities
    u = uh / h
    v = vh / h

    ## Build and compute mean vorticity discretely
    ω = Field(∂x(v) - ∂y(u))
    diver = Field(∂x(u) + ∂y(v))
    sp = @at (Center, Center, Center) sqrt(u^2 + v^2)
    compute!(ω)



    uh0_f, vh0_f, h0_f = create_initialization_functions(grid, parameters_dict)
    set!(model, uh = uh0_f, vh = vh0_f, h = h0_f)
    mean!(mean_h, h)

    #Create the simulation
    #simulation = Simulation(model, Δt = 1e-2, stop_time = 150)
    simulation = Simulation(model, Δt = timestep, stop_time = 86400.0 * simulation_length)

    function progress(sim)
        m = sim.model
        @info(
            @sprintf(
                "Iter: %d, time: %.1f, Δt: %.1f",
                m.clock.iteration,
                m.clock.time,
                sim.Δt
            )
        )

    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
    simulation.callbacks[:update_convective_helper_arrays] =
        Callback(update_convective_helper_arrays, IterationInterval(1); parameters)
    convec_heating_f = model -> convective_heating_output(model, parameters)
#    simulation.callbacks[:update_convec_heating] =
#        Callback(update_convec_heating, TimeInterval(save_every); parameters)
    #prepare output files
    outputfilename =
        haskey(parameters_dict, "output_filename") ? parameters_dict["output_filename"] :
        savename(
            shorten_names(parameters_dict, short_parameter_names),
            ignores = (
                "arch",
                "g",
                "Lx",
                "Nx",
                "output_interval_in_seconds",
                "simulation_length_in_days",
                "timestep_in_seconds",
                "output_filename",
                "initialization_amplitude",
                "initialization_style",
            ),
        )
    overwrite_existing_file = pickup ? false : true
    simulation.output_writers[:fields] = NetCDFOutputWriter(
        model,
        (
            h = h,
            v = v,
            u = u,
            ω = ω,
            sp = sp,
            convec_heating = convec_heating_f,
        ),
        dir = datadir(),
        filename = outputfilename * ".nc",
        schedule = TimeInterval(save_every),
        overwrite_existing = overwrite_existing_file,
        compression = 1,
        dimensions = Dict("convec_heating" => ("xC","yC","zC"))
    )
    #### CHECKPOINTERS
    simulation.output_writers[:checkpointer] = Checkpointer(model,dir = datadir(),
     schedule=TimeInterval(checkpoint_interval), prefix="checkpoint_"*outputfilename*"_",
     properties = [:architecture, :grid, :clock, :coriolis, :closure, :velocities, 
     :tracers, :timestepper])

    if pickup
        restore_helper_fields!(isconvecting,parameters.convection_triggered_time,"checkpoint_"*outputfilename*"_helper_arrays.jld2")
    end
    simulation.output_writers[:checkpointer_helpers] = JLD2OutputWriter(
        model,
        (;
            isconvecting = isconvecting,
            convection_triggered_time = parameters.convection_triggered_time,
        ),
        dir = datadir(),
        filename = "checkpoint_"*outputfilename*"_helper_arrays",
        schedule = TimeInterval(checkpoint_interval),
        overwrite_existing = true,
        with_halos = true,
        )


   

    run!(simulation; pickup)
    @info "Simulation finished in $(prettytime(simulation.run_wall_time))."
end #runsimulation
