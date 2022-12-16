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

    nghosts = compute_nghosts(parameters_dict)
    grid = RectilinearGrid(architecture,size = (parameters_dict["Nx"], parameters_dict["Ny"]),
                           x = (0, parameters_dict["Lx"]), y = (0, parameters_dict["Ly"]),
                           topology = (Periodic, Periodic, Flat), halo = (nghosts, nghosts))

#    @info "Building grid..."

    isconvecting, mean_h, convec_heating, parameters = build_convective_parameterization_tools(grid, parameters_dict)

    #build forcing
    convec_forcing = Forcing(model_forcing,discrete_form=true,parameters = parameters)
    u_forcing = Forcing(u_damping, parameters=parameters.relaxation_parameter, field_dependencies=:u)
    v_forcing = Forcing(v_damping, parameters=parameters.relaxation_parameter, field_dependencies=:v)

    ## Build the model

    model = ShallowWaterModel(;grid=grid,
                            timestepper=:RungeKutta3,
                            momentum_advection=WENO(grid=grid),
                            mass_advection=WENO(grid=grid),
                            tracer_advection=WENO(grid=grid),
                            gravitational_acceleration=parameters_dict["g"],
                            coriolis=FPlane(f=parameters_dict["f"]),
                            forcing=(h=convec_forcing,u = u_forcing, v = v_forcing)
                            )

    uhⁱ(x, y, z) = 0.0 
    h0_rand(x,y,z) = parameters_dict["Lz"] + parameters_dict["initialization_amplitude"]*rand()

    function h0_one_convecting_point(x,y,z) 
        if x == grid.xᶜᵃᵃ[parameters_dict["Nx"] ÷ 2] && y == grid.yᵃᶜᵃ[parameters_dict["Ny"] ÷ 2]
            parameters_dict["boundary_layer"] ? parameters_dict["convection_critical_height"] + parameters_dict["initialization_amplitude"] : parameters_dict["convection_critical_height"] - parameters_dict["initialization_amplitude"]
        else
            parameters_dict["boundary_layer"] ? parameters_dict["convection_critical_height"] - 1.0 : parameters_dict["convection_critical_height"] + 1.0
        end
    end

    function rankine_vortex(center, radius, amplitude, h0)
        
        function vortex_maker_u(x,y,z)
            amplitude2 = amplitude
            x_relative = (x - center[1])
            y_relative = (y - center[2])
            dist_from_center = hypot(x_relative,y_relative)
            theta = atan(y_relative,x_relative)
            if dist_from_center <= radius
                return -amplitude2 * h0 * dist_from_center/radius^2 * dist_from_center*sin(theta)
            else
                return -amplitude2 * h0 * sin(theta)
            end
        end


        function vortex_maker_v(x,y,z)
            amplitude2 = amplitude
            x_relative = (x - center[1])
            y_relative = (y - center[2])
            dist_from_center = hypot(x_relative,y_relative)
            theta = atan(y_relative,x_relative)
            if dist_from_center <= radius
                return amplitude2 * h0 * dist_from_center/radius^2 * dist_from_center*cos(theta)
            else
                return amplitude2 * h0 * cos(theta)
            end
        end

        return (vortex_maker_u,vortex_maker_v)
    end
    
    uh, vh, h = model.solution

    ## Build velocities
    u = uh / h
    v = vh / h

    ## Build and compute mean vorticity discretely
    ω = Field(∂x(v) - ∂y(u))
    diver = Field(∂x(u) + ∂y(v))
    sp = @at (Center,Center, Center) sqrt(u^2 + v^2)
    compute!(ω)


    # and finally set the "true" initial condition with noise,
    
    if parameters_dict["initialization_style"] == "rand"
        set!(model, uh = uhⁱ, h = h0_rand)
    elseif parameters_dict["initialization_style"] == "one_convecting_point"
        set!(model, uh = uhⁱ, h = h0_one_convecting_point)
    elseif parameters_dict["initialization_style"] == "rankine_vortex"
        u_maker,v_maker = rankine_vortex(parameters_dict["rankine_center"],parameters_dict["rankine_radius"],parameters_dict["rankine_amplitude"],parameters_dict["Lz"])
        set!(model, uh = u_maker, vh = v_maker, h = h0_rand)
    else 
        error("Intialization style must be either \"rand\" or \"one_convecting_point\"")
    end
    mean!(mean_h,h)


    #Create the simulation
    #simulation = Simulation(model, Δt = 1e-2, stop_time = 150)
    simulation = Simulation(model, Δt = timestep, stop_time = 86400.0*simulation_length)

    function progress(sim)
        m = sim.model
        @info(@sprintf("Iter: %d, time: %.1f, Δt: %.1f",
                       m.clock.iteration, m.clock.time,
                       sim.Δt))
        
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
    simulation.callbacks[:update_convective_helper_arrays] = Callback(update_convective_helper_arrays, IterationInterval(1); parameters)
    simulation.callbacks[:update_convec_heating] = Callback(update_convec_heating, TimeInterval(save_every); parameters)
    #prepare output files
    outputfilename = haskey(parameters_dict, "output_filename") ? parameters_dict["output_filename"] : savename(shorten_names(parameters_dict, short_parameter_names), ignores=("arch","g","Lx","Nx","output_interval_in_seconds", "simulation_length_in_days", "timestep_in_seconds", "output_filename", "initialization_amplitude", "initialization_style"))
    simulation.output_writers[:fields] =
        NetCDFOutputWriter(
            model,
            (h = h , v = v , u = u, isconvecting = isconvecting = isconvecting, ω = ω, sp = sp, diver = diver, convec_heating = convec_heating),
            dir = datadir(),
            filename = outputfilename*".nc",
            schedule = TimeInterval(save_every),
            overwrite_existing = true,
            compression = 1,
        )
    
    
    run!(simulation)
    @info "Simulation finished in $(prettytime(simulation.run_wall_time))."
end #runsimulation
