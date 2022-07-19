"""
    This is intended to be launched from scripts/read_parameter_file_and_launch_15d_simulation.jl
    """
function run_shallow_simulation_debug(arch; ultrashort = false)

    architecture = if arch == "CPU"
        CPU()
    elseif arch == "GPU" 
        GPU()
    end
    @info "Using architecture: ", architecture
    #Setup physicsl parameters
    parameters_dict = Dict{String, Union{Nothing, Real}}(
    "f" => 5e-4, #Coriolis
    "g" => 9.8, #Gravity
    "convection_timescale" => 10800.0, #Duration of convective events (in seconds)
    "convection_critical_height" => 130.0, #Critical height that triggers convection (in meters)
    "heating_amplitude" => 3e9, #The amplitude of the convective event (q0 in paper)
    "large_scale_forcing" => (1.12/3)*1.0e-8, #The amplitude of the large scale forcing
    "convective_radius"    => 30000.0, #Radius of convective event (in meters)
    "relaxation_parameter" => 1.0/(2*86400), # 1/τ where tau is the relaxation timescale for friction and h recovery
    "relaxation_height" => nothing, #131.0 #Target for the recovery of the h field if nothing it will relax to the mean
    "Lx" => 1.5e6, #Size of the domain (in meters)
    "Ly" => 1.5e6,
    "Lz" => 126.5, # A characteristic height
    "Nx" => 500, #Number of points
    "Ny" => 500,
    "boundary_layer" => true #If true, convection is a mass sink, otherwise is false
    )
    
    nghosts = compute_nghosts(parameters_dict)
    grid = RectilinearGrid(architecture,size = (parameters_dict["Nx"], parameters_dict["Ny"]),
                           x = (0, parameters_dict["Lx"]), y = (0, parameters_dict["Ly"]),
                           topology = (Periodic, Periodic, Flat), halo = (nghosts, nghosts))

    @info "Built grid successfully"

    isconvecting, mean_h, parameters = build_convective_parameterization_tools(grid, parameters_dict)

    #build forcing
    convec_forcing = Forcing(model_forcing,discrete_form=true,parameters = parameters)
    u_forcing = Forcing(u_damping, parameters=parameters.relaxation_parameter, field_dependencies=:u)
    v_forcing = Forcing(v_damping, parameters=parameters.relaxation_parameter, field_dependencies=:v)

    ## Build the model

    model = ShallowWaterModel(;grid=grid,
                            timestepper=:RungeKutta3,
                            momentum_advection=WENO5(grid=grid),
                            mass_advection=WENO5(grid=grid),
                            tracer_advection=WENO5(grid=grid),
                            gravitational_acceleration=parameters_dict["g"],
                            coriolis=FPlane(f=parameters_dict["f"]),
                            forcing=(h=convec_forcing,u = u_forcing, v = v_forcing)
                            )

    uhⁱ(x, y, z) = 0.0 #uⁱ(x, y, z) * hⁱ(x, y, z)
    h̄(x, y, z) = parameters_dict["Lz"] + 4rand()
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

    set!(model, uh = uhⁱ, h = h̄)
    set!(mean_h,mean(h))

    stop_time = if ultrashort 
        100
    else
        10_000
    end
    #Create the simulation
    #simulation = Simulation(model, Δt = 1e-2, stop_time = 150)
    simulation = Simulation(model; Δt = 5.0, stop_time)

    function progress(sim)
        m = sim.model
        @info(@sprintf("Iter: %d, time: %.1f, Δt: %.1f, max|h|: %.2f, min|h|: %.2f",
                       m.clock.iteration, m.clock.time,
                       sim.Δt, maximum(abs, m.solution.h),  minimum(abs, m.solution.h)))
        
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
    simulation.callbacks[:update_convective_helper_arrays] = Callback(update_convective_helper_arrays, IterationInterval(1); parameters)
    #prepare output files
    outputfilename = "debug_run_"*arch
    simulation.output_writers[:fields] =
        NetCDFOutputWriter(
            model,
            (;h ,v , u, isconvecting, ω, sp, diver),
            dir = datadir(),
            filename = outputfilename*".nc",
            schedule = IterationInterval(100),
            overwrite_existing = true)
 
#@profview run!(simulation)
run!(simulation)

end #runsimulation
