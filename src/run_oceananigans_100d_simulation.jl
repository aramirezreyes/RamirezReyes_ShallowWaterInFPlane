"""
    This is intended to be launched from scripts/read_parameter_file_and_launch_100d_simulation.jl
    """
function run_shallow_simulation_100d(params)

    architecture = if params["architecture"] == "CPU"
        CPU()
    elseif params["architecture"] == "GPU"
        GPU()
    end
    @info "Using architecture: ", architecture
    #Setup physicsl parameters
    f = params["f"]#5e-4 #5e-4
    g = params["g"]# 9.8
    τ_c = params["tauconvec"]
    h_c = params["hc"]
    heating_amplitude = params["q0"]
    radiative_cooling_rate = params["r"]
    convective_radius    = params["rconvec"]
    relaxation_parameter = params["taurelax"]
    relaxation_height = params["hrelax"]
    boundary_layer = params["boundarylayer"]
    Lx = params["Lx"]
    Ly = params["Ly"]
    Lz = params["Lz"]
    Nx = params["Nx"]
    Ny = params["Ny"]

    

    grid_spacing_x = Lx ÷ Nx #These two need to be equal for x and y!
    grid_spacing_y = Ly ÷ Ny

    numelements_to_traverse = Int(convective_radius ÷ grid_spacing_x)
    halo_indices = numelements_to_traverse-1

    grid = RectilinearGrid(architecture,size = (Nx, Ny),
                           x = (0, Lx), y = (0, Ly),
                           topology = (Periodic, Periodic, Flat), halo = (max(numelements_to_traverse,3), max(numelements_to_traverse,3)))

    @info "Built grid successfully"
    isconvecting  = CenterField(grid,Bool)
    convection_triggered_time  = CenterField(grid)
    ## Will create heating stencil with the spatial component
    q_stencil = CenterField(grid,Float64; indices=(-halo_indices:halo_indices,-halo_indices:halo_indices,:))

    fill_heating_stencil!(grid.architecture,q_stencil,heating_amplitude,grid_spacing_x,convective_radius^2)

    parameters = (; isconvecting = isconvecting, convection_triggered_time, τ_c, h_c, nghosts = numelements_to_traverse - 1, radiative_cooling_rate , q0 = heating_amplitude, R = convective_radius, relaxation_parameter, relaxation_height, Δx2 = grid_spacing_x^2, Δy2 = grid_spacing_y^2, heating_stencil = q_stencil, boundary_layer)


    #build forcing
    convec_forcing = Forcing(model_forcing,discrete_form=true,parameters = parameters)
    u_forcing = Forcing(u_damping, parameters=relaxation_parameter, field_dependencies=:u)
    v_forcing = Forcing(v_damping, parameters=relaxation_parameter, field_dependencies=:v)

    ## Build the model

    model = ShallowWaterModel(;timestepper=:RungeKutta3,
                              momentum_advection=WENO5(grid=grid),
                              grid=grid,
                              gravitational_acceleration=g,
                              coriolis=FPlane(f=f),
                              forcing=(h=convec_forcing,u = u_forcing, v = v_forcing)
                              )


    #Build background state and perturbations

    uhⁱ(x, y, z) = 0.0 #uⁱ(x, y, z) * hⁱ(x, y, z)
    h̄(x, y, z) = Lz + 4rand()
    uh, vh, h = model.solution

    ## Build velocities
    u = uh / h
    v = vh / h

    ## Build and compute mean vorticity discretely
    ω = Field(∂x(v) - ∂y(u))
    diver = Field(∂x(u) + ∂y(v))
    sp = @at (Center, Center, Center) sqrt(u^2 + v^2)
    compute!(ω)

    ## Copy mean vorticity to a new field
    ωⁱ = Field{Face, Face, Nothing}(model.grid)
    ωⁱ .= ω

    ## Use this new field to compute the perturbation vorticity
    ω′ = Field(ω - ωⁱ)

    # and finally set the "true" initial condition with noise,

    set!(model, uh = uhⁱ, h = h̄)


    #Create the simulation
    #simulation = Simulation(model, Δt = 1e-2, stop_time = 150)
    simulation = Simulation(model, Δt = 5.0, stop_time = 86400*100.0f0)

    function update_convective_helper_arrays(sim, parameters)
        p = parameters
        m = sim.model
        update_convective_events!(m.architecture,p.isconvecting,p.convection_triggered_time,m.solution.h,
                                  m.clock.time,p.τ_c,p.h_c,m.grid.Nx,m.grid.Ny, p.boundary_layer)
        Oceananigans.BoundaryConditions.fill_halo_regions!(p.isconvecting, m.architecture)
        Oceananigans.BoundaryConditions.fill_halo_regions!(p.convection_triggered_time, m.architecture)

    end


    function progress(sim)
        m = sim.model
        @info(@sprintf("Iter: %d, time: %.1f, Δt: %.1f",
                       m.clock.iteration, m.clock.time,
                       sim.Δt))
        
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
    simulation.callbacks[:update_convective_helper_arrays] = Callback(update_convective_helper_arrays, IterationInterval(1); parameters)
    #prepare output files
    outputfilename = savename(params, ignores=("architecture","g","Lx","Nx","Ny"))
    simulation.output_writers[:fields] =
        NetCDFOutputWriter(
            model,
            (h = h , v = v , u = u, isconvecting = isconvecting, ω, ω′, sp, diver),
            dir = joinpath(ENV["SCRATCH"],projectname(),"data"),
            filename = outputfilename*".nc",
            schedule = IterationInterval(1200),
            overwrite_existing = true)
    
    
    run!(simulation)
    
end #runsimulation
