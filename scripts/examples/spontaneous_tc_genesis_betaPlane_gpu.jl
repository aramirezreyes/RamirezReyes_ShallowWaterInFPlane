using RamirezReyes_ShallowWaterInFPlane


function run_shallow_simulation(arch = "CPU")
    architecture = if arch == "CPU"
        CPU()
    elseif arch == "GPU" 
        GPU()
    end
    @info "Using architecture: ", architecture
    #Setup physicsl parameters
    f = 3e-5
    β = 3e-5 / 2e7                                  #Coriolis
    g = 9.8                                   #Gravity
    τ_c = 28800.0                             #Duration of convective events [s]
    h_c = 130.0                               #Critical height that triggers convection [m]
    heating_amplitude = 1e9                   #Amplitude of the convective event (q0 in paper)
    radiative_cooling_rate = (1.12/3)*1.0e-9  #The amplitude of the large scale forcing
    convective_radius    = 30000.0            #Radius of convective event [m]
    relaxation_parameter = 1.0/(2*86400)      #1/τ where τ is the relaxation timescale for friction and h recovery
    relaxation_height = 131.0                 #Target for the recovery of the h field
    Lx = 8e6                                #Size of the domain [m]
    Ly = 4e6
    Lz = 126.5                                # Acharacteristic height used in initialization
    Nx = 1000                                  #Number of points
    Ny = 500
    boundary_layer = true                     #If true, convection is a mass sink, otherwise is false
    
    ## Setup other simulation prameters
    output_dir = datadir()
    output_filename = "spontaneous_tc_genesis_run_"*arch*".nc"
    save_every = 7200                          #in number of iterations
    simulation_length = 86400.0*20             #in seconds
    Δt = 40.0                                   #timestep in seconds

    ####################################
    ##### Simulation setup below ######
    ###################################
    grid_spacing_x = Lx ÷ Nx #These two need to be equal for x and y!
    grid_spacing_y = Ly ÷ Ny

    numelements_to_traverse = Int(convective_radius ÷ grid_spacing_x)
    halo_indices = numelements_to_traverse-1

    grid = RectilinearGrid(architecture,size = (Nx, Ny),
        x = (0, Lx), y = (0, Ly),
        topology = (Periodic, Bounded, Flat), 
        halo = (max(numelements_to_traverse,3), max(numelements_to_traverse,3)))

    @info "Built grid successfully"

    ## Fields needed for the convective parameterization
    isconvecting               = CenterField(grid,Bool)
    convection_triggered_time  = CenterField(grid)
    ## Will create heating stencil with the spatial component
    q_stencil = CenterField(grid,Float64; indices=(-halo_indices:halo_indices,-halo_indices:halo_indices,:))
    fill_heating_stencil!(grid.architecture,q_stencil,heating_amplitude,grid_spacing_x,convective_radius^2)


    parameters = (; isconvecting = isconvecting, convection_triggered_time, τ_c, h_c, 
    nghosts = numelements_to_traverse - 1, radiative_cooling_rate , q0 = heating_amplitude,
    R = convective_radius, relaxation_parameter, relaxation_height, Δx2 = grid_spacing_x^2,
    Δy2 = grid_spacing_y^2, heating_stencil = q_stencil, boundary_layer
    )

    #build forcing
    convec_forcing = Forcing(model_forcing,discrete_form=true,parameters = parameters)
    u_forcing      = Forcing(u_damping, parameters=relaxation_parameter, field_dependencies=:u)
    v_forcing      = Forcing(v_damping, parameters=relaxation_parameter, field_dependencies=:v)
    no_slip_bc = ValueBoundaryCondition(0.0)
    no_slip_field_bcs = FieldBoundaryConditions(no_slip_bc);
    # Build the model
    model = if iszero(f)
        ShallowWaterModel(;grid=grid,
        timestepper=:RungeKutta3,
        momentum_advection=WENO5(grid=grid),
        mass_advection=WENO5(grid=grid),
        tracer_advection=WENO5(grid=grid),
        gravitational_acceleration=g,
        forcing=(h=convec_forcing,u = u_forcing, v = v_forcing)
        )

    else
         ShallowWaterModel(;grid=grid,
        timestepper=:RungeKutta3,
        momentum_advection=WENO5(grid=grid),
        mass_advection=WENO5(grid=grid),
        tracer_advection=WENO5(grid=grid),
        gravitational_acceleration=g,
        coriolis=BetaPlane(f₀=f, β=β),
        forcing=(h=convec_forcing,u = u_forcing, v = v_forcing)
        )
end

#######################################
###########Initial conditions #########
#######################################

#### This shallow water model uses the fluxes as prognostic variables we will initialize with these functions
    uhⁱ(x, y, z) = 0.0 
    h̄(x, y, z)   = Lz + 4rand()
    set!(model, uh = uhⁱ, h = h̄)


    ### For output purposes we want the velocities and vorticity
    uh, vh, h = model.solution
    u         = uh / h
    v         = vh / h

    ## Build and compute mean vorticity discretely
    ω = Field(∂x(v) - ∂y(u))
    diver = Field(∂x(u) + ∂y(v))
    sp = @at (Center, Center, Center) sqrt(u^2 + v^2)
    compute!(ω)
    simulation = Simulation(model; Δt = Δt , stop_time = simulation_length)

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
    
    
    #### Prepare output files ####
    
    simulation.output_writers[:fields] =
        NetCDFOutputWriter(
            model,
            (;h ,v , u, isconvecting, sp, diver, ω),
            dir = output_dir,
            filename = output_filename,
            schedule = TimeInterval(save_every),
            overwrite_existing = true)
 
run!(simulation)

end 

run_shallow_simulation("GPU")
