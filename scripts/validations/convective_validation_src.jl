# Goal of this script is validate the convective parameterization that we will use.
# the first validation consists in measuring the mass added to the system by a convective event. This is based on Integrating equation 4 of
#Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
# In which the mass added to the system after a convective event is equal to
# ΔH = 1/3 q0 


function validate(arch)

    architecture = if arch == "CPU"
        CPU()
    elseif arch == "GPU"
        GPU()
    end
    g = 0.0
    τ_c = 3600.0
    h_c = 40
    heating_amplitude = 1e9
    radiative_cooling_rate = 0.0
    convective_radius    = 10_000.0
    relaxation_parameter = 0
    relaxation_height = 45.0
    boundary_layer = false
    Lx = 30_000
    Ly = 30_000
    Lz = 45
    Nx = 100
    Ny = 100
    simulation_length = 7200 / 86400
    timestep = 20.0

    grid_spacing_x = Lx ÷ Nx #These two need to be equal for x and y!
    grid_spacing_y = Ly ÷ Ny

    numelements_to_traverse = Int(convective_radius ÷ grid_spacing_x)
    halo_indices = numelements_to_traverse-1

    grid = RectilinearGrid(architecture,size = (Nx, Ny),
                            x = (0, Lx), y = (0, Ly),
                            topology = (Periodic, Periodic, Flat), halo = (max(numelements_to_traverse,3), max(numelements_to_traverse,3)))

    isconvecting  = CenterField(grid,Bool)
    convection_triggered_time  = CenterField(grid)
    q_stencil = CenterField(grid,Float64; indices=(-halo_indices:halo_indices,-halo_indices:halo_indices,:))
    fill_heating_stencil!(grid.architecture,q_stencil,heating_amplitude,grid_spacing_x,convective_radius^2)


    parameters = (; isconvecting = isconvecting, convection_triggered_time, τ_c, h_c, nghosts = numelements_to_traverse - 1, radiative_cooling_rate , q0 = heating_amplitude, R = convective_radius, relaxation_parameter, relaxation_height, Δx2 = grid_spacing_x^2, Δy2 = grid_spacing_y^2, heating_stencil = q_stencil, boundary_layer)

    convec_forcing = Forcing(model_forcing,discrete_form=true,parameters = parameters)

    model = ShallowWaterModel(;grid = grid,
    timestepper=:RungeKutta3,
    momentum_advection=WENO5(grid=grid),
    mass_advection=WENO5(grid=grid),
    tracer_advection=WENO5(grid=grid),
    gravitational_acceleration=g,
    forcing=(h=convec_forcing,)
    )


    uhⁱ(x, y, z) = 0.0 
    function h̄(x, y, z)
        if x == grid.xᶜᵃᵃ[Nx ÷ 2] && y == grid.yᵃᶜᵃ[Ny ÷ 2]
            return 39
        else
            return Lz 
        end
    end
    uh, vh, h = model.solution

   
    ## Build velocities
    u = uh / h
    v = vh / h

    set!(model, uh = uhⁱ, h = h̄)
    initial_state = deepcopy(model.solution)

    simulation = Simulation(model, Δt = timestep, stop_time = 86400.0*simulation_length)

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
        @info(@sprintf("Iter: %d, time: %.1f, Δt: %.1f, max|h|: %.2f, min|h|: %.2f",
                    m.clock.iteration, m.clock.time,
                    sim.Δt, maximum(m.solution.h), minimum(m.solution.h)))
        
    end

    simulation.output_writers[:fields] =
    NetCDFOutputWriter(
        model,
        (h = h , v = v , u = u, isconvecting = isconvecting ),
        dir = datadir(),
        filename = "validation.nc",
        schedule = IterationInterval(1),
        overwrite_existing = true)

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))
    simulation.callbacks[:update_convective_helper_arrays] = Callback(update_convective_helper_arrays, IterationInterval(1); parameters)

    run!(simulation)

    return initial_state, model.solution, heating_amplitude, grid.Δxᶠᵃᵃ, grid.Δyᵃᶠᵃ
end