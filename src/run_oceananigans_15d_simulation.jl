"""
This is intended to be launched from scripts/read_parameter_file_and_launch_15d_simulation.jl
"""
function run_shallow_simulation(params)

    architecture = if params["architecture"] == "CPU"
        CPU()
    elseif params["architecture"] == "GPU"
        GPU()
    end
    #Setup physicsl parameters
    f = params["coriolis_parameter"]#5e-4 #5e-4
    g = params["gravitational_acceleration"]# 9.8
    τ_c = params["convection_timescale"]
    h_c = params["critical_height"]
    heating_amplitude = params["heating_amplitude"]
    radiative_cooling_rate = params["radiative_cooling_rate"]
    convective_radius    = params["convective_radius"]
    relaxation_parameter = params["relaxation_timescale"]
    relaxation_height = params["relaxation_height"]
    Lx = params["Lx"]
    Ly = params["Ly"]
    Lz = params["Lz"]
    Nx = params["Nx"]
    Ny = params["Ny"]

    
    
    
    
    grid_spacing_x = Lx ÷ Nx
    grid_spacing_y = Ly ÷ Ny

    
    numelements_to_traverse_x = Int(convective_radius ÷ grid_spacing_x)
    numelements_to_traverse_y = Int(convective_radius ÷ grid_spacing_y)


    grid = RectilinearGrid(size = (Nx, Ny),
                           x = (0, Lx), y = (0, Ly),
                           topology = (Periodic, Periodic, Flat), halo = (max(numelements_to_traverse_x,3), max(numelements_to_traverse_y,3)))

    @info "Built grid successfully"
    isconvecting  = CenterField(grid,Bool)
    convection_triggered_time  = CenterField(grid)

    parameters = (; isconvecting = isconvecting, convection_triggered_time, τ_c, h_c, nghosts_x = numelements_to_traverse_x, nghosts_y = numelements_to_traverse_y, radiative_cooling_rate , q0 = heating_amplitude, R = convective_radius, relaxation_parameter, relaxation_height)


    #build forcing
    convec_forcing = Forcing(model_forcing,discrete_form=true,parameters = parameters)
    u_forcing = Forcing(u_damping, parameters=relaxation_parameter, field_dependencies=:u)
    v_forcing = Forcing(v_damping, parameters=relaxation_parameter, field_dependencies=:v)

    ## Build the model

    model = ShallowWaterModel(;timestepper=:RungeKutta3,
                              advection=WENO5(grid=grid),
                              grid=grid,
                              gravitational_acceleration=g,
                              coriolis=FPlane(f=f),
                              forcing=(h=convec_forcing,u = u_forcing, v = v_forcing)
                              )


    #Build background state and perturbations

    function h̄(x, y, z)
#        if x == grid.xᶜᵃᵃ[10] && y == grid.yᵃᶜᵃ[10]
#            return 39 
#        else
            return Lz  + 4rand()#-  Δη * exp(-((x-Lx/2)^2/(0.1*Lx)^2 + (y-Ly/2)^2/(0.1*Ly)^2 ))
#        end
    end
    ū(x, y, z) = 0 #U * sech(y)^2
    ω̄(x, y, z) = 0 #2 * U * sech(y)^2 * tanh(y)

    small_amplitude = 1e-4

    uⁱ(x, y, z) = ū(x, y, z) #+ small_amplitude * exp(-y^2) * randn()
    hⁱ(x, y, z) = h̄(x, y, z)
    uhⁱ(x, y, z) = 0.0 #uⁱ(x, y, z) * hⁱ(x, y, z)

    uh, vh, h = model.solution

    ## Build velocities
    u = uh / h
    v = vh / h

    #set initial condition

    set!(model, uh = uhⁱ, h = h̄)


    #Create the simulation
    #simulation = Simulation(model, Δt = 1e-2, stop_time = 150)
    simulation = Simulation(model, Δt = 20, stop_time = 86400*15)

    function update_convective_helper_arrays(sim, parameters = parameters)
        p = parameters
        m = sim.model
        update_convective_events!(m.architecture,p.isconvecting,p.convection_triggered_time,m.solution.h,
                                  m.clock.time,p.τ_c,p.h_c,m.grid.Nx,m.grid.Ny)
        Oceananigans.BoundaryConditions.fill_halo_regions!(p.isconvecting, m.architecture)
        Oceananigans.BoundaryConditions.fill_halo_regions!(p.convection_triggered_time, m.architecture)


    end

    function progress(sim)
        m = sim.model
        @info(@sprintf("Iter: %d, time: %.1f, Δt: %.3f, max|h|: %.2f, min|h|: %.2f",
                       m.clock.iteration, m.clock.time,
                       sim.Δt, maximum(abs, m.solution.h), minimumm(abs, m.solution.h) ))
        
    end
    filename = savename(params) ## produces name too long, do not use
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(500))
    simulation.callbacks[:update_convective_helper_arrays] = Callback(update_convective_helper_arrays, IterationInterval(1))
    simulation.output_writers[:fields] =
        NetCDFOutputWriter(
            model,
            (h = h , v = v , u = u, isconvecting = isconvecting),
            filepath = joinpath(ENV["SCRATCH"], "15daysim.nc"),
            schedule = TimeInterval(1),
            mode = "c")

    run!(simulation)
end #runsimulation
