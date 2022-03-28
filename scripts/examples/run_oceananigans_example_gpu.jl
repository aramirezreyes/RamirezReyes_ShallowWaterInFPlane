using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"

using RamirezReyes_ShallowWaterInFPlane: update_convective_events!, model_forcing, u_damping, v_damping, fill_heating_stencil!
using Oceananigans
using Printf
#using ProfileView

architecture = GPU()
#Setup physicsl parameters

f = 0.0#5e-4 #5e-4
g = 9.8# 9.8
U = 50.0
Δη =  5.5#f * U / g

τ_c = 50
h_c = 40
heating_amplitude    = 1.0e9#1.0e9 #originally 9 for heating, -8 for cooling
radiative_cooling_rate = 1.0e-8
convective_radius    = 20000.0
relaxation_parameter = 1.0/3600.0
relaxation_height = 38
#Setup domain

Lx, Ly = 1.0e5, 1.0e5
const Lz = 45
Nx, Ny = 128, 128

grid_spacing_x = Lx ÷ Nx #These two need to be equal for x and y!
grid_spacing_y = Ly ÷ Ny

numelements_to_traverse = Int(convective_radius ÷ grid_spacing_x)
halo_indices = numelements_to_traverse-1

grid = RectilinearGrid(size = (Nx, Ny),
                            x = (0, Lx), y = (0, Ly),
                              topology = (Periodic, Periodic, Flat), halo = (max(numelements_to_traverse,3), max(numelements_to_traverse,3)))

@info "Built grid successfully"
isconvecting  = CenterField(grid,Bool)
convection_triggered_time  = CenterField(grid)
## Will create heating stencil with the spatial component
q_stencil = CenterField(grid,Float64; indices=(-halo_indices:halo_indices,-halo_indices:halo_indices,:))

fill_heating_stencil!(q_stencil,heating_amplitude,grid_spacing_x,convective_radius^2)

parameters = (; isconvecting = isconvecting, convection_triggered_time, τ_c, h_c, nghosts = numelements_to_traverse - 1, radiative_cooling_rate , q0 = heating_amplitude, R = convective_radius, relaxation_parameter, relaxation_height, Δx2 = grid_spacing_x^2, Δy2 = grid_spacing_y^2, heating_stencil = q_stencil)


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

#h̄(x, y, z) = Lz +  Δη * tanh(y/Ly)*tanh(x/Lx)
#h̄(x, y, z) = model.grid.Lz -  Δη * exp(-((x-Lx/2)^2/(0.1*Lx)^2 + (y-Ly/2)^2/(0.1*Ly)^2 ))
function h̄(x, y, z)
    if x == grid.xᶜᵃᵃ[10] && y == grid.yᵃᶜᵃ[10]
        return 39
    else
        return Lz #-  Δη * exp(-((x-Lx/2)^2/(0.1*Lx)^2 + (y-Ly/2)^2/(0.1*Ly)^2 ))
    end
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

## Build and compute mean vorticity discretely
ω = Field(∂x(v) - ∂y(u))
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
simulation = Simulation(model, Δt = 1, stop_time = 600)

function update_convective_helper_arrays(sim, parameters = parameters)
    p = parameters
    #@info "Go run update_...!"
    m = sim.model
    update_convective_events!(m.architecture,p.isconvecting,p.convection_triggered_time,m.solution.h,
                              m.clock.time,p.τ_c,p.h_c,m.grid.Nx,m.grid.Ny)
    Oceananigans.BoundaryConditions.fill_halo_regions!(p.isconvecting, m.architecture)
    Oceananigans.BoundaryConditions.fill_halo_regions!(p.convection_triggered_time, m.architecture)
    #display(Array(p.convection_triggered_time.data.parent))
    #display(Array(p.isconvecting.data.parent))

end

function progress(sim)
    m = sim.model
    @info(@sprintf("Iter: %d, time: %.1f, Δt: %.3f, max|h|: %.2f",
                   m.clock.iteration, m.clock.time,
                   sim.Δt, maximum(abs, m.solution.h)))
    
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
simulation.callbacks[:update_convective_helper_arrays] = Callback(update_convective_helper_arrays, IterationInterval(1))
#prepare output files
 simulation.output_writers[:fields] =
     NetCDFOutputWriter(
         model,
         (h = h , v = v , u = u, isconvecting = isconvecting),
           filepath = joinpath(datadir(), "shallow_water_example_gpu.nc"),
           schedule = TimeInterval(1),
         mode = "c")

#@profview run!(simulation)
run!(simulation)
