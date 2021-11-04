using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"

using Oceananigans
using Oceananigans.Models: ShallowWaterModel
using Printf
using OffsetArrays
using CUDA
#using ProfileView
using ImageFiltering: padarray, Pad

include(joinpath(@__DIR__,"../src/convectiveparameterization.jl"))
include(joinpath(@__DIR__,"../src/arrayutils.jl"))


#Setup domain

Lx, Ly = 1.0e5, 1.0e5
const Lz = 45
Nx, Ny = 128, 128

grid = RegularRectilinearGrid(size = (Nx, Ny),
                            x = (0, Lx), y = (0, Ly),
                              topology = (Periodic, Periodic, Flat))


#Setup physicsl parameters

 f = 0 #5e-4
 g = 9.8# 9.8
 U = 50.0
Δη =  5.5#f * U / g

τ_c = 50
h_c = 40
heating_amplitude    = 1.0e9 #originally 9 for heating, -8 for cooling
convective_radius    = 20000.0
numelements_to_traverse_x = Int(convective_radius ÷ grid.Δx)
numelements_to_traverse_y = Int(convective_radius ÷ grid.Δy)

## convecting parameterization arrays
isconvecting = OffsetArray(CUDA.fill(false,Nx + 2numelements_to_traverse_x, Ny + 2numelements_to_traverse_y),
                           (-numelements_to_traverse_x + 1) : (Nx + numelements_to_traverse_x),
                           (-numelements_to_traverse_y + 1) : (Ny + numelements_to_traverse_y)) 
convection_triggered_time = OffsetArray(CUDA.fill(0.0f0,Nx + 2numelements_to_traverse_x, Ny + 2numelements_to_traverse_y),
                                        (-numelements_to_traverse_x + 1) : (Nx + numelements_to_traverse_x),
                                        (-numelements_to_traverse_y + 1) : (Ny + numelements_to_traverse_y)) 
#isconvecting[10,10] = true
#convection_triggered_time[10,10] = 0.0
#build forcing
convec_forcing = Forcing(model_forcing,discrete_form=true,parameters=(τ_c = τ_c,
                                                                      h_c = h_c,
                                                                      isconvecting = isconvecting,
                                                                      convection_triggered_time = convection_triggered_time,
                                                                      R = convective_radius,
                                                                      q0 = heating_amplitude,
                                                                      numelements_to_traverse_x = numelements_to_traverse_y,
                                                                      numelements_to_traverse_y = numelements_to_traverse_y
                                                                      ))


## Build the model

model = ShallowWaterModel(architecture = GPU(),
    timestepper=:RungeKutta3,
    advection=WENO5(),
    grid=grid,
    gravitational_acceleration=g,
    coriolis=FPlane(f=f),
    forcing=(h=convec_forcing,)
)


#Build background state and perturbations

#h̄(x, y, z) = model.grid.Lz +  Δη * tanh(y/Ly)*tanh(x/Lx)
#h̄(x, y, z) = model.grid.Lz -  Δη * exp(-((x-Lx/2)^2/(0.1*Lx)^2 + (y-Ly/2)^2/(0.1*Ly)^2 ))
function h̄(x, y, z)
    if x == grid.xC[10] && y == grid.yC[10]
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

set!(model, uh = uhⁱ, h = hⁱ)






function progress(sim)
    p = sim.parameters 
    update_convective_events_gpu!(p.isconvecting,p.convection_triggered_time,p.h.data,sim.model.clock.time,p.τ_c,p.h_c,
                              sim.model.grid.Nx,sim.model.grid.Ny, p.nghosts_x,p.nghosts_y)
    if sim.model.clock.iteration % 20 == 0
        @info(@sprintf("Iter: %d, time: %.1f, Δt: %.3f, max|h|: %.2f",
                       sim.model.clock.iteration, sim.model.clock.time,
                       sim.Δt, maximum(abs, sim.parameters.h.data.parent)))
    end
end
uh, vh, h = model.solution

        u = ComputedField(uh / h)
        v = ComputedField(vh / h)
        ω = ComputedField(∂x(v) - ∂y(u))
ω_pert = ComputedField(ω - ω̄)

#Create the simulation
#simulation = Simulation(model, Δt = 1e-2, stop_time = 150)
simulation = Simulation(model, Δt = 1, stop_time = 600, progress= progress, parameters = (h = h, isconvecting = isconvecting, convection_triggered_time = convection_triggered_time, τ_c, h_c, nghosts_x = numelements_to_traverse_x, nghosts_y = numelements_to_traverse_y ))


#prepare output files
using LinearAlgebra: norm

function perturbation_norm(model)
    compute!(h)
    return norm(interior(h))
end

outputs = (ω_total = ω, ω_pert = ω_pert)

simulation.output_writers[:fields] =
    NetCDFOutputWriter(
        model,
        (ω = ω, ω_pert = ω_pert, h = h , v = v , u = u),
          filepath = joinpath(@__DIR__, "../data/shallow_water_example.nc"),
          schedule = TimeInterval(1),
        mode = "c")

# simulation.output_writers[:growth] =
#     NetCDFOutputWriter(
#         model,
#         (perturbation_norm = perturbation_norm,),
#           filepath = joinpath(@__DIR__, "perturbation_norm_shallow_water.nc"),
#           schedule = IterationInterval(1),
#         dimensions = (perturbation_norm=(),),
#         mode = "c")

run!(simulation)
