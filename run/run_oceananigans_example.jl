using Oceananigans
using Oceananigans.Models: ShallowWaterModel
using Printf
include(joinpath(@__DIR__,"../src/convectiveparameterization.jl"))


#Setup domain

Lx, Ly, Lz = 1.0e6, 1.0e6, 45
Nx, Ny = 128, 128

grid = RegularRectilinearGrid(size = (Nx, Ny, 1),
                            x = (0, Lx), y = (0, Ly), z = (0, Lz),
                              topology = (Periodic, Periodic, Bounded))

isconvecting = Array{Bool,2}(undef,Nx,Ny)
timeconvecting = Array{Float32,2}(undef,Nx,Ny)


#Setup physicsl parameters

 f = 5e-4
 g = 9.8
 U = 50.0
Δη =  5.5#f * U / g


## Build the model

model = ShallowWaterModel(
    timestepper=:RungeKutta3,
    advection=WENO5(),
    grid=grid,
    gravitational_acceleration=g,
    coriolis=FPlane(f=f),
)

#Build background state and perturbations

#h̄(x, y, z) = model.grid.Lz +  Δη * tanh(y/Ly)*tanh(x/Lx)
h̄(x, y, z) = model.grid.Lz +  Δη * exp(-((x-Lx/2)^2/(0.1*Lx)^2 + (y-Ly/2)^2/(0.1*Ly)^2 ))
#h['g'] = H - amp*np.exp(- ((x-0.5*Lx)**2/(0.1*Lx)**2 + (y-0.5*Ly)**2/(0.1*Ly)**2 ))
ū(x, y, z) = 0 #U * sech(y)^2
ω̄(x, y, z) = 0 #2 * U * sech(y)^2 * tanh(y)

 small_amplitude = 1e-4

 uⁱ(x, y, z) = ū(x, y, z) #+ small_amplitude * exp(-y^2) * randn()
 hⁱ(x, y, z) = h̄(x, y, z)
uhⁱ(x, y, z) = 0.0 #uⁱ(x, y, z) * hⁱ(x, y, z)

set!(model, uh = uhⁱ, h = hⁱ)
progress(sim) = @info(@sprintf("Iter: %d, time: %.1f, Δt: %.3f, max|u|: %.2f",
                               sim.model.clock.iteration, sim.model.clock.time,
                                   sim.Δt, maximum(abs, u.data.parent)))
uh, vh, h = model.solution

        u = ComputedField(uh / h)
        v = ComputedField(vh / h)
        ω = ComputedField(∂x(v) - ∂y(u))
ω_pert = ComputedField(ω - ω̄)

#Create the simulation
#simulation = Simulation(model, Δt = 1e-2, stop_time = 150)
simulation = Simulation(model, Δt = 50, stop_time = 86400, progress= progress)


#prepare output files
using LinearAlgebra: norm

function perturbation_norm(model)
    compute!(v)
    return norm(interior(v))
end

outputs = (ω_total = ω, ω_pert = ω_pert)

simulation.output_writers[:fields] =
    NetCDFOutputWriter(
        model,
        (ω = ω, ω_pert = ω_pert, h = h , v = v , u = u),
          filepath = joinpath(@__DIR__, "shallow_water_Bickley_jet.nc"),
          schedule = TimeInterval(3600),
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
