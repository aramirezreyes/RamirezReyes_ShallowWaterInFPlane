using Oceananigans
using ProfileView
using BenchmarkTools

include("../src/convectiveparameterization.jl")

Lx, Ly, Lz = 1.0e2, 1.0e2, 40
Nx, Ny = 10, 10
grid = RegularRectilinearGrid(size = (Nx, Ny, 1),
                            x = (0, Lx), y = (0, Ly), z = (0, Lz),
                              topology = (Periodic, Periodic, Bounded))
h_c = 40.0
τ_c = 5.0
h = 40.5ones(10,10)
t1 = 10.
t2 = t1 + τ_c/2
t3 = t1 + τ_c + 0.5
isconvecting = zeros(Bool,10,10)
convection_triggered_time = zeros(10,10)
h[5,5] = 39.0
isconvecting[5,5] = true
convection_triggered_time[5,5] = t1
R2 = 20
q0 = 5
struct TestClock{T}
    time :: T
end
clock = TestClock(t2)

function profile_test(n)
    forcing = 0.00
    Lx, Ly, Lz = 1.0e2, 1.0e2, 40
    Nx, Ny = 10, 10
    grid = RegularRectilinearGrid(size = (Nx, Ny, 1),
                              x = (0, Lx), y = (0, Ly), z = (0, Lz),
                              topology = (Periodic, Periodic, Bounded))
    h_c = 40.0
    τ_c = 5.0
    h = 40.5ones(10,10)
    t1 = 10.
    t2 = t1 + τ_c/2
    t3 = t1 + τ_c + 0.5
    isconvecting = zeros(Bool,10,10)
    convection_triggered_time = zeros(10,10)
    h[5,5] = 39.0
    isconvecting[5,5] = true
    convection_triggered_time[5,5] = t1
    R2 = 200
    q0 = 5
    clock = TestClock(t2)
    for i in 1:n
        forcing += heat_at_point(5,5,i,clock,h,τ_c,h_c,isconvecting,convection_triggered_time,R2,q0,grid.xC.parent,grid.yC.parent,grid.Δx,grid.Δy,Nx,Ny)
    end
    forcing
end

