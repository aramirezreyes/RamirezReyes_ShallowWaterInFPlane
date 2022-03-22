Lx, Ly, Lz = 1.0e2, 1.0e2, 40
Nx, Ny = 10, 10

#grid = RegularRectilinearGrid(size = (Nx, Ny, 1),
#                            x = (0, Lx), y = (0, Ly), z = (0, Lz),
#                              topology = (Periodic, Periodic, Bounded))
h_c = 40.0
τ_c = 5.0
h = 40.5ones(10,10)
isconvecting = zeros(Bool,10,10)
h[5,5] = 39.0
convection_triggered_time = zeros(10,10)
t1 = 10.
t2 = t1 + τ_c/2
t3 = t1 + τ_c + 0.5
R2 = 20
q0 = 5
struct TestClock{T}
    time :: T
end
clock = TestClock(t2)
@testset "Detect convecting events" begin
    @test convection_triggered_time[5,5] == 0.0 #Make it not convecting
    update_convective_events!(CPU(),isconvecting,convection_triggered_time,h,t1,τ_c,h_c,0,0) #Should update 5,5
    @test only(findall(isconvecting)) == CartesianIndex(5,5)
    @test convection_triggered_time[5,5] == t1 # check that time got set
    update_convective_events!(CPU(),isconvecting,convection_triggered_time,h,t2,τ_c,h_c,0,0) #should keep everying the same
    @test only(findall(isconvecting)) == CartesianIndex(5,5) #should still convect
    @test convection_triggered_time[5,5] == t1 #time not modified
    update_convective_events!(CPU(),isconvecting,convection_triggered_time,h,t3,τ_c,h_c,0,0) #should start new convective event because height is still below the threshold
    @test only(findall(==(1.0),isconvecting)) == CartesianIndex(5, 5)
    @test t3 ≈ convection_triggered_time[5,5] #new convective event started
    h .= 40.5
    convection_triggered_time[5,5] = 0.0 #let's imagine it started convecting a long time ago
    update_convective_events!(CPU(),isconvecting,convection_triggered_time,h,t3,τ_c,h_c,0,0) #should convect because height is still below the threshold
    @test isempty(findall(isconvecting)) #nothing convecting
    @test t3 - convection_triggered_time[5,5] > τ_c #time since last convection start still greater than what is alloweble
end


#@testset "Heating function" begin
#    isconvecting = zeros(Bool,10,10)
#    h[5,5] = 39.0
#    convection_triggered_time = zeros(10,10)
#    @test RamirezReyes_ShallowWaterInFPlane.heat_at_point(0,0,1,grid,clock,h,τ_c,h_c,isconvecting,convection_triggered_time,R2,q0) == 0.0
#
#end

