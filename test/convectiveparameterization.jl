
h_c = 40.0
τ_c = 5.0
h = 40.5ones(10,10)
isconvecting = zeros(Bool,10,10)
h[5,5] = 39.0
convection_triggered_time = zeros(10,10)
t1 = 10.
t2 = t1 + τ_c/2
t3 = t1 + τ_c + 0.5

@testset "Detect convecting events" begin
    @test convection_triggered_time[5,5] == 0.0
    update_convective_events!(isconvecting,convection_triggered_time,h,t1,τ_c,h_c)
    @test only(findall(isconvecting)) == CartesianIndex(5,5)
    @test convection_triggered_time[5,5] == t1
    update_convective_events!(isconvecting,convection_triggered_time,h,t2,τ_c,h_c)
    @test only(findall(isconvecting)) == CartesianIndex(5,5)
    @test convection_triggered_time[5,5] == t1
    update_convective_events!(isconvecting,convection_triggered_time,h,t3,τ_c,h_c)
    @test isempty(findall(isconvecting)) 
    @test convection_triggered_time[5,5] == 0
end


