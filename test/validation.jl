include("convective_validation_src.jl")

@testset "Check amounts of mass inyected/substracted" begin

initial_state, final_state, heating_amplitude, Δx, Δy = validate("CPU", boundary_layer = false)
dif = sum(interior(final_state.h)) - sum(interior(initial_state.h))
total_mass_inyected  = 3 * dif * Δx * Δy
@test isapprox(total_mass_inyected, heating_amplitude, rtol=0.01)

initial_state, final_state, heating_amplitude, Δx, Δy = validate("CPU", boundary_layer = true)
dif = sum(interior(initial_state.h)) - sum(interior(final_state.h))
total_mass_inyected  = 3 * dif * Δx * Δy
@test isapprox(total_mass_inyected, heating_amplitude, rtol=0.01)

end

if CUDA.functional()

    @testset "Check amounts of mass inyected/substracted on the GPU" begin

        initial_state, final_state, heating_amplitude, Δx, Δy = validate("GPU", boundary_layer = false)
        dif = sum(interior(final_state.h)) - sum(interior(initial_state.h))
        total_mass_inyected  = 3 * dif * Δx * Δy
        @test isapprox(total_mass_inyected, heating_amplitude, rtol=0.01)
        
        initial_state, final_state, heating_amplitude, Δx, Δy = validate("GPU", boundary_layer = true)
        dif = sum(interior(initial_state.h)) - sum(interior(final_state.h))
        total_mass_inyected  = 3 * dif * Δx * Δy
        @test isapprox(total_mass_inyected, heating_amplitude, rtol=0.01)
        
        end

end