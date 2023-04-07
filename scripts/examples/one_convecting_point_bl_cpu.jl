# Goal of this script is validate the convective parameterization. Here, for the boundary layer formulation.
# the first validation consists in measuring the mass added to the system by a convective event. This is based on Integrating equation 4 of
#Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
# In which the mass added to the system after a convective event is equal to
# ΔH = 1/3 q0 

parameters = Dict(
    "architecture" => "CPU",
    "output_filename" => "one_convecting_point_bl_cpu",
    "f" => 0.0,
    "g" => 0.0,
    "convection_timescale" => 1000.0,
    "convection_critical_height" => 40.0,
    "heating_amplitude" => 1.0e9,
    "large_scale_forcing" => 0.0,
    "convective_radius" => 10_000,
    "relaxation_parameter" => 0.0,
    "relaxation_height" => nothing,
    "Lx" => 30000,
    "Ly" => 30000,
    "Lz" => 40.0,
    "Nx" => 100,
    "Ny" => 100,
    "boundary_layer" => true,
    "initialization_style" => "one_convecting_point",
    "initialization_amplitude" => 1,
    "simulation_length_in_days" => 2000 / 86400,
    "output_interval_in_seconds" => 20,
    "timestep_in_seconds" => 20,
    "restart" => false,
    "checkpoint_interval_in_seconds" => Inf,
)

using RamirezReyes_ShallowWaterInFPlane

run_shallow_simulation(parameters)
