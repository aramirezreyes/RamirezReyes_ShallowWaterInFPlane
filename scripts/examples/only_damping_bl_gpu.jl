# Goal of this script is validate the convective parameterization that we will use.
# the first validation consists in measuring the mass added to the system by a convective event. This is based on Integrating equation 4 of
#Yang, D., and A. P. Ingersoll, 2013: Triggered Convection, Gravity Waves, and the MJO: A Shallow-Water Model. J. Atmos. Sci., 70, 2476–2486, https://doi.org/10.1175/JAS-D-12-0255.1.
# In which the mass added to the system after a convective event is equal to
# ΔH = 1/3 q0 

parameters = Dict(
    "architecture" => "GPU",
    "output_filename" => "only_damping_bl_GPU",
    "f" => 0.0,# coriolis parameter5e-4 #5e-4
    "g" => 0.0,#gravitational acceleration 9.8
    "convection_timescale" => 1000.0, #convective time scale
    "convection_critical_height" => 40.0, #convection_triggering height
    "heating_amplitude" => 0.0,#convective heating 1.0e9 #originally 9 for heating, -8 for cooling convective heating
    "large_scale_forcing" => 0.0,#radiative cooling rate
    "convective_radius" => 10_000, #convective radius
    "relaxation_parameter" => 1/1000.0, #relaxation timescale
    "relaxation_height" => nothing, #height to relax to, if nothing, h relax to its mean
    "Lx" => 30000, #domain
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
