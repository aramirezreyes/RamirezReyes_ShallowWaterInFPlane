#!/bin/bash -l

/usr/bin/time /global/u2/a/aramreye/Software/julia-1.7.1/bin/julia -t 1 --project=@. -e 'using RamirezReyes_ShallowWaterInFPlane; run_shallow_simulation_debug("CPU")'


/usr/bin/time /global/u2/a/aramreye/Software/julia-1.7.1/bin/julia -t 2 --project=@. -e 'using RamirezReyes_ShallowWaterInFPlane; run_shallow_simulation_debug("CPU")'


/usr/bin/time /global/u2/a/aramreye/Software/julia-1.7.1/bin/julia -t 4 --project=@. -e 'using RamirezReyes_ShallowWaterInFPlane; run_shallow_simulation_debug("CPU")'


/usr/bin/time /global/u2/a/aramreye/Software/julia-1.7.1/bin/julia -t 8 --project=@. -e 'using RamirezReyes_ShallowWaterInFPlane; run_shallow_simulation_debug("CPU")'


/usr/bin/time /global/u2/a/aramreye/Software/julia-1.7.1/bin/julia -t 16 --project=@. -e 'using RamirezReyes_ShallowWaterInFPlane; run_shallow_simulation_debug("CPU")'


/usr/bin/time /global/u2/a/aramreye/Software/julia-1.7.1/bin/julia -t 32 --project=@. -e 'using RamirezReyes_ShallowWaterInFPlane; run_shallow_simulation_debug("CPU")'


/usr/bin/time /global/u2/a/aramreye/Software/julia-1.7.1/bin/julia -t 64 --project=@. -e 'using RamirezReyes_ShallowWaterInFPlane; run_shallow_simulation_debug("CPU")'




