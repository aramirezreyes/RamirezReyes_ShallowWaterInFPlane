### This file is temporary, I will define constructors for checkpointers that work on the ShallowWaterModel.
# This is type piracy and must go into a PR in Oceananigans
import Oceananigans.Fields: set!

using Oceananigans: fields, prognostic_fields
import Oceananigans.OutputWriters: has_reference, set_time_stepper!

"""
    Checkpointer(model; schedule,
                 dir = ".",
                 prefix = "checkpoint",
                 overwrite_existing = false,
                 cleanup = false,
                 additional_kwargs...)

Construct a `Checkpointer` that checkpoints the model to a JLD2 file on `schedule.`
The `model.clock.iteration` is included in the filename to distinguish between multiple checkpoint files.

To restart or "pickup" a model from a checkpoint, specify `pickup=true` when calling `run!`, ensuring
that the checkpoint file is the current working directory. See 

```julia
help> run!
```
for more details.

Note that extra model `properties` can be safely specified, but removing crucial properties
such as `:velocities` will make restoring from the checkpoint impossible.

The checkpointer attempts to serialize as much of the model to disk as possible,
but functions or objects containing functions cannot be serialized at this time.

Keyword arguments
=================
- `schedule` (required): Schedule that determines when to checkpoint.

- `dir`: Directory to save output to. Default: "." (current working directory).

- `prefix`: Descriptive filename prefixed to all output files. Default: "checkpoint".

- `overwrite_existing`: Remove existing files if their filenames conflict. Default: `false`.

- `verbose`: Log what the output writer is doing with statistics on compute/write times
             and file sizes. Default: `false`.

- `cleanup`: Previous checkpoint files will be deleted once a new checkpoint file is written.
             Default: `false`.

- `properties`: List of model properties to checkpoint. This list must contain
                `[:grid, :architecture, :timestepper, :particles]`.
                Default: [:architecture, :grid, :clock, :coriolis, :buoyancy, :closure,
                          :velocities, :tracers, :timestepper, :particles]
"""
function Oceananigans.Checkpointer(model :: ShallowWaterModel; schedule,
                      dir = ".",
                      prefix = "checkpoint",
                      overwrite_existing = false,
                      verbose = false,
                      cleanup = false,
                      properties = [:architecture, :grid, :clock, :coriolis,
                                    :buoyancy, :closure, :timestepper, :particles])

    # Certain properties are required for `restore_from_checkpoint` to work.
    required_properties = (:grid, :architecture, :timestepper, :particles)

    for rp in required_properties
        if rp ∉ properties && rp ∈ propertynames(model)
            @warn "$rp is required for checkpointing. It will be added to checkpointed properties"
            push!(properties, rp)
        end
    end

    for p in properties
        p isa Symbol || error("Property $p to be checkpointed must be a Symbol.")
        p ∉ propertynames(model) && error("Cannot checkpoint $p, it is not a model property!")

        if (p ∉ required_properties) && has_reference(Function, getproperty(model, p))
            @warn "model.$p contains a function somewhere in its hierarchy and will not be checkpointed."
            filter!(e -> e != p, properties)
        end
    end

    mkpath(dir)

    return Checkpointer(schedule, dir, prefix, properties, overwrite_existing, verbose, cleanup)
end

"""
    set!(model, filepath::AbstractString)

Set data in `model.velocities`, `model.tracers`, `model.timestepper.Gⁿ`, and
`model.timestepper.G⁻` to checkpointed data stored at `filepath`.
"""
function Oceananigans.set!(model :: ShallowWaterModel, filepath::AbstractString)

    jldopen(filepath, "r") do file

        # Validate the grid
        checkpointed_grid = file["grid"]

        model.grid == checkpointed_grid ||
             error("The grid associated with $filepath and model.grid are not the same!")

        model_fields = prognostic_fields(model)

        for name in propertynames(model_fields)
            try
                parent_data = file["$name/data"]
                model_field = model_fields[name]
                copyto!(model_field.data.parent, parent_data)
            catch
                @warn "Could not retore $name from checkpoint."
            end
        end

        set_time_stepper!(model.timestepper, file, model_fields)

        if hasproperty(model, :particles)
            if !isnothing(model.particles)
                copyto!(model.particles.properties, file["particles"])
            end
        end

        checkpointed_clock = file["clock"]

        # Update model clock
        model.clock.iteration = checkpointed_clock.iteration
        model.clock.time = checkpointed_clock.time
    end

    return nothing
end