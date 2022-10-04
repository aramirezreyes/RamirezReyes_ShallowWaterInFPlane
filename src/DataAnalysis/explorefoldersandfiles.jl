function list_ncfiles(data_dir; filter_length = false, minimum_length = 0, filter_nconvectivepoints = false, min_convectivepoints = 90)
    file_list = readdir(data_dir)
    filter!(filename ->  endswith(filename,".nc"),file_list)
    
    n_original_files = length(file_list)
    n_discarded_by_length = 0
    n_discarded_by_nconvectivepoints = 0
    if filter_length
        filter!(file -> (read_last_time(joinpath(data_dir,file)) > minimum_length), file_list)
        n_discarded_by_length = n_original_files - length(file_list)
    end
    if filter_nconvectivepoints
        filter!(file -> (count_last_convecting(joinpath(data_dir,file)) > 90), file_list)
        n_discarded_by_nconvectivepoints = n_original_files - n_discarded_by_length - length(file_list)
    end
    @info "I found: $n_original_files files, i discarded: $n_discarded_by_length,  by length and, $n_discarded_by_nconvectivepoints  by number of convective points"
    return file_list
end


function read_last_time(file)
    try
        NetCDF.open(file) do ds
            maximum(ds["time"])./86400
        end
    catch e
        @info "Failed reading: $file"
        0.0
    end
end
function read_last_frame_vort(file)
    NetCDF.open(file) do ds
        ω_mean = mean(ds["ω"][:,:,1,end])
        ω = ds["ω"][:,:,1,end]
        ω_perturb = (ω .- ω_mean)
        ω_perturb, ds["xC"][:],ds["yC"][:],ds["time"][end]
    end
end
function read_last_frame_generic(file, var)
    NetCDF.open(file) do ds
        ds[var][:,:,1,end], ds["xC"][:],ds["yC"][:], ds["time"][end]
    end
end

function count_last_convecting(file)
    try
        NetCDF.open(file) do ds
            count(==(1),ds["isconvecting"][:,:,1,end])
        end
    catch
        @info "Failed reading :$file"
        0
    end
end

function read_sim_length(file)
    NetCDF.open(file) do ds
        length(ds["time"])
    end
end

function maximum_abs_vorticity_perturb(file)
    try
        NetCDF.open(file) do ds
            ω_mean = mean(ds["ω"][:,:,1,end])
            ω = ds["ω"][:,:,1,end]
            filtered = imfilter(ω .- ω_mean, Kernel.gaussian(4));    
            ω_perturb = maximum(abs,filtered)
        end
    catch e
        @info "Failed reading :$file"
        rethrow(e)
    end
end

function maximum_sp(file)
    try
        NetCDF.open(file) do ds
            sp = ds["sp"][:,:,1,end]  
            max_sp = maximum(sp)
        end
    catch e
        @info "Failed reading :$file"
        rethrow(e)
    end
end

function parse_params_from_filenames(file_list)
    nfiles = length(file_list)
    without_ext = splitext(file_list[1])[1]
    matches = split(without_ext,'_')
    param_names = [split.(match,'=')[1] for match in matches]
    pushfirst!(param_names,"Index")
    push!(param_names,"filename")
    ncols = length(matches)+2
    param_table = Array{Union{Float64,String},2}(undef,nfiles,ncols)
    for i in 1:nfiles
        param_table[i,1] = Float64(i)
        filename = file_list[i]
        without_ext = splitext(filename)[1]
        matches = split(without_ext,'_')
        param_values = [split.(match,'=')[2] for match in matches]
        j = 2
        for value in param_values
            parsedvalue = tryparse(Float64,value)
            if isnothing(parsedvalue)
                parsedvalue = tryparse(Bool,value)
            end
            param_table[i,j] = Float64(parsedvalue)
            j = j + 1
        end
        param_table[i,end] = filename
    end
    param_table,param_names
end
