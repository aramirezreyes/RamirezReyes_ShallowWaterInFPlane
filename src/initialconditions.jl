function create_initialization_functions(grid, parameters_dict)
    p = parameters_dict
    initialization_style = p["initialization_style"]
    if initialization_style == "rand"
        return uhⁱ, uhⁱ, h0_rand(p)
    elseif initialization_style == "one_convecting_point"
        return uhⁱ, uhⁱ, h0_one_convecting_point(grid, p)
    elseif initialization_style == "rankine_vortex"
        uh_f, vh_f = rankine_vortex(p)
        return uh_f, vh_f, h0_rand(p)
    elseif initialization_style == "gaussian"
        return uhⁱ, uhⁱ, h0_gaussian(grid, p)
    else
        error("Intialization style must be either \"rand\" or \"one_convecting_point\"")
    end
end


uhⁱ(x, y, z) = 0.0

function h0_rand(p)
    function h0(x, y, z)
        return p["Lz"] + p["initialization_amplitude"] * rand()
    end
end

"""
    h0_one_convecting_point(x,y,z)
"""
function h0_one_convecting_point(grid, parameters_dict)
    p = parameters_dict
    isbl = p["boundary_layer"]
    τc = p["convection_critical_height"]
    Δh = p["initialization_amplitude"]
    function h0(x, y, z)
        if x == grid.xᶜᵃᵃ[p["Nx"]÷2] && y == grid.yᵃᶜᵃ[p["Ny"]÷2]
            isbl ? τc + Δh : τc - Δh
        else
            isbl ? τc - 1.0 : τc + 1.0
        end
    end
end

"""
    rankine_vortex(center::CartesianIndex{2}, radius::Number, amplitude::Number, h0::Number)
Return two functions of x,y,z for the u and v components of the mass flux `uh` and `vh`. The vortex is centered in `center` and has a radius of `radius`.
"""
function rankine_vortex(parameters_dict)
    p = parameters_dict
    center = p["rankine_center"]
    radius = p["rankine_radius"]
    amplitude = p["rankine_amplitude"]
    h0 = p["Lz"]

    function vortex_maker_u(x, y, z)
        amplitude2 = amplitude
        x_relative = (x - center[1])
        y_relative = (y - center[2])
        dist_from_center = hypot(x_relative, y_relative)
        theta = atan(y_relative, x_relative)
        if dist_from_center <= radius
            return -amplitude2 * h0 * dist_from_center / radius^2 *
                   dist_from_center *
                   sin(theta)
        else
            return -amplitude2 * h0 * sin(theta) * radius / dist_from_center
        end
    end


    function vortex_maker_v(x, y, z)
        amplitude2 = amplitude
        x_relative = (x - center[1])
        y_relative = (y - center[2])
        dist_from_center = hypot(x_relative, y_relative)
        theta = atan(y_relative, x_relative)
        if dist_from_center <= radius
            return amplitude2 * h0 * dist_from_center / radius^2 *
                   dist_from_center *
                   cos(theta)
        else
            return amplitude2 * h0 * cos(theta) * radius / dist_from_center
        end
    end

    return (vortex_maker_u, vortex_maker_v)
end

"""
following https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function

"""
function h0_gaussian(grid, parameters_dict)
    
    p = parameters_dict
    isbl = p["boundary_layer"]
    Δh = p["initialization_amplitude"]
    sigmax = p["gaussian_sigma_x"]
    sigmay = p["gaussian_sigma_y"] 
    rot = p["gaussian_rotation"]
    h0 = p["Lz"]

    a = ((cos(rot)^2) / (2*sigmax^2)) + ((sin(rot)^2) / (2*sigmay^2))
    b = -((sin(2*rot)) / (4*sigmax^2)) + ((sin(2*rot)) / (4*sigmay^2))
    c = ((sin(rot)^2) / (2*sigmax^2)) + ((cos(rot)^2) / (2*sigmay^2))

    function h(x, y, z)
        x0 = grid.xᶜᵃᵃ[p["Nx"]÷2] 
        y0 = grid.yᵃᶜᵃ[p["Ny"]÷2]
        h0 +Δh*exp(-(a*(x - x0)^2 - 2*b*(x - x0)*(y - y0) + c*(y - y0)^2))    
    end
end
