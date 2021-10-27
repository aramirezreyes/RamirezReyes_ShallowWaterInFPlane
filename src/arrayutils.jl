@inline mirror_index(i,n,N) = mod1(i + n,N)
"""
    fill_ghosts!(array,halo_x,halo_y)
This function will fill ghosts cells on the sides of an array to represent periodic boundary conditions
TODO: Needs to fill the corners as well.
"""
function fill_ghosts!(array,Nx,Ny,halo_x,halo_y)
    for j in (-halo_y + 1) : 0
        mirror_index = mod1(j,Ny)
        for i in 1:Nx
            array[i,j] = array[i,mirror_index]
        end
    end
    for j in (Ny + 1) : (Ny + halo_y)
        mirror_index = mod1(j,Ny)
        for i in 1:Nx
            array[i,j] = array[i,mirror_index]
        end
    end

    for i in (-halo_x + 1) : 0
        mirror_index = mod1(i,Nx)
        for j in 1:Ny
            array[i,j] = array[mirror_index,j]
        end
    end

    for i in (Nx + 1) : (Nx + halo_x)
        mirror_index = mod1(i,Nx)
        for j in 1:Ny
            array[i,j] = array[mirror_index,j]
        end
    end
end
