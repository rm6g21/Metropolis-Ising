using Random


grid_size = 50

grid_randoms = rand(Float64,(grid_size, grid_size))

lattice = zeros(Int,grid_size,grid_size)

lattice[grid_randoms .> 0.5] .= 1
lattice[grid_randoms .<= 0.5] .= -1



function relative_energy(lattice)
    grid_size = size(lattice, 1)
    energy = 0.0

    for i in 1:grid_size
        for j in 1:grid_size
            # Calculate the indices of the neighboring sites
            left = i == 1 ? grid_size : i - 1
            right = i == grid_size ? 1 : i + 1
            top = j == 1 ? grid_size : j - 1
            bottom = j == grid_size ? 1 : j + 1

            # Sum the values of the neighboring sites
            energy += lattice[i, j] * (lattice[left, j] + lattice[right, j] + lattice[i, top] + lattice[i, bottom])
        end
    end

    return -energy
    end

print(relative_energy(lattice))




