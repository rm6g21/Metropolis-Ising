

using Random


grid_size = 50

grid_randoms = rand(Float64,(grid_size, grid_size))

lattice = zeros(Int,grid_size,grid_size)

lattice[grid_randoms .>= 0.25] .= 1
lattice[grid_randoms .< 0.25] .= -1



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

function metropolis(lattice, temperature,iterations,energy)
    lattice = copy(lattice)
    num_spins = zeros(iterations-1)
    energy_list = zeros(iterations-1)

    for t in 1:iterations-1
        i = rand(1:grid_size)
        j = rand(1:grid_size)

        spin_old = lattice[i, j]
        spin_new = -spin_old

        E_old = 0
        E_new = 0

        if i > 1
            E_old += lattice[i - 1, j]
            E_new += spin_new * lattice[i - 1, j]
        end

        if i < grid_size
            E_old += lattice[i + 1, j]
            E_new += spin_new * lattice[i + 1, j]
        end

        if j > 1
            E_old += lattice[i, j - 1]
            E_new += spin_new * lattice[i, j - 1]
        end

        if j < grid_size
            E_old += lattice[i, j + 1]
            E_new += spin_new * lattice[i, j + 1]
        end



        delta_E = E_new - E_old

        if delta_E > 0 && rand() < exp(-delta_E / temperature)
            lattice[i, j] = spin_new
            energy += delta_E
        
        elseif delta_E <= 0
            lattice[i, j] = spin_new
            energy += delta_E
        end

        num_spins[t] = sum(lattice)
        energy_list[t] = energy
    end

    return num_spins, energy_list
    end

spins, energies = metropolis(lattice, 1.0, 100000, relative_energy(lattice))

using Plots
#plotlyjs()


display(plot(1:99999, spins./(grid_size^2), label="Energy", xlabel="Iteration", ylabel="Energy"))





