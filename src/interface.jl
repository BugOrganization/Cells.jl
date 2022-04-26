module Interface

using Printf, LinearAlgebra

using ..Cells: Cell

export read_poscar, print_poscar, read_fdf, print_fdf, print_openmx

####################################################################################################
# VASP Interface
####################################################################################################

function read_poscar(filename::String)
    return open(filename) do f
        readline(f)
        lattice_scale = parse(Float64, readline(f))
        lattice_vectors = lattice_scale * hcat([map(x -> parse(Float64, x), split(readline(f))) for _ in 1:3]...)
        atom_types = split(readline(f))
        atom_counts = map(x -> parse(Int64, x), split(readline(f)))
        atoms = [atom_types[i] for i in 1:length(atom_types) for _ in 1:atom_counts[i]]
        coordinates_specification = readline(f)
        if startswith(coordinates_specification, r"D|d")
            coordinates = hcat([map(x -> parse(Float64, x), split(readline(f))[1:3]) for _ in 1:sum(atom_counts)]...)
        elseif startswith(coordinates_specification, r"C|c")
            coordinates = inv(lattice_vectors) * hcat([map(x -> parse(Float64, x), split(readline(f))[1:3]) for _ in 1:sum(atom_counts)]...)
        else
            error("Unknown specification of coordinates.")
        end
        Cell(lattice_vectors, atoms, coordinates)
    end
end


function print_poscar(
    filename::String,
    cell::Cell;
    title="Hello",
    selective_dynamics::Bool=false,
    selective_dynamics_string::String=" F F T"
)
    open(filename, "w") do f
        write(f, title)
        @printf f "\n"
        @printf f "%6.1f\n" 1.0
        for i in 1:3
            for j in 1:3
                @printf f "%20.14f" cell.lattice_vectors[j, i]
            end
            @printf f "\n"
        end
        atom_types = unique(cell.atoms)
        atom_numbers = [count(isequal(atom_type), cell.atoms) for atom_type in atom_types]
        for atom_type in atom_types
            @printf f "    %s" atom_type
        end
        @printf f "\n"
        for atom_number in atom_numbers
            @printf f "%6d" atom_number
        end
        @printf f "\n"
        if selective_dynamics
            @printf f "Selective Dynamics\n"
        end
        @printf f "Direct\n"
        for i in 1:size(cell.coordinates, 2)
            for j in 1:3
                @printf f "%20.14f" cell.coordinates[j, i]
            end
            if selective_dynamics
                write(f, selective_dynamics_string)
            end
            @printf f "\n"
        end
    end
    return nothing
end

####################################################################################################
# SIESTA Interface
####################################################################################################

function _startswith_ignore_case(x, y)
    return startswith(strip(uppercase(x)), uppercase(y))
end

function _endswith_ignore_case(x, y)
    return endswith(strip(uppercase(x)), uppercase(y))
end

function _get_scale!(f, scale)
    mark(f)
    foo = readline(f)
    if !_startswith_ignore_case(foo, "LatticeConstant")
        reset(f)
        return false
    else
        unmark(f)
    end
    @assert _endswith_ignore_case(foo, "Ang")
    scale[1] = parse(Float64, split(foo)[2])
    return true
end

function _get_species!(f, species_labels, species_numbers, species_names)
    mark(f)
    if !_endswith_ignore_case(readline(f), "ChemicalSpeciesLabel")
        reset(f)
        return false
    else
        unmark(f)
    end

    while true
        foo = readline(f)
        _endswith_ignore_case(foo, "ChemicalSpeciesLabel") && break
        push!(species_labels, parse(Int64, split(foo)[1]))
        push!(species_numbers, parse(Int64, split(foo)[2]))
        push!(species_names, split(foo)[3])
    end
    return true
end

function _get_lattice!(f, lattice_vectors)
    mark(f)
    if !_endswith_ignore_case(readline(f), "LatticeVectors")
        reset(f)
        return false
    else
        unmark(f)
    end
    for i in 1:3
        lattice_vectors[:, i] = map(x -> parse(Float64, x), split(readline(f)))
    end
    @assert _endswith_ignore_case(readline(f), "LatticeVectors")
    return true
end

function _get_coordinates!(f, coordinates, atom_labels)
    mark(f)
    if !_endswith_ignore_case(readline(f), "AtomicCoordinatesAndAtomicSpecies")
        reset(f)
        return false
    else
        unmark(f)
    end
    while true
        foo = readline(f)
        _endswith_ignore_case(foo, "AtomicCoordinatesAndAtomicSpecies") && break
        push!(coordinates, map(x -> parse(Float64, x), split(foo)[1:3]))
        push!(atom_labels, parse(Int64, split(foo)[4]))
    end
    return true
end

function read_fdf(filename::String)
    return open(filename) do f
        scale = [1.0]
        species_labels = Vector{Int64}()
        species_numbers = Vector{Int64}()
        species_names = Vector{String}()
        lattice_vectors = zeros(3, 3)
        coordinates = Vector{Vector{Float64}}()
        atom_labels = Vector{Int64}()
        while !eof(f)
            matched = false
            matched = _get_scale!(f, scale) || matched
            matched = _get_species!(f, species_labels, species_numbers, species_names) || matched
            matched = _get_lattice!(f, lattice_vectors) || matched
            matched = _get_coordinates!(f, coordinates, atom_labels) || matched
            matched || readline(f)
        end
        atom_names = ["" for _ in atom_labels]
        for (i, species_name) in enumerate(species_names)
            atom_names[atom_labels .== species_labels[i]] .= species_name
        end
        Cell(lattice_vectors * scale[1], atom_names, reduce(hcat, coordinates))
    end
end

function print_fdf(filename::String, cell::Cell, species_numbers::Dict{String, Int64})
    open(filename, "w") do f
        species_names = unique(cell.atoms)
        species_labels = collect(1:length(species_names))
        atom_labels = zeros(Int64, length(cell.atoms))
        for (ispecies, species_name) in enumerate(species_names)
            atom_labels[cell.atoms .== species_name] .= species_labels[ispecies]
        end
        @printf f "NumberOfSpecies%8d\n" length(species_names)
        @printf f "NumberOfAtoms%8d\n" length(cell.atoms)
        @printf f "%%block ChemicalSpeciesLabel\n"
        for (ispecies, species_name) in enumerate(species_names)
            @printf f "%8d%8d%8s\n" species_labels[ispecies] species_numbers[species_name] species_name
        end
        @printf f "%%endblock ChemicalSpeciesLabel\n"
        @printf f "LatticeConstant    1.0 Ang\n"
        @printf f "%%block LatticeVectors\n"
        for i in 1:3
            for j in 1:3
                @printf f "%20.14f" cell.lattice_vectors[j, i]
            end
            @printf f "\n"
        end
        @printf f "%%endblock LatticeVectors\n"
        @printf f "AtomicCoordinatesFormat    Fractional\n"
        @printf f "%%block AtomicCoordinatesAndAtomicSpecies\n"
        for i in 1:size(cell.coordinates, 2)
            for j in 1:3
                @printf f "%20.14f" cell.coordinates[j, i]
            end
            @printf f "%8d" atom_labels[i]
            @printf f "\n"
        end
        @printf f "%%endblock AtomicCoordinatesAndAtomicSpecies\n"
    end
    return nothing
end

####################################################################################################
# OpenMX Interface
####################################################################################################

function print_openmx(
    filename::String,
    cell::Cell,
    species_settings::Dict{String,String},
    species_nelectrons::Dict{String,Int64}
)
    open(filename, "w") do f
        species = unique(cell.atoms)
        @printf f "Species.Number%8d\n" length(species)
        @printf f "<Definition.of.Atomic.Species\n"
        for s in species
            @printf f "%8s  %s\n" s species_settings[s]
        end
        @printf f "Definition.of.Atomic.Species>\n\n"
        
        @printf f "Atoms.UnitVectors.Unit    Ang\n"
        @printf f "<Atoms.UnitVectors\n"
        for i in 1:3
            for j in 1:3
                @printf f "%20.14f" cell.lattice_vectors[j, i]
            end
            @printf f "\n"
        end
        @printf f "Atoms.UnitVectors>\n\n"
        
        @printf f "Atoms.Number%8d\n" length(cell.atoms)
        @printf f "Atoms.SpeciesAndCoordinates.Unit    Frac\n"
        @printf f "<Atoms.SpeciesAndCoordinates\n"
        for i in 1:length(cell.atoms)
            @printf(
                f, "%8d%4s%20.14f%20.14f%20.14f%4.1f%4.1f%s\n", i, cell.atoms[i],
                cell.coordinates[1, i], cell.coordinates[2, i], cell.coordinates[3, i],
                species_nelectrons[cell.atoms[i]] / 2, species_nelectrons[cell.atoms[i]] / 2,
                " 0.0 0.0 0.0 0.0 0 off"
            )
        end
        @printf f "Atoms.SpeciesAndCoordinates>\n\n"
    end
    return nothing
end

end
