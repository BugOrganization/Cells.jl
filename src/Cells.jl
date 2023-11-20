module Cells

using StaticArrays

export Cell

struct Cell
    # lattice vectors in Angstrom stored in columns
    lattice_vectors::SMatrix{3,3,Float64,9}
    # atomic types
    atoms::Vector{String}
    # atom fractional coordinates stored in columns
    coordinates::Matrix{Float64}
    # number of atoms
    natoms::Int64
end

function Cell(
    lattice_vectors::AbstractMatrix{<:Real},
    atoms::Vector{<:AbstractString},
    coordinates::AbstractMatrix{<:Real}
)
    return Cell(lattice_vectors, atoms, coordinates, length(atoms))
end

function Base.show(io::IO, cell::Cell)
    print(io, "Cell: $(length(cell.atoms)) atoms")
end


function sort_atoms!(cell::Cell)
    perm = sortperm(cell.atoms)
    cell.atoms = cell.atoms[perm]
    cell.coordinates = cell.coordinates[:, perm]
    return cell
end


function move_to_home_unit_cell!(cell::Cell)
    cell.coordinates = mod.(cell.coordinates, 1.0)
    return cell
end


include("interface.jl")

end # module
