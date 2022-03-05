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
    perm = sortperm(atoms)
    return Cell(lattice_vectors, atoms[perm], mod.(coordinates[:, perm], 1.0), length(atoms))
end

function Base.show(io::IO, cell::Cell)
    print(io, "Cell: $(length(cell.atoms)) atoms")
end

include("interface.jl")

end # module
