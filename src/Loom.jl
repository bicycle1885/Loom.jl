module Loom

using HDF5: HDF5

struct Connection
    file::HDF5.HDF5File
end

function Base.getproperty(conn::Connection, name::Symbol)
    if name == :ca
        return ColumnAttributes(conn)
    elseif name == :ra
        return RowAttributes(conn)
    elseif name == :matrix
        return Dataset(conn, "/matrix")
    else
        return getfield(conn, name)
    end
end

function Base.show(io::IO, conn::Connection)
    print(io, summary(conn))
    if get(io, :compact, false)
        print(io, "(<filename=", repr(basename(HDF5.filename(conn.file))), ">)")
    else
        m, n = matrixsize(conn)
        print(
            io,
            summary(conn), ':', '\n',
            "  filename: ", basename(HDF5.filename(conn.file)), '\n',
            "  #genes: ", m, '\n',
            "  #cells: ", n, '\n',
            "  column attributes: ", join(colattrs(conn), ", "), '\n',
            "  row attributes: ", join(rowattrs(conn), ", "),
        )
    end
end

function matrixsize(conn::Connection)
    n, m = size(conn.file["/matrix"])
    return m, n
end

ngenes(conn::Connection) = matrixsize(conn)[1]
ncells(conn::Connection) = matrixsize(conn)[2]

colattrs(conn::Connection) = names(conn.file["/col_attrs"])
rowattrs(conn::Connection) = names(conn.file["/row_attrs"])

open(filepath::AbstractString; mode::AbstractString="r") = Connection(HDF5.h5open(filepath, mode))
Base.close(conn::Connection) = close(conn.file)

readmatrix(conn::Connection) = Matrix(read(conn.file["/matrix"])')

# Dataset
# -------

struct Dataset
    conn::Connection
    path::String
end

const Slice = Union{Colon,AbstractRange{<:Integer}}

Base.getindex(dataset::Dataset, i::Integer, j::Integer) = dataset.conn.file[dataset.path][j,i]
Base.getindex(dataset::Dataset, i::Integer, j::Slice) = vec(dataset.conn.file[dataset.path][j,i])
Base.getindex(dataset::Dataset, i::Slice, j::Integer) = vec(dataset.conn.file[dataset.path][j,i])
Base.getindex(dataset::Dataset, i::Slice, j::Slice) = Matrix(dataset.conn.file[dataset.path][j,i]')

# Column Attributes
# -----------------

struct ColumnAttributes
    conn::Connection
end

function Base.show(io::IO, ca::ColumnAttributes)
    print(io, summary(ca), ": ", join(colattrs(ca.conn), ", "))
end

Base.getindex(ca::ColumnAttributes, name::AbstractString) = read(ca.conn.file["/col_attrs/$(name)"])


# Row Attributes
# --------------

struct RowAttributes
    conn::Connection
end

function Base.show(io::IO, ra::RowAttributes)
    print(io, summary(ra), ": ", join(rowattrs(ra.conn), ", "))
end

Base.getindex(ra::RowAttributes, name::AbstractString) = read(ra.conn.file["/row_attrs/$(name)"])

end # module
