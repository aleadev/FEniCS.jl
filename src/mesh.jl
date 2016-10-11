## Mesh class
@fenicsclass Mesh
hmin(mesh::Mesh) = fenicscall(mesh, :hmin)
hmax(mesh::Mesh) = fenicscall(mesh, :hmax)
rmin(mesh::Mesh) = fenicscall(mesh, :rmin)
rmax(mesh::Mesh) = fenicscall(mesh, :rmax)
coordinates(mesh::Mesh) = fenicscall(mesh, :coordinates)
size(mesh::Mesh, dim::Int) = fenicscall(mesh, :size, dim)
# TODO: add more mesh functions

export
  Mesh, hmin, hmax, rmin, rmax, coordinates

## Mesh constructors
UnitTriangleMesh() = Mesh(dolfin.UnitTriangleMesh())
UnitSquareMesh(comm::MPICommunicator, nx::Int, ny::Int; diagonal::String="right") = Mesh(dolfin.UnitSquareMesh(comm, nx, ny, diagonal))
UnitSquareMesh(nx::Int, ny::Int; diagonal::String="right") = Mesh(dolfin.UnitSquareMesh(nx, ny, diagonal))
# TODO: add more mesh constructors (Box, Rectangle)
# TODO: add mesh loading functions

export
  UnitTriangleMesh,
  UnitSquareMesh
