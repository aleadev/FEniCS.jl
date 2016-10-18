## Mesh class
@fenicsclass Mesh
hmin(mesh::Mesh) = fenicscall(mesh, :hmin)
hmax(mesh::Mesh) = fenicscall(mesh, :hmax)
rmin(mesh::Mesh) = fenicscall(mesh, :rmin)
rmax(mesh::Mesh) = fenicscall(mesh, :rmax)
cells(mesh::Mesh) = fenicscall(mesh, :cells)
coordinates(mesh::Mesh) = fenicscall(mesh, :coordinates)
init(mesh::Mesh) = fenicscall(mesh, :init)
init(mesh::Mesh, dim::Int) = fenicscall(mesh, :init, dim)
size(mesh::Mesh, dim::Int) = fenicscall(mesh, :size, dim)
size(mesh::Mesh) = [size(mesh, dim) for dim = 0:3]
num_cells(mesh::Mesh) = fenicscall(mesh, :num_cells)
num_edges(mesh::Mesh) = fenicscall(mesh, :num_edges)
num_entities(mesh::Mesh, dim::Int) = fenicscall(mesh, :num_entities, dim)
num_faces(mesh::Mesh) = fenicscall(mesh, :num_faces)
num_facets(mesh::Mesh) = fenicscall(mesh, :num_facets)
num_vertices(mesh::Mesh) = fenicscall(mesh, :num_vertices)
# TODO: add more mesh functions
# closest_point(point::Point)::Point, closest_cell(point)::Int

export
  Mesh, hmin, hmax, rmin, rmax, init, coordinates, cells,
  num_cells, num_edges, num_entities, num_faces, num_facets, num_vertices



## Mesh constructors
UnitTriangleMesh() = Mesh(dolfin.UnitTriangleMesh())
UnitSquareMesh(comm::MPI_Comm, nx::Int, ny::Int; diagonal::String="right") = Mesh(dolfin.UnitSquareMesh(comm.pyobject, nx, ny, diagonal))
UnitSquareMesh(nx::Int, ny::Int; diagonal::String="right") = Mesh(dolfin.UnitSquareMesh(nx, ny, diagonal))
# TODO: add more mesh constructors (Box, Rectangle)
# TODO: add mesh loading functions

# TODO: sometime later: add MeshEditor (MeshTopology, MeshGeometry,...), Point
export
  UnitTriangleMesh,
  UnitSquareMesh
