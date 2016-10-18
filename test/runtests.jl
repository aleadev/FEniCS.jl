using FEniCS
using Base.Test

# write your own tests here
# @test 1 == 2

@testset "Mesh Tests" begin
  mesh = UnitSquareMesh(2, 2)
  coords = [(j==0?mod(i,3):div(i,3))*0.5 for i=0:8, j=0:1]
  @test coordinates(mesh) == coords
  @test hmin(mesh) ≈ sqrt(0.5)
  @test hmax(mesh) ≈ sqrt(0.5)
  @test rmin(mesh) ≈ 0.5*(1-sqrt(0.5))
  @test rmax(mesh) ≈ 0.5*(1-sqrt(0.5))
  @test size(mesh, 0) == 9 # 9 nodes
  @test size(mesh, 1) == 2*2*3+2*2 # edges
  @test size(mesh, 2) == 8 # triangles (cells)
  @test size(mesh, 3) == 0 # no 3d elements here
  @test size(mesh) == [9, 16, 8, 0]

  mesh = UnitSquareMesh(1, 1)
  init(mesh, 1)
  @test size(mesh) == [4, 5, 2, 0]
  @test cells(mesh) == [0 1 3; 0 2 3]

  mesh = UnitTriangleMesh()
  init(mesh)
  @test size(mesh) == [3, 3, 1, 0]
  @test cells(mesh) == [0 1 2]
end

@testset "MPI Tests" begin
  @test isa(mpi_comm_world(), MPI_Comm)
  @test isa(mpi_comm_self(), MPI_Comm)
  @test isa(UnitSquareMesh(mpi_comm_world(), 4, 4), Mesh) # just test that we have the right type as first argument
  @test isa(UnitSquareMesh(mpi_comm_self(), 4, 4), Mesh) # just test that we have the right type as first argument
end
