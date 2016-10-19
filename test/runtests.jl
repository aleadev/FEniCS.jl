using FEniCS
using Base.Test

include("helpers.jl")

@testset "Mesh Tests" begin
  mesh = @inferred_isa UnitSquareMesh(2, 2)
  coords = [(j==0?mod(i,3):div(i,3))*0.5 for i=0:8, j=0:1]
  #@inferred coordinates(mesh)
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
  @test num_cells(mesh) == 8
  @test num_edges(mesh) == 2*2*3+2*2 # edges
  @test num_entities(mesh, 0) == 9
  @test num_faces(mesh) == 8
  @test num_facets(mesh) == 2*2*3+2*2 # edges
  @test num_vertices(mesh) == 9



  mesh = UnitSquareMesh(1, 1)
  init(mesh, 1)
  #@inferred size(mesh)
  @test size(mesh) == [4, 5, 2, 0]
  #@inferred cells(mesh)
  @test cells(mesh) == [0 1 3; 0 2 3]

  mesh = UnitTriangleMesh()
  init(mesh)
  @test size(mesh) == [3, 3, 1, 0]
  @test cells(mesh) == [0 1 2]
end

@testset "MPI Tests" begin
  @inferred mpi_comm_world()
  @test isa(mpi_comm_world(), MPI_Comm)
  @inferred mpi_comm_self()
  @test isa(mpi_comm_self(), MPI_Comm)
  @test isa(UnitSquareMesh(mpi_comm_world(), 4, 4), Mesh) # just test that we have the right type as first argument
  @test isa(UnitSquareMesh(mpi_comm_self(), 4, 4), Mesh) # just test that we have the right type as first argument
end

@testset "FEM Tests" begin
  @testset "FunctionSpace Tests" begin
    mesh = UnitSquareMesh(2, 2)
    V = FunctionSpace(mesh, FEniCS.FAMILY_LAGRANGE, 3)

    @test (@inferred size(FEniCS.mesh(V),0)) == size(mesh, 0)
    @test isa( ufl_element(V), FiniteElement)
    @inferred dim(V)
    @test dim(V) == (3*2+1)^2
  end


  @testset "Expression Tests" begin
    mesh = UnitSquareMesh(1, 1)
    V = FunctionSpace(mesh, FEniCS.FAMILY_LAGRANGE, 1)
    result = [3.0,1.0,4.0,2.0] # ordering may be arbitrary
    u0 = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]")
    @test array(interpolate(u0, V)) == result
    u0 = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", 1)
    @test array(interpolate(u0, V)) == result
    u0 = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", V)
    @test array(interpolate(u0, V)) == result
  end

  @testset "Solve Laplace" begin
    mesh = UnitSquareMesh(6, 4)
    V = FunctionSpace(mesh, "Lagrange", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(6.0)
    a = inner(nabla_grad(u), nabla_grad(v)) * dx
    L = f * v * dx
    unum = FEniCS.Function(V)
    problem = LinearVariationalProblem(a, L, unum)
    solve(problem)
    @test 0*array(unum) ≈ zeros(dim(V))
  end
end

# using FEniCS
# reload("FEniCS")
# mesh = UnitSquareMesh(1, 1)
# V = FunctionSpace(mesh, FEniCS.FAMILY_LAGRANGE, 1)
# v=FEniCS.Vector(dim(V))
# u0 = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]")
# #dolfin.interpolate()
#
