using FEniCS
using Base.Test

# write your own tests here
# @test 1 == 2

@testset "Mesh Tests" begin
  mesh = UnitSquareMesh(2, 2)
  coords = coordinates(mesh)
  coordsex = [(j==0?mod(i,3):div(i,3))*0.5 for i=0:8, j=0:1]
  @test coords == coordsex
end
