## The Equation class
@fenicsclass Equation

#compare(form1::Form, form2::Form) = Equation(form1.pyobject[:__eq__](form2.pyobject) )
#Equation(form1::Form, form2::Form) = Equation(dolfin.Equation(form1.pyobject, form2.pyobject, 23) )
export Equation


@fenicsclass Vector
Vector(comm::MPI_Comm, n::Int) = Vector(dolfin.Vector(comm.pyobject, n))
Vector(n::Int) = Vector(mpi_comm_world(), n)
array(v::Vector) = fenicscall(v, :array)::Array{Float64}
export array

## The Function class
# (still problems with the import/export clash with Base.Function or Core.Function)
#@fenicsclass CoefficientFunction
@fenicsclass Function
Function(V::FunctionSpace) = Function(dolfin.Function(V.pyobject))

vector(f::Function) = Vector(fenicscall(f, :vector))
array(f::Function) = array(vector(f))
interpolate(ex::Expression, V::FunctionSpace) = Function(dolfin.interpolate(ex.pyobject, V.pyobject))

export Function, vector, interpolate



## LinearVariationalProblem and Solver
@fenicsclass LinearVariationalProblem
LinearVariationalProblem(a::Form, L::Form, u::Function) =
  LinearVariationalProblem(dolfin.LinearVariationalProblem(a.pyobject, L.pyobject, u.pyobject))
# TODO: add boundary conditions (class is also missing)
# TODO: add parameters

@fenicsclass LinearVariationalSolver
LinearVariationalSolver(problem::LinearVariationalProblem) =
  LinearVariationalSolver(dolfin.LinearVariationalSolver(problem.pyobject))

solve(solver::LinearVariationalSolver) = fenicscall(solver, :solve)
solve(problem::LinearVariationalProblem) = solve(LinearVariationalSolver(problem))
export
  LinearVariationalProblem,
  LinearVariationalSolver, solve
