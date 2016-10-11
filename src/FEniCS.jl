module FEniCS

# Import Python packages
using PyCall
@pyimport dolfin
@pyimport ufl

# Explicitly import overridden symbols
import Base: size, length, *, +, -, Function
import Core.Function

################3
## Base FEniCSObject
abstract FEniCSObject
@inline fenicscall(object::FEniCSObject, func::Union{Symbol,String}, args...) = object.pyobject[func](args...)

macro fenicsclass(name::Symbol, base1::Symbol=:FEniCSObject, base2::Symbol=:Any, base3::Symbol=:Any)
  impl = Symbol(name, "Impl")
  quote
    abstract $name <: $base1, $base2, $base3
    immutable $impl <: $name
      pyobject::PyObject
    end
    $(esc(name))(pyobject::PyObject) = $impl(pyobject)
  end
end

# A macro that formats a Python error if one is thrown
macro formatpyerror(expr)
  return quote
    try
      $(esc(expr))
    catch e
      @show e
      if isa(e, PyCall.PyError)
        println( e.val[:message])
      elseif isa(e, MethodError)
        throw(e)
      end
    end
  end
end

include("mpi.jl")
include("mesh.jl")
include("fem.jl")
include("solver.jl")

include("convenience.jl")

# some general todos
# TODO: use the name= argument in fenics classes and supply something sensible as defaults
# such that str() returns a readable representation
# TODO: Can we import the help? Maybe some function pythonhelp(obj)

end # module





module FenicsTest
using FEniCS, PyCall


@show mesh = UnitSquareMesh(6, 4)
@show V = FunctionSpace(mesh, "Lagrange", 1)

@show u = TrialFunction(V)
@show v = TestFunction(V)
@show f = Constant(6.0)
@show a = inner(nabla_grad(u), nabla_grad(v)) * dx
@show L = f * v * dx

@show unum = Function(V)

problem = LinearVariationalProblem(a, L, unum)
solve(problem)
#TODO: create vector function in CoefficientFunction returning a GenericVector
@show u_nodal_values = unum.pyobject[:vector]()
#TODO: create array function in GenericVector and CoefficientFunction
@show u_array = u_nodal_values[:array]()

end
