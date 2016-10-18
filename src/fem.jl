## Boundary conditions
# TODO: add class
# TODO: add constructors and other methods


## FunctionSpace class
@fenicsclass FiniteElement

## FunctionSpace class
@fenicsclass FunctionSpace
FunctionSpace(mesh::Mesh, family::String, degree::Int) = FunctionSpace(dolfin.FunctionSpace(mesh.pyobject, family, degree))
export FunctionSpace
# TODO: import operator(s), at least * (shall that be mapped to ⊗?)
const FAMILY_ARGYRIS = "ARG"
const FAMILY_ARNOLD_WINTHER = "AW"
const FAMILY_BREZZI_DOUGLAS_FORTIN_MARINI = "BDFM"
const FAMILY_BREZZI_DOUGLAS_MARINI = "BDM"
const FAMILY_BUBBLE = "B"
const FAMILY_CROUZEIX_RAVIART = "CR"
const FAMILY_DISCONTINUOUS_LAGRANGE = "DG"
const FAMILY_DISCONTINUOUS_RAVIART_THOMAS = "DRT"
const FAMILY_HERMITE = "HER"
const FAMILY_LAGRANGE = "CG"
const FAMILY_MARDAL_TAI_WINTHER = "MTW"
const FAMILY_MORLEY = "MOR"
const FAMILY_NEDELEC_1ST_KIND = "N1curl"
const FAMILY_NEDELEC_2ND_KIND = "N2curl"
const FAMILY_QUADRATURE = "Quadrature"
const FAMILY_RAVIART_THOMAS = "RT"
const FAMILY_REAL = "R"

mesh(V::FunctionSpace) = Mesh(fenicscall(V, :mesh))
ufl_element(V::FunctionSpace) = FiniteElement(fenicscall(V, :ufl_element))
family(V::FunctionSpace) = fenicscall(V, :family)::AbstractString
dim(V::FunctionSpace) = fenicscall(V, :dim)

export mesh, ufl_element, family, dim


## Expression class
@fenicsclass Expression
Expression(expr::String) = Expression(dolfin.Expression(expr))
Expression(expr::String, degree::Int) = Expression(dolfin.Expression(expr, degree=degree))
Expression(expr::String, element::FiniteElement) = Expression(dolfin.Expression(expr, element=element.pyobject))
Expression(expr::String, V::FunctionSpace) = Expression(expr, ufl_element(V))
export Expression
# TODO: add degree or element argument (make one necessary) (either degree=2, or V.ufl_element() or V directly)
# TODO: add optional domain arguments

## Expr class hierarchy (export only classes that are really needed in user code)
@fenicsclass Expr
@fenicsclass Argument Expr
@fenicsclass Constant Expr
TrialFunction(V::FunctionSpace) = Argument(dolfin.TrialFunction(V.pyobject))
TestFunction(V::FunctionSpace) = Argument(dolfin.TestFunction(V.pyobject))
Constant(x::Real) = Constant(dolfin.Constant(x, name="Constant($x)"))

export
  TrialFunction,
  TestFunction,
  Constant

## functions that operate on and return Expr's
nabla_grad(u::Expr) = Expr(dolfin.nabla_grad(u.pyobject))
inner(u::Expr, v::Expr) = Expr(dolfin.inner(u.pyobject, v.pyobject))
# TODO: add more functions (e.g. nabla_div, ...)
export
  nabla_grad, inner

## Measure class (export only concrete instances)
@fenicsclass Measure
dx = Measure(dolfin.dx)
ds = Measure(dolfin.ds)
export dx, ds

## Form class and operators that generate forms
@fenicsclass Form
*(expr::Expr, measure::Measure) = Form(measure.pyobject[:__rmul__](expr.pyobject) )
*(expr::Expr, expr2::Expr) = Expr(expr.pyobject[:__mul__](expr2.pyobject) )
# TODO: import more operators
