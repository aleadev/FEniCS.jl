## MPICommunicator class
immutable MPI_Comm <: FEniCSObject
  pyobject::PyObject
end
mpi_comm_world() = MPI_Comm(dolfin.mpi_comm_world())
mpi_comm_self() = MPI_Comm(dolfin.mpi_comm_world())

export MPI_Comm, mpi_comm_world, mpi_comm_self
