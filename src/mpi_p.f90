module mpi_info
  use mpi
  implicit none

  ! include "mpif.h"
  integer:: n                   ! the number of rank that i have(for each processor!
  integer:: isize               !the total number of ranks(size of)
  ! integer:: icommunicator     ! the communicator of comm_world
  integer::ierror,provided
  integer::status(mpi_status_size)
end module
