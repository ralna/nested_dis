program main
  use hsl_mc70_double
  implicit none

  ! Local variables
  integer :: n, ne

  integer, dimension (:), allocatable :: row, ptr, perm !!,weight

  type (mc70_control) :: nested_control
  type (mc70_info) :: info

  ! Read in the order n of the matrix and the number
  ! of non-zeros in its lower triangular part.
  read (5,*) n, ne

  ! Allocate arrays
  allocate (row(ne),ptr(n+1),perm(n))

  ! Read in pointers
  read (5,*) ptr(1:n+1)

  ! Read in row indices
  read (5,*) row(1:ne)

  ! Call nested dissection and switch to approximate minimum degree when
  !  sub matrix has order less than 8
  nested_control%nd_switch = 4
  call mc70_order(n,ptr,row,perm,nested_control,info)

  ! Print out nested dissection ordering
  write (6,'(a)') ' perm : '
  write (6,'(8i8)') perm
  write (6,'(a)') ' '

  ! Deallocate all arrays
  deallocate (row,ptr,perm)

end program main
