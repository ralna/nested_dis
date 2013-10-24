PROGRAM main
  USE hsl_mc70_double
      USE hsl_mc69_double
  IMPLICIT NONE


!--------------------------------------
! Parameters
!--------------------------------------
      INTEGER, PARAMETER :: myint1 = kind(1), wp = kind(1.0D+0), &
        myint = kind(1), leniw = 1000000, lena = 1000000, maxn = 10000, &
        maxnz = 1000000, lenirn = 1000000, maxrhs = 5
      REAL (wp), PARAMETER :: zero = 0.0_wp
!--------------------------------------
! Local variables
!--------------------------------------
  INTEGER :: i, n, ne, st, test,jj,k,netot,nf,mc69flag

      INTEGER :: j, kase, lcase, lfact, liw, lkeep, lp, &
        lrhs, ncase, ne1, nres, numed0, numed1, tour, &
        ix, iy

      INTEGER :: irn1(lenirn), iw(leniw), jcn1(lenirn), jcolst(maxn+1), &
        keep(5*maxn+2*maxnz+42), ym11_icntl(10)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ptr, row, perm, seps
      INTEGER, DIMENSION (:), ALLOCATABLE :: irn, row_out, ptr_out, &
        jcn
      REAL (wp) :: eps, resmax
      REAL (wp) :: aa1(maxnz), w(maxn,maxrhs)
      LOGICAL :: def, llcase, ltest, undef
      CHARACTER(len=8) :: key1
  TYPE (mc70_control) :: control 
  TYPE (mc70_info) :: info

! .. External Functions ..
      REAL (wp) fa14ad, fd15ad
      EXTERNAL fa14ad, fd15ad

  control%ml = 0

!--------------------------------------
! Test error flag n<1
!--------------------------------------
  test = 1
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 0
  ne = 0
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1) = 1

  control%nd_switch = 2
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    CALL mc70_order_full(n,ptr,row,perm ,control,info,seps)     
  END DO
  control%print_level = 0
  control%nd_switch = 50

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test diagonal matrix
!--------------------------------------
  test = 2
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 10
  ne = 0
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = 1

  control%nd_switch = 2
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    CALL mc70_order_full(n,ptr,row,perm ,control,info,seps)     
  END DO
  control%print_level = 0
  control%nd_switch = 50

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test diagonal submatrix
!--------------------------------------
  test = 3
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 10
  ne = 10
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n-2) = (/ (i,i=1,n-2) /)
  ptr(n-1) = ptr(n-2)+2
  ptr(n:n+1) = ne+1
  row(1:ne-2) = n
  row(ne-1) = n-1
  row(ne) = n

  control%nd_switch = 2
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)       
  END DO
  control%print_level = 0
  control%nd_switch = 50

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test expansion of matrix from lower triangle input
!--------------------------------------
  test = 4
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 10
  ne = 9
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n-1) = (/ (i,i=1,n-1) /)
  ptr(n:n+1) = n
  row(1:ne) = (/ (i,i=2,n) /)

  control%nd_switch = 2
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test expansion of matrix from lower triangle format and that diags are removed
!--------------------------------------

  test = 5
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 10
  ne = 19
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n) = (/ (2*i-1,i=1,n) /)
  ptr(n+1) = ne+1
  row(1:ne) = (/ 1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10 /)

  control%nd_switch = 2
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test expansion of matrix from lower and upper triangle input with no diags
!--------------------------------------
  test = 6
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 10
  ne = 18
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1) = 1
  ptr(2:n) = (/ (2*(i-1),i=2,n) /)
  ptr(n+1) = ne+1
  row(1) = 2
  DO i = 2,n-1
    row(ptr(i)) = i-1
    row(ptr(i)+1) = i+1
  END DO
  row(ptr(n)) = n-1

  control%nd_switch = 2
  DO i = 2,2
    control%print_level = i
    CALL mc70_order_full(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test expansion of matrix from lower and upper triangle input with no diags
!--------------------------------------
  test = 7
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 10
  ne = 19
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1) = 1
  ptr(2:n) = (/ (2*(i-1),i=2,n) /)
  ptr(n+1) = ne+1
  row(1) = 2
  DO i = 2,n-1
    row(ptr(i)) = i-1
    row(ptr(i)+1) = i+1
  END DO
  row(ptr(n)) = n-1
  row(ne) = n

  control%nd_switch = 2
  DO i = 0,2
    control%print_level = i
    CALL mc70_order_full(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test dense row removal - row 1 dense
!--------------------------------------
test = 8
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 800
  ne = n
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1) = 1
  ptr(2:n+1) = ne+1
  row(1:ne) = (/ (i,i=1,n) /)

  control%nd_switch = 2
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test dense row removal - first row dense
!--------------------------------------
test = 9
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 800
  ne = 2*n-3
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1) = 1
  ptr(2:n) = (/ (n+i-2,i=2,n) /)
  ptr(n+1) = ne+1
  row(1:ptr(2)-1) = (/ (i,i=2,n) /)
  row(ptr(2):ptr(n)-1) = (/ (i+1,i=2,n-1)  /)

  control%nd_switch = 2
  DO i = -1,1
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test dense row removal - last row dense
!--------------------------------------
test = 10
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 800
  ne = n-1
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n-1) = (/ (i,i=1,n-1)  /)
  ptr(n:n+1) = ne+1
  row(1:ne) = n

  control%nd_switch = 2
  DO i = -1,0
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test dense row removal - first two rows dense but row 2 has max degree
!--------------------------------------
  test = 11
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 800
  ne = n-3 + n - 2
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1) = 1
  ptr(2) = n-2
  ptr(3:n+1) = ne+1
  row(ptr(1):ptr(2)-1) = (/ (i+3,i=1,n-3)  /)
  row(ptr(2):ne) = (/ (i+2,i=1,n-2)  /)

  control%nd_switch = 2
  DO i = -1,0
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test independent component detection
!--------------------------------------
  test = 12
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 5
  ne = 3
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,2,2,3,4,4  /)
  row(1:ne) = (/ 2,4,5 /)

  control%nd_switch = 2
  DO i = -1,2
    write(*,*) 'i',i
    control%print_level = i
    control%partition_method = 2
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
    control%partition_method = 1
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20
!--------------------------------------
! Test Ashcraft method
!--------------------------------------
 
  test = 13
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 10
  ne = 14
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,3,5,7,9,11,12,13,14,15,15  /)
  row(1:ne) = (/ 2,7,7,3,4,8,5,9,6,10,10,8,9,10 /)

  control%nd_switch = 5
  control%partition_method = 1
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test Ashcraft method
!--------------------------------------
 
  test = 14
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 10
  ne = 15
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,4,6,8,9,10,11,13,15,16,16  /)
  row(1:ne) = (/ 2,3,4,3,5,4,8,6,7,9,8,10,9,10,10 /)

  control%nd_switch = 2
  control%partition_method = 1
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test Ashcraft method
!--------------------------------------
 
  test = 15
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 13
  ne = 16
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,2,4,6,8,10,12,13,13,14,14,16,17,17 /)
  row(1:ne) = (/ 2,3,4,4,11,5,6,7,8,9,10,8,10,12,13,13 /)

  control%nd_switch = 2
  control%partition_method = 1
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test Ashcraft method
!--------------------------------------
 
  test = 16
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 13
  ne = 16
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,2,4,6,8,10,12,13,13,14,14,16,17,17 /)
  row(1:ne) = (/ 2,3,4,4,11,5,6,7,8,9,10,8,10,12,13,13 /)

  control%nd_switch = 2
  control%partition_method = 1
  control%ratio = 2.0/3.0
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1
  control%ratio = 0.0

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20
!--------------------------------------
! Test one-sided levelset method
!--------------------------------------
 
  test = 17
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 10
  ne = 14
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,3,5,7,9,11,12,13,14,15,15  /)
  row(1:ne) = (/ 2,7,7,3,4,8,5,9,6,10,10,8,9,10 /)

  control%nd_switch = 5
  control%partition_method = 2
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test one-sided levelset method
!--------------------------------------
 
  test = 18
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 10
  ne = 14
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,3,5,7,9,11,12,13,14,15,15  /)
  row(1:ne) = (/ 2,7,7,3,4,8,5,9,6,10,10,8,9,10 /)

  control%nd_switch = 5
  control%partition_method = 2
  control%ratio = 2.0/3.0
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1
  control%ratio = 0.0

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20


!--------------------------------------
! Test DM refinement with 2-sided partition
!--------------------------------------
 
  test = 19
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 8
  ne = 11
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,4,5,6,8,9,11,12,12  /)
  row(1:ne) = (/ 2,3,4,4,4,5,6,6,7,8,8 /)

  control%nd_switch = 5
  control%partition_method = 1
  control%refinement = 2
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1
  control%refinement = 1

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test DM refinement with 2-sided partition
!--------------------------------------
 
  test = 20
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 21
  ne = 37
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/1,3,5,7,9,11,13,15,16,20,22,24,26,28,29,29,32,33,35,37,38,38 /)
  row(1:ne) = (/ 2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,11,16,17,11,13,12,13,13,15,&
        14,15,15,17,18,19,18,19,20,20,21,21    /)

  control%nd_switch = 5
  control%partition_method = 1
  control%refinement = 2
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1
  control%refinement = 1

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test DM refinement with 1-sided partition
!--------------------------------------
 
  test = 21
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 31
  ne = 49
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,4,6,8,10,12,13,15,16,18,20,21,22,25,26,27,29,31,32,34,36,&
                  38,39,40,42,44,45,47,48,49,50,50       /)
  row(1:ne) = (/ 2,3,4,4,5,6,7,7,8,8,9,10,10,11,11,11,12,13,14,14,14,15,17,&
                18,18,16,17,19,20,21,21,20,22,23,24,24,25,23,26,26,27,27,28,&
                29,29,30,30,31,31 /)

  control%nd_switch = 5
  control%partition_method = 2
  control%refinement = 2
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1
  control%refinement = 1

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test DM refinement with 2-sided partition
!--------------------------------------
 
  test = 22
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 7
  ne = 6
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,2,3,4,5,6,7,7      /)
  row(1:ne) = (/ 7,7,7,7,7,7 /)

  control%nd_switch = 5
  control%partition_method = 1
  control%refinement = 2
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1
  control%refinement = 1

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test DM refinement with 1-sided partition
!--------------------------------------
 
  test = 23
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 7
  ne = 6
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,2,3,4,5,6,7,7      /)
  row(1:ne) = (/ 7,7,7,7,7,7 /)

  control%nd_switch = 5
  control%partition_method = 2
  control%refinement = 2
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1
  control%refinement = 1

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20


!--------------------------------------
! Test DM refinement with 1-sided partition with balanced partition
!--------------------------------------
 
  test = 24
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 31
  ne = 49
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,4,6,8,10,12,13,15,16,18,20,21,22,25,26,27,29,31,32,34,36,&
                  38,39,40,42,44,45,47,48,49,50,50       /)
  row(1:ne) = (/ 2,3,4,4,5,6,7,7,8,8,9,10,10,11,11,11,12,13,14,14,14,15,17,&
                18,18,16,17,19,20,21,21,20,22,23,24,24,25,23,26,26,27,27,28,&
                29,29,30,30,31,31 /)

  control%nd_switch = 5
  control%partition_method = 2
  control%refinement = 2
  control%ratio = 0.5
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1
  control%refinement = 1
  control%ratio = 0.0

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test DM refinement with 2-sided partition with balanced partition
!--------------------------------------
 
  test = 25
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 31
  ne = 49
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,4,6,8,10,12,13,15,16,18,20,21,22,25,26,27,29,31,32,34,36,&
                  38,39,40,42,44,45,47,48,49,50,50       /)
  row(1:ne) = (/ 2,3,4,4,5,6,7,7,8,8,9,10,10,11,11,11,12,13,14,14,14,15,17,&
                18,18,16,17,19,20,21,21,20,22,23,24,24,25,23,26,26,27,27,28,&
                29,29,30,30,31,31 /)

  control%nd_switch = 5
  control%partition_method = 1
  control%refinement = 2
  control%ratio = 0.5
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1
  control%refinement = 1
  control%ratio = 0.0

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Test DM refinement with 2-sided partition with balanced partition
!--------------------------------------
 
  test = 26
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 31
  ne = 49
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,3,5,7,8,10,12,13,15,17,18,20,22,23,25,27,28,29,32,33,34,37,&
                 39,40,42,44,45,46,48,49,50,50 /)
  row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9,11,11,12,10,12,13,14,15,13,15,&
                 16,18,19,16,19,17,19,20,21,22,22,23,23,24,25,25,26,27,27,28,&
                 28,29,29,30,30,31,31,31 /)

  control%nd_switch = 5
  control%partition_method = 1
  control%refinement = 2
  control%ratio = 0.95
  control%refinement_band = n
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1
  control%refinement = 1
  control%ratio = 0.0
  control%refinement_band = -3

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20


!--------------------------------------
! Test multigrid, 2-sided partitioning
!--------------------------------------
 
  test = 27
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 31
  ne = 49
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,3,5,7,8,10,12,13,15,17,18,20,22,23,25,27,28,29,32,33,34,37,&
                 39,40,42,44,45,46,48,49,50,50 /)
  row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9,11,11,12,10,12,13,14,15,13,15,&
                 16,18,19,16,19,17,19,20,21,22,22,23,23,24,25,25,26,27,27,28,&
                 28,29,29,30,30,31,31,31 /)

  control%nd_switch = 5
  control%partition_method = 1
  control%refinement = 2
  control%ratio = 0.95
  control%refinement_band = n
  control%ml = 1
  control%ml_switch = 15
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1
  control%refinement = 1
  control%ratio = 0.0
  control%refinement_band = -3
  control%ml = 0
  control%ml_switch = 100

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20


!--------------------------------------
! Test multigrid, 1-sided partitioning
!--------------------------------------
 
  test = 28
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
  n = 31
  ne = 49
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,3,5,7,8,10,12,13,15,17,18,20,22,23,25,27,28,29,32,33,34,37,&
                 39,40,42,44,45,46,48,49,50,50 /)
  row(1:ne) = (/ 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9,11,11,12,10,12,13,14,15,13,15,&
                 16,18,19,16,19,17,19,20,21,22,22,23,23,24,25,25,26,27,27,28,&
                 28,29,29,30,30,31,31,31 /)

  control%nd_switch = 5
  control%partition_method = 2
  control%refinement = 2
  control%ratio = 0.95
  control%refinement_band = n
  control%ml = 1
  control%ml_switch = 15
  DO i = -1,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1
  control%refinement = 1
  control%ratio = 0.0
  control%refinement_band = -3
  control%ml = 0
  control%ml_switch = 100

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20
!--------------------------------------
! Call AMD with no nested
!--------------------------------------
 
  test = 29
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
   n = 31
  ne = 49
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1:n+1) = (/ 1,4,6,8,10,12,13,15,16,18,20,21,22,25,26,27,29,31,32,34,36,&
                  38,39,40,42,44,45,47,48,49,50,50       /)
  row(1:ne) = (/ 2,3,4,4,5,6,7,7,8,8,9,10,10,11,11,11,12,13,14,14,14,15,17,&
                18,18,16,17,19,20,21,21,20,22,23,24,24,25,23,26,26,27,27,28,&
                29,29,30,30,31,31 /)


  control%nd_switch = 32
  control%partition_method = 2
  DO i = 2,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Call AMD with no nested
!--------------------------------------
 
  test = 30
  WRITE(6,'(a1)') ' '
  WRITE(6,'(a20,i5,a4)') '****Test section = ', test,'****'
   n = 31
  ne = 1
  ALLOCATE(ptr(n+1),row(ne),perm(n),seps(n),stat=st)
  IF (st .ne. 0) GOTO 10
  
  ptr(1) = 1
  ptr(2:n+1) = 2
  row(1:ne) = (/ 2 /)


  control%nd_switch = 32
  control%partition_method = 2
  DO i = 2,2
    control%print_level = i
    CALL mc70_order(n,ptr,row,perm ,control,info,seps)     
    WRITE(6,'(a9)') 'perm = '
    WRITE(6,'(10i5)') perm(:)
  END DO
  control%print_level = 0
  control%nd_switch = 50
  control%partition_method = 1

  DEALLOCATE(ptr,row,perm,seps,stat=st)
  IF (st .ne. 0) GOTO 20

!--------------------------------------
! Tests for AMD coverage
!--------------------------------------

      ALLOCATE (irn(500000),jcn(500000),ptr(50000))
      DO n = 10, 50
! N is the matrix dimension
        DO nf = 1, n/5
! NF is the number of full rows
          netot = 0
          DO jj = 1, nf
! Add full rows
            DO k = jj + 1, n
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
          END DO
          DO jj = nf + 1, n
            netot = netot + 1
            iw(netot) = jj
            jcn(netot) = jj - 1
          END DO

! Copy column indices into IW
          DO jj = 1, netot
            iw(netot+jj) = jcn(jj)
          END DO
          ne = netot
          irn(1:ne) = iw(1:ne)
          IF (allocated(ptr_out)) DEALLOCATE (ptr_out)
          ALLOCATE (ptr_out(n+1))

          CALL mc69_coord_convert(4,n,n,ne,irn,jcn,ptr_out,row_out,mc69flag)
          IF (mc69flag<0) THEN
            GO TO 70
          END IF

          IF (allocated(perm)) DEALLOCATE (perm)
          IF (allocated(seps)) DEALLOCATE (seps)
          ALLOCATE (perm(n),seps(n))
            control%nd_switch = n-1
            CALL mc70_order(n,ptr_out,row_out(1:ptr_out(n+ &
              1)-1),perm ,control,info,seps)     
            control%nd_switch = 50
            CALL testpr(n,perm)
            WRITE (6,'(a10)') 'perm='
            WRITE (6,'(10I6)') (perm(1:min(20,n)))
        END DO
      END DO

      DO n = 10, 100
! N is the matrix dimension
        DO nf = 1, n/5
! NF is the number of full rows
          netot = 0
          DO jj = 1, nf
! Add quasi-dense rows
            DO k = jj + 7, n
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
          END DO
! Add a full row
          jj = nf + 1
          DO k = jj + 1, n
            netot = netot + 1
            iw(netot) = k
            jcn(netot) = jj
          END DO
          DO k = 1, jj - 1
            netot = netot + 1
            jcn(netot) = k
            iw(netot) = jj
          END DO
! Add sparse rows
          DO jj = nf + 2, n
            netot = netot + 1
            iw(netot) = jj
            jcn(netot) = jj - 1
          END DO

! Copy column indices into IW
          DO jj = 1, netot
            iw(netot+jj) = jcn(jj)
          END DO
          ne = netot
          irn(1:ne) = iw(1:ne)
          IF (allocated(ptr_out)) DEALLOCATE (ptr_out)
          ALLOCATE (ptr_out(n+1))

          CALL mc69_coord_convert(4,n,n,ne,irn,jcn,ptr_out,row_out,mc69flag)
          IF (mc69flag<0) THEN
            GO TO 70
          END IF
! call mc68_setup(a,n,ne,jcn,control,iinfo,irn=irn)
          IF (allocated(perm)) DEALLOCATE (perm)
          IF (allocated(seps)) DEALLOCATE (seps)
          ALLOCATE (perm(n),seps(n))
            control%nd_switch = n-1
            CALL mc70_order(n,ptr_out,row_out(1:ptr_out(n+ &
              1)-1),perm ,control,info,seps)     
            control%nd_switch = 50
            CALL testpr(n,perm)
            WRITE (6,'(a10)') 'perm='
            WRITE (6,'(10I6)') (perm(1:min(20,n)))
        END DO
      END DO

      DO n = 10, 100
! N is the matrix dimension
        DO nf = 1, n/5
! NF is the number of full rows
          netot = 0
          DO jj = 1, nf
! Add quasi-dense rows
            DO k = jj + 1, n - n/4
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
          END DO

! Add a full rows
          DO jj = nf + 1, nf + 1
            DO k = jj + 1, n
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
            DO k = 1, jj - 1
              netot = netot + 1
              jcn(netot) = k
              iw(netot) = jj
            END DO
          END DO
          DO jj = nf + 2, n
            netot = netot + 1
            iw(netot) = jj
            jcn(netot) = jj - 1
          END DO

! Copy column indices into IW
          DO jj = 1, netot
            iw(netot+jj) = jcn(jj)
          END DO
          ne = netot
          irn(1:ne) = iw(1:ne)
          IF (allocated(ptr_out)) DEALLOCATE (ptr_out)
          ALLOCATE (ptr_out(n+1))

          CALL mc69_coord_convert(4,n,n,ne,irn,jcn,ptr_out,row_out,mc69flag)
          IF (mc69flag<0) THEN
            GO TO 70
          END IF
!  call mc68_setup(a,n,ne,jcn,control,iinfo,irn=irn)
          IF (allocated(perm)) DEALLOCATE (perm)
          IF (allocated(seps)) DEALLOCATE (seps)
          ALLOCATE (perm(n),seps(n))
            control%nd_switch = n-1
            CALL mc70_order(n,ptr_out,row_out(1:ptr_out(n+ &
              1)-1),perm ,control,info,seps)     
            control%nd_switch = 50
            CALL testpr(n,perm)
            WRITE (6,'(a10)') 'perm='
            WRITE (6,'(10I6)') (perm(1:min(20,n)))
        END DO
      END DO


      n = 1000
! N is the matrix dimension
      DO nf = 1, n/4 + 1, n/8
! NF is the number of full rows
        netot = 0
        DO jj = 1, nf
! Add quasi-dense rows
          DO k = jj + 1, n - 2*jj
            netot = netot + 1
            iw(netot) = k
            jcn(netot) = jj
          END DO
        END DO
        DO jj = nf + 2, n
          netot = netot + 1
          iw(netot) = jj
          jcn(netot) = jj - 2
        END DO

! Copy column indices into IW
        DO jj = 1, netot
          iw(netot+jj) = jcn(jj)
        END DO
        ne = netot
        irn(1:ne) = iw(1:ne)
        IF (allocated(ptr_out)) DEALLOCATE (ptr_out)
        ALLOCATE (ptr_out(n+1))

        CALL mc69_coord_convert(4,n,n,ne,irn,jcn,ptr_out,row_out,mc69flag)
        IF (mc69flag<0) THEN
          GO TO 70
        END IF
!call mc68_setup(a,n,ne,jcn,control,iinfo,irn=irn)
        IF (allocated(perm)) DEALLOCATE (perm)
          IF (allocated(seps)) DEALLOCATE (seps)
        ALLOCATE (perm(n),seps(n))
            control%nd_switch = n-1
            CALL mc70_order(n,ptr_out,row_out(1:ptr_out(n+ &
              1)-1),perm ,control,info,seps)     
            control%nd_switch = 50
          CALL testpr(n,perm)
          WRITE (6,'(a10)') 'perm='
          WRITE (6,'(10I6)') (perm(1:min(20,n)))
      END DO

      DO n = 1000, 1000, 1000
! N is the matrix dimension
        DO nf = 1, 1
! NF is the number of full rows
          netot = 0
          DO jj = 1, nf
! Add quasi-dense rows
            DO k = jj + 10, jj + 100
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
            DO k = jj + 110, n
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
          END DO
! Add sparse rows
          DO jj = 2, 40
            DO k = 98 + jj, 98 + jj
              netot = netot + 1
              iw(netot) = k
              jcn(netot) = jj
            END DO
          END DO
! Copy column indices into IW
          DO jj = 1, netot
            iw(netot+jj) = jcn(jj)
          END DO
          ne = netot
          irn(1:ne) = iw(1:ne)
          IF (allocated(ptr_out)) DEALLOCATE (ptr_out)
          ALLOCATE (ptr_out(n+1))

          CALL mc69_coord_convert(4,n,n,ne,irn,jcn,ptr_out,row_out,mc69flag)
          IF (mc69flag<0) THEN
            GO TO 70
          END IF
! call mc68_setup(a,n,ne,jcn,control,iinfo,irn=irn)
          IF (allocated(perm)) DEALLOCATE (perm)
          IF (allocated(seps)) DEALLOCATE (seps)
          ALLOCATE (perm(n),seps(n))


            control%nd_switch = n-1
            CALL mc70_order(n,ptr_out,row_out(1:ptr_out(n+ &
              1)-1),perm ,control,info,seps)     
            control%nd_switch = 50
            CALL testpr(n,perm)
            WRITE (6,'(a10)') 'perm='
            WRITE (6,'(10I6)') (perm(1:min(20,n)))
        END DO
      END DO

      ltest = .FALSE.

      resmax = zero
      nres = 0

      numed0 = 0
      numed1 = 0

      lrhs = maxn
      lkeep = 5*maxn + 2*maxnz + 42
      keep(:) = 0

      lp = 6

! Initialize LLCASE.   Again only used for serious debug.
      llcase = .FALSE.

! Generate value of eps
      eps = fd15ad('E')

      CALL fa14id(ix)
      CALL ym11id(ym11_icntl,iy)

      tour = 1

           ! DO ncase = -5, 24
            DO ncase = -5, 24
              lcase = 20
              IF (ncase<0) lcase = 10
              IF (ncase==0) lcase = 11
              IF (ncase==-5) lcase = 6
              IF (ncase==21) lcase = 10
              IF (ncase==22) lcase = 2
              IF (ncase==23) lcase = 1
              IF (ncase==24) lcase = 1
              DO kase = 1, lcase
                control%print_level = 2
                WRITE (lp,'(//A,2I8)') '***** Case *****', ncase, kase

                IF (llcase) THEN
                  GO TO 40
                ELSE
! Set LLCASE
! LTEST  = TOUR.EQ.1 .AND. NCASE.EQ.-4.AND. KASE.EQ.4
                  llcase = ltest

! Set LIW and LFACT
                  liw = leniw
                  lfact = lena

! jump to reading special cases

! Set N, IN
                    n = 50*ncase
                    IF (ncase<0) THEN
                      n = 100
                    END IF
                    IF (ncase==0) n = 12

! Set NE
                    CALL fa14bd(ix,n,ne)
                    CALL fa14bd(ix,n,ne1)
                    ne1 = ne*ne1
                    IF (fa14ad(ix,1)>=0.01) ne1 = max(ne1,1+n/2)
! IF (fa04ad(1)>=0.01) ne1 = max(ne1,1+n/2)
! Generate test matrix
                    ym11_icntl(2) = 0
                    ym11_icntl(3) = 0
                    key1 = 'file'
                    CALL ym11ad(n,n,ne1,ne,irn1,aa1,jcolst,iw,ym11_icntl,key1, &
                      iy)
                    DO j = 1, ne
                      aa1(j) = float(j)
                    END DO

                    def = .FALSE.
                    undef = .FALSE.
                    IF (ncase>0 .AND. (kase==3 .OR. kase==5)) undef = .TRUE.
                    IF (ncase>0 .AND. (kase==6 .OR. kase==8)) def = .TRUE.
                    IF (undef) CALL ndef(n,ne,aa1,jcolst,irn1,w,.TRUE.)
                    IF (def) CALL pdef(n,ne,aa1,jcolst,irn1,w,.TRUE.)
! Input negative definite matrix
                    IF (ncase>0 .AND. kase==8) THEN
                      DO k = 1, ne
                        aa1(k) = -aa1(k)
                      END DO
                    END IF
                    IF (ne==-1) THEN
                      GO TO 100
                    END IF

                  IF (llcase) THEN
                    control%print_level = 2
                  END IF

                  IF (ncase==-3) THEN
                    IF (kase==1) control%print_level = 0
                  END IF
                  IF (ncase==-2) THEN
                    jcn1(1) = irn1(2)
                    IF (kase==1) control%print_level = 0
                  END IF
                  IF (ncase==-1) THEN
                    jcn1(1) = irn1(2)
                  END IF

                  IF (ncase==-5) THEN
                    IF (kase==1) n = 1
                    IF (kase==1 .OR. kase==2) control%print_level = 0
                  END IF
                  IF (ncase==22 .AND. kase==1) lkeep = 2
                  IF (ltest) THEN
                    control%print_level = 1
                    WRITE (lp,'(A)') 'Testing test ... MC68 (first)'
                  END IF

                  IF (allocated(perm)) DEALLOCATE (perm)
                  IF (allocated(seps)) DEALLOCATE (seps)
                  ALLOCATE (perm(n),seps(n))
                  IF (ncase.GT.10) control%nd_switch = max(50,n/4)
                      control%print_level = 0
            CALL mc70_order(n,jcolst(1:n+1),irn1(1:jcolst(n+ &
              1)-1),perm ,control,info,seps)     
                      control%print_level = 0
                  IF (ncase.GT.20) control%nd_switch = 50
                      CALL testpr(n,perm)

                  IF (ncase==-5 .AND. (kase==1 .OR. kase==2)) &
                    control%print_level = -1
                  IF (ncase==22 .AND. kase==1) lkeep = 5*maxn + 2*maxnz + 42
                  IF (ncase/=22 .OR. kase/=1) THEN
                    IF (.TRUE.) THEN
                      IF (ncase<=20 .AND. kase<=4) THEN
                        IF (ncase==-5) THEN
                         ! n = 10
                          IF (kase==3) keep(1) = keep(2)
                          IF (kase==4) keep(1) = 0
                        END IF
                        IF (ncase==-4) THEN
                          DO i = 1, n
                            keep(i) = i
                          END DO
                        END IF
                        IF (ltest) THEN
                          control%print_level = 1
                          WRITE (lp,'(A)') &
                            'Testing test ... MC68 (second)'
                        END IF

                        IF (allocated(perm)) DEALLOCATE (perm)
                        IF (allocated(seps)) DEALLOCATE (seps)
                      ALLOCATE (perm(n),seps(n))
                      control%print_level = 0
            CALL mc70_order(n,jcolst(1:n+1),irn1(1:jcolst(n+ &
              1)-1),perm ,control,info,seps)     
                      control%print_level = 0
                            CALL testpr(n,perm)

                      END IF


                    ELSE IF (ncase>0) THEN
                      GO TO 30
                    END IF
                  END IF
                END IF
100              CONTINUE
              END DO
            END DO
! End of extra loop

      GO TO 50

30    CONTINUE
      STOP
40    STOP
50    CONTINUE


70    IF (allocated(irn)) DEALLOCATE (irn)
      IF (allocated(ptr)) DEALLOCATE (ptr)
      IF (allocated(perm)) DEALLOCATE (perm)
      IF (allocated(seps)) DEALLOCATE (seps)
      IF (allocated(ptr_out)) DEALLOCATE (ptr_out)
      IF (allocated(row_out)) DEALLOCATE (row_out)

   GOTO 80

10 WRITE(*,'(a,i4)') 'Allocation error during test section ', test
   GOTO 80

20 WRITE(*,'(a,i4)') 'Deallocation error during test section ', test

80 CONTINUE

END PROGRAM main

    SUBROUTINE testpr(n,perm)
!     .. Scalar Arguments ..
      INTEGER n
!     ..
!     .. Array Arguments ..
      INTEGER perm(n)
!     ..
!     .. Local Scalars ..
      INTEGER i, ip
!     .. Local Arrays
      INTEGER, ALLOCATABLE, DIMENSION (:) :: iw
!     ..
! Initialize array used to test for valid permutation.
      ALLOCATE (iw(n))
      iw(1:n) = 1
      DO i = 1, n
        ip = abs(perm(i))
        IF (iw(ip)==0) THEN
          WRITE (68,*) i, ip
          GO TO 10
        ELSE
          iw(ip) = 0
        END IF
      END DO
      DEALLOCATE (iw)
      RETURN
10    WRITE (6,'(A)') '**** Error in permutation'
      DEALLOCATE (iw)
    END SUBROUTINE testpr

    SUBROUTINE ndef(n,nz,a,ip,ind,w,sym)
! This subroutine augments the diagonal to make the matrix pos def
! Then large entries put in off-diagonal to make it a little not so.
! .. Parameters ..
      DOUBLE PRECISION zero, one
      PARAMETER (zero=0.0D0,one=1.0D0)
! ..
! .. Scalar Arguments ..
      INTEGER n, nz
      LOGICAL sym
! ..
! .. Array Arguments ..
      DOUBLE PRECISION a(nz), w(n)
      INTEGER ind(nz), ip(*)
! ..
! .. Local Scalars ..
      DOUBLE PRECISION rmax
      INTEGER i, i1, i2, idiag, ii, ioff, j, maxoff, mprint, numoff
! ..
! .. Intrinsic Functions ..
      INTRINSIC abs, max
! ..
      numoff = 0
! !! MAXOFF was 1 .. now 10 to try to have more 2 x 2 pivots
      maxoff = 10
      mprint = 6
      DO i = 1, n
        w(i) = zero
      END DO
      DO j = 1, n
        rmax = zero
        IF (sym) rmax = w(j)
        idiag = 0
        ioff = 0
        i1 = ip(j)
        i2 = ip(j+1) - 1
        IF (i2<i1) THEN
          GO TO 30
        ELSE
          DO ii = i1, i2
            i = ind(ii)
            rmax = rmax + abs(a(ii))
            w(i) = w(i) + abs(a(ii))
            IF (i==j) idiag = ii
            IF (i/=j) ioff = ii
          END DO
          IF (idiag==0) THEN
            GO TO 30
          ELSE
            a(idiag) = rmax + one
            IF (ioff/=0 .AND. numoff<maxoff) THEN
              a(ioff) = 1.1*a(idiag)
! !! added
              a(idiag) = zero
              numoff = numoff + 1
            END IF
          END IF
        END IF
      END DO
      IF ( .NOT. sym) THEN
        DO j = 1, n
          i1 = ip(j)
          i2 = ip(j+1) - 1
          DO ii = i1, i2
            i = ind(ii)
            IF (i==j) GO TO 10
          END DO
          GO TO 20
10        a(ii) = max(a(ii),w(i)+one)
20        CONTINUE
        END DO
      END IF
      GO TO 40
30    WRITE (mprint,fmt=90000) j
      nz = -1
40    RETURN
90000 FORMAT (' NO DIAGONAL ENTRY IN COLUMN ',I8,/, &
        ' SO CANNOT MAKE MATRIX DIAGONALLY DOMINANT')
    END SUBROUTINE ndef

    SUBROUTINE pdef(n,nz,a,ip,ind,w,sym)
! THIS SUBROUTINE AUGMENTS THE DIAGONAL TO MAKE THE MATRIX POS DEF.
! .. Parameters ..
      DOUBLE PRECISION zero, one
      PARAMETER (zero=0.0D0,one=1.0D0)
! ..
! .. Scalar Arguments ..
      INTEGER n, nz
      LOGICAL sym
! ..
! .. Array Arguments ..
      DOUBLE PRECISION a(nz), w(n)
      INTEGER ind(nz), ip(*)
! ..
! .. Local Scalars ..
      DOUBLE PRECISION rmax
      INTEGER i, i1, i2, idiag, ii, j, mprint
! ..
! .. Intrinsic Functions ..
      INTRINSIC abs, max
! ..
      mprint = 6
      DO i = 1, n
        w(i) = zero
      END DO
      DO j = 1, n
        rmax = zero
        IF (sym) rmax = w(j)
        idiag = 0
        i1 = ip(j)
        i2 = ip(j+1) - 1
        IF (i2<i1) THEN
          GO TO 30
        ELSE
          DO ii = i1, i2
            i = ind(ii)
            IF (i/=j) THEN
              rmax = rmax + abs(a(ii))
              w(i) = w(i) + abs(a(ii))
            END IF
            IF (i==j) idiag = ii
          END DO
          IF (idiag==0) THEN
            GO TO 30
          ELSE
            a(idiag) = rmax + one
          END IF
        END IF
      END DO
      IF ( .NOT. sym) THEN
        DO j = 1, n
          i1 = ip(j)
          i2 = ip(j+1) - 1
          DO ii = i1, i2
            i = ind(ii)
            IF (i==j) GO TO 10
          END DO
          GO TO 20
10        a(ii) = max(a(ii),w(i)+one)
20        CONTINUE
        END DO
      END IF
      GO TO 40
30    WRITE (mprint,fmt=90000) j
      nz = -1
40    RETURN
90000 FORMAT (' NO DIAGONAL ENTRY IN COLUMN ',I8,/, &
        ' SO CANNOT MAKE MATRIX DIAGONALLY DOMINANT')
    END SUBROUTINE pdef

