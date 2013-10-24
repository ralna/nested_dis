MODULE hsl_mc70_double
  USE hsl_mc79_integer
  USE hsl_mc78_integer
  USE hsl_zd11_double
  USE hsl_mc65_double
  USE hsl_fa14_double

      IMPLICIT NONE
      PRIVATE

! ---------------------------------------------------
! Precision
! ---------------------------------------------------
      INTEGER, PARAMETER :: myreal_mc70 = kind(1.0D0)

! ---------------------------------------------------
! Error flags
! ---------------------------------------------------
      INTEGER, PARAMETER :: mc70_err_memory_alloc = -1, & ! memory
! alloc error
        mc70_err_memory_dealloc = -2, & ! memory dealloc error
        mc70_err_n = -3 ! n<1

! ---------------------------------------------------
! Partition flags
! ---------------------------------------------------
      INTEGER, PARAMETER :: mc70_part1_flag = 0, & ! node in partition 1
        mc70_part2_flag = 2, & ! node in partition 2
        mc70_sep_flag = 1 ! node in separator
! ---------------------------------------------------
! Real version of partition flags
! ---------------------------------------------------
      REAL (myreal_mc70), PARAMETER :: mc70_part1_flag_real = 0.0, & ! node in partition 1
        mc70_part2_flag_real = 2.0, & ! node in partition 2
        mc70_sep_flag_real = 1.0

! ---------------------------------------------------
! Derived type definitions
! ---------------------------------------------------

! *****************************************************************
! ------ control ------
      TYPE, PUBLIC :: mc70_control
        INTEGER :: print_level = 0 ! amount of informational output required
        INTEGER :: unit_diagnostics = 6 ! stream number for diagnostic messages
        INTEGER :: unit_error = 6 ! stream number for error messages
        INTEGER :: ml = 2 ! Are we allowed to use a multilevel 
           ! strategy
           ! <= 0 : do not use multilevel
           ! == 1 : use multilevel
           ! >= 2 : automatic choice based on size of levelsets
        INTEGER :: nd_max_levels = 20 ! maximum number of ND levels to be performed
        INTEGER :: nd_switch = 50 ! switch to (halo)AMD if matrix size
                 ! is less than or equal to nd_switch 
        INTEGER :: partition_method = 1 ! Which partition method to use at coarsest level
           ! <=1 : Ashcraft method (half-level set)
           ! >=2 : Level-set method
        INTEGER :: refinement = 2 ! Which sort of refinement to use
           ! <1 : trim + fm
           ! =1 : DM + fm
           ! >1 : automatic choice
        INTEGER :: refinement_band = 4 ! band width for FM refinement. Values less than 1 
         ! mean that full FM refinement is done
        LOGICAL :: remove_dense = .TRUE. ! test the input for dense rows and place
                 ! them at the end of the ordering
        INTEGER :: ml_max_levels = 20 ! Max number of levels in the multilevel grid
        INTEGER :: ml_switch = 50 ! Stop coarsening once matrix has order at most ml_switch

! minimum and maximum grid reduction factors
! that must be achieved during coarsening.
! If cgrid%size is greater than ml_max_reduction*grid%size
! or cgrid%size is less than ml_min_reduction*grid%size
! then carry on coarsening
        REAL (myreal_mc70) :: ml_min_reduction = 0.01 ! size of next multigrid 
   !matrix must be greater than ml_min_reduction*(size of current multigrid matrix)
        REAL (myreal_mc70) :: ml_max_reduction = 0.9 ! size of next multigrid 
   !matrix must be less than ml_max_reduction*(size of current multigrid matrix)
        REAL (myreal_mc70) :: ratio = huge(mc70_sep_flag_real) ! Try to make sure that 
                                           ! max(P1,P2)/min(P1/P2) <= ratio

        REAL (myreal_mc70) :: ml_bandwidth = 1.0 ! Let B be matrix A after 
   ! dense rows removed and compressed. If ml>1 and bandwidth of B after RCM 
   ! ordering is larger than ml_bandwidth*order(B), continue as if ml=1; 
   ! otherwise continue as if ml=0; If B is separable, perform test on each 
   ! component separately use accordingly with each component. Note: RCM 
   ! ordering is not computed.
        LOGICAL :: block = .true.
        INTEGER :: expand = 0
        LOGICAL :: expand_scale = .true.
        INTEGER :: nstrt = 1  ! 1 = min degree
                               ! 2 = max degree
        INTEGER :: ww = 2
      END TYPE mc70_control

! *****************************************************************
! ------ info ------
      TYPE, PUBLIC :: mc70_info
        INTEGER :: flag = 0 ! error/warning flag
        INTEGER :: ndense = 0 ! holds number of dense rows 
        INTEGER :: stat = 0 ! holds Fortran stat parameter
        INTEGER :: nsuper = 0 ! holds number of supervariables + number of zero rows
        REAL (myreal_mc70) :: band = -1 ! holds L, where L is the size 
           ! of the largest level set at the top level of nested dissection. If 
           ! the matrix is reducible, then it holds the maximum value over all 
           ! of the irreducible components. 
           ! Not returned if control%ml==1.
      END TYPE mc70_info

! *****************************************************************

      TYPE mc70_multigrid

        INTEGER :: size ! size of this level (number of rows)
        TYPE (zd11_type), POINTER :: graph ! this level of matrix
        INTEGER, POINTER, DIMENSION (:) :: where ! where each row of this 
! level of matrix will go (ie ordering for this level)
        INTEGER, POINTER, DIMENSION (:) :: row_wgt !number of vertices 
! this vertex of the coarse graph matrix represents
        INTEGER :: level ! the level
        INTEGER :: part_div(2) ! number of vertices in each part
        TYPE (mc70_multigrid), POINTER :: coarse ! pointer to the coarse grid
        TYPE (mc70_multigrid), POINTER :: fine ! pointer to the fine grid
        TYPE (zd11_type), POINTER :: p ! the prolongation operator

      END TYPE mc70_multigrid
! *****************************************************************

      TYPE list_node_type
        INTEGER :: id
        TYPE (list_node_type), POINTER :: next
        TYPE (list_node_type), POINTER :: prev
      END TYPE list_node_type
! *****************************************************************

! this is the first node in the bucket, any new node will be added
! before this node and become the first node
! thus this is a last-in-first-out scheme
      TYPE list_node_pointer
        TYPE (list_node_type), POINTER :: current_node
      END TYPE list_node_pointer
! *****************************************************************

      TYPE queue_type
        INTEGER :: low_bound, up_bound
! maxnodes = max no. of nodes allowed (as dimension of mynode)
! maxgain = max gain so far in the buckets (gives the highest non-empty bucket)
! num_nodes = total number of nodes in the buckets
        INTEGER :: maxnodes, maxgain, num_nodes
! the starting node of each bucket (point to NULL if empty)
        TYPE (list_node_pointer), POINTER, DIMENSION (:) :: buckets
! the location of each node in the buckets
        TYPE (list_node_type), POINTER, DIMENSION (:) :: mynode
      END TYPE queue_type

! ---------------------------------------------------
! Interfaces
! ---------------------------------------------------
      INTERFACE mc70_order
        MODULE PROCEDURE mc70_nested_lower_double
      END INTERFACE

      INTERFACE mc70_order_full
        MODULE PROCEDURE mc70_nested_lower_upper_double
      END INTERFACE

   !   INTERFACE mc70_print_message
   !     MODULE PROCEDURE mc70_print_message_double
   !   END INTERFACE

      PUBLIC mc70_order, mc70_order_full

! ---------------------------------------------------
! The main code
! ---------------------------------------------------

    CONTAINS

! ---------------------------------------------------
! mc70_print_message
! ---------------------------------------------------
      SUBROUTINE mc70_print_message(flag,unit,context)
! Prints out errors and warnings according to value of flag

! flag: is an integer scaler of intent(in). It is the information flag
! whose corresponding error message is printed
        INTEGER, INTENT (IN) :: flag

! unit: is an optional integer scaler of intent(in). It is the unit number the
! error/warning message should be printed on
        INTEGER, INTENT (IN) :: unit

! context: is an optional assumed size character array of intent(in).
! It describes the context under which the error occured
        CHARACTER (len=*), INTENT (IN) :: context
        INTEGER :: length, p_unit

        p_unit = unit

        IF (p_unit<=0) RETURN

        IF (flag<0) THEN
          WRITE (p_unit,advance='yes',fmt='('' ERROR: '')')
        END IF

        length = len_trim(context)
        WRITE (p_unit,advance='no',fmt='('' '', a,'':'')') &
            context(1:length)

        SELECT CASE (flag)
        CASE (0)
          WRITE (p_unit,'(A)') ' successful completion'

        CASE (mc70_err_memory_alloc)
          WRITE (p_unit,'(A)') ' memory allocation failure'

        CASE (mc70_err_memory_dealloc)
          WRITE (p_unit,'(A)') ' memory deallocation failure'

        CASE (mc70_err_n)
          WRITE (p_unit,'(A)') ' n<1'
        END SELECT
        RETURN

      END SUBROUTINE mc70_print_message

! ---------------------------------------------------
! mc70_nested_lower
! ---------------------------------------------------
      SUBROUTINE mc70_nested_lower_double(n,ptr,row,perm,control,info,seps)
        INTEGER, INTENT(IN) :: n  ! number of rows in the matrix
        INTEGER, INTENT(IN) :: ptr(n+1) ! column pointers
        INTEGER, INTENT(IN) :: row(ptr(n+1)-1) ! row indices (lower triangle)
        INTEGER, INTENT(OUT) :: perm(n) ! permutation: row i becomes row perm(i)
        TYPE (mc70_control), INTENT(IN) :: control 
        TYPE (mc70_info), INTENT(INOUT) :: info
        INTEGER, INTENT(OUT),OPTIONAL :: seps(n)
                 ! seps(i) is -1 if vertex is not in a separator; otherwise it
                 ! is equal to l, where l is the nested dissection level at
                 ! which it became part of the separator

! ---------------------------------------------------
! LOCAL VARIABLES
        INTEGER :: a_n ! dimension of the expanded matrix
        INTEGER :: a_ne ! number off-diagonal entries stored in expanded matrix
        INTEGER :: i, j, k
        INTEGER :: unit_error ! unit on which to print errors
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER, ALLOCATABLE, DIMENSION(:) :: a_ptr ! a_ptr(i) will contain the 
!                 position in a_row that column i ends for the expanded matrix
        INTEGER, ALLOCATABLE, DIMENSION(:) :: a_row ! contains for row indices 
!                 for the pattern of the expanded matrix. The row indices of 
!                 column j are stored before those of column j+1, j=1,..,n-1.
        LOGICAL :: printe, printi, printd

! ---------------------------------------------
! Printing levels
        unit_diagnostics = control%unit_diagnostics
        unit_error = control%unit_error
        printe = (control%print_level>=0 .AND. unit_error>=0)
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
! ---------------------------------------------------
 
! Expand the matrix to hold all of its pattern
! Throughout the code we will use a modified compressed sparse column/row
! format a_ptr(i) will contain the position in a_row that column i begins. This 
! is to avoid lots of allocations during the nested part of the algorithm. a_ne
! contains the number of entries stored for the expanded matrix.

        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'MC70_nested_lower:'
        END IF

 

! Set the dimension of the expanded matrix
        a_n = n
        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a,i10)') 'n = ', n
        END IF
        IF (n .LT. 1) THEN
           info%flag = mc70_err_n
           IF (printe) CALL mc70_print_message(info%flag,unit_error, &
            'mc70_nested')
           RETURN
        END IF


! Allocate space to store pointers for expanded matrix
        ALLOCATE (a_ptr(a_n),STAT=info%stat)
        IF (info%stat/=0)  GOTO 10

! Fill a_col and a_ptr removing any diagonal entries
        a_ptr(:) = 0

! Set a_ptr(j) to hold no. nonzeros in column j
        DO j = 1,n
          DO k = ptr(j), ptr(j+1) - 1
            i = row(k)
            IF (j/=i) THEN
              a_ptr(i) = a_ptr(i) + 1
              a_ptr(j) = a_ptr(j) + 1
            END IF
          END DO
        END DO

! Set a_ptr(j) to point to where row indices will end in a_row
        DO j = 2,n
          a_ptr(j) = a_ptr(j-1) + a_ptr(j)
        END DO
        a_ne = a_ptr(a_n)
        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a,i10)') &
            'entries in expanded matrix with diags removed = ', a_ne
        END IF

! Allocate space to store row indices of expanded matrix
        ALLOCATE (a_row(a_ne),STAT=info%stat)
        IF (info%stat/=0)  GOTO 10

! Initialise all of a_row to 0
        a_row(:) = 0
           
! Fill a_row and a_ptr
        DO j = 1,n
            DO k = ptr(j), ptr(j+1) - 1
              i = row(k)
              IF (j/=i) THEN
                a_row(a_ptr(i)) = j
                a_row(a_ptr(j)) = i
                a_ptr(i) = a_ptr(i) - 1
                a_ptr(j) = a_ptr(j) - 1
              END IF
            END DO
        END DO

! Reset a_ptr to point to where column starts
        DO j = 1,a_n
           a_ptr(j) = a_ptr(j)+1
        END DO

        IF (printd) THEN
! Print out a_ptr and a_row
          WRITE (unit_diagnostics,'(a8)') 'a_ptr = '
          WRITE (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,a_n)
          WRITE (unit_diagnostics,'(a8)') 'a_row = '
          WRITE (unit_diagnostics,'(5i15)') (a_row(i),i=1,a_ne)
        ELSE IF (printi) THEN
! Print out first few entries of a_ptr and a_row
          WRITE (unit_diagnostics,'(a21)') 'a_ptr(1:min(5,a_n)) = '
          WRITE (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,min(5,a_n))
          WRITE (unit_diagnostics,'(a21)') 'a_row(1:min(5,a_ne)) = '
          WRITE (unit_diagnostics,'(5i15)') (a_row(i),i=1,min(5,a_ne))
        END IF
        IF (present(seps)) THEN
          CALL mc70_nested_both(a_n,a_ne,a_ptr,a_row,perm,control,info,seps)
        ELSE
          CALL mc70_nested_both(a_n,a_ne,a_ptr,a_row,perm,control,info)
        END IF

        IF (printd) THEN
! Print out perm
          WRITE (unit_diagnostics,'(a8)') 'perm = '
          WRITE (unit_diagnostics,'(5i15)') (perm(i),i=1,a_n)
        ELSE IF (printi) THEN
! Print out first few entries of perm
          WRITE (unit_diagnostics,'(a21)') 'perm(1:min(5,a_n)) = '
          WRITE (unit_diagnostics,'(5i15)') (perm(i),i=1,min(5,a_n))
        END IF

! Deallocate arrays
        DEALLOCATE (a_ptr,STAT=info%stat)
        IF (info%stat/=0)  GOTO 20
        DEALLOCATE (a_row,STAT=info%stat)
        IF (info%stat/=0)  GOTO 20

        info%flag = 0
        IF (printi .or. printd) THEN
          CALL mc70_print_message(info%flag,unit_diagnostics, &
            'mc70_nested')
        END IF
        RETURN
    
10      info%flag = mc70_err_memory_alloc
        IF (printe) CALL mc70_print_message(info%flag,unit_error, &
            'mc70_nested')
        RETURN
    
20      info%flag = mc70_err_memory_dealloc
        IF (printe) CALL mc70_print_message(info%flag,unit_error, &
            'mc70_nested')
        RETURN

      END SUBROUTINE mc70_nested_lower_double


! ---------------------------------------------------
! mc70_nested
! ---------------------------------------------------
      SUBROUTINE mc70_nested_lower_upper_double(n,ptr,row,perm,control,info,seps)
        INTEGER, INTENT(IN) :: n  ! number of rows in the matrix
        INTEGER, INTENT(IN) :: ptr(n+1) ! column pointers
        INTEGER, INTENT(IN) :: row(ptr(n+1)-1) ! row indices (lower triangle)
        INTEGER, INTENT(OUT) :: perm(n) ! permutation: row i becomes row perm(i)
        TYPE (mc70_control), INTENT(IN) :: control 
        TYPE (mc70_info), INTENT(INOUT) :: info
        INTEGER, INTENT(OUT), OPTIONAL :: seps(n)
                 ! seps(i) is -1 if vertex is not in a separator; otherwise it
                 ! is equal to l, where l is the nested dissection level at
                 ! which it became part of the separator

! ---------------------------------------------------
! LOCAL VARIABLES
        INTEGER :: a_n ! dimension of the expanded matrix
        INTEGER :: a_ne ! number off-diagonal entries stored in expanded matrix
        INTEGER :: i, j, k, l, ndiags
        INTEGER :: unit_error ! unit on which to print errors
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER, ALLOCATABLE, DIMENSION(:) :: a_ptr ! a_ptr(i) will contain the 
!                 position in a_row that column i ends for the expanded matrix
        INTEGER, ALLOCATABLE, DIMENSION(:) :: a_row ! contains for row indices 
!                 for the pattern of the expanded matrix. The row indices of 
!                 column j are stored before those of column j+1, j=1,..,n-1.
        LOGICAL :: printe, printi, printd

! ---------------------------------------------
! Printing levels
        unit_diagnostics = control%unit_diagnostics
        unit_error = control%unit_error
        printe = (control%print_level>=0 .AND. unit_error>=0)
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
! ---------------------------------------------------
 
! Expand the matrix to hold all of its pattern
! Throughout the code we will use a modified compressed sparse column/row
! format a_ptr(i) will contain the position in a_row that column i begins. This 
! is to avoid lots of allocations during the nested part of the algorithm. a_ne
! contains the number of entries stored for the expanded matrix.

        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'MC70_nested_lower_upper:'
        END IF

! Set the dimension of the expanded matrix
        a_n = n
        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a,i10)') 'n = ', n
        END IF
        IF (n .LT. 1) THEN
           info%flag = mc70_err_n
           IF (printe) CALL mc70_print_message(info%flag,unit_error, &
            'mc70_nested')
           RETURN
        END IF


! Work out how many diagonal entries need removing
        ndiags = 0
        DO j=1,n
          DO l = ptr(j),ptr(j+1)-1
            i = row(l)
            IF (i .EQ. j) ndiags = ndiags + 1
          END DO
        END DO
        a_ne = ptr(n+1)-1-ndiags

! Allocate space to store pointers and rows for expanded matrix
        ALLOCATE (a_ptr(a_n),STAT=info%stat)
        IF (info%stat/=0)  GOTO 10
        ALLOCATE (a_row(a_ne),STAT=info%stat)
        IF (info%stat/=0)  GOTO 10

        IF (ndiags .EQ. 0) THEN
          ! No diagonal entries so do direct copy
          a_ptr(1:a_n) = ptr(1:n)
          a_row(1:a_ne) = row(1:a_ne)

        ELSE
          ! Diagonal entries present
          k = 1
          DO i=1,n
            a_ptr(i) = k
            DO l = ptr(i),ptr(i+1)-1
              j = row(l)
              IF (i.NE.j) THEN
                 a_row(k) = j
                 k = k+1
              END IF
            END DO
          END DO
        END IF


        IF (printd) THEN
! Print out a_ptr and a_row
          WRITE (unit_diagnostics,'(a8)') 'a_ptr = '
          WRITE (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,a_n)
          WRITE (unit_diagnostics,'(a8)') 'a_row = '
          WRITE (unit_diagnostics,'(5i15)') (a_row(i),i=1,a_ne)
        ELSE IF (printi) THEN
! Print out first few entries of a_ptr and a_row
          WRITE (unit_diagnostics,'(a21)') 'a_ptr(1:min(5,a_n)) = '
          WRITE (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,min(5,a_n))
          WRITE (unit_diagnostics,'(a21)') 'a_row(1:min(5,a_ne)) = '
          WRITE (unit_diagnostics,'(5i15)') (a_row(i),i=1,min(5,a_ne))
        END IF
        IF (present(seps)) THEN
          CALL mc70_nested_both(a_n,a_ne,a_ptr,a_row,perm,control,info,seps)
        ELSE
          CALL mc70_nested_both(a_n,a_ne,a_ptr,a_row,perm,control,info)
        END IF

        IF (printd) THEN
! Print out perm
          WRITE (unit_diagnostics,'(a8)') 'perm = '
          WRITE (unit_diagnostics,'(5i15)') (perm(i),i=1,a_n)
        ELSE IF (printi) THEN
! Print out first few entries of perm
          WRITE (unit_diagnostics,'(a21)') 'perm(1:min(5,a_n)) = '
          WRITE (unit_diagnostics,'(5i15)') (perm(i),i=1,min(5,a_n))
        END IF

! Deallocate arrays
        DEALLOCATE (a_ptr,STAT=info%stat)
        IF (info%stat/=0)  GOTO 20
        DEALLOCATE (a_row,STAT=info%stat)
        IF (info%stat/=0)  GOTO 20

        info%flag = 0
        IF (printi .or. printd) THEN
          CALL mc70_print_message(info%flag,unit_diagnostics, &
            'mc70_nested_lower_upper')
        END IF
        RETURN
    
10      info%flag = mc70_err_memory_alloc
        IF (printe) CALL mc70_print_message(info%flag,unit_error, &
            'mc70_nested_lower_upper')
        RETURN
    
20      info%flag = mc70_err_memory_dealloc
        IF (printe) CALL mc70_print_message(info%flag,unit_error, &
            'mc70_nested_lower_upper')
        RETURN

      END SUBROUTINE mc70_nested_lower_upper_double

! ---------------------------------------------------
! mc70_nested
! ---------------------------------------------------
      SUBROUTINE mc70_nested_both(a_n,a_ne,a_ptr,a_row,perm,control,info,seps)
        INTEGER, INTENT(IN) :: a_n  ! number of rows in the matrix
        INTEGER, INTENT(IN) :: a_ne  ! number of entries in the matrix
        INTEGER, INTENT(INOUT) :: a_ptr(a_n) ! column pointers
        INTEGER, INTENT(INOUT) :: a_row(a_ne) ! row indices (lower and upper)
        INTEGER, INTENT(OUT) :: perm(a_n) ! permutation: row i becomes row perm(i)
        TYPE (mc70_control), INTENT(IN) :: control 
        TYPE (mc70_info), INTENT(INOUT) :: info
        INTEGER, INTENT(OUT),OPTIONAL :: seps(a_n) 
                 ! seps(i) is -1 if vertex is not in a separator; otherwise it
                 ! is equal to l, where l is the nested dissection level at
                 ! which it became part of the separator

! ---------------------------------------------------
! LOCAL VARIABLES
        INTEGER :: i, j, k, l,ll, a_ne_new, a_n_new, lwork, lirn
        INTEGER :: a_n_curr,a_ne_curr, num_zero_row
        INTEGER :: hamd_irn,hamd_ip,hamd_sep,hamd_perm,hamd_work,hamd_iperm
        INTEGER :: work_iperm, work_seps
        INTEGER :: unit_error ! unit on which to print errors
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER :: st, nsvar,svinfo
        INTEGER :: sv_ptr,sv_perm, sv_invp, sv_svar, sv_ptr2, sv_row2,sumweight
        INTEGER, ALLOCATABLE, DIMENSION(:) :: a_weight ! a_weight(i) will 
!                 contain the weight of variable (column) i ends for the 
!                 expanded matrix
        INTEGER, ALLOCATABLE, DIMENSION(:) :: iperm ! row iperm(i) will 
!                 become row i when matrix reordered
        INTEGER, ALLOCATABLE, DIMENSION(:) :: work ! space for doing work
        INTEGER, ALLOCATABLE, DIMENSION(:) :: svwork ! supervariable work space
        LOGICAL :: printe, printi, printd
        LOGICAL :: use_multilevel

! ---------------------------------------------
! Printing levels
        unit_diagnostics = control%unit_diagnostics
        unit_error = control%unit_error
        printe = (control%print_level>=0 .AND. unit_error>=0)
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
! ---------------------------------------------------
 
! Expand the matrix to hold all of its pattern
! Throughout the code we will use a modified compressed sparse column/row
! format a_ptr(i) will contain the position in a_row that column i begins. This 
! is to avoid lots of allocations during the nested part of the algorithm. a_ne
! contains the number of entries stored for the expanded matrix.

        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'MC70_nested:'
        END IF
        work_seps = 0
        work_iperm = 0

! Allocate iperm to have length n
        ALLOCATE (iperm(a_n),STAT=info%stat)
        IF (info%stat/=0)  GOTO 10

! Initialise iperm to identity permutation
        iperm(:) = (/ (i,i=1,a_n) /)
        IF (present(seps)) THEN
           seps(:) = -1
        END IF

! Check for dense rows
        IF (control%remove_dense) THEN
! Initially allocate work to have length 4*a_n
        ALLOCATE (work(4*a_n),STAT=info%stat)
        IF (info%stat/=0)  GOTO 10
        
        CALL mc70_dense_rows(a_n,a_ne,a_ptr,a_row,i,j,&
         iperm,work(1:4*a_n),control,info)
        a_n_new = i
        a_ne_new = j

        DEALLOCATE (work,STAT=info%stat)
        IF (info%stat/=0)  GOTO 20

        ELSE
        a_n_new = a_n
        a_ne_new = a_ne

        END IF
        
! Check whether matrix is diagonal
        IF (a_ne_new .eq. 0) THEN
          
          info%nsuper = a_n_new
          ! Create perm from iperm
          DO i = 1,a_n
           j = iperm(i)
           perm(j) = i
          END DO
          RETURN
        END IF


! Check for supervariables
        ALLOCATE (svwork(5*a_n_new+a_ne_new+2),STAT=info%stat)
        IF (info%stat/=0)  GOTO 10
        sv_ptr = 0
        sv_perm = sv_ptr + a_n_new+1
        sv_invp = sv_perm + a_n_new
        sv_svar = sv_invp + a_n_new
        sv_ptr2 = sv_svar + a_n_new
        sv_row2 = sv_ptr2 + a_n_new + 1
        

        svwork(sv_ptr+1:sv_ptr+a_n_new) = a_ptr(1:a_n_new)
        svwork(sv_ptr+a_n_new+1) = a_ne_new+1
        svwork(sv_perm+1:sv_perm+a_n_new) = (/ (i,i=1,a_n_new) /)
        svwork(sv_invp+1:sv_invp+a_n_new) = (/ (i,i=1,a_n_new) /)
        i = a_n_new
        CALL mc78_supervars(i,svwork(sv_ptr+1:sv_ptr+a_n_new+1),&
           a_row(1:a_ne_new),svwork(sv_perm+1:sv_perm+a_n_new),&
           svwork(sv_invp+1:sv_invp+a_n_new), nsvar,&
           svwork(sv_svar+1:sv_svar+a_n_new),st )
        num_zero_row = a_n_new - i
        info%nsuper = nsvar+num_zero_row
        IF (st/=0)  GOTO 10
         
        IF (nsvar+num_zero_row .EQ. a_n_new) THEN
          DEALLOCATE (svwork,STAT=info%stat)
          IF (info%stat/=0)  GOTO 20
          a_n_curr = a_n_new
          a_ne_curr = a_ne_new
        ALLOCATE (a_weight(a_n_curr),STAT=info%stat)
        IF (info%stat/=0)  GOTO 10

! Initialise a_weight          
        a_weight(:) = 1

        ELSE
          CALL mc78_compress_by_svar(a_n_new,svwork(sv_ptr+1:sv_ptr+a_n_new+1),&
           a_row(1:a_ne_new),svwork(sv_invp+1:sv_invp+a_n_new), nsvar,&
           svwork(sv_svar+1:sv_svar+a_n_new),&
           svwork(sv_ptr2+1:sv_ptr2+a_n_new+1),&
           a_ne_new,svwork(sv_row2+1:sv_row2+a_ne_new),svinfo, st )

             a_n_curr = nsvar
             ! Fill a_ptr removing any diagonal entries
             a_ptr(:) = 0

             ! Set a_ptr(j) to hold no. nonzeros in column j
             DO j = 1,a_n_curr
               DO k =  svwork(sv_ptr2+j),  svwork(sv_ptr2+j+1) - 1
                i = svwork(sv_row2+k)
                IF (j.lt.i) THEN
                 a_ptr(i) = a_ptr(i) + 1
                 a_ptr(j) = a_ptr(j) + 1
                END IF
               END DO
             END DO

             ! Set a_ptr(j) to point to where row indices will end in a_row
             DO j = 2,a_n_curr
               a_ptr(j) = a_ptr(j-1) + a_ptr(j)
             END DO
             a_ne_curr = a_ptr(a_n_curr)
             ! Initialise all of a_row to 0
             a_row(1:a_ne_curr) = 0
           
             ! Fill a_row and a_ptr
             DO j = 1,a_n_curr
              DO k =  svwork(sv_ptr2+j),  svwork(sv_ptr2+j+1) - 1
               i = svwork(sv_row2+k)
               IF (j.lt.i) THEN
                a_row(a_ptr(i)) = j
                a_row(a_ptr(j)) = i
                a_ptr(i) = a_ptr(i) - 1
                a_ptr(j) = a_ptr(j) - 1
               END IF
              END DO
             END DO

             ! Reset a_ptr to point to where column starts
             DO j = 1,a_n_curr
              a_ptr(j) = a_ptr(j)+1
             END DO

            ALLOCATE (a_weight(a_n_curr+num_zero_row),STAT=info%stat)
            IF (info%stat/=0)  GOTO 10

            ! Initialise a_weight          
            a_weight(1:a_n_curr) = svwork(sv_svar+1:sv_svar+a_n_curr)
            a_weight(a_n_curr+1:a_n_curr+num_zero_row) = 1

            ! Add zero rows/cols to matrix
            a_ptr(a_n_curr+1:a_n_curr+num_zero_row) = a_ne_curr+1
            a_n_curr = a_n_curr+num_zero_row

           ! set svwork(sv_svar+1:sv_svar+a_n_new) such that svwork(sv_svar+i) 
           ! points to the end of the list of variables in sv_invp for 
           ! supervariable i
           DO i = 2,nsvar
             svwork(sv_svar+i) = svwork(sv_svar+i)+svwork(sv_svar+i-1)
           END DO
           j = svwork(sv_svar+nsvar)
           DO i = 1,num_zero_row
             svwork(sv_svar+nsvar+i) = j+1
             j = j+1
           END DO

          

        END IF

! Carryout nested dissection on matrix once dense rows removed

        IF (control%nd_max_levels .LE. 0 .or. a_n_curr .LE. &
               max(2,control%nd_switch)) THEN
          ! Apply AMD to matrix
! Allocate work to have length 5*a_n_curr+a_ne_curr
        ALLOCATE (work(11*a_n_curr+a_ne_curr+a_n_new),STAT=info%stat)
        IF (info%stat/=0)  GOTO 10
          lirn = a_ne_curr + a_n_curr
          hamd_irn = 0
          hamd_ip = hamd_irn + lirn
          hamd_sep = hamd_ip + a_n_curr
          hamd_perm = hamd_sep + a_n_curr
          hamd_work = hamd_perm + a_n_curr
          hamd_iperm = hamd_work + 7*a_n_curr

          work(hamd_irn+1:hamd_irn+a_ne_curr) = a_row(1:a_ne_curr)          
          work(hamd_irn+a_ne_curr+1:hamd_irn+lirn) = 0
          work(hamd_ip+1:hamd_ip+a_n_curr) = a_ptr(1:a_n_curr)
          work(hamd_sep+1:hamd_sep+a_n_curr) = mc70_sep_flag-1
          call hamd(a_n_curr,a_ne_curr,lirn,work(hamd_irn+1:hamd_irn+lirn),&
             work(hamd_ip+1:hamd_ip+a_n_curr),work(hamd_sep+1:hamd_sep+a_n_curr),&
             work(hamd_perm+1:hamd_perm+a_n_curr),&
             work(hamd_work+1:hamd_work+7*a_n_curr))

          ! Extract perm from hamd to apply to iperm
          IF (nsvar+num_zero_row .EQ. a_n_new) THEN
            DO i = 1,a_n_curr
             j = work(hamd_perm+i)
             work(hamd_work+i) = iperm(j)            
            END DO
            iperm(1:a_n_curr) = work(hamd_work+1:hamd_work+a_n_curr)
          ELSE
            DO i = 1,a_n_curr
             j = work(hamd_perm+i)
             work(hamd_work+i) = j           
            END DO
            ! Expand to matrix before supervariables detected
            k = 1
            DO i = 1,a_n_curr
              j = work(hamd_work+i)
              IF (j .EQ. 1) THEN
               ll = 1
              ELSE
               ll = svwork(sv_svar+j-1)+1
              END IF
              DO l=ll,svwork(sv_svar+j)
                 work(hamd_iperm+k) = iperm(svwork(sv_invp+l))
                 k = k+1
              END DO
            END DO
            
            iperm(1:a_n_new) = work(hamd_iperm+1:hamd_iperm+a_n_new)
          END IF

        ELSE
          ! Apply ND to matrix


          IF (nsvar+num_zero_row .EQ. a_n_new) THEN
           ! Allocate work to have length 14*a_n_curr+a_ne_curr
           ALLOCATE (work(a_n_new+14*a_n_curr+a_ne_curr),STAT=info%stat)
          ELSE
           ! Allocate work to have length a_n_new+14*a_n_curr+a_ne_curr
           ALLOCATE (work(3*a_n_new+14*a_n_curr+a_ne_curr),STAT=info%stat)
           work_iperm = 14*a_n_curr+a_ne_curr + a_n_new
           work(work_iperm+1:work_iperm+a_n_curr) = (/ (i,i=1,a_n_curr) /)
           IF (present(seps)) THEN
             work_seps = work_iperm+a_n_curr
             work(work_seps+1:work_seps+a_n_curr) = -1
           END IF
          END IF
          IF (info%stat/=0)  GOTO 10

          use_multilevel = .true.
          IF (nsvar+num_zero_row .EQ. a_n_new) THEN
             sumweight = sum(a_weight(1:a_n_curr))
             lwork = 12*a_n_curr+sumweight+a_ne_curr
             IF (present(seps)) THEN
               call mc70_nested_internal(a_n_curr,a_ne_curr,a_ptr(1:a_n_curr),&
                a_row(1:a_ne_curr),a_weight(1:a_n_curr),sumweight,iperm(1:a_n_curr),&
                work(1:lwork),work(lwork+1:lwork+a_n_curr),&
                work(lwork+a_n_curr+1:lwork+2*a_n_curr),0,control,info,&
                use_multilevel,seps(1:a_n_curr))
             ELSE
               call mc70_nested_internal(a_n_curr,a_ne_curr,a_ptr(1:a_n_curr),&
                a_row(1:a_ne_curr),a_weight(1:a_n_curr),sumweight,iperm(1:a_n_curr),&
                work(1:lwork),work(lwork+1:lwork+a_n_curr),&
                work(lwork+a_n_curr+1:lwork+2*a_n_curr),0,control,info,&
                use_multilevel)

             END IF
          ELSE
             sumweight = sum(a_weight(1:a_n_curr))
             lwork = 12*a_n_curr+sumweight+a_ne_curr
             IF (present(seps)) THEN
               call mc70_nested_internal(a_n_curr,a_ne_curr,a_ptr(1:a_n_curr),&
                a_row(1:a_ne_curr),a_weight(1:a_n_curr),sumweight,&
                work(work_iperm+1:work_iperm+a_n_curr),&
                work(1:lwork),work(lwork+1:lwork+a_n_curr),&
                work(lwork+a_n_curr+1:lwork+2*a_n_curr),0,control,info,&
                use_multilevel,work(work_seps+1:work_seps+a_n_curr))
             ELSE
               call mc70_nested_internal(a_n_curr,a_ne_curr,a_ptr(1:a_n_curr),&
                a_row(1:a_ne_curr),a_weight(1:a_n_curr),sumweight,&
                work(work_iperm+1:work_iperm+a_n_curr),&
                work(1:lwork),work(lwork+1:lwork+a_n_curr),&
                work(lwork+a_n_curr+1:lwork+2*a_n_curr),0,control,info,&
                use_multilevel)

             END IF
          END IF

          IF (nsvar+num_zero_row .LT. a_n_new) THEN
           IF (present(seps)) THEN
             ! Expand and reorder seps
            DO i =1,a_n_curr
             j = work(work_iperm+i)
              IF (j .EQ. 1) THEN
               ll = 1
              ELSE
               ll = svwork(sv_svar+j-1)+1
              END IF
              DO l=ll,svwork(sv_svar+j)
                 seps(svwork(sv_invp+l)) = work(work_seps+i)
              END DO
            END DO
           END IF
            
            ! Expand iperm to matrix before supervariables detected
            k = a_n_new
            DO i = a_n_curr,1,-1
              j = work(work_iperm+i)
              IF (j .EQ. 1) THEN
               ll = 1
              ELSE
               ll = svwork(sv_svar+j-1)+1
              END IF
              DO l=ll,svwork(sv_svar+j)
                 work(work_iperm+k) = iperm(svwork(sv_invp+l))
                 k = k-1
              END DO
            END DO
            
            iperm(1:a_n_new) = work(work_iperm+1:work_iperm+a_n_new)
            
          ELSE
           IF (present(seps)) THEN
            ! reorder seps! 
            DO i =1,a_n_curr
             j = iperm(i)
             work(j) = seps(i)
            END DO
            seps(1:a_n_curr) = work(1:a_n_curr)
           END IF
          END IF
        END IF
! Create perm from iperm
        DO i = 1,a_n
           j = iperm(i)
           perm(j) = i
        END DO

! Deallocate arrays
        IF (nsvar+num_zero_row .LT. a_n_new) THEN
          DEALLOCATE (svwork,STAT=info%stat)
          IF (info%stat/=0)  GOTO 20
        END IF
        DEALLOCATE (a_weight,STAT=info%stat)
        IF (info%stat/=0)  GOTO 20
        DEALLOCATE (iperm,STAT=info%stat)
        IF (info%stat/=0)  GOTO 20
        DEALLOCATE (work,STAT=info%stat)
        IF (info%stat/=0)  GOTO 20

        info%flag = 0
        IF (printi .or. printd) THEN
          CALL mc70_print_message(info%flag,unit_diagnostics, &
            'mc70_nested_both')
        END IF
        RETURN
    
10      info%flag = mc70_err_memory_alloc
        IF (printe) CALL mc70_print_message(info%flag,unit_error, &
            'mc70_nested')
        RETURN
    
20      info%flag = mc70_err_memory_dealloc
        IF (printe) CALL mc70_print_message(info%flag,unit_error, &
            'mc70_nested')
        RETURN

      END SUBROUTINE mc70_nested_both

! ---------------------------------------------------
! mc70_dense_rows
! ---------------------------------------------------
! Identifies and removes dense rows
      SUBROUTINE mc70_dense_rows(a_n_in,a_ne_in,a_ptr,a_row,a_n_out,a_ne_out,&
         iperm,work,control,info)
        INTEGER, INTENT(IN) :: a_n_in ! dimension of subproblem before dense 
             ! rows removed
        INTEGER, INTENT(IN) :: a_ne_in ! no. nonzeros of subproblem before 
             ! dense rows removed
        INTEGER, INTENT(INOUT) :: a_ptr(a_n_in) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. This is then 
             ! used to hold positions for submatrices after dense row removed
        INTEGER, INTENT(INOUT) :: a_row(a_ne_in) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.This is then used to hold row indices for
             ! submatrices after partitioning
        INTEGER, INTENT(INOUT) :: iperm(a_n_in) ! On input, iperm(i) contains 
             ! the row in the original matrix (when mc70_nested was called) that 
             ! row i in this sub problem maps to. On output, this is updated to
             ! reflect the computed permutation.
        INTEGER, INTENT(OUT) :: a_n_out ! dimension of subproblem after dense 
             ! rows removed
        INTEGER, INTENT(OUT) :: a_ne_out ! no. nonzeros of subproblem after 
             ! dense rows removed
        INTEGER, INTENT(OUT) :: work(4*a_n_in) ! Used during the algorithm to 
             ! reduce need for allocations. The output is garbage.
        TYPE (mc70_control), INTENT(IN) :: control 
        TYPE (mc70_info), INTENT(INOUT) :: info

! ---------------------------------------------
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER :: deg, prev, next, dense ! pointers into work array
        INTEGER :: ndense ! number of dense rows found
        INTEGER :: max_deg ! maximum degree
        INTEGER :: degree,i,j,k,l,inext,ilast,l1,l2,m,m1
        LOGICAL :: printi, printd

! ---------------------------------------------
! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
! ---------------------------------------------------
        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Find and remove dense rows'
        END IF

        ! Set pointers into work array
        deg = 0
        prev = deg+a_n_in
        next = prev + a_n_in
        dense = next + a_n_in

        ! By the end of this loop work(dense+i) will be
        !       0 if row is not dense
        !       <0 otherwise. The larger the number, the earlier the row was 
        !          was determined to be dense.
        ndense = 0
        max_deg = 0
        work(deg+1:deg+a_n_in) = 0

        ! Calculate degree of each row before anything removed
        DO i = 1, a_n_in
          k = a_ptr(i)
          IF (i<a_n_in) THEN
           degree = a_ptr(i+1) - k
          ELSE
           degree = a_ne_in - a_ptr(a_n_in)+1
          END IF
          work(dense+i) = degree
          IF (degree .ne. 0) THEN
            max_deg = max(max_deg,degree)
            CALL add_to_list(i,degree)
          END IF
        END DO
        degree = max_deg
        a_n_out = a_n_in
        a_ne_out = a_ne_in

        DO WHILE (real(degree) - real(a_ne_out)/real(a_n_out)>= &
          20*(real(a_n_out-1)/real(a_n_out))*log(real(a_n_out)) &
          .AND. degree>0)
   !     DO WHILE (real(degree) - real(a_ne_out)/real(a_n_out)>= &
   !       300*(real(a_n_out-1)/real(a_n_out)) &
   !       .AND. degree>0)
          i = work(deg+degree)
          ndense = ndense+1
          work(dense+i) = - ndense
          CALL remove_from_list(i,degree)
          ! update degrees of adjacent vertices
          IF (i<a_n_in) THEN
            l = a_ptr(i+1) - 1
          ELSE
            l = a_ne_in
          END IF
          DO k = a_ptr(i), l
              j = a_row(k)
              IF (work(dense+j)>0) THEN
                CALL remove_from_list(j,work(dense+j))
                work(dense+j) = work(dense+j) - 1
                IF (work(dense+j)>0) THEN
                  CALL add_to_list(j,work(dense+j))
                END IF
              END IF
          END DO
          a_n_out = a_n_out-1
          a_ne_out = a_ne_out - 2*degree
          IF (work(deg+degree)==0) THEN
          ! Find next largest degree
              degree = degree - 1
              DO
                IF (degree==0 ) EXIT
                IF (work(deg+degree)>0) EXIT
                degree = degree - 1
              END DO
            END IF
        END DO

        ! By the end of this loop work(dense+i) will be
        !       >=0 if row is not dense
        !       <0 otherwise. The larger the number, the earlier the row was 
        !          was determined to be dense.
        !!!!!

        IF (ndense>0) THEN
        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(i10,a)') ndense, ' dense rows detected'
        END IF
        info%ndense = ndense

            a_n_out = 0
            l = a_n_in+1
            DO i = 1, a_n_in
              k = work(dense+i)
              IF (k>=0) THEN
                a_n_out = a_n_out + 1
                work(dense+i) = a_n_out
                work(next+a_n_out) = i
              ELSE
                work(next+l+k) = i
              END IF
            END DO
            
            k=1
            j=1
            
            DO i=1,a_n_in
              l1 = a_ptr(i)
              IF (i<a_n_in) THEN
              l2 = a_ptr(i+1)-1
              ELSE
                l2 = a_ne_in
              END IF
              if (work(dense+i) >= 0) THEN
                a_ptr(j) = k
                DO l=l1,l2
                  m = a_row(l)
                  m1 = work(dense+m)
                  IF (m1 >= 0) THEN
                    a_row(k) = m1
                    k = k+1
                  END IF                
                END DO
                j = j+1              
              END IF     
            END DO
            a_ptr(j) = k
        IF (printd) THEN
! Print out a_ptr and a_row
          WRITE (unit_diagnostics,'(a11)') 'a_n_out = '
          WRITE (unit_diagnostics,'(i15)') a_n_out
          WRITE (unit_diagnostics,'(a11)') 'a_ne_out = '
          WRITE (unit_diagnostics,'(i15)') a_ne_out
          WRITE (unit_diagnostics,'(a8)') 'a_ptr = '
          WRITE (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,a_n_out)
          WRITE (unit_diagnostics,'(a8)') 'a_row = '
          WRITE (unit_diagnostics,'(5i15)') (a_row(i),i=1,a_ne_out)
        ELSE IF (printi) THEN
! Print out first few entries of a_ptr and a_row
          WRITE (unit_diagnostics,'(a11)') 'a_n_out = '
          WRITE (unit_diagnostics,'(i15)') a_n_out
          WRITE (unit_diagnostics,'(a11)') 'a_ne_out = '
          WRITE (unit_diagnostics,'(i15)') a_ne_out
          WRITE (unit_diagnostics,'(a21)') 'a_ptr(1:min(5,a_n_out)) = '
          WRITE (unit_diagnostics,'(5i15)') (a_ptr(i),i=1,min(5,a_n_out))
          WRITE (unit_diagnostics,'(a21)') 'a_row(1:min(5,a_ne_out)) = '
          WRITE (unit_diagnostics,'(5i15)') (a_row(i),i=1,min(5,a_ne_out))
        END IF
        ELSE
          
            a_n_out = a_n_in
            a_ne_out = a_ne_in
            work(next+1:next+a_n_in) = (/ (i,i=1,a_n_in)/)
        END IF

        DO i=1,a_n_in
           j = work(next+i)
           work(next+i) = iperm(j)
        END DO
        
        DO i=1,a_n_in
           iperm(i) = work(next+i)
        END DO
        
        info%flag = 0
        IF (printi .or. printd) THEN
          CALL mc70_print_message(info%flag,unit_diagnostics, &
            'mc70_dense_rows')
        END IF

      CONTAINS
        SUBROUTINE remove_from_list(irm,ig)
          INTEGER :: irm, ig

          inext = work(next+irm)
          ilast = work(prev+irm)
          IF (ilast==0) THEN
            work(deg+ig) = inext
            IF (inext/=0) work(prev+inext) = 0
          ELSE
            work(next+ilast) = inext
            IF (inext/=0) work(prev+inext) = ilast
          END IF
        END SUBROUTINE remove_from_list

        SUBROUTINE add_to_list(irm,ig)
          INTEGER :: irm, ig

          inext = work(deg+ig)
          work(deg+ig) = irm
          work(next+irm) = inext
          IF (inext/=0) work(prev+inext) = irm
          work(prev+irm) = 0
        END SUBROUTINE add_to_list

      END SUBROUTINE mc70_dense_rows

! ---------------------------------------------------
! mc70_nested_internal
! ---------------------------------------------------
! Does the recursive nested dissection

      RECURSIVE SUBROUTINE mc70_nested_internal(a_n,a_ne,a_ptr,a_row,a_weight,&
         sumweight,iperm,work,work_comp_n,work_comp_nz,level,&
         control,info,use_multilevel,seps)

        INTEGER, INTENT(IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT(IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT(INOUT) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. This is then 
             ! used to hold positions for submatrices after partitioning
        INTEGER, INTENT(INOUT) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.This is then used to hold row indices for
             ! submatrices after partitioning
        INTEGER, INTENT(INOUT) :: a_weight(a_n) ! On input a_weight(i) contains 
             ! weight of column i. This is then 
             ! used to hold the weights for submatrices after partitioning
        INTEGER, INTENT(IN) :: sumweight ! sum entries in a_weight (unchanged)
        INTEGER, INTENT(INOUT) :: iperm(a_n) ! On input, iperm(i) contains the
             ! row in the original matrix (when mc70_nested was called) that 
             ! row i in this sub problem maps to. On output, this is updated to
             ! reflect the computed permutation.
        INTEGER, INTENT(OUT) :: work_comp_n(a_n)
        INTEGER, INTENT(OUT) :: work_comp_nz(a_n)
        INTEGER, INTENT(OUT) :: work(12*a_n+sumweight+a_ne) ! Used during the 
             ! algorithm to reduce need for allocations. The output is garbage.
        INTEGER, INTENT(IN) :: level ! which level of nested dissection is this
        TYPE (mc70_control), INTENT(IN) :: control 
        TYPE (mc70_info), INTENT(INOUT) :: info
        LOGICAL, INTENT(INOUT) :: use_multilevel
        INTEGER, INTENT(INOUT),OPTIONAL :: seps(a_n)
                 ! seps(i) is -1 if vertex i of permuted submatrix is not in a 
                 ! separator; otherwise it is equal to l, where l is the nested
                 ! dissection level at which it became part of the separator

! ---------------------------------------------
! Local variables
        INTEGER :: i,j,k,l,m,s
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER :: num_components ! Number of independent components found
        INTEGER :: compwork
        INTEGER :: sumweight_sub
        INTEGER :: offset_ptr, offset_row
        INTEGER :: a_n1, a_n2, a_ne1, a_ne2
        INTEGER :: hamd_irn,hamd_ip,hamd_sep,hamd_perm,hamd_work,lirn
        LOGICAL :: printi, printd

! ---------------------------------------------
! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
 
        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a,i6)') 'Nested dissection level ', level
        END IF

! Check whether matrix is diagonal and act accordingly
        IF (a_ne .EQ. 0) THEN
          IF (printi .or. printd) THEN
            WRITE (unit_diagnostics,'(a)') ' '
            WRITE (unit_diagnostics,'(a)') 'Submatrix is diagonal'
          END IF
          RETURN
        END IF


! Check whether max number of levels has been reached or if matrix size is below
! nd_switch
        IF (level .GE. control%nd_max_levels .OR. a_n .LE. max(2,control%nd_switch)) &
           GOTO 30
          CALL mc70_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,level,a_n1,a_n2,&
            a_ne1,a_ne2,iperm,work(1:12*a_n+sumweight+a_ne),control,info,&
            use_multilevel)
        
        IF (a_n1 .EQ. a_n) THEN
           GOTO 30
        END IF
        

        IF (a_n1 .ne. 0 .and. a_n2.ne. 0 .and. a_n1+a_n2.eq.a_n) THEN
          ! matrix is reducible
          compwork = 0
          ! work array needs to be total length 5*a_n+a_ne
          call mc70_find_indep_comps(a_n,a_ne,a_ptr,a_row,a_weight,iperm,&
             num_components,work_comp_n(1:a_n),work_comp_nz(1:a_n),&
             work(compwork+1:compwork+3*a_n+a_ne),control,info)
          IF (num_components.EQ.1) THEN
             k = control%nd_max_levels
          ELSE
             k = level
          END IF

            ! Apply the ND to each component - do not test for indep comps 
            offset_ptr = a_n+1
            offset_row = a_ne+1
            i = num_components
            ! work from last component so that space at end of work_comp_n and 
            ! work_comp_nz can be reused without worrying about overwriting 
            ! important data
           ! write(*,*) work_comp_n(1:i)
           ! write(*,*) work_comp_nz(1:i)
            j = a_n
            DO WHILE (i .GE. 1)
              l = work_comp_n(i)
              m = work_comp_nz(i)
              s = sum(a_weight(offset_ptr-l:offset_ptr-1))
              IF (m .GT. 0) THEN
               ! Matrix not diagonal
               IF (present(seps)) THEN
                CALL mc70_nested_internal(l,m,&
                 a_ptr(offset_ptr-l:offset_ptr-1),&
                 a_row(offset_row-m:offset_row-1),&
                 a_weight(offset_ptr-l:offset_ptr-1),&
                 s,&
                 iperm(offset_ptr-l:offset_ptr-1),&
                 work(compwork+1:compwork+12*l+s+m),&
                 work_comp_n(j-l+1:j),&
                 work_comp_nz(j-l+1:j),k,control,info,&
                 use_multilevel,seps(offset_ptr-l:offset_ptr-1))
               ELSE
                CALL mc70_nested_internal(l,m,&
                 a_ptr(offset_ptr-l:offset_ptr-1),&
                 a_row(offset_row-m:offset_row-1),&
                 a_weight(offset_ptr-l:offset_ptr-1),&
                 s,&
                 iperm(offset_ptr-l:offset_ptr-1),&
                 work(compwork+1:compwork+12*l+s+m),&
                 work_comp_n(j-l+1:j),&
                 work_comp_nz(j-l+1:j),k,control,info,&
                 use_multilevel)

               END IF
              END IF
              offset_ptr = offset_ptr-l
              offset_row = offset_row-m
              j = j-l
              i = i-1
            END DO
          RETURN
        END IF  
        IF (present(seps)) THEN
          seps(a_n1+a_n2+1:a_n) = level
        END IF
        IF (a_n1 .GT. max(2,control%nd_switch)) THEN
           sumweight_sub = sum(a_weight(1:a_n1))
           IF (present(seps)) THEN
            CALL mc70_nested_internal(a_n1,a_ne1,a_ptr(1:a_n1),a_row(1:a_ne1),&
             a_weight(1:a_n1),sumweight_sub,iperm(1:a_n1),&
             work(1:12*a_n1+sumweight_sub+a_ne1),work_comp_n(1:a_n1),&
             work_comp_nz(1:a_n1),level+1,control,info,use_multilevel,&
             seps(1:a_n1))
           ELSE
            CALL mc70_nested_internal(a_n1,a_ne1,a_ptr(1:a_n1),a_row(1:a_ne1),&
             a_weight(1:a_n1),sumweight_sub,iperm(1:a_n1),&
             work(1:12*a_n1+sumweight_sub+a_ne1),work_comp_n(1:a_n1),&
             work_comp_nz(1:a_n1),level+1,control,info,use_multilevel)
           END IF

        END IF

        IF (a_n2 .GT. max(2,control%nd_switch)) THEN
           IF (a_n1 .GT. max(2,control%nd_switch)) THEN
            sumweight_sub = sum(a_weight(a_n1+1:a_n1+a_n2))
            IF (present(seps)) THEN
             CALL mc70_nested_internal(a_n2,a_ne2,a_ptr(a_n1+1:a_n1+a_n2),&
             a_row(a_ne1+1:a_ne1+a_ne2),&
             a_weight(a_n1+1:a_n1+a_n2),sumweight_sub,iperm(a_n1+1:a_n1+a_n2),&
             work(1:12*a_n2+sumweight_sub+a_ne2),work_comp_n(1:a_n2),&
             work_comp_nz(1:a_n2),level+1,control,info,use_multilevel,&
             seps(a_n1+1:a_n1+a_n2))
            ELSE
             CALL mc70_nested_internal(a_n2,a_ne2,a_ptr(a_n1+1:a_n1+a_n2),&
             a_row(a_ne1+1:a_ne1+a_ne2),&
             a_weight(a_n1+1:a_n1+a_n2),sumweight_sub,iperm(a_n1+1:a_n1+a_n2),&
             work(1:12*a_n2+sumweight_sub+a_ne2),work_comp_n(1:a_n2),&
             work_comp_nz(1:a_n2),level+1,control,info,use_multilevel)
            END IF
           ELSE
            sumweight_sub = sum(a_weight(a_n1+1:a_n1+a_n2))
            IF (present(seps)) THEN
             CALL mc70_nested_internal(a_n2,a_ne2,a_ptr(1:a_n2),&
             a_row(1:a_ne2),&
             a_weight(a_n1+1:a_n1+a_n2),sumweight_sub,iperm(a_n1+1:a_n1+a_n2),&
             work(1:12*a_n2+sumweight_sub+a_ne2),work_comp_n(1:a_n2),&
             work_comp_nz(1:a_n2),level+1,control,info,use_multilevel,&
             seps(a_n1+1:a_n1+a_n2))
            ELSE
             CALL mc70_nested_internal(a_n2,a_ne2,a_ptr(1:a_n2),&
             a_row(1:a_ne2),&
             a_weight(a_n1+1:a_n1+a_n2),sumweight_sub,iperm(a_n1+1:a_n1+a_n2),&
             work(1:12*a_n2+sumweight_sub+a_ne2),work_comp_n(1:a_n2),&
             work_comp_nz(1:a_n2),level+1,control,info,use_multilevel)

            END IF

           END IF

        END IF
        GOTO 20
        RETURN

        ! No partition found or max number of levels have been reached
30      CONTINUE
          ! Apply AMD to matrix
         IF (a_n.GT.1) THEN
          lirn = a_ne + a_n
          hamd_irn = 0
          hamd_ip = hamd_irn + lirn
          hamd_sep = hamd_ip + a_n
          hamd_perm = hamd_sep + a_n
          hamd_work = hamd_perm + a_n

          work(hamd_irn+1:hamd_irn+a_ne) = a_row(1:a_ne)          
          work(hamd_irn+a_ne+1:hamd_irn+lirn) = 0
          work(hamd_ip+1:hamd_ip+a_n) = a_ptr(1:a_n)
          work(hamd_sep+1:hamd_sep+a_n) = mc70_sep_flag-1
         ! write(*,*) work(hamd_ip+1:hamd_ip+a_n),a_ne+1
         ! write(*,*) work(hamd_irn+1:hamd_irn+a_ne)
          call hamd(a_n,a_ne,lirn,work(hamd_irn+1:hamd_irn+lirn),&
             work(hamd_ip+1:hamd_ip+a_n),work(hamd_sep+1:hamd_sep+a_n),&
             work(hamd_perm+1:hamd_perm+a_n),&
             work(hamd_work+1:hamd_work+7*a_n))

          ! Extract perm from hamd to apply to iperm
          DO i = 1,a_n
             j = work(hamd_perm+i)
             work(hamd_work+i) = iperm(j)            
          END DO
             iperm(1:a_n) = work(hamd_work+1:hamd_work+a_n)
         ELSE
             iperm(1) = 1
         END IF
        
20      info%flag = 0
        IF (printi .or. printd) THEN
          CALL mc70_print_message(info%flag,unit_diagnostics, &
            'mc70_nested_internal')
        END IF
        RETURN

      END SUBROUTINE mc70_nested_internal

! ---------------------------------------------------
! mc70_find_indep_comps
! ---------------------------------------------------
! Finds and forms independent components in a matrix
      SUBROUTINE mc70_find_indep_comps(a_n,a_ne,a_ptr,a_row,a_weight,iperm,&
             comp_num,compsizes,compnzs,work,control,info)
      INTEGER, INTENT(IN) :: a_n ! size of matrix
      INTEGER, INTENT(IN) :: a_ne ! no. nonzeros in matrix
      INTEGER, INTENT(INOUT) :: a_ptr(a_n) ! On entry, column ptrs for input 
        ! matrix.
        ! On exit, a_ptr(1:compsizes(1)) contains column ptrs for compontent 1;
        ! a_ptr(compsizes(1)+1:compsizes(1)+compsizes(2)) contains column ptrs 
        ! for compontent 2; etc.
      INTEGER, INTENT(INOUT) :: a_row(a_ne) ! On entry, row indices for input 
        ! matrix.
        ! On exit, a_row(1:compnzs(1)) contains row indices for compontent 1;
        ! a_ptr(compnzs(1)+1:compnzs(1)+compnzs(2)) contains row indices 
        ! for compontent 2; etc.
      INTEGER, INTENT(INOUT) :: a_weight(a_n) ! On entry, a_weight(i) contains
        ! weight of column i for input matrix.
        ! On exit, a_weight(1:compsizes(1)) contains column weights for 
        ! compontent 1;
        ! a_weight(compsizes(1)+1:compsizes(1)+compsizes(2)) contains column 
        ! weights or compontent 2; etc.
        INTEGER, INTENT(INOUT) :: iperm(a_n) ! On input, iperm(i) contains the
             ! row in the original matrix (when mc70_nested was called) that 
             ! row i in this sub problem maps to. On output, this is updated to
             ! reflect the computed permutation.
      INTEGER, INTENT(OUT) :: comp_num ! number independent components found
      INTEGER, INTENT(OUT) :: compsizes(a_n) ! compsizes(i) will contain the 
        ! size of compontent i
      INTEGER, INTENT(OUT) :: compnzs(a_n) ! compnzs(i) will contain the 
        ! number of nonzeros in compontent i
      INTEGER, INTENT(OUT) :: work(3*a_n+a_ne) ! used as work arrays during 
        ! computation
        TYPE (mc70_control), INTENT(IN) :: control 
        TYPE (mc70_info), INTENT(INOUT) :: info

! ---------------------------------------------
! Local variables
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER :: mask, ptr_temp, row_temp, front ! pointers into work array
        INTEGER :: root
        INTEGER :: front_sta, front_sto ! start and end of front list
        INTEGER :: num_assigned ! number of cols that have been assigned to a 
          ! component
        INTEGER :: i, j,l,u, v,p
        INTEGER :: ptr_temp_pos, row_temp_pos, front_pos
        LOGICAL :: printi, printd
       
! ---------------------------------------------
! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)

        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Start independent component search'
        END IF

        mask = 0
        front = mask + a_n
        ptr_temp = front + a_n
        row_temp = ptr_temp + a_n  ! total length of work array 3*a_n+a_ne

        ! Initialise arrays
        work(1:3*a_n) = 0
        compsizes(1:a_n) = 0
        compnzs(1:a_n) = 0

        ! Check from all possible roots
        num_assigned = 0
        comp_num = 0
        DO root = 1,a_n
          IF (work(mask+root) .EQ. 0) THEN
            comp_num = comp_num+1
            front_sta = num_assigned+1
            front_sto = front_sta
            work(front+front_sta) = root
            compsizes(comp_num) = 1
            compnzs(comp_num) = 0
            work(mask+root) = compsizes(comp_num)
            num_assigned = num_assigned+1

            DO WHILE (front_sto-front_sta .GE. 0)
              DO i = front_sta,front_sto
                ! pick vertex from front
                v = work(front+i)
                ! update compnzs
                IF (v .LT. a_n) THEN
                   l = a_ptr(v+1)
                ELSE
                   l = a_ne + 1
                END IF
                   compnzs(comp_num) = compnzs(comp_num) + l - a_ptr(v)
                DO j = a_ptr(v),l-1
                   ! pick a neighbour
                   u = a_row(j)
                   IF (work(mask+u) .NE. 0) cycle
                   ! found unmasked vertex
                   compsizes(comp_num) = compsizes(comp_num)+1
                   num_assigned = num_assigned + 1
                   ! mask this vertex
                   work(mask+u) = compsizes(comp_num)
                   ! add vertex to component
                   work(front+num_assigned) = u
                END DO
              END DO
              front_sta = front_sto+1
              front_sto = num_assigned
            END DO
          END IF
        END DO

        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(i5,a)') comp_num, ' components found'
        END IF

       IF (comp_num .GT. 1) THEN
        ! Reorder matrix into block diagonal form
        ! Indexing will be local to each block
        ptr_temp_pos = 1
        row_temp_pos = 1
        front_pos = 1
        DO l=1,comp_num  ! for each component
          work(ptr_temp+ptr_temp_pos) = 1
          ptr_temp_pos = ptr_temp_pos+1
          DO u = 1,compsizes(l)
            v = work(front+front_pos) ! for each column in the component
            front_pos = front_pos+1
            IF (v .EQ. a_n) THEN
              p = a_ne
            ELSE
              p = a_ptr(v+1)-1
            END IF
            DO i = a_ptr(v),p ! for each nonzero in the column
              work(row_temp+row_temp_pos) = work(mask+a_row(i))
              row_temp_pos = row_temp_pos + 1
            END DO
            IF (u .LT. compsizes(l)) THEN
              work(ptr_temp+ptr_temp_pos) = work(ptr_temp+ptr_temp_pos-1) + &
                    p+1 - a_ptr(v)
              ptr_temp_pos = ptr_temp_pos+1
            END IF
          END DO
        END DO

        a_ptr(1:a_n) = work(ptr_temp+1:ptr_temp+a_n)       
        a_row(1:a_ne) = work(row_temp+1:row_temp+a_ne)

        ! Reorder iperm and a_weight
        DO front_pos = 1,a_n
          j = work(front+front_pos)
          work(front+front_pos) = iperm(j)
          work(mask+front_pos) = a_weight(j)
        END DO
        iperm(1:a_n) = work(front+1:front+a_n)
        a_weight(1:a_n) = work(mask+1:mask+a_n)
       END IF
        
        info%flag = 0
        IF (printi .or. printd) THEN
          CALL mc70_print_message(info%flag,unit_diagnostics, &
            'mc70_find_indep_comps')
        END IF
        RETURN

      END SUBROUTINE mc70_find_indep_comps

! ---------------------------------------------------
! mc70_partition matrix
! ---------------------------------------------------
! Partition the matrix and if one (or more) of the generated submatrices is 
! small enough, apply halo amd

      SUBROUTINE mc70_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,level,a_n1,&
         a_n2,a_ne1,a_ne2,iperm,work,control,info,use_multilevel)

        INTEGER, INTENT(IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT(IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT(INOUT) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. This is then 
             ! used to hold positions for submatrices after partitioning
        INTEGER, INTENT(INOUT) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.This is then used to hold row indices for
             ! submatrices after partitioning
        INTEGER, INTENT(INOUT) :: a_weight(a_n) ! On input a_weight(i) contains 
             ! the weight of column i. This is then used to hold the weights for
             ! the submatrices after partitioning.
        INTEGER, INTENT(IN) :: sumweight ! Sum entries in a_weight. Unchanged.
        INTEGER, INTENT(IN) :: level ! Current nested dissection level
        INTEGER, INTENT(OUT) :: a_n1, a_n2 ! size of the two submatrices
        INTEGER, INTENT(OUT) :: a_ne1, a_ne2 ! no. nonzeros in two submatrices
        INTEGER, INTENT(INOUT) :: iperm(a_n) ! On input, iperm(i) contains the
             ! row in the original matrix (when mc70_nested was called) that 
             ! row i in this sub problem maps to. On output, this is updated to
             ! reflect the computed permutation.
        LOGICAL, INTENT(INOUT) :: use_multilevel
        INTEGER, INTENT(OUT) :: work(12*a_n+sumweight+a_ne)
        TYPE (mc70_control), INTENT(IN) :: control 
        TYPE (mc70_info), INTENT(INOUT) :: info

! ---------------------------------------------
! Local variables
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        LOGICAL :: printi, printd
        INTEGER :: partition_ptr ! pointer into work array
        INTEGER :: work_ptr ! pointer into work array
        INTEGER :: part_ptr ! pointer into work array
        INTEGER :: a_ptr_sub_ptr ! pointer into work array
        INTEGER :: a_row_sub_ptr ! pointer into work array
        INTEGER :: partition_method
        INTEGER :: i,j,k,p,q,l
        INTEGER :: a_weight_1,a_weight_2,a_weight_sep,aw1,aw2,aws,an1,an2,ans
        INTEGER ::aw1o,aw2o,awso,an1o,an2o,anso,aw1t,aw2t,awst,an1t,an2t,anst
        INTEGER :: ref_method
        INTEGER :: a_n1_new,a_n2_new,a_weight_1_new,a_weight_2_new,a_weight_sep_new
        INTEGER :: nstrt,nend
        REAL(myreal_mc70) :: ratio,tau_best,tau
        LOGICAL :: imbal
        

! ---------------------------------------------
! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
 
        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Start finding a partition'
        END IF

        ! If matrix is full, then don't partition
        IF (a_ne .EQ. a_n*(a_n-1)) THEN
           a_n1 = a_n
           a_n2 = 0
           a_ne1 = a_ne
           a_ne2 = 0
           GOTO 10
        END IF

        ! Find the partition
        IF (control%partition_method .LE. 1) THEN
          partition_method = 1
        ELSE 
          partition_method = 2
        END IF


        partition_ptr = 0 ! length a_n
        work_ptr = partition_ptr + a_n ! max length 9*a_n+sumweight+a_ne/2
     !  DO p = 1,a_n
     !   DO q = p+1,a_n
        nstrt = -1
        nend = -1
        IF (partition_method.EQ.1) THEN
           ! Ashcraft method
           CALL mc70_ashcraft(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,&
               level,nstrt,nend,a_n1,&
               a_n2,a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work(work_ptr+1:work_ptr+9*a_n+sumweight+a_ne/2),control,&
               info%flag,info%band,use_multilevel)
        ELSE
           ! Level set method
           CALL mc70_level_set(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,&
               level,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n), &
               work(work_ptr+1:work_ptr+4*a_n),control,info%flag,info%band,&
               use_multilevel)

        END IF
        IF (a_n1+a_n2.EQ. a_n) THEN
             RETURN
        END IF

        IF (a_n1 .NE. 0 .AND. a_n2 .NE. 0 .AND. a_n .GT. 3 ) THEN
           IF (.not. use_multilevel) THEN
          ! Refine the partition
          !  write(*,*) '******************'
          ! write(*,*) 'input',a_n1,a_n2,a_n-a_n1-a_n2,a_weight_1,&
          !      a_weight_2,a_weight_sep, &
          !      (real(a_weight_sep)/real(a_weight_1))/real(a_weight_2),&
          !      real(a_weight_sep)*(1.0+0.3*abs(real(a_weight_1-a_weight_2))/real(sumweight))
           IF (control%refinement .GE. 2) THEN
              IF (min(a_weight_1,a_weight_2) + a_weight_sep .LT. &
               max(a_weight_1,a_weight_2)) THEN
                 ref_method = 2
              ELSE
                 ref_method = 1
              END IF
           ELSE
              IF (control%refinement .LT. 1) THEN
                 ref_method = 1
              ELSE
                 ref_method = 2
              END IF
           END IF
            
            IF (ref_method .EQ. 1) THEN
             IF (control%block) THEN
              CALL mc70_refine_block_trim(a_n,a_ne,a_ptr,a_row,a_weight,&
               sumweight,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control) 
             ELSE
              CALL mc70_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control) 
             END IF
         
             CALL mc70_refine(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control,ref_method) 
            ELSE
             CALL mc70_refine(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control,ref_method) 
             CALL mc70_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control) 

            END IF
          ! Refine the partition
         !   write(*,*) '******************'
         !  write(*,*) 'output',a_n1,a_n2,a_n-a_n1-a_n2,a_weight_1,&
         !       a_weight_2,a_weight_sep, &
         !       (real(a_weight_sep)/real(a_weight_1))/real(a_weight_2),&
         !       real(a_weight_sep)*(1.0+0.3*abs(real(a_weight_1-a_weight_2))/real(sumweight))
           IF (control%expand.GT.0) THEN
               ratio = control%ratio
               IF (ratio .GT. real(sumweight-2)) THEN
                 imbal = .FALSE.
               ELSE
                 imbal = .TRUE.
               END IF
             CALL cost_function(a_weight_1,a_weight_2,a_weight_sep,sumweight,&
              ratio,imbal,tau_best)
             a_n1_new = a_n1
             a_n2_new = a_n2
             a_weight_1_new = a_weight_1
             a_weight_2_new = a_weight_2
             a_weight_sep_new = a_weight_sep
           END IF
           IF (control%expand_scale) THEN
             k = control%expand-level
           ELSE
             k = control%expand
           END IF
            part_ptr = work_ptr+8*a_n+sumweight+a_ne/2+1
            work(part_ptr+1:part_ptr+a_n) = work(partition_ptr+1:partition_ptr+a_n)
           DO i=1,k
          !  CALL expand_partition_simple(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,&
          !     a_n1_new,a_n2_new,a_weight_1_new,a_weight_2_new,&
          !     a_weight_sep_new,work(part_ptr+1:part_ptr+a_n),&
          !     work(work_ptr+1:work_ptr+a_n),control)
          ! write(*,*) 'zzz',a_weight_1_new,a_weight_2_new,&
          !     a_weight_sep_new

            CALL expand_partition_kinks(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,&
               3,10.0_myreal_mc70,1.0_myreal_mc70,a_n1_new,a_n2_new,&
                a_weight_1_new,a_weight_2_new,&
               a_weight_sep_new,work(part_ptr+1:part_ptr+a_n),&
               work(work_ptr+1:work_ptr+5*a_n),control)

            IF (ref_method .EQ. 1) THEN
             IF (control%block) THEN
              CALL mc70_refine_block_trim(a_n,a_ne,a_ptr,a_row,a_weight,&
               sumweight,a_n1_new,a_n2_new,a_weight_1_new,a_weight_2_new,&
               a_weight_sep_new,work(part_ptr+1:part_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control) 
             ELSE
              CALL mc70_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,&
               a_n1_new,a_n2_new,&
               a_weight_1_new,a_weight_2_new,a_weight_sep_new,&
               work(part_ptr+1:part_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control) 
             END IF
             !IF (i.EQ.k) THEN
             CALL mc70_refine(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1_new,a_n2_new,&
               a_weight_1_new,a_weight_2_new,a_weight_sep_new,&
               work(part_ptr+1:part_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control,ref_method) 
             !END IF
            ELSE
             CALL mc70_refine(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,&
               a_n1_new,a_n2_new,&
               a_weight_1_new,a_weight_2_new,a_weight_sep_new,&
               work(part_ptr+1:part_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control,ref_method) 
             CALL mc70_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,&
               a_n1_new,a_n2_new,&
               a_weight_1_new,a_weight_2_new,a_weight_sep_new,&
               work(part_ptr+1:part_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control) 

            END IF

            CALL cost_function(a_weight_1_new,a_weight_2_new,a_weight_sep_new,&
              sumweight,ratio,imbal,tau)
            IF (tau .LT. tau_best) THEN
               tau_best = tau
               work(partition_ptr+1:partition_ptr+a_n) = &
                   work(part_ptr+1:part_ptr+a_n)
               a_n1 = a_n1_new
               a_n2 = a_n2_new
               a_weight_1 = a_weight_1_new
               a_weight_2 = a_weight_2_new
               a_weight_sep = a_weight_sep_new
            END IF
           END DO


           END IF  
        ELSE
          GOTO 10
        END IF
     !   write(*,*) a_weight_1,a_weight_2,a_weight_sep
     !  END DO
     !  END DO

        IF ((a_n1 .le. max(2,control%nd_switch) .and. &
            a_n2 .le. max(2,control%nd_switch))) THEN
          !  apply halo amd to submatrics
           CALL hamd_both_old1(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,&
            work(partition_ptr+1:partition_ptr+a_n),iperm,a_weight,&
            work(work_ptr+1:work_ptr+12*a_n+a_ne))
           a_ne1 = 0
           a_ne2 = 0

        ELSE IF (a_n1 .le. max(2,control%nd_switch)) THEN
          !  apply halo amd to [A1, B1'; B1, I] using two levels
          !  return A2 and apply ND to it 
          CALL hamd_one(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,&
            work(partition_ptr+1:partition_ptr+a_n),iperm,&
            a_weight,a_ne2,work(work_ptr+1:work_ptr+11*a_n+a_ne))
          a_ne1 = 0


        ELSE IF (a_n2 .le. max(2,control%nd_switch)) THEN
          !  apply halo amd to [A2, B2'; B2, I] using two levels
          !  return A1 and apply ND to it 
          CALL hamd_one(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,&
            work(partition_ptr+1:partition_ptr+a_n),iperm,a_weight,&
            a_ne1,work(work_ptr+1:work_ptr+11*a_n+a_ne))
          a_ne2 = 0

         ELSE 
          !  return A1 and A2 and apply ND to them
            a_ptr_sub_ptr = work_ptr + a_n ! length a_n
            a_row_sub_ptr = a_ptr_sub_ptr + a_n ! length a_ne
          
            CALL extract_both_matrices(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,&
             work(partition_ptr+1:partition_ptr+a_n1+a_n2),a_ne1,a_ne2,&
             work(a_ptr_sub_ptr+1:a_ptr_sub_ptr+a_n1+a_n2),a_ne,&
             work(a_row_sub_ptr+1:a_row_sub_ptr+a_ne),&
             work(work_ptr+1:work_ptr+a_n))

            ! Copy extracted matrices
            a_ptr(1:a_n1+a_n2) = work(a_ptr_sub_ptr+1:a_ptr_sub_ptr+a_n1+a_n2)
            a_row(1:a_ne1+a_ne2) = &
              work(a_row_sub_ptr+1:a_row_sub_ptr+a_ne1+a_ne2)

            ! Update iperm and a_weight
            DO i = 1,a_n
              j = work(partition_ptr+i)
              work(a_ptr_sub_ptr+i) = iperm(j)
              work(work_ptr+i) = a_weight(j)
            END DO 
            DO i = 1,a_n
              iperm(i) = work(a_ptr_sub_ptr+i)
              a_weight(i) = work(work_ptr+i)   
            END DO 
            
         END IF

        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Partition found'
        END IF
        GOTO 20

10      IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'No partition found'
        END IF

20        info%flag = 0
        IF (printi .or. printd) THEN
          CALL mc70_print_message(info%flag,unit_diagnostics, &
            'mc70_partition')
        END IF
        RETURN
      
      END SUBROUTINE mc70_partition

! ---------------------------------------------------
! mc70_ashcraft
! ---------------------------------------------------
! Partition the matrix using the Ashcraft method      
    RECURSIVE  SUBROUTINE mc70_ashcraft(a_n,a_ne,a_ptr,a_row,a_weight,&
        sumweight,level,nstrt,nend,a_n1,&
        a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,control,info,&
        band,use_multilevel)

        INTEGER, INTENT(IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT(IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT(IN) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. 
        INTEGER, INTENT(IN) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.
        INTEGER, INTENT(IN) :: a_weight(a_n) ! On input a_weight(i) contains 
             ! the weight of column i.
        INTEGER, INTENT(IN) :: sumweight ! sum of entries in a_weight
        INTEGER, INTENT(IN) :: level ! current level of nested dissection
        INTEGER, INTENT(INOUT) :: nstrt ! Index of suggest node to start 
           ! search. If negative or greater than a_n, code determines where to 
           ! start search and where end node is
        INTEGER, INTENT(INOUT) :: nend! Index of suggested end node. If 
           ! either nend or nstrt are negative or greater than a_n, code 
           ! determines nend.
        INTEGER, INTENT(OUT) :: a_n1, a_n2 ! size of the two submatrices
        INTEGER, INTENT(OUT) :: a_weight_1,a_weight_2,a_weight_sep ! Weighted 
             ! size of partitions and separator
        INTEGER, INTENT(OUT) :: partition(a_n) ! First a_n1 entries will contain
             ! list of (local) indices in partition 1; next a_n2 entries will 
             ! contain list of (local) entries in partition 2; entries in 
             ! separator are listed at the end
        INTEGER, INTENT(OUT) :: work(9*a_n+sumweight+a_ne/2) ! used as work array
        TYPE (mc70_control), INTENT(IN) :: control 
        INTEGER, INTENT(INOUT) :: info
        REAL (myreal_mc70), INTENT(INOUT) :: band ! If level = 0, then on 
          ! output band = max(band,100*L/a_n), where L is the size of the 
          ! largest levelset
        LOGICAL,INTENT(INOUT) :: use_multilevel ! are we allowed to use a multilevel
                ! partitioning strategy

! ---------------------------------------------
! Local variables
        INTEGER :: i,j,dptr,p1sz,p2sz,sepsz,lwork,k
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER :: mask_p,level_p,level_ptr_p,level2_p,level2_ptr_p,work_p
        INTEGER :: num_levels_nend ! no. levels in structure rooted at nend
        INTEGER :: num_levels_nstrt ! no. levels in structure rooted at nstrt
        INTEGER :: num_entries ! no. entries in level set structure
        INTEGER :: best_sep_start
        INTEGER :: distance
        INTEGER :: distance_ptr
        INTEGER :: lwidth,mindeg,degree,max_search
        INTEGER :: ww
        INTEGER :: ml_max_levels ! max no. multigrid levels
        REAL(myreal_mc70) :: bestval
        REAL(myreal_mc70) :: val
        REAL(myreal_mc70) :: ratio
        LOGICAL :: printi, printd
        LOGICAL :: imbal
        TYPE(fa14_seed) :: seed

! ---------------------------------------------
! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
 
        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Use two-sided level set method'
        END IF
        ratio = control%ratio
        IF (ratio .GT. real(sumweight-2)) THEN
          imbal = .FALSE.
        ELSE
          imbal = .TRUE.
        END IF
        p2sz=0
        sepsz=0

        IF (control%ml==1 .AND. use_multilevel) THEN
           use_multilevel = .TRUE.
            GOTO 222
        END IF

        ! Find pseudoperipheral nodes nstart and nend, and the level structure 
        ! rooted at nend
        mask_p = 0 ! size a_n
        level_ptr_p = mask_p + a_n ! size a_n
        level_p = level_ptr_p + a_n ! size a_n
        level2_ptr_p = level_p + a_n ! size a_n
        level2_p = level2_ptr_p + a_n ! size a_n
        work_p = level2_p + a_n ! size 2*a_n 

       IF (nstrt .LT. 1 .OR. nstrt.GT.a_n) THEN
        nend = -1

        ! Choose nstrt
        SELECT CASE (control%nstrt)

        CASE (1)
         ! node with minimum degree
         mindeg = a_n+1
         DO i = 1,a_n
          IF (i .LT. a_n) THEN
            degree = a_ptr(i+1)-a_ptr(i)
          ELSE
            degree = a_ne + 1 - a_ptr(i)
          END IF
          IF (degree .LT. mindeg) THEN
            mindeg = degree
            nstrt = i
          END IF
         END DO

        CASE (2)
         ! node with maximum degree
         mindeg = -1
         DO i = 1,a_n
          IF (i .LT. a_n) THEN
            degree = a_ptr(i+1)-a_ptr(i)
          ELSE
            degree = a_ne + 1 - a_ptr(i)
          END IF
          IF (degree .GT. mindeg) THEN
            mindeg = degree
            nstrt = i
          END IF
         END DO
 
        CASE (3)
         nstrt = 1

        CASE (4)
         call fa14_random_integer(seed,a_n,nstrt)
        
        CASE (5)
         call fa14_random_integer(seed,a_n,nstrt)
       !  call fa14_random_integer(seed,a_n,nstrt)
       !  call fa14_random_integer(seed,a_n,nstrt)
         DO 
            IF (nend .GT.0 .AND. nend .LT. a_n+1 .AND. nend .NE. nstrt) exit           
            call fa14_random_integer(seed,a_n,nend)       
            call fa14_random_integer(seed,a_n,nend)
         END DO
        END SELECT
       END IF
       max_search = 1
        
       IF (nend .LT. 1 .OR. nend.GT. a_n) THEN
        CALL mc70_find_pseudo(a_n,a_ne,a_ptr,a_row,&
           work(level_ptr_p+1:level_ptr_p+a_n),work(level_p+1:level_p+a_n),&
           nstrt,nend,max_search,work(work_p+1:work_p+2*a_n),num_levels_nend,&
           num_entries,lwidth)
       ELSE
        ! Find level structure rooted at nend
        work(mask_p+1:mask_p+a_n) = 1
        CALL mc70_level_struct(nend,a_n,a_ne,a_ptr,a_row,a_n,&
          work(mask_p+1:mask_p+a_n),&
          work(level_ptr_p+1:level_ptr_p+a_n),work(level_p+1:level_p+a_n),&
          num_levels_nend,lwidth,num_entries)
       ! write(*,*) work(level_ptr_p+1:level_ptr_p+num_levels_nend)
       END IF

        IF (num_entries .lt. a_n) THEN
           ! matrix is separable
           a_n1 = num_entries
           a_n2 = a_n - a_n1
           a_weight_sep = 0
           a_weight_1 = 0
           work(work_p+1:work_p+a_n) = 0
           DO i = 1,num_entries
             j = work(level_p+i)
             partition(i) = j
             work(work_p+j) = 1
             a_weight_1 = a_weight_1 + a_weight(j)
           END DO
           a_weight_2 = sumweight - a_weight_1
           j = num_entries+1
           DO i = 1,a_n
             IF (work(work_p+i) .eq. 0) THEN
               partition(j) = i
               j = j+1
             END IF

           END DO                     

           RETURN
        END IF
   
        ! ***********************************************************************  
        IF (level .EQ. 0) THEN
             band = max(band,100.0*real(lwidth,kind(1.0D0))/real(a_n,kind(1.0D0)))
        END IF
        IF (control%ml_max_levels .LE. 0 .OR. control%ml .LT. 1) THEN
            use_multilevel = .false.
        END IF
        IF (control%ml>=2 .AND. level .EQ. 0) THEN
         IF (real(lwidth,kind(1.0D0))/real(a_n,kind(1.0D0)) .LE. &
               control%ml_bandwidth   ) THEN
            use_multilevel = .false.
         ELSE
            use_multilevel = .true.
         END IF
        END IF

222     CONTINUE

        IF (use_multilevel) THEN
           ml_max_levels = control%ml_max_levels
          lwork = 9*a_n+sumweight+a_ne/2
          CALL multilevel_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,&
          partition,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
          control,info,lwork,work,ml_max_levels)
          RETURN
        END IF

        ! Find level structure rooted at nstrt
        work(mask_p+1:mask_p+a_n) = 1
        CALL mc70_level_struct(nstrt,a_n,a_ne,a_ptr,a_row,a_n,&
          work(mask_p+1:mask_p+a_n),&
          work(level2_ptr_p+1:level2_ptr_p+a_n),work(level2_p+1:level2_p+a_n),&
          num_levels_nstrt,lwidth,num_entries)
        !write(*,*) work(level2_ptr_p+1:level2_ptr_p+num_levels_nstrt)
      !  write(*,*) nend, nstrt, work(level_p+1), work(level2_p+1)


        ! Calculate difference in distances from nstart and nend, and set up 
        ! lists D_i of nodes with same distance
          distance_ptr = work_p
          distance = mask_p
        CALL mc70_distance(a_n,num_levels_nend,&
          work(level_ptr_p+1:level_ptr_p+a_n),&
          work(level_p+1:level_p+a_n),num_levels_nstrt,&
          work(level2_ptr_p+1:level2_ptr_p+a_n),work(level2_p+1:level2_p+a_n),&
          work(distance_ptr+1:distance_ptr+2*a_n-1), &
          work(distance+1:distance+a_n))
        ! write(*,*)  work(distance_ptr+1:distance_ptr+2*a_n-1)

        ! Do not need the information in work(level_ptr_p+1:level2_ptr_p)
        ! Calculate total weight in each distance level
         dptr = level_ptr_p + a_n
         work(dptr+1-a_n:dptr+a_n-1) = 0
         DO i = 1-num_levels_nstrt,num_levels_nend-2
           DO j = work(distance_ptr+a_n+i),work(distance_ptr+a_n+i+1)-1
             work(dptr+i) = work(dptr+i) + a_weight(work(distance+j))
           END DO
        !   write(*,*) i,work(dptr+i)
         END DO
         i = num_levels_nend-1
         DO j = work(distance_ptr+a_n+i),a_n
             work(dptr+i) = work(dptr+i) + a_weight(work(distance+j))
         END DO
         !  write(*,*) i,work(dptr+i)

         ww = max(2,control%ww)
        !write(*,*) ' '
        
         ! Find first possible separator
         p1sz = 0
         DO i = 1-num_levels_nstrt,num_levels_nend-ww-2
           p1sz = p1sz + work(dptr+i)
           sepsz = work(dptr+i+1)
           DO j = 1,ww-1
             sepsz = sepsz + work(dptr+i+1+j)
           END DO
           p2sz = sumweight - p1sz-sepsz
           !  write(*,*) p1sz,sepsz,p2sz
           IF (p1sz .GT. 0 .AND. sepsz .GT. 0) THEN
           !  write(*,*) p1sz,sepsz,p2sz
             EXIT
           END IF
         END DO

         IF (i+ww .GE. num_levels_nend-1 .OR. p2sz.EQ.0) THEN
           ! Not possible to find separator
           ! This can only be invoked for a fully connected graph. The partition
           ! subroutine checks for this case so it should never be called
           a_n1 = a_n
           a_n2 = 0
           partition(1:a_n) = (/ (i,i=1,a_n) /)
           RETURN
         END IF

         ! Entries in level i will form first possible partition
         best_sep_start = i+1
         a_n1 = work(distance_ptr+a_n+i+1) - 1
         a_n2 = a_n - work(distance_ptr+a_n+i+1+ww) +1
         a_weight_1 = p1sz
         a_weight_2 = p2sz
         a_weight_sep = sepsz
           CALL cost_function(p1sz,p2sz,&
              sepsz,sumweight,ratio,imbal,bestval)

         ! Search for best separator using tau
         
         DO j = i+1,num_levels_nend-ww-2
           p1sz = p1sz+work(dptr+j)
           sepsz = work(dptr+j+1)
           DO k = 1,ww-1
             sepsz = sepsz + work(dptr+j+1+k)
           END DO
           p2sz = sumweight-sepsz-p1sz
           IF (p2sz.EQ.0) exit
          !   write(*,*) j,num_levels_nend-1,p1sz,sepsz,p2sz, work(dptr+j+1:dptr+j+ww+1),sumweight,&
          !        sum(work(dptr+1-num_levels_nstrt:dptr+num_levels_nend-1))
           CALL cost_function(p1sz,p2sz,&
              sepsz,sumweight,ratio,imbal,val)
           IF (val .LT. bestval ) THEN
             bestval = val
             best_sep_start = j+1
             a_n1 = work(distance_ptr+a_n+j+1) - 1
             a_n2 = a_n - work(distance_ptr+a_n+j+1+ww) +1     
             a_weight_1 = p1sz
             a_weight_2 = p2sz
             a_weight_sep = sepsz
           END IF
         END DO
        ! write(*,*) a_weight_1,a_weight_2,a_weight_sep,sumweight - a_weight_1-a_weight_2-a_weight_sep

         ! Rearrange partition
         ! Entries in partition 1
         j = 1
         DO i = 1,work(distance_ptr+a_n+best_sep_start)-1
           partition(j) = work(distance+i)
           j = j+1
         END DO
        !  write(*,*) 'a_n1',a_n1,j-1

         ! Entries in partition 2
         DO i = work(distance_ptr+a_n+best_sep_start+ww),a_n
           partition(j) = work(distance+i)
           j = j+1
         END DO
        !  write(*,*) 'a_n1+a_n2',a_n1+a_n2,j-1

         ! Entries in separator
         DO i = work(distance_ptr+a_n+best_sep_start),&
          work(distance_ptr+a_n+best_sep_start+ww)-1
           partition(j) = work(distance+i)
           j = j+1
         END DO
        !  write(*,*) 'a_n',a_n,j-1
         
        info = 0
        IF (printi .or. printd) THEN
          CALL mc70_print_message(info,unit_diagnostics, &
            'mc70_ashcraft')
        END IF
        RETURN

      END SUBROUTINE mc70_ashcraft


! ---------------------------------------------------
! mc70_level_set
! ---------------------------------------------------
! Partition the matrix using the level set method
      SUBROUTINE mc70_level_set(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,level,a_n1,a_n2,&
        a_weight_1,a_weight_2,a_weight_sep,partition,work,control,info,band,use_multilevel)

        INTEGER, INTENT(IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT(IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT(IN) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. 
        INTEGER, INTENT(IN) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.
        INTEGER, INTENT(IN) :: a_weight(a_n) ! On input a_weight(i) contains 
             ! the weight of column i.
        INTEGER, INTENT(IN) :: sumweight ! sum of entries in a_weight
        INTEGER, INTENT(IN) :: level ! current nested dissection level
        INTEGER, INTENT(OUT) :: a_n1, a_n2 ! size of the two submatrices
        INTEGER, INTENT(OUT) :: a_weight_1,a_weight_2,a_weight_sep ! Weighted 
             ! size of partitions and separator
        INTEGER, INTENT(OUT) :: partition(a_n) ! First a_n1 entries will contain
             ! list of (local) indices in partition 1; next a_n2 entries will 
             ! contain list of (local) entries in partition 2; entries in 
             ! separator are listed at the end
        INTEGER, INTENT(OUT) :: work(4*a_n)
        TYPE (mc70_control), INTENT(IN) :: control 
        INTEGER, INTENT(INOUT) :: info
        REAL (myreal_mc70), INTENT(INOUT) :: band ! If level = 1, then on 
          ! output band = max(band,100*L/a_n), where L is the size of the 
          ! largest levelset
        LOGICAL,INTENT(INOUT) :: use_multilevel ! are we allowed to use a multilevel
                ! partitioning strategy

! ---------------------------------------------
! Local variables
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        INTEGER :: nstrt,nend
        LOGICAL :: printi, printd
        INTEGER :: level_p,level_ptr_p,work_p
        INTEGER :: num_levels_nend ! no. levels in structure rooted at nend
        INTEGER :: num_entries ! no. entries in level structure rooted at nend
        INTEGER :: best_sep_start
        INTEGER :: i,j,p1sz,p2sz,sepsz,lwidth
        INTEGER :: ml_max_levels, lwork
        INTEGER :: mindeg, degree, max_search
        REAL(myreal_mc70) :: bestval
        REAL(myreal_mc70) :: val
        REAL(myreal_mc70) :: ratio
        LOGICAL :: imbal

! ---------------------------------------------
! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
 
        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Use one-sided level set method'
        END IF

        ratio = control%ratio
        IF (ratio.GT.real(sumweight-2)) THEN
            imbal = .FALSE.
        ELSE
            imbal = .TRUE.
        END IF



        IF (control%ml==1 .AND. use_multilevel) THEN
           use_multilevel = .TRUE.
            GOTO 222
        END IF


        ! Find pseudoperipheral nodes nstart and nend, and the level structure 
        ! rooted at nend
        level_ptr_p = 0 ! size a_n
        level_p = level_ptr_p + a_n ! size a_n
        work_p = level_p + a_n ! size 2*a_n 
        mindeg = a_n+1
        DO i = 1,a_n
          IF (i .LT. a_n) THEN
            degree = a_ptr(i+1)-a_ptr(i)
          ELSE
            degree = a_ne + 1 - a_ptr(i)
          END IF
          IF (degree .LT. mindeg) THEN
            mindeg = degree
            nstrt = i
          END IF
        END DO
        max_search = a_n

        CALL mc70_find_pseudo(a_n,a_ne,a_ptr,a_row,&
           work(level_ptr_p+1:level_ptr_p+a_n),work(level_p+1:level_p+a_n),&
           nstrt,nend,max_search,work(work_p+1:work_p+2*a_n),num_levels_nend,&
           num_entries,lwidth)

        IF (num_entries .lt. a_n) THEN
           ! matrix is separable
           a_n1 = num_entries
           a_n2 = a_n - a_n1
           work(work_p+1:work_p+a_n) = 0
           DO i = 1,num_entries
             j = work(level_p+i)
             partition(i) = j
             work(work_p+j) = 1
           END DO
           j = num_entries+1
           DO i = 1,a_n
             IF (work(work_p+i) .eq. 0) THEN
               partition(j) = i
               j = j+1
             END IF
           END DO
           RETURN
        END IF

        IF (level .EQ. 0) THEN
             band = max(band,100.0*real(lwidth,kind(1.0D0))/real(a_n,kind(1.0D0)))
        END IF

        IF ((control%ml<=0) .OR. (use_multilevel &
           .AND.control%ml_max_levels .LE. 0 ) ) THEN
            use_multilevel = .false.
        END IF
        IF (control%ml>=2 .AND. level .EQ. 0) THEN
         IF (real(lwidth,kind(1.0D0))/real(a_n,kind(1.0D0)) .LE. &
               control%ml_bandwidth   ) THEN
            use_multilevel = .false.
         ELSE
            use_multilevel = .true.
         END IF
        END IF

222     CONTINUE

        IF (use_multilevel) THEN
             ml_max_levels = control%ml_max_levels
          lwork = 9*a_n+sumweight+a_ne/2
          CALL multilevel_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,&
          partition,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
          control,info,lwork,work,ml_max_levels)
          RETURN
        END IF

         IF (num_levels_nend .LE. 2) THEN
           ! Not possible to find separator
           ! This can only be invoked for a full connected graph. The partition
           ! subroutine checks for this case so it should never be called
           a_n1 = a_n
           a_n2 = 0
           partition(1:a_n) = (/ (i,i=1,a_n) /)
           RETURN
         END IF

        ! Calculate total weight in each level set
         work(work_p+1:work_p+num_levels_nend) = 0
         DO i = 1,num_levels_nend-1
           DO j = work(level_ptr_p+i),work(level_ptr_p+i+1)-1
             work(work_p+i) = work(work_p+i) + a_weight(work(level_p+j))
           END DO
         END DO
         i = num_levels_nend
         DO j = work(level_ptr_p+i),a_n
             work(work_p+i) = work(work_p+i) + a_weight(work(level_p+j))
         END DO


         ! First possible separator contains all of the nodes in level 2
         p1sz = work(work_p+1)
         sepsz = work(work_p+2)
         p2sz = sum(work(work_p+3:work_p+num_levels_nend))
         a_weight_1 = p1sz
         a_weight_2 = p2sz
         a_weight_sep = sepsz
         best_sep_start = 2
         a_n1 = work(level_ptr_p+2) - 1
         a_n2 = a_n - work(level_ptr_p+3) +1
           CALL cost_function(p1sz,p2sz,&
              sepsz,sumweight,ratio,imbal,bestval)

         ! Search for best separator using tau
         DO j = 2,num_levels_nend-4
           p1sz = p1sz+work(work_p+j)
           sepsz = work(work_p+j+1) 
           p2sz = p2sz - work(work_p+j+1)
           CALL cost_function(p1sz,p2sz,&
              sepsz,sumweight,ratio,imbal,val)
           IF (val .LT. bestval) THEN
             bestval = val
             best_sep_start = j+1
             a_n1 = work(level_ptr_p+j+1) - 1
             a_n2 = a_n - work(level_ptr_p+j+2) +1  
             a_weight_1 = p1sz
             a_weight_2 = p2sz
             a_weight_sep = sepsz           
           END IF
         END DO


         ! Rearrange partition
         ! Entries in partition 1
         j = 1
         DO i = 1,work(level_ptr_p+best_sep_start)-1
           partition(j) = work(level_p+i)
           j = j+1
         END DO

         ! Entries in partition 2
         DO i = work(level_ptr_p+best_sep_start+1),a_n
           partition(j) = work(level_p+i)
           j = j+1
         END DO

         ! Entries in separator
         DO i = work(level_ptr_p+best_sep_start),&
          work(level_ptr_p+best_sep_start+1)-1
           partition(j) = work(level_p+i)
           j = j+1
         END DO

        info = 0
        IF (printi .or. printd) THEN
          CALL mc70_print_message(info,unit_diagnostics, &
            'mc70_level_set')
        END IF
        RETURN

      END SUBROUTINE mc70_level_set


! ---------------------------------------------------
! mc70_find_pseudo
! ---------------------------------------------------
! Find pseudoperipheral pairs of nodes for irreducible graph
      SUBROUTINE mc70_find_pseudo(a_n,a_ne,a_ptr,a_row,level_ptr,level,&
           nstrt,nend,max_search,work,num_levels,num_entries,lwidth)
        INTEGER, INTENT(IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT(IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT(IN) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. 
        INTEGER, INTENT(IN) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.
        INTEGER, INTENT(OUT) :: level_ptr(a_n) ! On output level_ptr(i) contains 
             ! position in level that entries for level i start. 
        INTEGER, INTENT(OUT) :: level(a_n) ! On output level contains lists of 
             ! rows according to the level set that they are in
        INTEGER, INTENT(INOUT) :: nstrt ! Starting pseudoperipheral node
        INTEGER, INTENT(OUT) :: nend ! End pseudoperipheral node
        INTEGER, INTENT(IN) :: max_search
        INTEGER, INTENT(OUT) :: work(2*a_n)
        INTEGER, INTENT(OUT) :: num_levels !
        INTEGER, INTENT(OUT) :: num_entries ! number of entries in level structure
        INTEGER, INTENT(OUT) :: lwidth
        ! Based on MC60HD

! ---------------------------------------------
! Local variables
        INTEGER :: i, j, k
        INTEGER :: mindeg, degree,maxdep,main,lwidth1
        INTEGER :: mask,list
        INTEGER :: nstop ! ending pseudoperipheral node
        INTEGER :: node ! Index of graph node
        INTEGER :: lsize ! size of final levelset
        INTEGER :: nlsize ! no. nodes of differing degrees in final level set
        INTEGER :: minwid ! Minimum levelset width
        INTEGER :: nlvl ! number of levels in level structure
        INTEGER :: minwid1
        j = 0
        mask = 0
        list = mask+a_n
! Initialise work(mask+1:mask+a_n) and work(list+1:list+a_n)
        work(mask+1:mask+a_n) = 1
        work(list+1:list+a_n) = 0
        level_ptr(:) = 0

! First guess for starting node is input nstrt

! Generate level structure for node with min supervariable degree
        CALL mc70_level_struct(nstrt,a_n,a_ne,a_ptr,a_row,a_n, &
          work(mask+1:mask+a_n),level_ptr,level,maxdep,lwidth,num_entries)
        IF (num_entries .lt. a_n) THEN
          ! matrix is separable
          RETURN

        END IF

        nstop = level(a_n)
        DO main = 1,a_n
        ! Store nodes in the final level set and their degrees
          lsize = 0
          DO i=level_ptr(maxdep),a_n
            node = level(i)
            lsize = lsize+1
            work(list+lsize) = node
            IF (node .EQ. a_n) THEN
              level_ptr(node) = a_ne+1 - a_ptr(node)
            ELSE
              level_ptr(node) = a_ptr(node+1)-a_ptr(node)
            END IF
          END DO

         ! Choose at most max_search nodes
          DO nlsize = 1,max_search
            ! Look for candiate with least degree
            mindeg = a_n+1
            !mindeg = -1
            DO i = nlsize,lsize
              IF (level_ptr(work(list+i)) .LT. mindeg) THEN
                j = i
                mindeg = level_ptr(work(list+i))
              END IF
            END DO
            ! Jump out of loop if no candidates left
            IF (mindeg .EQ. a_n+1) GO TO 55
          !  IF (mindeg .EQ. -1) GO TO 55
            ! Swap chose candidate to next position
            node = work(list+j)
            work(list+j) = work(list+nlsize)
            work(list+nlsize) = node
            ! Rule out the neighbours of the chosen node
            IF (node .EQ. a_n) THEN
              k = a_ne
            ELSE
              k = a_ptr(node+1)-1
            END IF
            DO i = a_ptr(node),k
              level_ptr(a_row(i)) = a_n+1
            END DO
          END DO
55        nlsize = nlsize-1

          ! Loop over nodes in list
          minwid = huge(a_n)
          minwid1 = huge(a_n)

          DO i = 1,nlsize
            node = work(list+i)

            ! Form rooted level structures for node
            CALL mc70_level_struct(node,a_n,a_ne,a_ptr,a_row,minwid, &
              work(mask+1:mask+a_n),level_ptr,level,nlvl,lwidth,num_entries)
          ! IF (lwidth .LE. minwid) THEN
            lwidth1 = 0
            DO k = 1,nlvl-1
              lwidth1 = lwidth1 + k*(level_ptr(i+1)-level_ptr(i))
            END DO 
            lwidth1 = lwidth1 + nlvl*(a_ne-level_ptr(nlvl))

              IF (nlvl .GT. maxdep) THEN
                ! Level structure of greater depth. Begin a new iteration.
                 nstrt = node
                 maxdep = nlvl
   
                 GOTO 70
              ELSE
                 IF (lwidth1 .LT.minwid1) THEN
                  
                  nstop = node
                  minwid = lwidth
                  minwid1 = lwidth1
                 END IF
              END IF
         !  END IF
          END DO
          GOTO 80
70        CONTINUE
        END DO
80      IF (nstop .NE. node) THEN
            CALL mc70_level_struct(node,a_n,a_ne,a_ptr,a_row,a_n, &
              work(mask+1:mask+a_n),level_ptr,level,nlvl,lwidth,num_entries)
        END IF
        num_levels = maxdep
        nend = nstop

        RETURN
      END SUBROUTINE mc70_find_pseudo


! ---------------------------------------------------
! mc70_level_struct
! ---------------------------------------------------
! Given a root, calculate the level structure of a given graph from that root

      SUBROUTINE mc70_level_struct(root,a_n,a_ne,a_ptr,a_row,maxwid,mask,&
          level_ptr,level,num_levels,lwidth,num_entries)

        INTEGER, INTENT(IN) :: root ! Root node for level structure
        INTEGER, INTENT(IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT(IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT(IN) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. 
        INTEGER, INTENT(IN) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.
        INTEGER, INTENT(IN) :: maxwid ! Max permissible width of level structure
        INTEGER, INTENT(INOUT) :: mask(a_n) ! Always restored to input value at
             ! end of the call. mask(node) > 0 for all visible nodes 
        INTEGER, INTENT(OUT) :: level_ptr(a_n) ! On output level_ptr(i) contains 
             ! position in level that entries for level i start. 
        INTEGER, INTENT(OUT) :: level(a_n) ! On output level contains lists of 
             ! rows according to the level set that they are in
        INTEGER, INTENT(OUT) :: num_levels ! On output num_levels contains the 
             ! number of levels
        INTEGER, INTENT(OUT) :: lwidth ! On output, contains the width of the 
             ! structure
        INTEGER, INTENT(OUT) :: num_entries ! On output, contains number of 
             ! entries in the tree structure containing root
        
! ---------------------------------------------
! Local variables
        INTEGER :: lvlend ! End of previous level in level
        INTEGER :: lnbr ! Next position in level
        INTEGER :: lbegin ! Beginning of present level in level
        INTEGER :: lw ! Level width
        INTEGER :: nbr ! neighbour node
        INTEGER :: node ! current node
        INTEGER :: nvars
        INTEGER :: i,j,k

! Based on mc60ld
! Establish level 1
        level_ptr(:) = 0
        mask(root) = - mask(root)
        level(1) = root
        lvlend = 0
        nvars = 0
        lnbr = 1
        lwidth = 1
        DO num_levels = 1,a_n
          ! Generate next level by finding all unmasked neighbours of nodes in
          ! present level
          lbegin = lvlend+1
          lvlend = lnbr
          level_ptr(num_levels) = lbegin
          lw = 0
          DO i = lbegin,lvlend
            node = level(i)
            IF (node .EQ. a_n) THEN
              k = a_ne
            ELSE
              k = a_ptr(node+1)-1
            END IF
            DO j = a_ptr(node),k
              nbr = a_row(j)
              IF (mask(nbr) .GT. 0) THEN
                lnbr = lnbr+1
                level(lnbr) = nbr
                mask(nbr) = - mask(nbr)
                lw = lw + 1
              END IF
            END DO
          END DO
          lwidth = MAX(lwidth,lw)
          nvars = nvars + lw
          ! If no neighbours found, we are done
          IF (lnbr .EQ. lvlend) EXIT
          ! Abort construnction if level structure too wide
         ! IF (lwidth .GT. maxwid) EXIT
        END DO
        ! Reset mask
        DO i = 1,lnbr
          mask(level(i)) = abs(mask(level(i)))
        END DO
        num_entries = lnbr        

      END SUBROUTINE mc70_level_struct


! ---------------------------------------------------
! mc70_distance
! ---------------------------------------------------
! Given two level structures, calculate the difference in each nodes distance 
! from the start and end node
      SUBROUTINE mc70_distance(a_n,num_levels_nend,level_ptr_nend,level_nend,&
          num_levels_nstrt,level_ptr_nstrt,level_nstrt,distance_ptr,distance)
        INTEGER, INTENT(IN) :: a_n ! order of matrix
        INTEGER, INTENT(IN) :: num_levels_nend ! number of levels with root nend
        INTEGER, INTENT(IN) :: level_ptr_nend(a_n) ! level_ptr(i) contains 
             ! position in level that entries for level i start (root = nend)
        INTEGER, INTENT(IN) :: level_nend(a_n) ! Contains lists of rows
             ! according to the level set that they are in (root = nend)
        INTEGER, INTENT(IN) :: num_levels_nstrt ! no. of levels with root nstrt
        INTEGER, INTENT(INOUT) :: level_ptr_nstrt(a_n) ! level_ptr(i) contains 
             ! position in level that entries for level i start (root = nstrt)
             ! Reused during subroutine
        INTEGER, INTENT(IN) :: level_nstrt(a_n) ! Contains lists of rows
             ! according to the level set that they are in (root = nstrt)
        INTEGER, INTENT(OUT) :: distance_ptr(2*a_n-1) ! distance(i) contains
             ! position in distance where entries with distance i-a_n  
        INTEGER, INTENT(OUT) :: distance(a_n) ! Contains lists of rows ordered 
             ! according to their distance
 
! ---------------------------------------------
! Local variables
        INTEGER :: j,k
        INTEGER :: dptr ! used as offset into distance_ptr
        INTEGER :: lev ! stores current level
     
        dptr = a_n
        distance_ptr(dptr+1-a_n:dptr+a_n-1) = 0
        distance(1:a_n) = 0

        ! set distance(i) to hold the level that i belongs to in the structure 
        ! rooted at nend (= distance from end node + 1)
        DO lev = 1,num_levels_nend-1
          DO j = level_ptr_nend(lev),level_ptr_nend(lev+1)-1
            distance(level_nend(j)) = lev
          END DO
        END DO
        DO j = level_ptr_nend(num_levels_nend),a_n
            distance(level_nend(j)) = num_levels_nend
        END DO
        
        ! now consider level structure rooted at start node
        DO lev = 1,num_levels_nstrt-1
          DO j = level_ptr_nstrt(lev),level_ptr_nstrt(lev+1)-1
            distance(level_nstrt(j)) = distance(level_nstrt(j))-lev
            k = distance(level_nstrt(j))
            distance_ptr(dptr+k) = distance_ptr(dptr+k) + 1
          END DO
        END DO
        DO j = level_ptr_nstrt(num_levels_nstrt),a_n
            distance(level_nstrt(j)) = distance(level_nstrt(j))-lev
            k = distance(level_nstrt(j))
            distance_ptr(dptr+k) = distance_ptr(dptr+k) + 1
        END DO

        ! Copy distance into level_ptr_nstrt to save memory
        level_ptr_nstrt(1:a_n) = distance(1:a_n)

        ! Set distance_ptr(i) to point one place to the right of where entries 
        ! with distance i will be
        DO j = 1-a_n,a_n-1
          IF (distance_ptr(dptr+j) .GT. 0) THEN
            EXIT
          ELSE
            distance_ptr(dptr+j) = 1
          END IF
        END DO
        distance_ptr(dptr+j) = distance_ptr(dptr+j) + 1
        DO k = j+1,a_n-1
          distance_ptr(dptr+k) = distance_ptr(dptr+k) + distance_ptr(dptr+k-1)
        END DO

        ! Set up lists of rows with same value of distance
        DO j = 1,a_n
          k = level_ptr_nstrt(j)
          distance_ptr(dptr+k) = distance_ptr(dptr+k) -1
          distance(distance_ptr(dptr+k)) = j
        END DO


      END SUBROUTINE mc70_distance

! ---------------------------------------------------
! mc70_refine
! ---------------------------------------------------
! Given a partition, refine the partition to improve the (weighted) value of
! a_n1*a_n2/(a_n-a_n1-a_n2). Dulmage-Mendelsohn is used 
      SUBROUTINE mc70_refine(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,partition,work,control,ref_method)
        INTEGER, INTENT(IN) :: a_n ! order of matrix
        INTEGER, INTENT(IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT(IN) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. 
        INTEGER, INTENT(IN) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.
        INTEGER, INTENT(IN) :: a_weight(a_n) ! On input a_weight(i) contains 
             ! the weight of column i 
        INTEGER, INTENT(IN) :: sumweight ! Sum of weights in a_weight
        INTEGER, INTENT(INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT(INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT(INOUT) :: a_weight_1,a_weight_2,a_weight_sep ! Weighted 
             ! size of partitions and separator
        INTEGER, INTENT(INOUT) :: partition(a_n) !First a_n1 entries contain
             ! list of (local) indices in partition 1; next a_n2 entries  
             ! contain list of (local) entries in partition 2; entries in 
             ! separator are listed at the end. This is updated to the new 
             ! partition
        INTEGER, INTENT(OUT) :: work(8*a_n+sumweight+a_ne/2) ! Work array
        TYPE(mc70_control), INTENT(IN) :: control
        INTEGER, INTENT(IN) :: ref_method ! <= 1 just do FM refinement,
                                            ! >= 2 do DM and them FM

! ---------------------------------------------
! Local variables
        INTEGER :: dm_m ! number of rows in part of matrix DM applied to
        INTEGER :: dm_n ! number of columns in part of matrix DM applied to
        INTEGER :: dm_ne ! number of entries in part of matrix DM applied to
        INTEGER :: dm_ptr ! pointer into work array for start of col pointers 
                          !for submatrix
        INTEGER :: dm_row ! pointer into work array for start of row indices 
                          !for submatrix
        INTEGER :: dm_rowperm ! pointer into work array for start of rowperm 
                          ! from DM
        INTEGER :: dm_colperm ! pointer into work array for start of colperm 
                          ! from DM
        INTEGER :: dm_rowptr ! pointer into work array for start of rowptr 
                          ! from DM
        INTEGER :: dm_colptr ! pointer into work array for start of colptr 
                          ! from DM
        INTEGER :: fm_flags ! pointer into work array for start of flags from FM
        INTEGER :: fm_ipart ! pointer into work array for start of ipart from FM
        INTEGER :: fm_next ! pointer into work array for start of next from FM
        INTEGER :: fm_last ! pointer into work array for start of last from FM
        INTEGER :: fm_gain1 ! pointer into work array for start of gain1 from FM
        INTEGER :: fm_gain2 ! pointer into work array for start of gain2 from FM
        INTEGER :: fm_done ! pointer into work array for start of done from FM
        INTEGER :: fm_head ! pointer into work array for start of head from FM
        INTEGER :: fm_distance ! pointer into work array for start of head from FM
        INTEGER :: imap ! pointer into work array for start of imap
        INTEGER :: part_no ! partition number used with DM
        INTEGER :: ibest,a_weight_1_best,a_weight_2_best,a_weight_sep_best,band 
                  ! Used in finding best partition from DM
        INTEGER :: wrow, wcol ! Used in finding best partition from DM
        INTEGER :: icut,mult ! Used within FM refinement
        INTEGER :: i,j,k,l,p,q,r
        REAL(myreal_mc70) :: tau ! Current value of tau
        REAL(myreal_mc70) :: tau_best ! Best value of tau so far
        REAL(myreal_mc70) :: ratio
        LOGICAL :: imbal ! Should we check for imbalance?

        TYPE(mc79_control) :: dm_control ! control for DM call
        TYPE(mc79_info) :: dm_info ! info from DM call
        ratio = control%ratio
        IF (ratio .GT. real(sumweight-2)) THEN
          imbal = .FALSE.
        ELSE
          imbal = .TRUE.
        END IF
        IF ( ref_method .LE. 1) THEN
          GOTO 111
        END IF
        ! Decide whether to apply DM to the sub-matrix formed from the 
        ! intersection of the separator with partition 1 or the intersection of 
        ! the separator with partition 2
        IF (a_weight_1 .GT. a_weight_2) THEN
          ! Apply to partition 1
          part_no = 1
          dm_m = a_n1
        ELSE
          ! Apply to partition 2
          part_no = 2
          dm_m = a_n2
        END IF
        dm_n = a_n - a_n1 - a_n2

        ! Setup pointers into work array
        dm_ptr = 0 ! length dm_n+1
        dm_row = dm_ptr + dm_n + 1 ! length bounded by a_ne/2 -1 (irreduc. mtrx)
        dm_rowperm = dm_row + a_ne/2 -1 ! length dm_m
        dm_colperm = dm_rowperm + dm_m ! length dm_n
        dm_rowptr = dm_colperm + dm_n ! length dm_m+2
        dm_colptr = dm_rowptr + dm_m + 2 ! length dm_n+2
        imap = dm_colptr + dm_n + 2 ! length a_n
        ! dm_n + dm_m + 1 .LE. a_n =>
        ! 3*dm_n +4 + 2*dm_m + a_ne/2 + a_n 
        !    .LE. 3*dm_n +3 + 3*dm_m + a_ne/2 + a_n 
        !    .LE. 4*a_n + a_ne/2

        ! Setup imap such that imap(i) if i is a row in A, then it becomes row 
        ! imap(i) in the matrix we apply DM to
        work(imap+1:imap+a_n) = 0
        IF (part_no .EQ. 1) THEN
          j = 0
        ELSE
          j = a_n1
        END IF
        k = 1
        DO i = 1,dm_m
          work(imap+partition(j+k)) = k
          k = k+1
        END DO
        ! Setup imap such that imap(i) if i is a column in A, then it becomes column 
        ! -imap(i) in the matrix we apply DM to
        j = a_n1 + a_n2
        k = 1
        DO i = 1,dm_n
          work(imap+partition(j+k)) = -k
          k = k+1
        END DO

        ! Set dm_ptr(i) to hold the number of entries in column i of submatrix
        work(dm_ptr+1:dm_ptr+dm_n+1) = 0
        j = a_n1 + a_n2
        DO i = 1,dm_n
          l = partition(i+j)
          IF (l .EQ. a_n) THEN
            k = a_ne
          ELSE 
            k = a_ptr(l+1)-1
          END IF
          DO p = a_ptr(l),k
            IF (work(imap+a_row(p)) .GT. 0 ) work(dm_ptr+i) = work(dm_ptr+i) + 1
          END DO 
        END DO
  
        ! Set dm_ptr(i) to point one place to the right of where entries for 
        !column i are stored
        work(dm_ptr+1) =  work(dm_ptr+1) + 1
        DO i = 2,dm_n
          work(dm_ptr+i) =  work(dm_ptr+i) + work(dm_ptr+i-1)
        END DO
        work(dm_ptr+dm_n+1) =  work(dm_ptr+dm_n)
        dm_ne = work(dm_ptr+dm_n+1)-1

        ! Store row entries
        j = a_n1 + a_n2
        DO i = 1,dm_n
          l = partition(i+j)
          IF (l .EQ. a_n) THEN
            k = a_ne
          ELSE 
            k = a_ptr(l+1)-1
          END IF
          DO p = a_ptr(l),k
            IF (work(imap+a_row(p)) .GT. 0 ) THEN
              q = work(imap+a_row(p))
              work(dm_ptr+i) =  work(dm_ptr+i) - 1
              r = work(dm_ptr+i)
              work(dm_row+r) = q
            END IF
          END DO 
        END DO
        
        ! Compute fine Dulmage-Mendelsohn decomposition
        CALL mc79_fine(dm_m,dm_n,work(dm_ptr+1:dm_ptr+dm_n+1),&
          work(dm_row+1:dm_row+dm_ne),work(dm_rowperm+1:dm_rowperm+dm_m),&
          work(dm_colperm+1:dm_colperm+dm_n),&
          work(dm_rowptr+1:dm_rowptr+dm_m+2),&
          work(dm_colptr+1:dm_colptr+dm_n+2),dm_control,dm_info)

        ! Map local indices to those of partitioned matrix
        ! First do colperm
        j = a_n1 + a_n2
        DO i = 1,dm_n
          k = work(dm_colperm+i)
          work(dm_colperm+i) = partition(k+j)
        END DO
        ! Now do rowperm
        IF (part_no .EQ. 1) THEN
          j = 0
        ELSE
          j = a_n1
        END IF
        DO i = 1,dm_m
          k = work(dm_rowperm+i)
          work(dm_rowperm+i) = partition(k+j)
        END DO 

        ! Find best partition
        ibest = -1 ! Original partition
        CALL cost_function(a_weight_1,a_weight_2,&
              a_weight_sep,sumweight,ratio,imbal,tau_best)
        a_weight_1_best = a_weight_1
        a_weight_2_best = a_weight_2
        a_weight_sep_best = a_weight_sep
        wrow = 0
        wcol = 0
        DO i = 1, work(dm_rowptr+dm_info%hz_comps+1)-1
          k = work(dm_rowperm+i) ! Index of row wrt matrix being partitioned
          wrow = wrow + a_weight(k)
        END DO
        DO i = work(dm_colptr+dm_info%hz_comps+1),dm_n
          k = work(dm_colperm+i) ! Index of col wrt matrix being partitioned
          wcol = wcol + a_weight(k)
        END DO

        IF (part_no .EQ. 1) THEN
           CALL cost_function(a_weight_1-wrow,a_weight_2+a_weight_sep-wcol,&
              wrow+wcol,sumweight,ratio,imbal,tau)
        ELSE
           CALL cost_function(a_weight_1+a_weight_sep-wcol,a_weight_2-wrow,&
              wrow+wcol,sumweight,ratio,imbal,tau)
        END IF
        IF (tau .LT. tau_best) THEN
          ibest = dm_info%hz_comps
          tau_best = tau
          a_weight_sep_best = wrow+wcol

          IF (part_no .EQ. 1) THEN
           a_weight_1_best = a_weight_1-wrow
           a_weight_2_best = a_weight_2+a_weight_sep-wcol
          ELSE
           a_weight_2_best = a_weight_2-wrow
           a_weight_1_best = a_weight_1+a_weight_sep-wcol
          END IF

        END IF

        ! Check all possible separators
        DO i = dm_info%hz_comps+1,dm_info%hz_comps+dm_info%sq_comps
          DO j = work(dm_rowptr+i),work(dm_rowptr+i+1)-1
            k = work(dm_rowperm+j) ! Index of row wrt matrix being partitioned
            wrow = wrow + a_weight(k)
          END DO
          DO j = work(dm_colptr+i),work(dm_colptr+i+1)-1
            k = work(dm_colperm+j) ! Index of row wrt matrix being partitioned
            wcol = wcol - a_weight(k)
          END DO

         IF (part_no .EQ. 1) THEN
           CALL cost_function(a_weight_1-wrow,a_weight_2+a_weight_sep-wcol,&
              wrow+wcol,sumweight,ratio,imbal,tau)
         ELSE
          
           CALL cost_function(a_weight_1+a_weight_sep-wcol,a_weight_2-wrow,&
              wrow+wcol,sumweight,ratio,imbal,tau)

         END IF
         IF (tau .LT. tau_best) THEN
           ibest = i
           tau_best = tau
           a_weight_sep_best = wrow+wcol

           IF (part_no .EQ. 1) THEN
            a_weight_1_best = a_weight_1-wrow
            a_weight_2_best = a_weight_2+a_weight_sep-wcol
           ELSE
            a_weight_2_best = a_weight_2-wrow
            a_weight_1_best = a_weight_1+a_weight_sep-wcol
           END IF

          END IF
        END DO

        IF (ibest .EQ. -1) THEN
          ! No separator has been found that lowers tau
          GOTO 111
        END IF

        IF (part_no .EQ. 1) THEN
           ! Entries in new partition separator lie in union of 
           ! rowperm(1:rowptr(ibest+1)-1) and
           ! colperm(colptr(ibest+1),dm_n)
           ! Entries in new partition 1 lie in rowperm(rowptr(ibest+1):dm_m)
           ! Entries in new partition 2 lie in union of
           ! partition(a_n1+1:a_n1+a_n2) and
           ! colperm(1:colptr(ibest+1)-1)

           ! Move entries already in partition 2 to correct place in partition
           j = dm_m-(work(dm_rowptr+ibest+1)-1)+1
           k = a_n1+a_n2
           DO WHILE ((j .LE. a_n1) .AND. &
                (j .LE. a_n2+dm_m-(work(dm_rowptr+ibest+1)-1)+1) )
              partition(j) = partition(k)
              k = k-1
              j = j+1
           END DO
           
           ! Copy entries in partition 1 into correct place
           DO j = 1,dm_m-work(dm_rowptr+ibest+1)+1
             k = work(dm_rowptr+ibest+1)+j-1
             partition(j) = work(dm_rowperm+k)
           END DO
           a_n1 = dm_m-work(dm_rowptr+ibest+1)+1

           ! Copy entries moving from sep into partition 2 to correct place
           DO j = a_n1+a_n2+1,a_n1+a_n2+work(dm_colptr+ibest+1)-1
             partition(j) = work(dm_colperm+j-a_n1-a_n2)
           END DO
           a_n2 = a_n2 + work(dm_colptr+ibest+1)-1

           ! Copy entries moving from partition 1 to separator
           DO j = a_n1 + a_n2 +1, a_n1+a_n2+work(dm_rowptr+ibest+1)-1
             k = j - a_n1-a_n2
             partition(j) = work(dm_rowperm+k)
           END DO

           ! Copy entries remaining in separator
           DO j = a_n1+a_n2+work(dm_rowptr+ibest+1),a_n
             k = j-a_n1-a_n2-work(dm_rowptr+ibest+1)+1
             l = work(dm_colptr+ibest+1)+k-1
             partition(j) = work(dm_colperm+l)
           END DO
           a_weight_1 = a_weight_1_best
           a_weight_2 = a_weight_2_best
           a_weight_sep = a_weight_sep_best
        ELSE

           ! Entries in new partition separator lie in union of 
           ! rowperm(1:rowptr(ibest+1)-1) and
           ! colperm(colptr(ibest+1),dm_n)
           ! Entries in new partition 2 lie in rowperm(rowptr(ibest+1):dm_m)
           ! Entries in new partition 1 lie in union of
           ! partition(1:a_n1) and
           ! colperm(1:colptr(ibest+1)-1)

           ! Entries already in partition 1 do not need moving

           ! Copy entries moving from sep to partition 1 to correct place
           DO j = a_n1+1,a_n1+work(dm_colptr+ibest+1)-1
              k = j-a_n1
              partition(j) = work(dm_colperm+k)
           END DO
           a_n1 = a_n1 + work(dm_colptr+ibest+1)-1

           ! Copy entries into partition 2
           DO j = a_n1 + 1, a_n1+dm_m+1-work(dm_rowptr+ibest+1)
              k = j-a_n1-1+work(dm_rowptr+ibest+1)
              partition(j) = work(dm_rowperm+k)
           END DO
           a_n2 = dm_m+1-work(dm_rowptr+ibest+1)

           ! Copy entries from partition 2 into separator
           DO j = a_n1 + a_n2 +1, a_n1+a_n2+work(dm_rowptr+ibest+1)-1
             k = j - a_n1-a_n2
             partition(j) = work(dm_rowperm+k)
           END DO

           ! Copy entries remaining in separator
           DO j = a_n1+a_n2+work(dm_rowptr+ibest+1),a_n
             k = j-a_n1-a_n2-work(dm_rowptr+ibest+1)+1
             l = work(dm_colptr+ibest+1)+k-1
             partition(j) = work(dm_colperm+l)
           END DO
           a_weight_1 = a_weight_1_best
           a_weight_2 = a_weight_2_best
           a_weight_sep = a_weight_sep_best
        END IF

111     CONTINUE      
        
        fm_flags = 0 ! length a_n
        fm_ipart = fm_flags + a_n ! length a_n
        fm_next  = fm_ipart + a_n ! length a_n
        fm_last  = fm_next + a_n ! length a_n
        fm_gain1 = fm_last + a_n ! length a_n
        fm_gain2 = fm_gain1 + a_n ! length a_n
        fm_done  = fm_gain2 + a_n ! length a_n
        fm_head  = fm_done + a_n ! length icut+mult+1
        icut = min(sumweight-1,min(sumweight/5,100*sumweight/a_n))
       ! icut = sumweight/2
       ! mult = min(sumweight/20,10*sumweight/a_n) - 1
        mult = sumweight -icut-1
       ! mult = sumweight/2-1
        fm_distance = fm_head+icut+mult+1
        
        ! Initialise work(fm_flags+1:fm_flags+a_n)
        DO j = 1,a_n1
          k = partition(j)
          work(fm_flags+k) = mc70_part1_flag
        END DO
        DO j = a_n1+1,a_n1+a_n2
          k = partition(j)
          work(fm_flags+k) = mc70_part2_flag
        END DO
        DO j = a_n1+a_n2+1,a_n
          k = partition(j)
          work(fm_flags+k) = mc70_sep_flag
        END DO

        IF (control%refinement_band.gt. 0) THEN
         band = min(control%refinement_band,a_n)
        ELSE
         band = a_n
        END IF
        IF (.true.) THEN
         CALL mc70_fm_refinement_band_balance(a_n,a_ne,a_ptr,a_row,a_weight,&
          sumweight,&
          icut,mult,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,&
          band,ratio,work(fm_flags+1:fm_flags+a_n),&
          work(fm_ipart+1:fm_ipart+a_n),work(fm_next+1:fm_next+a_n),&
          work(fm_last+1:fm_last+a_n),work(fm_gain1+1:fm_gain1+a_n),&
          work(fm_gain2+1:fm_gain2+a_n),work(fm_done+1:fm_done+a_n),&
          work(fm_head+1:fm_head+icut+mult+1),&
          work(fm_distance+1:fm_distance+a_n))
        ELSE
         CALL mc70_fm_refinement_band(a_n,a_ne,a_ptr,a_row,a_weight,&
          sumweight,&
          icut,mult,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,&
          band,work(fm_flags+1:fm_flags+a_n),&
          work(fm_ipart+1:fm_ipart+a_n),work(fm_next+1:fm_next+a_n),&
          work(fm_last+1:fm_last+a_n),work(fm_gain1+1:fm_gain1+a_n),&
          work(fm_gain2+1:fm_gain2+a_n),work(fm_done+1:fm_done+a_n),&
          work(fm_head+1:fm_head+icut+mult+1),&
          work(fm_distance+1:fm_distance+a_n))
        END IF

        ! Update partition
        i = 1
        j = a_n1+1
        k = a_n1+a_n2+1
        DO l = 1,a_n
          SELECT CASE (work(fm_flags+l))
           CASE (mc70_part1_flag)
            partition(i) = l
            i = i+1
           CASE (mc70_part2_flag)
            partition(j) = l
            j = j+1
           CASE (mc70_sep_flag)
            partition(k) = l
            k = k+1
          END SELECT
        END DO

      END SUBROUTINE mc70_refine



      ! The subroutine mc70_fm_refinement_band uses a version of the Fiduccia-
! Mattheyses refinement algorithm on a tripartite partitioing of the nodes of a 
! graph where a node in the first partition is not connected to any node  
! in the second partition and any path between nodes in partition 1 and 
! partition 2 must go through a node in the cutset (partition 3). 
! The intention of the algorithm is to reduce the value of |P3|/(|P2||P3|). 
! This is a banded version so only nodes a distance of at most band from the 
! input separator can be moved into the new separator

      SUBROUTINE mc70_fm_refinement_band(n,a_ne,ptr,col,weight,&
          sumweight,icut,mult,a_n1,a_n2,wnv1,wnv2,wns,band,flags,&
          ipart,next,last,gain1,gain2,done,head,dist)

! Matrix is held in matrix using compressed column scheme 
        INTEGER, INTENT(IN) :: n  ! size of matrix
        INTEGER, INTENT(IN) :: a_ne  ! no. nonzeros in matrix
        INTEGER, INTENT(IN) :: ptr(n) ! row pointers
        INTEGER, INTENT(IN) :: col(a_ne) ! column indices
       ! TYPE (zd11_type), INTENT (INOUT) :: matrix
! The array weight is used to hold a weight on the vertices indicating 
! how many vertices from the finer graphs have been combined into the 
! current coarse graph vertex. 
        INTEGER, INTENT (IN) :: weight(n)
        INTEGER, INTENT(IN) :: sumweight
        INTEGER, INTENT(IN) :: icut ! Used to limit search
        INTEGER, INTENT(IN) :: mult ! Used to bound search
        INTEGER, INTENT (INOUT) :: a_n1 ! No. vertices partition 1
        INTEGER, INTENT (INOUT) :: a_n2 ! No. vertices partition 2
        INTEGER, INTENT (INOUT) :: wnv1 ! Weighted sum of vertices partition 1
        INTEGER, INTENT (INOUT) :: wnv2 ! Weighted sum of vertices partition 2
        INTEGER, INTENT (INOUT) :: wns ! Weighted sum of vertices separator
        INTEGER, INTENT (IN) :: band ! width of band around initial separator
! that the separator can lie in
!
! flags holds a list of nodes stating which partition node i is in.
! The whole point of this routine is to return a revised partition 
! with better properties.  Normally less nodes in the cutset while maintaining 
! a balance between the number of nodes in the two components. 
! flags(i) == mc70_part1_flag : i is in partition 1
! flags(i) == mc70_part2_flag : i is in partition 2
! flags(i) == mc70_sep_flag   : i is in separator/cutset
        INTEGER, INTENT (INOUT) :: flags(n)
! info holds parameters giving information about the performance of the  
! subroutine 
        INTEGER, INTENT(OUT) :: ipart(n), next(n), last(n)
        INTEGER, INTENT(OUT) :: gain1(n),gain2(n),done(n),head(-mult:icut)
        INTEGER, INTENT(OUT) :: dist(n)

! Number nodes in each partition 
        INTEGER :: nv1, ns, nv2, inv1, inv2, ins
! Weighted nodes in each partition 
        INTEGER :: winv1, winv2, wins
        INTEGER :: i, j, ii, jj, eye, k,l
        INTEGER :: inn, outer
        INTEGER :: move, ming, gain, old_gain, inext, ilast, idummy
        INTEGER :: first,tail
        REAL (myreal_mc70) :: eval, evalc, evalo,eval1, eval2

! Set-up distance to hold (min) distance of node from separator. Only find 
! distance for nodes within band distance from separator
        first = 0
        tail = 0
        DO i = 1,n
          IF (flags(i) .EQ. mc70_sep_flag) THEN
           dist(i) = 0
           IF (first .EQ. 0) THEN
             first = i
             tail = i
           ELSE
             next(tail) = i
             tail = i
           END IF
          ELSE
           dist(i) = -2
          END IF
        END DO
        
         DO WHILE (first .NE. 0)
          j = dist(first)
          IF (j .EQ. band-1) THEN
            IF (first .EQ. n) THEN
              l = a_ne
            ELSE
              l = ptr(first+1)-1
            END IF
            DO i = ptr(first),l
              IF (dist(col(i)) .EQ. -2) THEN
                k = col(i)
                dist(k) = j+1
              END IF
            END DO            

          ELSE
            IF (first .EQ. n) THEN
              l = a_ne
            ELSE
              l = ptr(first+1)-1
            END IF
            DO i = ptr(first),l
              IF (dist(col(i)) .EQ. -2) THEN
                k = col(i)
                dist(k) = j+1
                next(tail) = k
                tail = k
              END IF
            END DO

          END IF
          IF (first .EQ. tail) THEN
            first = 0
          ELSE
            k = next(first)
            first = k
          END IF
         END DO
        next(1:n) = 0

! nv1,nv2,ns are the number of nodes in partitions 1, 2 and the cutset in 
! the current partition 
! inv1,inv2,ins,ipart are the equivalent quantities within the inner loop 
! The same identifiers prefixed by w refer to weighted counts 
! inner and outer are the two main loop indices 

! Initialize nv1,nv2,ns 
        nv1 = a_n1
        nv2 = a_n2
        ns = n - (a_n1+a_n2)
        ii = 1
        jj = ii + nv1

!Initialize ipart 
        ipart(1:n) = flags(1:n)

! Initialize array done that flags that a node has been considered in 
! an inner loop pass 
        done = 0

! Compute evaluation function for current partitioning 
           evalc = (real(wns)/real(wnv1+1))/real(wnv2+1)


! icut is set to limit search in inner loop .. may later be a parameter 
! we allow gains of up to max(weight)*5 

       ! icut = (sumweight/n)*5
        !  ALLOCATE (head(-mult:icut),STAT=stat)
        !  IF (stat/=0) GO TO 100
          head(-mult:icut) = 0
        !  head(0-mult:icut) = 0

! Set up doubly linked list linking nodes with same gain and headers 
! to starts (up to cut off value icut) 

! Compute gains for nodes in cutset 
        ming = sumweight
        DO i = 1, n

         IF (flags(i).eq.mc70_sep_flag) THEN
! Node i is in cutset 
! gain1(i) is change to cutset size if node i is moved to partition 1. 
! gain2(i) is change to cutset size if node i is moved to partition 2. 
! Run through all neighbours of node i to see what loss/gain is if node 
! i is moved 

          CALL compute_gain(i,flags)
          gain = max(-mult,min(gain1(i),gain2(i)))
          
          IF (gain<ming) ming = gain
          IF (gain>icut) CYCLE
! New node is put at head of list 
          CALL add_to_list(i,gain)
         END IF
        END DO
       
! Initilialization finished.  Now perform F-M algorithm in two loops. 
! In each inner loop we choose the best obtained and if the evaluation 
! function for this is better than previous best we perform another inner 
! loop; otherwise we terminate. 
        evalo = evalc
        DO outer = 1, n
! Set partition that will be altered in inner loop 
          inv1 = nv1
          inv2 = nv2
          ins = ns
          winv1 = wnv1
          winv2 = wnv2
          wins = wns
          ipart(1:n) = flags(1:n)
INNER:    DO inn = 1, n

! Choose best eligible move 
            DO idummy = 1, n

               DO gain = ming, icut
                IF (head(gain)/=0)  EXIT
               END DO
              IF (gain>icut) EXIT INNER

! Now cycle through nodes of least gain 
! Currently inefficient because of re-searching linked list 
              inext = head(gain)
              k = 0
111           i = inext
              IF (i==0) CYCLE
! Use node if it has not been considered already 
              IF (done(i)<outer) GO TO 222
              inext = next(i)
!!! Extra statements to trap infinite loop 
              k = k + 1
              IF (k>ins) THEN
                !WRITE (*,*) 'Bug in code because of infinite loop'
!!! You may wish to change this to a stop 
                EXIT
              END IF
              GO TO 111
            END DO
            EXIT INNER
! Node i has been selected as the best eligible node 
! Set flag so only considered once in this pass 
222         done(i) = outer
! As i will not be chosen again in this pass, remove from list 
            CALL remove_from_list(i,gain)
! Move the node to the appropriate partition and reset partition information 
! We will try both weighted and unweighted 

            IF (wnv1 .eq. 0 .and. wnv2 .gt. 0) THEN
            
! Move node i to partition 1 
              move = mc70_part1_flag
              inv1 = inv1 + 1
              winv1 = winv1 + weight(i)
              ins = ins - 1
              wins = wins - weight(i)
            ELSE IF (wnv2 .eq. 0 .and. wnv1 .gt. 0) THEN
! Move node i to partition 2 
              move = mc70_part2_flag
              inv2 = inv2 + 1
              winv2 = winv2 + weight(i)
              ins = ins - 1
              wins = wins - weight(i)
            ELSE
              eval1 = ((real(wins+gain1(i)-1))/real(winv1+weight(i)))&
                      /real(max(1,winv2+1-gain1(i)-weight(i)))
              eval2 = ((real(wins+gain2(i)-1))/real(winv2+weight(i)))&
                      /real(max(1,winv1+1-gain2(i)-weight(i)))
            IF ((eval1<eval2) .OR. ( (eval1==eval2) .AND. (wnv1<  wnv2) )) THEN
! Move node i to partition 1 
              move = mc70_part1_flag
              inv1 = inv1 + 1
              winv1 = winv1 + weight(i)
              ins = ins - 1
              wins = wins - weight(i)
            ELSE
! Move node i to partition 2 
              move = mc70_part2_flag
              inv2 = inv2 + 1
              winv2 = winv2 + weight(i)
              ins = ins - 1
              wins = wins - weight(i)
            END IF
            END IF
! Set new partition for node i 
            ipart(i) = move
! Run through neigbours of node i to update data 
            IF (i.EQ.n) THEN
              l = a_ne
            ELSE
              l = ptr(i+1) - 1
            END IF
            DO jj = ptr(i), l
              j = col(jj)
! Check which partition node j is in and take appropriate action 
              IF (ipart(j)==move) CYCLE
! If node j is in cutset, update its gain value 
              IF (ipart(j)==mc70_sep_flag) THEN
! If it has already been chosen in this pass just skip it 
                IF (done(j)==outer .or. dist(j) .eq. -2) CYCLE
! old_gain is present gain 
             
                old_gain = max(-mult,min(gain1(j),gain2(j)))
               ! old_gain = min(gain1(j),gain2(j))

                IF (move==mc70_part1_flag) gain2(j) = gain2(j) + weight(i)
                IF (move==mc70_part2_flag) gain1(j) = gain1(j) + weight(i)
                gain = max(-mult,min(gain1(j),gain2(j)))
               ! gain = min(gain1(j),gain2(j))

                IF (old_gain==gain) CYCLE
! Remove from old list 
                IF (old_gain<=icut) THEN
                  CALL remove_from_list(j,old_gain)
                END IF
! gain has changed so move to new linked list if less than icut 
                IF (gain<=icut) THEN
! Reset ming if necessary 
                  IF (gain<ming) ming = gain
                  CALL add_to_list(j,gain)
                END IF
              END IF
              IF (ipart(j)==2-move) THEN
! We have a new node in the cutset. 
                ipart(j) = mc70_sep_flag
! Compute gains for this new node in the cutset and place in linked list 
! We intentionally did not do this earlier but we do now 
! [maybe not since won't access this node again in this pass] 
! We use done array to record this but probably not necessary as not put 
! in head linked list so won't be accessed 
! First check that it was not earlier moved from cutset 
                IF (done(j)/=outer .and. dist(j) .ne. -2) THEN
! Compute gain 
                  CALL compute_gain(j,ipart)
                  gain = max(-mult,min(gain1(j),gain2(j)))
                !  gain = min(gain1(j),gain2(j))
!!! Just added this 
                  IF (gain<ming) ming = gain
! Add to  list 
                  IF (gain<=icut) THEN
                      CALL add_to_list(j,gain)
                  END IF
                END IF
! Update partition and gain of any nodes in cutset connected to node j 
                ins = ins + 1
                wins = wins + weight(j)
                IF (move==mc70_part1_flag) THEN
                  inv2 = inv2 - 1
                  winv2 = winv2 - weight(j)
                END IF
                IF (move==mc70_part2_flag) THEN
                  inv1 = inv1 - 1
                  winv1 = winv1 - weight(j)
                END IF
! Check neighbours of j since any in cut set will have gain changed 
                IF (j.EQ.n) THEN
                  l = a_ne
                ELSE
                  l = ptr(j+1) - 1
                END IF
                DO ii = ptr(j), l
                  eye = col(ii)
                  IF (ipart(eye)/=mc70_sep_flag) CYCLE
                  IF (dist(eye) .eq. -2) CYCLE
! Neighbour is in cutset. Recompute gain and insert in linked list. 
                  IF (done(eye)==outer) CYCLE
! old_gain is present gain 
                  old_gain = max(-mult,min(gain1(eye),gain2(eye)))
                 ! old_gain = min(gain1(eye),gain2(eye))
                  
                  
                  IF (move==mc70_part1_flag) THEN
                     gain1(eye) = gain1(eye) - weight(j)
                  END IF
                  IF (move==mc70_part2_flag) THEN
                     gain2(eye) = gain2(eye) - weight(j)
                  END IF
! gain is new gain 
                  gain = max(-mult,min(gain1(eye),gain2(eye)))
                 ! gain = min(gain1(eye),gain2(eye))
                  IF (old_gain==gain) CYCLE
! Remove from old list 
                  IF (old_gain<=icut) THEN
                      CALL remove_from_list(eye,old_gain)
                  END IF
! gain has changed so move to new linked list if less than icut 
                  IF (gain<=icut) THEN
! Reset ming if necessary 
                    IF (gain<ming) ming = gain
                    CALL add_to_list(eye,gain)
                  END IF
                END DO
              END IF
! end of neighbours loop 
            END DO

!         ii = 0 
!         do i = 1,n 
!           if (ipart(i) == 2) ii = ii + 1 
!         enddo 
!         if (ii .ne. inv2) write(6,*) 'problem in partition',ii,inv2 

! Evaluate new partition 
           eval = (real(wins)/real(winv1+1))/real(winv2+1)
! Compare this with best so far in inner loop and store partition 
! information if it is the best 
            IF (inv1*inv2>0 .AND. nv1*nv2.eq.0) THEN
! Might have to store gains and who is in the cutset 
              evalc = eval
              nv1 = inv1
              nv2 = inv2
              ns = ins
              wnv1 = winv1
              wnv2 = winv2
              wns = wins
              flags = ipart

            ELSE IF (eval<evalc .AND. (inv1*inv2>0 )) THEN
! Might have to store gains and who is in the cutset 
              evalc = eval
              nv1 = inv1
              nv2 = inv2
              ns = ins
              wnv1 = winv1
              wnv2 = winv2
              wns = wins
              flags(1:n) = ipart(1:n)
            END IF
! End inner loop 
          END DO INNER
! Leave loop if inner loop has not found better partition 
          IF (evalc>=(1.0-1.0/(log(real(sumweight))**2.3))*evalo) EXIT
! Otherwise we reset evalo and go back to inner loop 
         ! write(*,*) 'e',evalo,evalc,evalc/evalo
          evalo = evalc
! Recompute gains for this new partition 
! Compute gains for nodes in cutset 
! This is very inefficient but is in now to test functionality 
          head(-mult:icut) = 0
          ming = icut + 1
          DO i = 1, n
            IF (flags(i)/=mc70_sep_flag) CYCLE
            IF (dist(i) .eq. -2) CYCLE
! Node i is in cutset 
! gain1(i) is change to cutset size if node i is moved to partition 1. 
! gain2(i) is change to cutset size if node i is moved to partition 2. 
! Run through all neighbours of node i to see what loss/gain is if node 
! i is moved 
            CALL compute_gain(i,flags)
! Recalculate doubly linked list linking nodes with same gain and headers 
! to starts (up to cut off value icut) 
! Initialize array done that flags that a node has been considered in 
! an inner loop pass 
            gain = max(-mult,min(gain1(i),gain2(i)))
           ! gain = min(gain1(i),gain2(i))
            IF (gain>icut) CYCLE
            IF (gain<ming) ming = gain
! New node is put at head of list 
            CALL add_to_list(i,gain)
          END DO
! End of outer loop 
        END DO

        a_n1 = nv1
        a_n2 = nv2
 
        RETURN

      CONTAINS
        SUBROUTINE remove_from_list(irm,ig)
          INTEGER :: irm, ig

          inext = next(irm)
          ilast = last(irm)
          IF (ilast==0) THEN
            head(ig) = inext
            IF (inext/=0) last(inext) = 0
          ELSE
            next(ilast) = inext
            IF (inext/=0) last(inext) = ilast
          END IF
        END SUBROUTINE remove_from_list


        SUBROUTINE add_to_list(irm,ig)
          INTEGER :: irm, ig

          inext = head(ig)
          head(ig) = irm
          next(irm) = inext
          IF (inext/=0) last(inext) = irm
          last(irm) = 0
        END SUBROUTINE add_to_list

        SUBROUTINE compute_gain(i,partit)
          INTEGER :: i, partit(:)
          INTEGER :: j, jj,l
! Initialize gain ... knowing node i will be removed from cutset 
! The +1 is to give identical result to previous code when unit weights 
          gain1(i) = -weight(i) + 1
          gain2(i) = -weight(i) + 1
          IF (i.EQ.n) THEN
              l = a_ne
          ELSE
              l = ptr(i+1) - 1
          END IF
          DO jj = ptr(i), l
            j = col(jj)
! Check which partition node j is in and adjust gain array appropriately 
            IF (partit(j)==mc70_part1_flag) THEN
              gain2(i) = gain2(i) + weight(j)
            END IF
            IF (partit(j)==mc70_part2_flag) THEN
              gain1(i) = gain1(i) + weight(j)
            END IF
          END DO
        END SUBROUTINE compute_gain
      END SUBROUTINE mc70_fm_refinement_band

      ! The subroutine mc70_fm_refinement uses a version of the Fiduccia-
! Mattheyses refinement algorithm on a tripartite partitioing of the nodes of a 
! graph where a node in the first partition is not connected to any node  
! in the second partition and any path between nodes in partition 1 and 
! partition 2 must go through a node in the cutset (partition 3). 
! The intention of the algorithm is to reduce f(P1,P2,P3), where
! f(P1,P2,P3) =   |P3|/(|P1||P2|) if min(|P1|,|P2|)/max(|P1|,|P2|) >= ratio
!             =  sumweight - 2 + max(|P1|,|P2|)/min(|P1|,|P2|), otherwise
! This is a banded version so only nodes a distance of at most band from the 
! input separator can be moved into the new separator

      SUBROUTINE mc70_fm_refinement_band_balance(n,a_ne,ptr,col,weight,&
          sumweight,icut,mult,a_n1,a_n2,wnv1,wnv2,wns,band,ratio, &
          flags,ipart,next,last,gain1,gain2,done,head,dist)

! Matrix is held in matrix using compressed column scheme 
        INTEGER, INTENT(IN) :: n  ! size of matrix
        INTEGER, INTENT(IN) :: a_ne  ! no. nonzeros in matrix
        INTEGER, INTENT(IN) :: ptr(n) ! row pointers
        INTEGER, INTENT(IN) :: col(a_ne) ! column indices
       ! TYPE (zd11_type), INTENT (INOUT) :: matrix
! The array weight is used to hold a weight on the vertices indicating 
! how many vertices from the finer graphs have been combined into the 
! current coarse graph vertex. 
        INTEGER, INTENT (IN) :: weight(n)
        INTEGER, INTENT(IN) :: sumweight
        INTEGER, INTENT(IN) :: icut ! Used to limit search
        INTEGER, INTENT(IN) :: mult ! Used to bound search
        INTEGER, INTENT (INOUT) :: a_n1 ! No. vertices partition 1
        INTEGER, INTENT (INOUT) :: a_n2 ! No. vertices partition 2
        INTEGER, INTENT (INOUT) :: wnv1 ! Weighted sum of vertices partition 1
        INTEGER, INTENT (INOUT) :: wnv2 ! Weighted sum of vertices partition 2
        INTEGER, INTENT (INOUT) :: wns ! Weighted sum of vertices separator
        INTEGER, INTENT (IN) :: band ! width of band around initial separator
! that the separator can lie in
        REAL (myreal_mc70), INTENT(IN) :: ratio ! ratio to determine whether 
!partition is balanced
!
! flags holds a list of nodes stating which partition node i is in.
! The whole point of this routine is to return a revised partition 
! with better properties.  Normally less nodes in the cutset while maintaining 
! a balance between the number of nodes in the two components. 
! flags(i) == mc70_part1_flag : i is in partition 1
! flags(i) == mc70_part2_flag : i is in partition 2
! flags(i) == mc70_sep_flag   : i is in separator/cutset
        INTEGER, INTENT (INOUT) :: flags(n)
! info holds parameters giving information about the performance of the  
! subroutine 
        INTEGER, INTENT(OUT) :: ipart(n), next(n), last(n)
        INTEGER, INTENT(OUT) :: gain1(n),gain2(n),done(n),head(-mult:icut)
        INTEGER, INTENT(OUT) :: dist(n)

! Number nodes in each partition 
        INTEGER :: nv1, ns, nv2, inv1, inv2, ins
! Weighted nodes in each partition 
        INTEGER :: winv1, winv2, wins
        INTEGER :: i, j, ii, jj, eye, k,l
        INTEGER :: inn, outer
        INTEGER :: move, ming, gain, old_gain, inext, ilast, idummy
        INTEGER :: first,tail
        REAL (myreal_mc70) :: eval, evalc, evalo,eval1, eval2
        LOGICAL :: imbal


        IF (ratio.GT.real(sumweight-2)) THEN
           imbal = .FALSE.
        ELSE
           imbal = .TRUE.
        END IF
        
! Set-up distance to hold (min) distance of node from separator. Only find 
! distance for nodes within band distance from separator
        first = 0
        tail = 0
        DO i = 1,n
          IF (flags(i) .EQ. mc70_sep_flag) THEN
           dist(i) = 0
           IF (first .EQ. 0) THEN
             first = i
             tail = i
           ELSE
             next(tail) = i
             tail = i
           END IF
          ELSE
           dist(i) = -2
          END IF
        END DO
        
         DO WHILE (first .NE. 0)
          j = dist(first)
          IF (j .EQ. band-1) THEN
            IF (first .EQ. n) THEN
              l = a_ne
            ELSE
              l = ptr(first+1)-1
            END IF
            DO i = ptr(first),l
              IF (dist(col(i)) .EQ. -2) THEN
                k = col(i)
                dist(k) = j+1
              END IF
            END DO            

          ELSE
            IF (first .EQ. n) THEN
              l = a_ne
            ELSE
              l = ptr(first+1)-1
            END IF
            DO i = ptr(first),l
              IF (dist(col(i)) .EQ. -2) THEN
                k = col(i)
                dist(k) = j+1
                next(tail) = k
                tail = k
              END IF
            END DO

          END IF
          IF (first .EQ. tail) THEN
            first = 0
          ELSE
            k = next(first)
            first = k
          END IF
         END DO
        next(1:n) = 0
        
! nv1,nv2,ns are the number of nodes in partitions 1, 2 and the cutset in 
! the current partition 
! inv1,inv2,ins,ipart are the equivalent quantities within the inner loop 
! The same identifiers prefixed by w refer to weighted counts 
! inner and outer are the two main loop indices 

! Initialize nv1,nv2,ns 
        nv1 = a_n1
        nv2 = a_n2
        ns = n - (a_n1+a_n2)
        ii = 1
        jj = ii + nv1

!Initialize ipart 
        ipart(1:n) = flags(1:n)

! Initialize array done that flags that a node has been considered in 
! an inner loop pass 
        done = 0

! Compute evaluation function for current partitioning 

        CALL cost_function(wnv1+1,wnv2+1,wns,&
              sumweight,ratio,imbal,evalc)

! icut is set to limit search in inner loop .. may later be a parameter 
! we allow gains of up to max(weight)*5 

          head(-mult:icut) = 0

! Set up doubly linked list linking nodes with same gain and headers 
! to starts (up to cut off value icut) 

! Compute gains for nodes in cutset 
        ming = sumweight
        DO i = 1, n

         IF (flags(i).eq.mc70_sep_flag) THEN
! Node i is in cutset 
! gain1(i) is change to cutset size if node i is moved to partition 1. 
! gain2(i) is change to cutset size if node i is moved to partition 2. 
! Run through all neighbours of node i to see what loss/gain is if node 
! i is moved 

          CALL compute_gain(i,flags)
          gain = max(-mult,min(gain1(i),gain2(i)))
          
          IF (gain<ming) ming = gain
          IF (gain>icut) CYCLE
! New node is put at head of list 
          CALL add_to_list(i,gain)
         END IF
        END DO
       
! Initilialization finished.  Now perform F-M algorithm in two loops. 
! In each inner loop we choose the best obtained and if the evaluation 
! function for this is better than previous best we perform another inner 
! loop; otherwise we terminate. 
        evalo = evalc
        DO outer = 1, n
! Set partition that will be altered in inner loop 
          inv1 = nv1
          inv2 = nv2
          ins = ns
          winv1 = wnv1
          winv2 = wnv2
          wins = wns
          ipart(1:n) = flags(1:n)
INNER:    DO inn = 1, n

! Choose best eligible move 
            DO idummy = 1, n

               DO gain = ming, icut
                IF (head(gain)/=0)  EXIT
               END DO
              IF (gain>icut) EXIT INNER

! Now cycle through nodes of least gain 
! Currently inefficient because of re-searching linked list 
              inext = head(gain)
              k = 0
111           i = inext
              IF (i==0) CYCLE
! Use node if it has not been considered already 
              IF (done(i)<outer) GO TO 222
              inext = next(i)
!!! Extra statements to trap infinite loop 
              k = k + 1
              IF (k>ins) THEN
                !WRITE (*,*) 'Bug in code because of infinite loop'
!!! You may wish to change this to a stop 
                EXIT
              END IF
              GO TO 111
            END DO
            EXIT INNER
! Node i has been selected as the best eligible node 
! Set flag so only considered once in this pass 
222         done(i) = outer
! As i will not be chosen again in this pass, remove from list 
            CALL remove_from_list(i,gain)
! Move the node to the appropriate partition and reset partition information 
! We will try both weighted and unweighted 

            IF (wnv1 .eq. 0 .and. wnv2 .gt. 0) THEN
            
! Move node i to partition 1 
              move = mc70_part1_flag
              inv1 = inv1 + 1
              winv1 = winv1 + weight(i)
              ins = ins - 1
              wins = wins - weight(i)
            ELSE IF (wnv2 .eq. 0 .and. wnv1 .gt. 0) THEN
! Move node i to partition 2 
              move = mc70_part2_flag
              inv2 = inv2 + 1
              winv2 = winv2 + weight(i)
              ins = ins - 1
              wins = wins - weight(i)
            
            ELSE
             CALL cost_function(winv1+weight(i),winv2+1-gain1(i)-weight(i),&
              wins+gain1(i)-1,&
              sumweight,ratio,imbal,eval1)
             
             CALL cost_function(winv1+1-gain2(i)-weight(i),&
              winv2+weight(i),wins+gain2(i)-1,&
              sumweight,ratio,imbal,eval2)
            IF ((eval1<eval2) .OR. ( (eval1==eval2) .AND. (wnv1<  wnv2) )) THEN
! Move node i to partition 1 
              move = mc70_part1_flag
              inv1 = inv1 + 1
              winv1 = winv1 + weight(i)
              ins = ins - 1
              wins = wins - weight(i)
            ELSE
! Move node i to partition 2 
              move = mc70_part2_flag
              inv2 = inv2 + 1
              winv2 = winv2 + weight(i)
              ins = ins - 1
              wins = wins - weight(i)
            END IF
            END IF
! Set new partition for node i 
            ipart(i) = move
! Run through neigbours of node i to update data 
            IF (i.EQ.n) THEN
              l = a_ne
            ELSE
              l = ptr(i+1) - 1
            END IF
            DO jj = ptr(i), l
              j = col(jj)
! Check which partition node j is in and take appropriate action 
              IF (ipart(j)==move) CYCLE
! If node j is in cutset, update its gain value 
              IF (ipart(j)==mc70_sep_flag) THEN
! If it has already been chosen in this pass just skip it 
                IF (done(j)==outer .or. dist(j) .eq. -2) CYCLE
! old_gain is present gain 
             
                old_gain = max(-mult,min(gain1(j),gain2(j)))
               ! old_gain = min(gain1(j),gain2(j))

                IF (move==mc70_part1_flag) gain2(j) = gain2(j) + weight(i)
                IF (move==mc70_part2_flag) gain1(j) = gain1(j) + weight(i)
                gain = max(-mult,min(gain1(j),gain2(j)))
               ! gain = min(gain1(j),gain2(j))

                IF (old_gain==gain) CYCLE
! Remove from old list 
                IF (old_gain<=icut) THEN
                  CALL remove_from_list(j,old_gain)
                END IF
! gain has changed so move to new linked list if less than icut 
                IF (gain<=icut) THEN
! Reset ming if necessary 
                  IF (gain<ming) ming = gain
                  CALL add_to_list(j,gain)
                END IF
              END IF
              IF (ipart(j)==2-move) THEN
! We have a new node in the cutset. 
                ipart(j) = mc70_sep_flag
! Compute gains for this new node in the cutset and place in linked list 
! We intentionally did not do this earlier but we do now 
! [maybe not since won't access this node again in this pass] 
! We use done array to record this but probably not necessary as not put 
! in head linked list so won't be accessed 
! First check that it was not earlier moved from cutset 
                IF (done(j)/=outer .and. dist(j) .ne. -2) THEN
! Compute gain 
                  CALL compute_gain(j,ipart)
                  gain = max(-mult,min(gain1(j),gain2(j)))
                !  gain = min(gain1(j),gain2(j))
!!! Just added this 
                  IF (gain<ming) ming = gain
! Add to  list 
                  IF (gain<=icut) THEN
                      CALL add_to_list(j,gain)
                  END IF
                END IF
! Update partition and gain of any nodes in cutset connected to node j 
                ins = ins + 1
                wins = wins + weight(j)
                IF (move==mc70_part1_flag) THEN
                  inv2 = inv2 - 1
                  winv2 = winv2 - weight(j)
                END IF
                IF (move==mc70_part2_flag) THEN
                  inv1 = inv1 - 1
                  winv1 = winv1 - weight(j)
                END IF
! Check neighbours of j since any in cut set will have gain changed 
                IF (j.EQ.n) THEN
                  l = a_ne
                ELSE
                  l = ptr(j+1) - 1
                END IF
                DO ii = ptr(j), l
                  eye = col(ii)
                  IF (ipart(eye)/=mc70_sep_flag) CYCLE
                  IF (dist(eye) .eq. -2) CYCLE
! Neighbour is in cutset. Recompute gain and insert in linked list. 
                  IF (done(eye)==outer) CYCLE
! old_gain is present gain 
                  old_gain = max(-mult,min(gain1(eye),gain2(eye)))
                 ! old_gain = min(gain1(eye),gain2(eye))
                  
                  
                  IF (move==mc70_part1_flag) THEN
                     gain1(eye) = gain1(eye) - weight(j)
                  END IF
                  IF (move==mc70_part2_flag) THEN
                     gain2(eye) = gain2(eye) - weight(j)
                  END IF
! gain is new gain 
                  gain = max(-mult,min(gain1(eye),gain2(eye)))
                 ! gain = min(gain1(eye),gain2(eye))
                  IF (old_gain==gain) CYCLE
! Remove from old list 
                  IF (old_gain<=icut) THEN
                      CALL remove_from_list(eye,old_gain)
                  END IF
! gain has changed so move to new linked list if less than icut 
                  IF (gain<=icut) THEN
! Reset ming if necessary 
                    IF (gain<ming) ming = gain
                    CALL add_to_list(eye,gain)
                  END IF
                END DO
              END IF
! end of neighbours loop 
            END DO

!         ii = 0 
!         do i = 1,n 
!           if (ipart(i) == 2) ii = ii + 1 
!         enddo 
!         if (ii .ne. inv2) write(6,*) 'problem in partition',ii,inv2 

! Evaluate new partition 
             CALL cost_function(winv1+1,winv2+1,wins,&
              sumweight,ratio,imbal,eval)
! Compare this with best so far in inner loop and store partition 
! information if it is the best 
            IF (inv1*inv2>0 .AND. nv1*nv2.eq.0) THEN
! Might have to store gains and who is in the cutset 
              evalc = eval
              nv1 = inv1
              nv2 = inv2
              ns = ins
              wnv1 = winv1
              wnv2 = winv2
              wns = wins
              flags = ipart

            ELSE IF (eval<evalc .AND. (inv1*inv2>0 )) THEN
! Might have to store gains and who is in the cutset 
              evalc = eval
              nv1 = inv1
              nv2 = inv2
              ns = ins
              wnv1 = winv1
              wnv2 = winv2
              wns = wins
              flags(1:n) = ipart(1:n)
            END IF
! End inner loop 
          END DO INNER
! Leave loop if inner loop has not found better partition 
       !   write(*,*) evalc,evalo
          IF (evalc>=(1.0-1.0/(log(real(sumweight))**2.3))*evalo) EXIT
! Otherwise we reset evalo and go back to inner loop 
         ! write(*,*) 'e',evalo,evalc,evalc/evalo
          evalo = evalc
! Recompute gains for this new partition 
! Compute gains for nodes in cutset 
! This is very inefficient but is in now to test functionality 
          head(-mult:icut) = 0
          ming = icut + 1
          DO i = 1, n
            IF (flags(i)/=mc70_sep_flag) CYCLE
            IF (dist(i) .eq. -2) CYCLE
! Node i is in cutset 
! gain1(i) is change to cutset size if node i is moved to partition 1. 
! gain2(i) is change to cutset size if node i is moved to partition 2. 
! Run through all neighbours of node i to see what loss/gain is if node 
! i is moved 
            CALL compute_gain(i,flags)
! Recalculate doubly linked list linking nodes with same gain and headers 
! to starts (up to cut off value icut) 
! Initialize array done that flags that a node has been considered in 
! an inner loop pass 
            gain = max(-mult,min(gain1(i),gain2(i)))
           ! gain = min(gain1(i),gain2(i))
            IF (gain>icut) CYCLE
            IF (gain<ming) ming = gain
! New node is put at head of list 
            CALL add_to_list(i,gain)
          END DO
! End of outer loop 
        END DO
        a_n1 = nv1
        a_n2 = nv2 
        RETURN

      CONTAINS
        SUBROUTINE remove_from_list(irm,ig)
          INTEGER :: irm, ig

          inext = next(irm)
          ilast = last(irm)
          IF (ilast==0) THEN
            head(ig) = inext
            IF (inext/=0) last(inext) = 0
          ELSE
            next(ilast) = inext
            IF (inext/=0) last(inext) = ilast
          END IF
        END SUBROUTINE remove_from_list

        SUBROUTINE add_to_list(irm,ig)
          INTEGER :: irm, ig

          inext = head(ig)
          head(ig) = irm
          next(irm) = inext
          IF (inext/=0) last(inext) = irm
          last(irm) = 0
        END SUBROUTINE add_to_list

        SUBROUTINE compute_gain(i,partit)
          INTEGER :: i, partit(:)
          INTEGER :: j, jj,l
! Initialize gain ... knowing node i will be removed from cutset 
! The +1 is to give identical result to previous code when unit weights 
          gain1(i) = -weight(i) + 1
          gain2(i) = -weight(i) + 1
          IF (i.EQ.n) THEN
              l = a_ne
          ELSE
              l = ptr(i+1) - 1
          END IF
          DO jj = ptr(i), l
            j = col(jj)
! Check which partition node j is in and adjust gain array appropriately 
            IF (partit(j)==mc70_part1_flag) THEN
              gain2(i) = gain2(i) + weight(j)
            END IF
            IF (partit(j)==mc70_part2_flag) THEN
              gain1(i) = gain1(i) + weight(j)
            END IF
          END DO
        END SUBROUTINE compute_gain
      END SUBROUTINE mc70_fm_refinement_band_balance

! -------------------------------------------------------------------
! hamd is an implementation of the halo_amd method by Pellegrini, Roman and 
! Amestoy. 

! We use the term Le to denote the set of all supervariables in element
! E.
! -------------------------------------------------------------------
      SUBROUTINE hamd(n,ne,lirn,irn,ip,sep,perm,work)
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: lirn, n,ne
! ..
! .. Array Arguments ..
        INTEGER, INTENT (INOUT) :: irn(lirn), ip(n)
        INTEGER, OPTIONAL, INTENT (IN) :: sep(n)
        INTEGER, INTENT (OUT) :: perm(n)
        INTEGER, INTENT (OUT) ::  work(7*n)


! N must be set to the matrix order.
! Restriction:  N .ge. 1
! NE is set to number of non-zeros stored

! lirn must be set to the length of irn. It is not altered. On input,
! the matrix is stored in irn (1..ne).
! *** We do not recommend running this algorithm with ***
! ***      lirn .LT. NE+1 + N.                      ***
! *** Better performance will be obtained if          ***
! ***      lirn .GE. NE+1 + N                       ***
! *** or better yet                                   ***
! ***      lirn .GT. 1.2 *(NE+1)                    ***
! Restriction: lirn .GE. NE

! irn(1..NE) must be set to  hold the patterns of the rows of
! the matrix.  The matrix must be symmetric, and both upper and
! lower triangular parts must be present.  The diagonal must not be
! present.  Row I is held as follows:
! irn(ip(I)...ip(I) + work(len+I) - 1) must hold the list of
! column indices for entries in row I (simple
! supervariables), excluding the diagonal.  All
! supervariables start with one row/column each
! (supervariable I is just row I). If work(len+I) is zero on
! input, then ip(I) is ignored on input. Note that the
! rows need not be in any particular order, and there may
! be empty space between the rows.
! During execution, the supervariable I experiences fill-in. This
! is represented by constructing a list of the elements that cause
! fill-in in supervariable I:
! IE(ip(i)...ip(I) + perm(I) - 1) is the list of elements
! that contain I. This list is kept short by removing
! absorbed elements. irn(ip(I)+perm(I)...ip(I)+work(len+I)-1)
! is the list of supervariables in I. This list is kept
! short by removing nonprincipal variables, and any entry
! J that is also contained in at least one of the
! elements in the list for I.
! When supervariable I is selected as pivot, we create an element E
! of the same name (E=I):
! IE(ip(E)..ip(E)+work(len+E)-1) is the list of supervariables
! in element E.
! An element represents the fill-in that occurs when supervariable
! I is selected as pivot.
! CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
! The contents of irn are undefined on output.

! ip(i) must be set to the the index in irn of the start of row I, or be
! zero if row I has no off-diagonal entries. During execution,
! it is used for both supervariables and elements:
! * Principal supervariable I:  index into irn of the
! list of supervariable I.  A supervariable
! represents one or more rows of the matrix
! with identical pattern.
! * Non-principal supervariable I:  if I has been absorbed
! into another supervariable J, then ip(I) = -J.
! That is, J has the same pattern as I.
! Note that J might later be absorbed into another
! supervariable J2, in which case ip(I) is still -J,
! and ip(J) = -J2.
! * Unabsorbed element E:  the index into irn of the list
! of element E.  Element E is created when
! the supervariable of the same name is selected as
! the pivot.
! * Absorbed element E:  if element E is absorbed into element
! E2, then ip(E) = -E2.  This occurs when one of its
! variables is eliminated and when the pattern of
! E (that is, Le) is found to be a subset of the pattern
! of E2 (that is, Le2).  If element E is "null" (it has
! no entries outside its pivot block), then ip(E) = 0.

! On output, ip holds the assembly tree/forest, which implicitly
! represents a pivot order with identical fill-in as the actual
! order (via a depth-first search of the tree). If work(nv+I) .GT. 0,
! then I represents a node in the assembly tree, and the parent of
! I is -ip(I), or zero if I is a root. If work(nv+I)=0, then (I,-ip(I))
! represents an edge in a subtree, the root of which is a node in
! the assembly tree.

! sep is an array that is used indicate whether a row/column is in the seperator 
! or not. If row i is in the seperator, then sep(i) must equal mc70_sep_flag;
! otherwise, sep(i) most be set to a value that is not equal to mc70_sep_flag

! perm(I) need not be set. See the description of irn above. At the
! start of execution, perm(I) is set to zero. For a supervariable,
! perm(I) is the number of elements in the list for supervariable
! I. For an element, perm(E) is the negation of the position in the
! pivot sequence of the supervariable that generated it. perm(I)=0
! if I is nonprincipal.
! On output perm(1..N) holds the permutation. That is, if K = perm(I),
! then row K is the ith pivot row.  Row K of A appears as the
! I-th row in the permuted matrix, PAP^T.

! CONTROL is of type mc70_control and contains control
! parameters and must be set by the user. 


! Local arrays:
! ---------------

! WORK(NV+1:NV+N) During execution, ABS(WORK(NV+I)) is equal to the
! number of rows represented by the principal supervariable I. If I
! is a nonprincipal variable, then WORK(NV+I) = 0. Initially, WORK(NV+I) = 1
! for all I.  WORK(NV+I) .LT. 0 signifies that I is a principal variable
! in the pattern Lme of the current pivot element ME. On termination,
! WORK(NV+E) holds the true degree of element E at the time it was
! created (including the diagonal part).

! WORK(LAST+1:LAST+N) In a degree list, work(last+I) is the
! supervariable preceding I, or zero if I is the head of the list.
! In a hash bucket, work(last+I) is the hash key for I. work(last+work(head+HASH))
! is also used as the head of a hash bucket if work(head+HASH) contains
! a degree list (see HEAD, below).
! On output, work(last+1..last+N) holds the permutation (the same as the
! 'PERM' argument in Sparspak). That is, if I = work(last+K), then row I
! is the Kth pivot row.  Row work(last+K) of A is the K-th row in the
! permuted matrix, PAP^T.

! work(len+I) is initialised to hold the number of entries in row I of the
! matrix, excluding the diagonal.  The contents of work(len+1..N) are
! undefined on output.

! DEGREE If I is a supervariable and sparse,
! then work(degree+I) holds the current approximation of the external
! degree of row I (an upper bound). The external degree is the
! number of entries in row I, minus ABS(work(nv+I)) (the diagonal
! part). The bound is equal to the external degree if perm(I) is
! less than or equal to two. We also use the term "external degree"
! for elements E to refer to |Le \ Lme|. If I is full in the reduced
! matrix, then work(degree+I)=N+1. If I is dense in the reduced matrix,
! then work(degree+I)=N+1+last_approximate_external_deg of I.

! work(head+DEG) is used for degree lists.
! work(head+DEG) is the first supervariable in a degree list (all
! supervariables I in a degree list DEG have the same approximate
! degree, namely, DEG = work(degree+I)). If the list DEG is empty then
! work(head+DEG) = 0.
! During supervariable detection work(head+HASH) also serves as a
! pointer to a hash bucket.
! If work(head+HASH) .GT. 0, there is a degree list of degree HASH. The
! hash bucket head pointer is work(last+work(head+HASH)).
! If work(head+HASH) = 0, then the degree list and hash bucket are
! both empty.
! If work(head+HASH) .LT. 0, then the degree list is empty, and
! -work(head+HASH) is the head of the hash bucket.
! After supervariable detection is complete, all hash buckets are
! empty, and the (work(last+work(head+HASH)) = 0) condition is restored for
! the non-empty degree lists.

! work(denxt+I)  For supervariable I, work(denxt+I) is
! the supervariable following I in a link list, or zero if I is
! the last in the list. Used for two kinds of lists: degree lists
! and hash buckets (a supervariable can be in only one kind of
! list at a time). For element E, work(denxt+E) is the number of
! variables with dense or full rows in the element E.

! work(w+I) The flag array W determines the status
! of elements and variables, and the external degree of elements.
! For elements:
! if work(w+E) = 0, then the element E is absorbed.
! if work(w+E) .GE. WFLG, then work(w+E)-WFLG is the size of the set
! |Le \ Lme|, in terms of nonzeros (the sum of ABS(work(nv+I))
! for each principal variable I that is both in the
! pattern of element E and NOT in the pattern of the
! current pivot element, ME).
! if WFLG .GT. WE(E) .GT. 0, then E is not absorbed and has
! not yet been seen in the scan of the element lists in
! the computation of |Le\Lme| in loop 150 below.
! ***SD: change comment to remove reference to label***
! For variables:
! during supervariable detection, if work(w+J) .NE. WFLG then J is
! not in the pattern of variable I.
! The W array is initialized by setting work(w+I) = 1 for all I, and by
! setting WFLG = 2. It is reinitialized if WFLG becomes too large
! (to ensure that WFLG+N does not cause integer overflow).

! PFREE must be set to the position in irn of the first free variable.
! During execution, additional data is placed in irn, and PFREE is
! modified so that components  of irn from PFREE are free.
! On output, PFREE is set equal to the size of irn that would have
! caused no compressions to occur.  If NCMPA is zero, then
! PFREE (on output) is less than or equal to lirn, and the space
! irn(PFREE+1 ... lirn) was not used. Otherwise, PFREE (on output)
! is greater than lirn, and all the memory in irn was used.

! Local variables:
! ---------------

! DEG:        the degree of a variable or element
! DEGME:      size (no. of variables), |Lme|, of the current element,
! ME (= work(degree+ME))
! DEXT:       external degree, |Le \ Lme|, of some element E
! DMAX:       largest |Le| seen so far
! E:          an element
! permME:     the length, perm(ME), of element list of pivotal var.
! ELN:        the length, perm(...), of an element list
! HASH:       the computed value of the hash function
! HMOD:       the hash function is computed modulo HMOD = MAX(1,N-1)
! I:          a supervariable
! IDUMMY:     loop counter
! ILAST:      the entry in a link list preceding I
! INEXT:      the entry in a link list following I
! IOVFLO:     local copy of ICNTL(5)
! J:          a supervariable
! JDUMMY:     loop counter
! JLAST:      the entry in a link list preceding J
! JNEXT:      the entry in a link list, or path, following J
! K:          the pivot order of an element or variable
! KNT1:       loop counter used during element construction
! KNT2:       loop counter used during element construction
! KNT3:       loop counter used during element construction
! LENJ:       work(len+J)
! LN:         length of a supervariable list
! MAXMEM:     amount of memory needed for no compressions
! ME:         current supervariable being eliminated, and the
! current element created by eliminating that
! supervariable
! MEM:        memory in use assuming no compressions have occurred
! MINDEG:     current approximate minimum degree
! NCMPA:      counter for the number of times irn was compressed
! NEL:        number of pivots selected so far
! NEWMEM:     amount of new memory needed for current pivot element
! NLEFT:      N-NEL, the number of nonpivotal rows/columns remaining
! NRLADU:     counter for the forecast number of reals in matrix factor
! NVI:        the number of variables in a supervariable I (= work(nv+I))
! NVJ:        the number of variables in a supervariable J (= work(nv+J))
! NVPIV:      number of pivots in current element
! P:          pointer into lots of things
! P1:         ip (i) for some variable i (start of element list)
! P2:         ip (i) + perm (i) -  1 for some var. i (end of el. list)
! P3:         index of first supervariable in clean list
! PJ:         pointer into an element or variable
! PDST:       destination pointer, for compression
! PEND:       end of memory to compress
! PFREE is set to the position in irn of the first free variable.
! At the end, PFREE is set equal to the size of irn that would have
! caused no compressions to occur.  If NCMPA is zero, then
! PFREE (on output) is less than or equal to lirn, and the space
! irn(PFREE+1 ... lirn) was not used. Otherwise, PFREE (on output)
! is greater than lirn, and all the memory in irn was used.
! PME:        pointer into the current element (PME1...PME2)
! PME1:       the current element, ME, is stored in irn(PME1...PME2)
! PME2:       the end of the current element
! PN:         pointer into a "clean" variable, also used to compress
! PSRC:       source pointer, for compression
! SLENME:     number of variables in variable list of pivotal variable
! WE:         work(w+E)
! WFLG:       used for flagging the W array.  See description of W.
! WNVI:       WFLG-work(nv+I)
! X:          either a supervariable or an element

! OPS:        counter for forecast number of flops

! IDENSE is true if supervariable I is dense

! -------------------------------------------------------------------
! FUNCTIONS CALLED:
! -------------------------------------------------------------------

! ====================================================================
! INITIALIZATIONS
! ====================================================================

! ..
! .. Local Arrays ..

  
! ..
! .. Local Scalars ..
        INTEGER :: nv, last,degree,head,denxt,w,len
        INTEGER deg, degme, dext, dmax, e, permme, eln, hash, hmod, i, &
          idummy, ilast, inext, iovflo, j, jdummy, jlast, jnext, k, knt1, &
          knt2, knt3, lenj, ln, maxmem, me, mem, mindeg, ncmpa, &
          nel, newmem, nleft, nvi, nvj, nvpiv, p, p1, p2, p3, pdst, pend, pj, &
          pme, pme1, pme2, pn, psrc, slenme, we, wflg, wnvi, x, pfree, &
          nosep, l,temp
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, max, min, mod
! ..
       
        me = 0
        nosep = 0
        DO i = 1, n
          IF (sep(i)==mc70_sep_flag) nosep = nosep + 1
        END DO
        nv = 0
        last = nv + n
        degree = last +n
        head = degree+n
        denxt = head+n
        w = denxt+n
        len = w+n

        dmax = 0
        hmod = max(1,n-1)
        iovflo = huge(0)
        pfree = ne+1
        mem = pfree - 1
        maxmem = mem
        mindeg = 1
        ncmpa = 0
        nel = 0
        wflg = 2

! Assign len and pfree
        work(len+1:len+n-1) = ip(2:n) - ip(1:n-1)
        work(len+n) = ne+1 - ip(n)

! ----------------------------------------------------------
! initialize arrays and eliminate rows with no off-diag. nz.
! ----------------------------------------------------------
        work(last+1:last+n) = 0
        work(head+1:head+n) = 0
        work(nv+1:nv+n) = 1
        work(degree+1:degree+n) = work(len+1:len+n)
        DO i = 1, n
          IF (work(degree+i)==0 .AND. sep(i)/=mc70_sep_flag) THEN
            nel = nel + 1
            perm(i) = -nel
            ip(i) = 0
            work(w+i) = 0
          ELSE
            work(w+i) = 1
            perm(i) = 0
          END IF
        END DO
! ----------------------------------------------------------------
! initialize degree lists
! ----------------------------------------------------------------
        DO i = 1, n
          deg = work(degree+i)
          IF (deg>0 .AND. sep(i)/=mc70_sep_flag) THEN
! ----------------------------------------------------------
! place i in the degree list corresponding to its degree
! or in the dense row list if i is dense
! ----------------------------------------------------------
! place i in the degree list corresponding to its degree
            inext = work(head+deg)
            IF (inext/=0) work(last+inext) = i
            work(denxt+i) = inext
            work(head+deg) = i
          END IF
        END DO

        DO WHILE (nel<n-nosep)

! ==================================================================
! GET PIVOT OF MINIMUM APPROXIMATE DEGREE
! ==================================================================
! -------------------------------------------------------------
! find next supervariable for elimination
! -------------------------------------------------------------
          DO deg = mindeg, n
            me = work(head+deg)
            IF (me>0) GO TO 20
          END DO
20        mindeg = deg

! -------------------------------------------------------------
! remove chosen variable from linked list
! -------------------------------------------------------------
          inext = work(denxt+me)
          IF (inext/=0) work(last+inext) = 0
          work(head+deg) = inext
! -------------------------------------------------------------
! me represents the elimination of pivots nel+1 to nel+work(nv+me).
! place me itself as the first in this set.  It will be moved
! to the nel+work(nv+me) position when the permutation vectors are
! computed.
! -------------------------------------------------------------
          permme = perm(me)
          perm(me) = -(nel+1)
          nvpiv = work(nv+me)
          nel = nel + nvpiv
          work(denxt+me) = 0

! ====================================================================
! CONSTRUCT NEW ELEMENT
! ====================================================================

! -------------------------------------------------------------
! At this point, me is the pivotal supervariable.  It will be
! converted into the current element.  Scan list of the
! pivotal supervariable, me, setting tree pointers and
! constructing new list of supervariables for the new element,
! me.  p is a pointer to the current position in the old list.
! -------------------------------------------------------------

! flag the variable "me" as being in the front by negating work(nv+me)
          work(nv+me) = -nvpiv
          degme = 0
          IF (permme==0) THEN
! ----------------------------------------------------------
! There are no elements involved.
! Construct the new element in place.
! ----------------------------------------------------------
            pme1 = ip(me)
            pme2 = pme1 - 1
            DO p = pme1, pme1 + work(len+me) - 1
              i = irn(p)
              nvi = work(nv+i)
              IF (nvi>0) THEN
! ----------------------------------------------------
! i is a principal variable not yet placed in the
! generated element. Store i in new list
! ----------------------------------------------------
                degme = degme + nvi
! flag i as being in Lme by negating nv (i)
                work(nv+i) = -nvi
                pme2 = pme2 + 1
                irn(pme2) = i

! ----------------------------------------------------
! remove variable i from degree list.
! ----------------------------------------------------
                IF (sep(i)/=mc70_sep_flag) THEN
                  ilast = work(last+i)
                  inext = work(denxt+i)
                  IF (inext/=0) work(last+inext) = ilast
                  IF (ilast/=0) THEN
                    work(denxt+ilast) = inext
                  ELSE
! i is at the head of the degree list
                    temp = work(degree+i)
                    work(head+temp) = inext
                  END IF
                END IF
              END IF
            END DO
! this element takes no new memory in irn:
            newmem = 0
          ELSE
! ----------------------------------------------------------
! construct the new element in empty space, irn (pfree ...)
! ----------------------------------------------------------
            p = ip(me)
            pme1 = pfree
            slenme = work(len+me) - permme
            DO knt1 = 1, permme
! search the elements in me.
              e = irn(p)
              p = p + 1
              pj = ip(e)
              ln = work(len+e)
! -------------------------------------------------------
! search for different supervariables and add them to the
! new list, compressing when necessary.
! -------------------------------------------------------
              DO knt2 = 1, ln
                i = irn(pj)
                pj = pj + 1
                nvi = work(nv+i)
                IF (nvi>0) THEN
! -------------------------------------------------
! compress irn, if necessary
! -------------------------------------------------
                  IF (pfree>lirn) THEN
! prepare for compressing irn by adjusting
! pointers and lengths so that the lists being
! searched in the inner and outer loops contain
! only the remaining entries.
! ***** SD: Seperate compression subroutine tried
! but found to be inefficient in comparison ****
                    ip(me) = p
                    work(len+me) = work(len+me) - knt1
! Check if anything left in supervariable ME
                    IF (work(len+me)==0) ip(me) = 0
                    ip(e) = pj
                    work(len+e) = ln - knt2
! Check if anything left in element E
                    IF (work(len+e)==0) ip(e) = 0
                    ncmpa = ncmpa + 1
! store first item in ip
! set first entry to -item
                    DO j = 1, n
                      pn = ip(j)
                      IF (pn>0) THEN
                        ip(j) = irn(pn)
                        irn(pn) = -j
                      END IF
                    END DO

! psrc/pdst point to source/destination
                    pdst = 1
                    psrc = 1
                    pend = pme1 - 1

! while loop:
                    DO idummy = 1, lirn
                      IF (psrc>pend) THEN
                        GO TO 30
                      ELSE
! search for next negative entry
                        j = -irn(psrc)
                        psrc = psrc + 1
                        IF (j>0) THEN
                          irn(pdst) = ip(j)
                          ip(j) = pdst
                          pdst = pdst + 1
! copy from source to destination
                          lenj = work(len+j)
                          DO knt3 = 0, lenj - 2
                            irn(pdst+knt3) = irn(psrc+knt3)
                          END DO
                          pdst = pdst + lenj - 1
                          psrc = psrc + lenj - 1
                        END IF
                      END IF
                    END DO

! move the new partially-constructed element
30                  p1 = pdst
                    DO psrc = pme1, pfree - 1
                      irn(pdst) = irn(psrc)
                      pdst = pdst + 1
                    END DO
                    pme1 = p1
                    pfree = pdst
                    pj = ip(e)
                    p = ip(me)
                  END IF

! -------------------------------------------------
! i is a principal variable not yet placed in Lme
! store i in new list
! -------------------------------------------------
                  degme = degme + nvi
! flag i as being in Lme by negating nv (i)
                  work(nv+i) = -nvi
                  irn(pfree) = i
                  pfree = pfree + 1

! -------------------------------------------------
! remove variable i from degree link list
! -------------------------------------------------
                  IF (sep(i)/=mc70_sep_flag) THEN
                    ilast = work(last+i)
                    inext = work(denxt+i)
                    IF (inext/=0) work(last+inext) = ilast
                    IF (ilast/=0) THEN
                      work(denxt+ilast) = inext
                    ELSE
! i is at the head of the degree list
                      temp = work(degree+i)
                      work(head+temp) = inext
                    END IF
                  END IF
                END IF
              END DO

! set tree pointer and flag to indicate element e is
! absorbed into new element me (the parent of e is me)
              IF (e/=me) THEN
                IF (sep(e)/=mc70_sep_flag) THEN
                  ip(e) = -me
                  work(w+e) = 0
                ELSE
                  ip(e) = 0
                  work(w+e) = 0
                END IF
              END IF
            END DO

! search the supervariables in me.
            knt1 = permme + 1
            e = me
            pj = p
            ln = slenme

! -------------------------------------------------------
! search for different supervariables and add them to the
! new list, compressing when necessary.
! -------------------------------------------------------
            DO knt2 = 1, ln
              i = irn(pj)
              pj = pj + 1
              nvi = work(nv+i)
              IF (nvi>0) THEN
! -------------------------------------------------
! compress irn, if necessary
! -------------------------------------------------
                IF (pfree>lirn) THEN
! prepare for compressing irn by adjusting
! pointers and lengths so that the lists being
! searched in the inner and outer loops contain
! only the remaining entries.
                  ip(me) = p
                  work(len+me) = work(len+me) - knt1
! Check if anything left in supervariable ME
                  IF (work(len+me)==0) ip(me) = 0
                  ip(e) = pj
                  work(len+e) = ln - knt2
! Check if anything left in element E
                  IF (work(len+e)==0) ip(e) = 0
                  ncmpa = ncmpa + 1
! store first item in ip
! set first entry to -item
                  DO j = 1, n
                    pn = ip(j)
                    IF (pn>0) THEN
                      ip(j) = irn(pn)
                      irn(pn) = -j
                    END IF
                  END DO

! psrc/pdst point to source/destination
                  pdst = 1
                  psrc = 1
                  pend = pme1 - 1

! while loop:
! 122              CONTINUE
                  DO idummy = 1, lirn
                    IF (psrc>pend) THEN
                      GO TO 40
                    ELSE
! search for next negative entry
                      j = -irn(psrc)
                      psrc = psrc + 1
                      IF (j>0) THEN
                        irn(pdst) = ip(j)
                        ip(j) = pdst
                        pdst = pdst + 1
! copy from source to destination
                        lenj = work(len+j)
                        DO knt3 = 0, lenj - 2
                          irn(pdst+knt3) = irn(psrc+knt3)
                        END DO
                        pdst = pdst + lenj - 1
                        psrc = psrc + lenj - 1
                      END IF
                    END IF
                  END DO

! move the new partially-constructed element
40                p1 = pdst
                  DO psrc = pme1, pfree - 1
                    irn(pdst) = irn(psrc)
                    pdst = pdst + 1
                  END DO
                  pme1 = p1
                  pfree = pdst
                  pj = ip(e)
                  p = ip(me)
                END IF

! -------------------------------------------------
! i is a principal variable not yet placed in Lme
! store i in new list
! -------------------------------------------------
                degme = degme + nvi
! flag i as being in Lme by negating nv (i)
                work(nv+i) = -nvi
                irn(pfree) = i
                pfree = pfree + 1

! -------------------------------------------------
! remove variable i from degree link list
! -------------------------------------------------
                IF (sep(i)/=mc70_sep_flag) THEN
                  ilast = work(last+i)
                  inext = work(denxt+i)
                  IF (inext/=0) work(last+inext) = ilast
                  IF (ilast/=0) THEN
                    work(denxt+ilast) = inext
                  ELSE
! i is at the head of the degree list
                    temp = work(degree+i)
                    work(head+temp) = inext
                  END IF
                END IF
              END IF
            END DO

            pme2 = pfree - 1
! this element takes newmem new memory in irn (possibly zero)
            newmem = pfree - pme1
            mem = mem + newmem
            maxmem = max(maxmem,mem)
          END IF

! -------------------------------------------------------------
! me has now been converted into an element in irn (pme1..pme2)
! -------------------------------------------------------------
! degme holds the external degree of new element
          work(degree+me) = degme
          ip(me) = pme1
          work(len+me) = pme2 - pme1 + 1

! -------------------------------------------------------------
! make sure that wflg is not too large.  With the current
! value of wflg, wflg+n must not cause integer overflow
! -------------------------------------------------------------
          IF (wflg>iovflo-n) THEN
            DO x = 1, n
              IF (work(w+x)/=0) work(w+x) = 1
            END DO
            wflg = 2
          END IF

! ====================================================================
! COMPUTE (work(w+e) - wflg) = |Le(G')\Lme(G')| FOR ALL ELEMENTS
! where G' is the subgraph of G containing just the sparse rows)
! ====================================================================
! -------------------------------------------------------------
! Scan 1:  compute the external degrees of elements touched
! with respect to the current element.  That is:
! (w (e) - wflg) = |Le \ Lme|
! for each element e involving a supervariable in Lme.
! The notation Le refers to the pattern (list of
! supervariables) of a previous element e, where e is not yet
! absorbed, stored in irn (ip (e) + 1 ... ip (e) + irn (ip (e))).
! The notation Lme refers to the pattern of the current element
! (stored in irn (pme1..pme2)).
! -------------------------------------------------------------
          DO pme = pme1, pme2
            i = irn(pme)
            eln = perm(i)
            IF (eln>0) THEN
! note that nv (i) has been negated to denote i in Lme:
              nvi = -work(nv+i)
              wnvi = wflg - nvi
              DO p = ip(i), ip(i) + eln - 1
                e = irn(p)
                we = work(w+e)
                IF (we>=wflg) THEN
! unabsorbed element e has been seen in this loop
                  we = we - nvi
                ELSE IF (we/=0) THEN
! e is an unabsorbed element - this is
! the first we have seen e in all of Scan 1
                  we = work(degree+e) + wnvi
                END IF
                work(w+e) = we
              END DO
            END IF
          END DO

! ====================================================================
! DEGREE UPDATE AND ELEMENT ABSORPTION
! ====================================================================

! -------------------------------------------------------------
! Scan 2:  for each sparse i in Lme, sum up the external degrees
! of each Le for the elements e appearing within i, plus the
! supervariables in i.  Place i in hash list.
! -------------------------------------------------------------

          DO pme = pme1, pme2
            i = irn(pme)
! remove absorbed elements from the list for i
            p1 = ip(i)
            p2 = p1 + perm(i) - 1
            pn = p1
            hash = 0
            deg = 0

! -------------------------------------------------------
! scan the element list associated with supervariable i
! -------------------------------------------------------
            DO p = p1, p2
              e = irn(p)
! dext = | Le | - | (Le \cap Lme)\D | - work(denxt+e)
              IF (work(w+e)/=0) THEN
                dext = work(w+e) - wflg
                IF (dext>0) THEN
                  deg = deg + dext
                  irn(pn) = e
                  pn = pn + 1
                  hash = hash + e
                ELSE IF (dext==0) THEN
! aggressive absorption: e is not adjacent to me, but
! |Le(G') \ Lme(G')| is 0, so absorb it into me
                  IF (sep(e)/=mc70_sep_flag) THEN
                    ip(e) = -me
                    work(w+e) = 0
                  ELSE
                    ip(e) = 0
                    work(w+e) = 0
                  END IF
                END IF
              END IF
            END DO

! count the number of elements in i (including me):
! if (i.eq.78) !!! write(*,*)'5699',pn - p1 + 1
            perm(i) = pn - p1 + 1

! ----------------------------------------------------------
! scan the supervariables in the list associated with i
! ----------------------------------------------------------
            p3 = pn
            DO p = p2 + 1, p1 + work(len+i) - 1
              j = irn(p)
              nvj = work(nv+j)
              IF (nvj>0) THEN
! j is unabsorbed, and not in Lme.
! add to degree and add to new list
                deg = deg + nvj
                irn(pn) = j
                pn = pn + 1
                hash = hash + j
              END IF
            END DO

! ----------------------------------------------------------
! update the degree and check for mass elimination
! ----------------------------------------------------------
            IF (deg==0 .AND. sep(i)/=mc70_sep_flag) THEN
! -------------------------------------------------------
! mass elimination - supervariable i can be eliminated
! -------------------------------------------------------
              ip(i) = -me
              nvi = -work(nv+i)
              degme = degme - nvi
              nvpiv = nvpiv + nvi
              nel = nel + nvi
              work(nv+i) = 0
              perm(i) = 0
            ELSE
! -------------------------------------------------------
! update the upper-bound degree of i
! A bound for the new external degree is the old bound plus
! the size of the generated element
! -------------------------------------------------------

! the following degree does not yet include the size
! of the current element, which is added later:
              work(degree+i) = min(deg,work(degree+i))

! -------------------------------------------------------
! add me to the list for i
! -------------------------------------------------------
! move first supervariable to end of list
              irn(pn) = irn(p3)
! move first element to end of element part of list
              irn(p3) = irn(p1)
! add new element to front of list.
              irn(p1) = me
! store the new length of the list in len (i)
              work(len+i) = pn - p1 + 1

! -------------------------------------------------------
! place in hash bucket.  Save hash key of i in last (i).
! -------------------------------------------------------
              hash = abs(mod(hash,hmod)) + 1
              j = work(head+hash)
              IF (j<=0) THEN
! the degree list is empty, hash head is -j
                work(denxt+i) = -j
                work(head+hash) = -i
              ELSE
! degree list is not empty - has j as its head
! last is hash head
                work(denxt+i) = work(last+j)
                work(last+j) = i
              END IF
              work(last+i) = hash
            END IF
          END DO
          work(degree+me) = degme

! -------------------------------------------------------------
! Clear the counter array, w (...), by incrementing wflg.
! -------------------------------------------------------------
          dmax = max(dmax,degme)
          wflg = wflg + dmax

! make sure that wflg+n does not cause integer overflow
          IF (wflg>=iovflo-n) THEN
            DO x = 1, n
              IF (work(w+x)/=0) work(w+x) = 1
            END DO
            wflg = 2
          END IF
! at this point, w (1..n) .lt. wflg holds

! ====================================================================
! SUPERVARIABLE DETECTION
! ====================================================================
          DO pme = pme1, pme2
            i = irn(pme)
            IF ((work(nv+i)<0)) THEN
! replace i by head of its hash bucket, and set the hash
! bucket header to zero

! -------------------------------------------------------
! examine all hash buckets with 2 or more variables.  We
! do this by examing all unique hash keys for super-
! variables in the pattern Lme of the current element, me
! -------------------------------------------------------
              hash = work(last+i)
! let i = head of hash bucket, and empty the hash bucket
              j = work(head+hash)
              IF (j/=0) THEN
                IF (j<0) THEN
! degree list is empty
                  i = -j
                  work(head+hash) = 0
                ELSE
! degree list is not empty, restore last () of head
                  i = work(last+j)
                  work(last+j) = 0
                END IF
                IF (i/=0) THEN

! while loop:
                  DO jdummy = 1, n
                    IF (work(denxt+i)==0) THEN
                      GO TO 80
                    ELSE
! ----------------------------------------------------
! this bucket has one or more variables following i.
! scan all of them to see if i can absorb any entries
! that follow i in hash bucket.  Scatter i into w.
! ----------------------------------------------------
                      ln = work(len+i)
                      eln = perm(i)
! do not flag the first element in the list (me)
                      DO p = ip(i) + 1, ip(i) + ln - 1
                        work(w+irn(p)) = wflg
                      END DO

! ----------------------------------------------------
! scan every other entry j following i in bucket
! ----------------------------------------------------
                      jlast = i
                      j = work(denxt+i)

! while loop:
                      DO idummy = 1, n
                        IF (j==0) THEN
                          GO TO 70
                        ELSE

! -------------------------------------------------
! check if j and i have identical nonzero pattern
! -------------------------------------------------
! jump if i and j do not have same size data structure
! jump if i and j do not have same number adj elts
                          IF (work(len+j)==ln .AND. perm(j)==eln .AND. &
                              sep(i)==sep(j)) THEN
! do not flag the first element in the list (me)

                            DO p = ip(j) + 1, ip(j) + ln - 1
! jump if an entry (irn(p)) is in j but not in i
                              IF (work(w+irn(p))/=wflg) GO TO 50
                            END DO

! -------------------------------------------------
! found it!  j can be absorbed into i
! -------------------------------------------------
                            ip(j) = -i
! both nv (i) and nv (j) are negated since they
! are in Lme, and the absolute values of each
! are the number of variables in i and j:
                            work(nv+i) = work(nv+i) + work(nv+j)
                            work(nv+j) = 0
                            perm(j) = 0
! delete j from hash bucket
                            j = work(denxt+j)
                            work(denxt+jlast) = j
                            GO TO 60
                          END IF

! -------------------------------------------------
50                        CONTINUE
! j cannot be absorbed into i
! -------------------------------------------------
                          jlast = j
                          j = work(denxt+j)
                        END IF
60                      CONTINUE
                      END DO

! ----------------------------------------------------
! no more variables can be absorbed into i
! go to next i in bucket and clear flag array
! ----------------------------------------------------
70                    wflg = wflg + 1
                      i = work(denxt+i)
                      IF (i==0) GO TO 80
                    END IF
                  END DO
                END IF
              END IF
            END IF
80          CONTINUE
          END DO

! ====================================================================
! RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
! Squeeze out absorbed variables
! ====================================================================
          p = pme1
          nleft = n - nel
          DO pme = pme1, pme2
            i = irn(pme)
            nvi = -work(nv+i)
            IF (nvi>0) THEN
! i is a principal variable in Lme
! restore nv (i) to signify that i is principal
              work(nv+i) = nvi
! -------------------------------------------------------
! compute the external degree (add size of current elem)
! -------------------------------------------------------
              deg = min(work(degree+i)+degme-nvi,nleft-nvi)
              work(degree+i) = deg
! -------------------------------------------------------
! place the supervariable at the head of the degree list
! -------------------------------------------------------
              IF (sep(i)/=mc70_sep_flag) THEN
                inext = work(head+deg)
                IF (inext/=0) work(last+inext) = i
                work(denxt+i) = inext
                work(last+i) = 0
                work(head+deg) = i
! -------------------------------------------------------
! save the new degree, and find the minimum degree
! -------------------------------------------------------
                mindeg = min(mindeg,deg)
              END IF
! -------------------------------------------------------
! place the supervariable in the element pattern
! -------------------------------------------------------
              irn(p) = i
              p = p + 1
            END IF
          END DO

! =====================================================================
! FINALIZE THE NEW ELEMENT
! =====================================================================
          work(nv+me) = nvpiv + degme
! nv (me) is now the degree of pivot (including diagonal part)
! save the length of the list for the new element me
          work(len+me) = p - pme1
          IF (work(len+me)==0) THEN
! there is nothing left of the current pivot element
            ip(me) = 0
            work(w+me) = 0
          END IF
          IF (newmem/=0) THEN
! element was not constructed in place: deallocate part
! of it (final size is less than or equal to newmem,
! since newly nonprincipal variables have been removed).
            pfree = p
            mem = mem - newmem + work(len+me)
          END IF

! =====================================================================
! END WHILE (selecting pivots)

        END DO
! ===================================================================
! COMPUTE THE PERMUTATION VECTORS
! ===================================================================

! ----------------------------------------------------------------
! The time taken by the following code is O(n).  At this
! point, perm (e) = -k has been done for all elements e,
! and perm (i) = 0 has been done for all nonprincipal
! variables i.  At this point, there are no principal
! supervariables left, and all elements are absorbed.
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! compute the ordering of unordered nonprincipal variables
! ----------------------------------------------------------------
        l = n
        DO i = 1, n
          IF (perm(i)==0 .AND. sep(i)/=mc70_sep_flag) THEN
! ----------------------------------------------------------
! i is an un-ordered row.  Traverse the tree from i until
! reaching an element, e.  The element, e, was the
! principal supervariable of i and all nodes in the path
! from i to when e was selected as pivot.
! ----------------------------------------------------------
            j = -ip(i)
! while (j is a variable) do:
            DO jdummy = 1, n
              IF (perm(j)<0) THEN
                GO TO 140
              ELSE
                j = -ip(j)
              END IF
           ! if (n.lt.51) write(*,*) 'jj',j,i
            END DO
140         e = j
! ----------------------------------------------------------
! get the current pivot ordering of e
! ----------------------------------------------------------
            k = -perm(e)

! ----------------------------------------------------------
! traverse the path again from i to e, and compress the
! path (all nodes point to e).  Path compression allows
! this code to compute in O(n) time.  Order the unordered
! nodes in the path, and place the element e at the end.
! ----------------------------------------------------------
            j = i
! while (j is a variable) do:
            DO idummy = 1, n
              IF (perm(j)<0) THEN
                GO TO 150
              ELSE
                jnext = -ip(j)
                ip(j) = -e
                IF (perm(j)==0) THEN
! j is an unordered row
                  perm(j) = k
                  k = k + 1
                END IF
                j = jnext
              END IF
            END DO
! leave perm (e) negative, so we know it is an element
150         perm(e) = -k
          ELSE
            IF (sep(i)==mc70_sep_flag) THEN
              perm(i) = l
              l = l - 1
            END IF
          END IF
        END DO

! ----------------------------------------------------------------
! reset the permutation (perm (1..n)) to be positive and ignore the halo
! ----------------------------------------------------------------
        DO i = 1, n
          k = abs(perm(i))
          ip(k) = i
        END DO
        perm(1:n) = ip(1:n)

! ====================================================================
! RETURN THE MEMORY USAGE IN irn AND SET INFORMATION ARRAYS
! ====================================================================
! If maxmem is less than or equal to lirn, then no compressions
! occurred, and irn (maxmem+1 ... lirn) was unused.  Otherwise
! compressions did occur, and lirn would have had to have been
! greater than or equal to maxmem for no compressions to occur.
! Return the value of maxmem in the pfree argument.


        pfree = maxmem
      END SUBROUTINE hamd

      
     SUBROUTINE hamd_both_old(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition,iperm,&
         a_weight,work)
       INTEGER, INTENT(IN) :: a_n ! order of matrix being partitioned
       INTEGER, INTENT(IN) :: a_ne ! no. entries in matrix being partitioned
       INTEGER, INTENT(IN) :: a_ptr(a_n) ! col ptrs for matrix being partitioned
       INTEGER, INTENT(IN) :: a_row(a_ne) ! row indices for matrix
         ! being partitioned. 
       INTEGER, INTENT(IN) :: a_n1 ! no. rows in partition 1
       INTEGER, INTENT(IN) :: a_n2 ! no. rows in partition 2
       INTEGER, INTENT(IN) :: partition(a_n) ! the partitions
       INTEGER, INTENT(INOUT) :: iperm(a_n) ! maps current permuation to the 
         ! column indices of the matrix whose ordering is being computed
       INTEGER, INTENT(INOUT) :: a_weight(a_n) ! weights of vertices
       INTEGER, INTENT(OUT) :: work(11*a_n + a_ne)


       ! Local variables
       INTEGER :: i,j
       INTEGER :: extract_work ! pointers into work array for mask arrays
       INTEGER :: hamd_perm ! pointer into work array for perm array
       INTEGER :: hamd_work ! pointer into work array for hamd work array
       INTEGER :: a_ptr_sub ! pointer into work for col ptrs of submatrix
       INTEGER :: a_irn_sub ! pointer into work for irn array of submatrix
       INTEGER :: rows_sub ! pointer into work for rows_sub array
       INTEGER :: a_lirn_sub ! length of irn array of submatrix
       INTEGER :: a_n_1 ! order of submatrix 1
       INTEGER :: a_n_2 ! order of submatrix 2
       INTEGER :: a_n_sep ! number entries in separator
       INTEGER :: a_ne_sub ! number entries in submatrix
       INTEGER :: len_a_row_sub ! used when extracting submatrices


       ! Set orders of submatrices
       a_n_1 = a_n - a_n2
       a_n_2 = a_n - a_n1
       a_n_sep = a_n - a_n1 - a_n2
       

       ! Set pointers into work array
       hamd_perm = 0 ! length a_n
       a_ptr_sub = hamd_perm + a_n ! length a_n
       a_lirn_sub = a_ne + max(a_n_1,a_n_2) + 1
       ! max(a_n_1,a_n_2) + 1 .le. a_n
       a_irn_sub = a_ptr_sub + a_n ! length a_lirn_sub
       rows_sub = a_irn_sub + a_lirn_sub ! length a_n
       extract_work = rows_sub + a_n ! length a_n
       hamd_work = rows_sub + a_n ! length 7*a_n


       ! Form submatrix 1
     !  IF (a_n1 .NE. 1) THEN
       work(rows_sub+1:rows_sub+a_n1) = partition(1:a_n1)
       work(rows_sub+a_n1+1:rows_sub+a_n_1) = partition(a_n1+a_n2+1:a_n)
       len_a_row_sub = a_ne
       call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n1,a_n_sep,&
        work(rows_sub+1:rows_sub+a_n_1),&
        a_ne_sub,work(a_ptr_sub+1:a_ptr_sub+a_n_1),len_a_row_sub,&
        work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+a_n))

       ! Apply hamd
       work(rows_sub+1:rows_sub+a_n1) = mc70_part1_flag
       work(rows_sub+a_n1+1:rows_sub+a_n_1) = mc70_sep_flag
       call hamd(a_n_1,a_ne_sub,a_lirn_sub,&
         work(a_irn_sub+1:a_irn_sub+a_lirn_sub),&
         work(a_ptr_sub+1:a_ptr_sub+a_n_1),work(rows_sub+1:rows_sub+a_n_1),&
         work(hamd_perm+1:hamd_perm+a_n_1),work(hamd_work+1:hamd_work+7*a_n_1))
      ! ELSE
      !   work(hamd_perm+1) = 1
      ! END IF

       ! Overwrite first a_n1 entries of hamd_perm with first a_n1 entries 
       ! that will form new iperm. Similarly, overwrite first a_n1 entries of 
       !rows_sub with first a_n1 entries that will form new a_weight        
       ! no longer need info in a_ptr
       DO i = 1,a_n1
          j = work(hamd_perm+i)
            work(a_ptr_sub+i) = iperm(partition(j))
            work(rows_sub+i) = a_weight(partition(j))
       END DO
       work(hamd_perm+1:hamd_perm+a_n1) = work(a_ptr_sub+1:a_ptr_sub+a_n1)

       ! Build second submatrix
       !IF (a_n2 .NE. 1) THEN
       work(rows_sub+a_n1+1:rows_sub+a_n1+a_n_2) = partition(a_n1+1:a_n)
       len_a_row_sub = a_ne
       call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n2,a_n_sep,&
        work(rows_sub+a_n1+1:rows_sub+a_n1+a_n_2),&
        a_ne_sub,work(a_ptr_sub+1:a_ptr_sub+a_n_2),len_a_row_sub,&
        work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+a_n))
       ! Apply hamd
       work(rows_sub+a_n1+1:rows_sub+a_n1+a_n2) = mc70_part1_flag
       work(rows_sub+a_n1+a_n2+1:rows_sub+a_n1+a_n_2) = mc70_sep_flag
       call hamd(a_n_2,a_ne_sub,a_lirn_sub,&
         work(a_irn_sub+1:a_irn_sub+a_lirn_sub),&
         work(a_ptr_sub+1:a_ptr_sub+a_n_2),&
         work(rows_sub+a_n1+1:rows_sub+a_n1+a_n_2),&
         work(hamd_perm+a_n1+1:hamd_perm+a_n),&
         work(hamd_work+1:hamd_work+7*a_n_2))
      ! ELSE
      !   work(hamd_perm+a_n1+1) = 1
      ! END IF
       DO i = 1,a_n_2
          j = work(hamd_perm+a_n1+i)
            work(a_ptr_sub+i+a_n1) = iperm(partition(a_n1+j))
            work(rows_sub+i+a_n1) = a_weight(partition(a_n1+j))
       END DO
       iperm(a_n1+1:a_n) = work(a_ptr_sub+1+a_n1:a_ptr_sub+a_n)
       iperm(1:a_n1) = work(hamd_perm+1:hamd_perm+a_n1)
       a_weight(1:a_n) = work(rows_sub+1:rows_sub+a_n)

     END SUBROUTINE hamd_both_old

    SUBROUTINE hamd_both_old1(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition,iperm,&
         a_weight,work)
       INTEGER, INTENT(IN) :: a_n ! order of matrix being partitioned
       INTEGER, INTENT(IN) :: a_ne ! no. entries in matrix being partitioned
       INTEGER, INTENT(IN) :: a_ptr(a_n) ! col ptrs for matrix being partitioned
       INTEGER, INTENT(IN) :: a_row(a_ne) ! row indices for matrix
         ! being partitioned. 
       INTEGER, INTENT(IN) :: a_n1 ! no. rows in partition 1
       INTEGER, INTENT(IN) :: a_n2 ! no. rows in partition 2
       INTEGER, INTENT(IN) :: partition(a_n) ! the partitions
       INTEGER, INTENT(INOUT) :: iperm(a_n) ! maps current permutation to the 
         ! column indices of the matrix whose ordering is being computed
       INTEGER, INTENT(INOUT) :: a_weight(a_n) ! weights of vertices
       INTEGER, INTENT(OUT) :: work(11*a_n + a_ne)


       ! Local variables
       INTEGER :: i,j
       INTEGER :: extract_work ! pointers into work array for mask arrays
       INTEGER :: hamd_perm ! pointer into work array for perm array
       INTEGER :: hamd_work ! pointer into work array for hamd work array
       INTEGER :: a_ptr_sub ! pointer into work for col ptrs of submatrix
       INTEGER :: a_irn_sub ! pointer into work for irn array of submatrix
       INTEGER :: rows_sub ! pointer into work for rows_sub array
       INTEGER :: a_lirn_sub ! length of irn array of submatrix
       INTEGER :: a_n_1 ! order of submatrix 1
       INTEGER :: a_n_2 ! order of submatrix 2
       INTEGER :: a_n_sep ! number entries in separator
       INTEGER :: a_ne_sub ! number entries in submatrix
       INTEGER :: len_a_row_sub ! used when extracting submatrices


       ! Set orders of submatrices
       a_n_1 = a_n1
       a_n_2 = a_n2
       a_n_sep = 0
       

       ! Set pointers into work array
       hamd_perm = 0 ! length a_n
       a_ptr_sub = hamd_perm + a_n ! length a_n
       a_lirn_sub = a_ne + max(a_n_1,a_n_2) + 1
       ! max(a_n_1,a_n_2) + 1 .le. a_n
       a_irn_sub = a_ptr_sub + a_n ! length a_lirn_sub
       rows_sub = a_irn_sub + a_lirn_sub ! length a_n
       extract_work = rows_sub + a_n ! length a_n
       hamd_work = rows_sub + a_n ! length 7*a_n


       ! Form submatrix 1
     !  IF (a_n1 .NE. 1) THEN
       work(rows_sub+1:rows_sub+a_n1) = partition(1:a_n1)
       len_a_row_sub = a_ne
       call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n1,a_n_sep,&
        work(rows_sub+1:rows_sub+a_n_1),&
        a_ne_sub,work(a_ptr_sub+1:a_ptr_sub+a_n_1),len_a_row_sub,&
        work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+a_n))

       ! Apply hamd
       work(rows_sub+1:rows_sub+a_n1) = mc70_part1_flag
       call hamd(a_n_1,a_ne_sub,a_lirn_sub,&
         work(a_irn_sub+1:a_irn_sub+a_lirn_sub),&
         work(a_ptr_sub+1:a_ptr_sub+a_n_1),work(rows_sub+1:rows_sub+a_n_1),&
         work(hamd_perm+1:hamd_perm+a_n_1),work(hamd_work+1:hamd_work+7*a_n_1))

       ! Overwrite first a_n1 entries of hamd_perm with first a_n1 entries 
       ! that will form new iperm. Similarly, overwrite first a_n1 entries of 
       !rows_sub with first a_n1 entries that will form new a_weight        
       ! no longer need info in a_ptr
       DO i = 1,a_n1
          j = work(hamd_perm+i)
            work(a_ptr_sub+i) = iperm(partition(j))
            work(rows_sub+i) = a_weight(partition(j))
       END DO
       work(hamd_perm+1:hamd_perm+a_n1) = work(a_ptr_sub+1:a_ptr_sub+a_n1)

       ! Build second submatrix
       !IF (a_n2 .NE. 1) THEN
       work(rows_sub+a_n1+1:rows_sub+a_n1+a_n_2) = partition(a_n1+1:a_n1+a_n2)
       len_a_row_sub = a_ne
       call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n2,a_n_sep,&
        work(rows_sub+a_n1+1:rows_sub+a_n1+a_n_2),&
        a_ne_sub,work(a_ptr_sub+1:a_ptr_sub+a_n_2),len_a_row_sub,&
        work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+a_n))
       ! Apply hamd
       work(rows_sub+a_n1+1:rows_sub+a_n1+a_n2) = mc70_part1_flag
       call hamd(a_n_2,a_ne_sub,a_lirn_sub,&
         work(a_irn_sub+1:a_irn_sub+a_lirn_sub),&
         work(a_ptr_sub+1:a_ptr_sub+a_n_2),&
         work(rows_sub+a_n1+1:rows_sub+a_n1+a_n_2),&
         work(hamd_perm+a_n1+1:hamd_perm+a_n),&
         work(hamd_work+1:hamd_work+7*a_n_2))
      ! ELSE
      !   work(hamd_perm+a_n1+1) = 1
      ! END IF
       DO i = 1,a_n_2
          j = work(hamd_perm+a_n1+i)
            work(a_ptr_sub+i+a_n1) = iperm(partition(a_n1+j))
            work(rows_sub+i+a_n1) = a_weight(partition(a_n1+j))
       END DO
       DO i = a_n_2+1,a_n_2 + (a_n-a_n1-a_n2)
            j = i
            work(a_ptr_sub+i+a_n1) = iperm(partition(a_n1+j))
            work(rows_sub+i+a_n1) = a_weight(partition(a_n1+j))

       END DO
       iperm(a_n1+1:a_n) = work(a_ptr_sub+1+a_n1:a_ptr_sub+a_n)
       iperm(1:a_n1) = work(hamd_perm+1:hamd_perm+a_n1)
       a_weight(1:a_n) = work(rows_sub+1:rows_sub+a_n)

     END SUBROUTINE hamd_both_old1


     SUBROUTINE extract_matrix(a_n,a_ne,a_ptr,a_row,a_n_part,a_n_sep,rows_sub,&
        a_ne_sub,a_ptr_sub,len_a_row_sub,a_row_sub,work)
       INTEGER, INTENT(IN) :: a_n ! order of matrix being partitioned
       INTEGER, INTENT(IN) :: a_ne ! no. entries in matrix being partitioned
       INTEGER, INTENT(IN) :: a_ptr(a_n) ! col ptrs for matrix being partitioned
       INTEGER, INTENT(IN) :: a_row(a_ne) ! row indices for matrix
         ! being partitioned. 
       INTEGER, INTENT(IN) :: a_n_part ! no. rows in partition
       INTEGER, INTENT(IN) :: a_n_sep ! no. rows in partition
       INTEGER, INTENT(IN) :: rows_sub(a_n_part+a_n_sep) ! rows/cols of matrix 
         ! to be extracted. Intersecting rows/cols of separator will be replaced
         ! by matrix of all zeros
       INTEGER, INTENT(OUT) :: a_ne_sub ! no. entries stored in extracted matrix
       INTEGER, INTENT(OUT) :: a_ptr_sub(a_n_part+a_n_sep) ! col ptrs for
         ! extracted matrix
       INTEGER, INTENT(IN) :: len_a_row_sub ! length of a_row_sub
       INTEGER, INTENT(OUT) :: a_row_sub(len_a_row_sub) ! row indices for
         ! extracted matrix
       INTEGER, INTENT(OUT) :: work(a_n)

       ! Local variables
       INTEGER :: i,j,k,l,m,p
       INTEGER :: a_n_sub ! Order of extracted matrix
       INTEGER :: mask! pointer into work array for mask arrays

       ! Set pointers into work array
       mask = 0 ! length a_n
       
       ! Set mask
       a_n_sub = a_n_part+a_n_sep
       work(mask+1:mask+a_n) = 0
       DO i = 1,a_n_part
         j = rows_sub(i)
         work(mask+j) = i
       END DO
       DO i = a_n_part+1,a_n_sub
         j = rows_sub(i)
         work(mask+j) = -i
       END DO
       a_row_sub(:)=0
 
       ! Count number of entries in  submatrix and set-up column ptrs
       a_ptr_sub(1:a_n_sub) = 0
       DO j = 1,a_n_part
         a_ptr_sub(j) = 0
         i = rows_sub(j)
         IF (i .EQ. a_n) THEN
           l = a_ne
         ELSE
           l = a_ptr(i+1) - 1
         END IF
         DO k = a_ptr(i),l
           m = a_row(k)
           IF (work(mask+m) .NE. 0) THEN
              a_ptr_sub(j) = a_ptr_sub(j) + 1

           END IF
         END DO
       END DO
       DO j = a_n_part+1,a_n_sub
         a_ptr_sub(j) = 0
         i = rows_sub(j)
         IF (i .EQ. a_n) THEN
           l = a_ne
         ELSE
           l = a_ptr(i+1) - 1
         END IF
         DO k = a_ptr(i),l
           m = a_row(k)
           IF (work(mask+m) .GT. 0) a_ptr_sub(j) = a_ptr_sub(j) + 1
         END DO
       END DO
       a_ptr_sub(1) = a_ptr_sub(1) + 1
       DO j = 2,a_n_sub
         a_ptr_sub(j) = a_ptr_sub(j) + a_ptr_sub(j-1)
       END DO
       a_ne_sub = a_ptr_sub(a_n_sub)-1

       ! Form a_row_sub
       DO j = 1,a_n_part
         i = rows_sub(j)
         IF (i .EQ. a_n) THEN
           l = a_ne
         ELSE
           l = a_ptr(i+1) - 1
         END IF
         DO k = a_ptr(i),l
           m = a_row(k)
           IF (work(mask+m) .NE. 0) THEN
            a_ptr_sub(j) = a_ptr_sub(j) - 1
            p = a_ptr_sub(j)
            a_row_sub(p) = abs(work(mask+m))
           END IF
         END DO
       END DO
       DO j = a_n_part+1,a_n_sub
         i = rows_sub(j)
         IF (i .EQ. a_n) THEN
           l = a_ne
         ELSE
           l = a_ptr(i+1) - 1
         END IF
         DO k = a_ptr(i),l
           m = a_row(k)
           IF (work(mask+m) .GT. 0) THEN
            a_ptr_sub(j) = a_ptr_sub(j) - 1
            p = a_ptr_sub(j)
            a_row_sub(p) = work(mask+m)
           END IF
         END DO
       END DO
     END SUBROUTINE extract_matrix

     SUBROUTINE extract_both_matrices(a_n,a_ne,a_ptr,a_row,a_n_part1,a_n_part2,&
        rows_sub,a_ne_sub1,a_ne_sub2,a_ptr_sub,len_a_row_sub,a_row_sub,work)
       INTEGER, INTENT(IN) :: a_n ! order of matrix being partitioned
       INTEGER, INTENT(IN) :: a_ne ! no. entries in matrix being partitioned
       INTEGER, INTENT(IN) :: a_ptr(a_n) ! col ptrs for matrix being partitioned
       INTEGER, INTENT(IN) :: a_row(a_ne) ! row indices for matrix
         ! being partitioned. 
       INTEGER, INTENT(IN) :: a_n_part1 ! no. rows in partition 1
       INTEGER, INTENT(IN) :: a_n_part2 ! no. rows in partition 2
       INTEGER, INTENT(IN) :: rows_sub(a_n_part1+a_n_part2) ! rows/cols of  
         ! matrices to be extracted. First a_n_part1 entries contain the 
         ! rows/cols forming the first matrix to be extracted.
       INTEGER, INTENT(OUT) :: a_ne_sub1 ! no. entries in extracted matrix 1
       INTEGER, INTENT(OUT) :: a_ne_sub2 ! no. entries in extracted matrix 2
       INTEGER, INTENT(OUT) :: a_ptr_sub(a_n_part1+a_n_part2) ! col ptrs for
         ! extracted matrices. First a_n_part1 are for first matrix, etc
       INTEGER, INTENT(IN) :: len_a_row_sub ! length of a_row_sub
       INTEGER, INTENT(OUT) :: a_row_sub(len_a_row_sub) ! row indices for
         ! extracted matrices. First a_ne_part1 entries are for first matrix; 
         ! the immediately following a_ne_part2 entries are for second matrix
       INTEGER, INTENT(OUT) :: work(a_n)

       ! Local variables
       INTEGER :: i,j,k,l,m,p
       INTEGER :: mask! pointer into work array for mask arrays


       ! Set pointers into work array
       mask = 0 ! length a_n
       
       ! Set mask
       work(mask+1:mask+a_n) = 0
       DO i = 1,a_n_part1
         j = rows_sub(i)
         work(mask+j) = i
       END DO
       DO i = a_n_part1+1,a_n_part1+a_n_part2
         j = rows_sub(i)
         work(mask+j) = -i+a_n_part1
       END DO
       a_row_sub(:)=0
 
       ! Count number of entries in each submatrix and set-up column ptrs
       a_ptr_sub(1:a_n_part1+a_n_part2) = 0
       DO j = 1,a_n_part1
         a_ptr_sub(j) = 0
         i = rows_sub(j)
         IF (i .EQ. a_n) THEN
           l = a_ne
         ELSE
           l = a_ptr(i+1) - 1
         END IF
         DO k = a_ptr(i),l
           m = a_row(k)
           IF (work(mask+m) .GT. 0) THEN
              a_ptr_sub(j) = a_ptr_sub(j) + 1

           END IF
         END DO
       END DO

       DO j = a_n_part1+1,a_n_part1+a_n_part2
         a_ptr_sub(j) = 0
         i = rows_sub(j)
         IF (i .EQ. a_n) THEN
           l = a_ne
         ELSE
           l = a_ptr(i+1) - 1
         END IF
         DO k = a_ptr(i),l
           m = a_row(k)
           IF (work(mask+m) .LT. 0) a_ptr_sub(j) = a_ptr_sub(j) + 1
         END DO
       END DO

       a_ptr_sub(1) = a_ptr_sub(1) + 1
       DO j = 2,a_n_part1
         a_ptr_sub(j) = a_ptr_sub(j) + a_ptr_sub(j-1)
       END DO
       a_ne_sub1 = a_ptr_sub(a_n_part1)-1

       a_ptr_sub(a_n_part1+1) = a_ptr_sub(a_n_part1+1) + 1
       DO j = a_n_part1+2,a_n_part1+a_n_part2
         a_ptr_sub(j) = a_ptr_sub(j) + a_ptr_sub(j-1)
       END DO
       a_ne_sub2 = a_ptr_sub(a_n_part1+a_n_part2)-1

       ! Form a_row_sub
       DO j = 1,a_n_part1
         i = rows_sub(j)
         IF (i .EQ. a_n) THEN
           l = a_ne
         ELSE
           l = a_ptr(i+1) - 1
         END IF
         DO k = a_ptr(i),l
           m = a_row(k)
           IF (work(mask+m) .GT. 0) THEN
            a_ptr_sub(j) = a_ptr_sub(j) - 1
            p = a_ptr_sub(j)
            a_row_sub(p) = abs(work(mask+m))
           END IF
         END DO
       END DO

       DO j = a_n_part1+1,a_n_part1+a_n_part2
         i = rows_sub(j)
         IF (i .EQ. a_n) THEN
           l = a_ne
         ELSE
           l = a_ptr(i+1) - 1
         END IF
         DO k = a_ptr(i),l
           m = a_row(k)
           IF (work(mask+m) .LT. 0) THEN
            a_ptr_sub(j) = a_ptr_sub(j) - 1
            p = a_ptr_sub(j)
            a_row_sub(p+a_ne_sub1) = -work(mask+m)
           END IF
         END DO
       END DO
     END SUBROUTINE extract_both_matrices

     SUBROUTINE hamd_one(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,partition,iperm,&
            a_weight,a_ne_sub,work)
       ! Apply hamd to the smaller partition and extract the other matrix into 
       ! a_ptr and a_row
       INTEGER, INTENT(IN) :: a_n ! order of matrix being partitioned
       INTEGER, INTENT(IN) :: a_ne ! no. entries in matrix being partitioned
       INTEGER, INTENT(INOUT) :: a_ptr(a_n) ! col ptrs for matrix being 
          !partitioned and, on return, for the extracted submatrix
       INTEGER, INTENT(INOUT) :: a_row(a_ne) ! row indices for matrix
         ! being partitioned and, on return, for the extracted submatrix
       INTEGER, INTENT(IN) :: a_n1 ! no. rows in partition 1
       INTEGER, INTENT(IN) :: a_n2 ! no. rows in partition 2
       INTEGER, INTENT(IN) :: partition(a_n) ! the partitions
       INTEGER, INTENT(INOUT) :: iperm(a_n) ! maps current permuation to the 
         ! column indices of the matrix whose ordering is being computed
       INTEGER, INTENT(INOUT) :: a_weight(a_n) ! weights of vertices
       INTEGER, INTENT(OUT) :: a_ne_sub ! number entries in returned submatrix
       INTEGER, INTENT(OUT) :: work(11*a_n + a_ne)

       ! Local variables
       INTEGER :: i,j
       INTEGER :: extract_work ! pointers into work array for mask arrays
       INTEGER :: hamd_perm ! pointer into work array for perm array
       INTEGER :: hamd_work ! pointer into work array for hamd work array
       INTEGER :: a_ptr_sub ! pointer into work for col ptrs of submatrix
       INTEGER :: a_irn_sub ! pointer into work for irn array of submatrix
       INTEGER :: rows_sub ! pointer into work for rows_sub array
       INTEGER :: a_lirn_sub ! length of irn array of submatrix
       INTEGER :: a_n_1 ! order of submatrix 1
       INTEGER :: a_n_2 ! order of submatrix 2
       INTEGER :: a_n_sep ! number entries in separator
       INTEGER :: len_a_row_sub ! used when extracting submatrices

       ! Set orders of submatrices
       IF (a_n1 .LT. a_n2) THEN
        ! Applying hamd to first partition
        a_n_1 = a_n - a_n2
        a_n_2 = a_n2
        a_n_sep = a_n - a_n1 - a_n2
       ELSE
        ! Applying hamd to second partition
        a_n_2 = a_n - a_n1
        a_n_1 = a_n1
        a_n_sep = a_n - a_n1 - a_n2

       END IF
       

       ! Set pointers into work array
       hamd_perm = 0 ! length a_n
       a_ptr_sub = hamd_perm + a_n ! length a_n
       IF (a_n1 .LT. a_n2) THEN
         a_lirn_sub = a_ne + a_n_1 + 1
       ELSE
         a_lirn_sub = a_ne + a_n_2 + 1
       END IF
       ! max(a_n_1,a_n_2) + 1 .le. a_n
       a_irn_sub = a_ptr_sub + a_n ! length a_lirn_sub
       rows_sub = a_irn_sub + a_lirn_sub ! length a_n
       extract_work = rows_sub + a_n ! length a_n
       hamd_work = rows_sub + a_n ! length 7*a_n

       ! Extract matrix that hamd is applied to
       IF (a_n1 .LT. a_n2) THEN

       ! Start by updating iperm and a_weight
       DO i = 1,a_n
         work(hamd_perm+i) = iperm(partition(i))
            work(rows_sub+i) = a_weight(partition(i))
       END DO
       iperm(1:a_n) = work(hamd_perm+1:hamd_perm+a_n)
       a_weight(1:a_n) = work(rows_sub+1:rows_sub+a_n)
       
       ! Form submatrix 1
       work(rows_sub+1:rows_sub+a_n1) = partition(1:a_n1)
       work(rows_sub+a_n1+1:rows_sub+a_n_1) = partition(a_n1+a_n2+1:a_n)
       len_a_row_sub = a_ne
       call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n1,a_n_sep,&
        work(rows_sub+1:rows_sub+a_n_1),&
        a_ne_sub,work(a_ptr_sub+1:a_ptr_sub+a_n_1),len_a_row_sub,&
        work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+a_n))
       

       ! Apply hamd
       work(rows_sub+1:rows_sub+a_n1) = mc70_part1_flag
       work(rows_sub+a_n1+1:rows_sub+a_n_1) = mc70_sep_flag
       call hamd(a_n_1,a_ne_sub,a_lirn_sub,&
         work(a_irn_sub+1:a_irn_sub+a_lirn_sub),&
         work(a_ptr_sub+1:a_ptr_sub+a_n_1),work(rows_sub+1:rows_sub+a_n_1),&
         work(hamd_perm+1:hamd_perm+a_n_1),work(hamd_work+1:hamd_work+7*a_n_1))

       ! Overwrite first a_n1 entries of hamd_perm with first a_n1 entries 
       ! that will form new iperm. Similarly, overwrite first a_n1 entries of 
       !rows_sub with first a_n1 entries that will form new a_weight        
       ! no longer need info in a_ptr   
       ! no longer need info in a_ptr
       DO i = 1,a_n_1
          j = work(hamd_perm+i)
          IF (j.LE.a_n1) THEN
            work(a_ptr_sub+i) = iperm(j)
            work(rows_sub+i) = a_weight(j)
          END IF
       END DO
      ! work(hamd_perm+1:hamd_perm+a_n1) = work(a_ptr_sub+1:a_ptr_sub+a_n1)
       

        iperm(1:a_n1) = work(a_ptr_sub+1:a_ptr_sub+a_n1)
        a_weight(1:a_n1) = work(rows_sub+1:rows_sub+a_n1)

       ! Build second submatrix
       work(rows_sub+1:rows_sub+a_n_2) = partition(a_n1+1:a_n1+a_n2)
       len_a_row_sub = a_ne
       call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n2,0,&
        work(rows_sub+1:rows_sub+a_n_2),&
        a_ne_sub,work(a_ptr_sub+1:a_ptr_sub+a_n_2),len_a_row_sub,&
        work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+a_n))
       
       ! Copy matrix to a_ptr and a_row
        a_ptr(1:a_n_2) = work(a_ptr_sub+1:a_ptr_sub+a_n_2)
        a_row(1:a_ne_sub) = work(a_irn_sub+1:a_irn_sub+a_ne_sub)

       ELSE

        
        ! Form submatrix 2
       work(rows_sub+1:rows_sub+a_n_2) = partition(a_n1+1:a_n)
       len_a_row_sub = a_ne
       call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n2,a_n_sep,&
        work(rows_sub+1:rows_sub+a_n_2),&
        a_ne_sub,work(a_ptr_sub+1:a_ptr_sub+a_n_2),len_a_row_sub,&
        work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+a_n))
       
       ! Apply hamd
       work(rows_sub+1:rows_sub+a_n2) = mc70_part2_flag
       work(rows_sub+a_n2+1:rows_sub+a_n_2) = mc70_sep_flag
       call hamd(a_n_2,a_ne_sub,a_lirn_sub,&
         work(a_irn_sub+1:a_irn_sub+a_lirn_sub),&
         work(a_ptr_sub+1:a_ptr_sub+a_n_2),work(rows_sub+1:rows_sub+a_n_2),&
         work(hamd_perm+a_n1+1:hamd_perm+a_n),&
         work(hamd_work+1:hamd_work+7*a_n_2))

         DO i = 1,a_n_2
            j = work(hamd_perm+a_n1+i)
            work(a_ptr_sub+i+a_n1) = iperm(partition(a_n1+j))
            work(rows_sub+i+a_n1) = a_weight(partition(a_n1+j))
         END DO

        DO i = 1,a_n_1
         work(a_ptr_sub+i) = iperm(partition(i))
         work(rows_sub+i) = a_weight(partition(i))
        END DO
        iperm(1:a_n) = work(a_ptr_sub+1:a_ptr_sub+a_n)
        a_weight(1:a_n) = work(rows_sub+1:rows_sub+a_n)

       ! Form submatrix 1
       work(rows_sub+1:rows_sub+a_n1) = partition(1:a_n1)
       len_a_row_sub = a_ne
       call extract_matrix(a_n,a_ne,a_ptr,a_row,a_n1,0,&
        work(rows_sub+1:rows_sub+a_n_1),&
        a_ne_sub,work(a_ptr_sub+1:a_ptr_sub+a_n_1),len_a_row_sub,&
        work(a_irn_sub+1:a_irn_sub+a_ne),work(extract_work+1:extract_work+a_n))

       ! Copy matrix to a_ptr and a_row
        a_ptr(1:a_n_1) = work(a_ptr_sub+1:a_ptr_sub+a_n_1)
        a_row(1:a_ne_sub) = work(a_irn_sub+1:a_irn_sub+a_ne_sub)

       END IF

     END SUBROUTINE hamd_one

     SUBROUTINE multilevel_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,&
          partition,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep, &
          control,info1,lwork,work,ml_max_levels)

       INTEGER, INTENT(IN) :: a_n ! order of matrix being partitioned
       INTEGER, INTENT(IN) :: a_ne ! no. entries in matrix being partitioned
       INTEGER, INTENT(IN) :: a_ptr(a_n) ! col ptrs for matrix being 
          !partitioned and, on return, for the extracted submatrix
       INTEGER, INTENT(IN) :: a_row(a_ne) ! row indices for matrix
       INTEGER, INTENT (IN) :: a_weight(a_n) ! weights associated with rows
! of matrix (useful if matrix has already been compressed)
       INTEGER, INTENT(IN) :: sumweight ! sum of entries in a_weight
       INTEGER, INTENT (OUT) :: partition(a_n) ! computed partition
       INTEGER, INTENT (OUT) :: a_n1 ! number of entries in partition 1
       INTEGER, INTENT (OUT) :: a_n2 ! number of entries in partition 2
       INTEGER, INTENT(OUT) :: a_weight_1,a_weight_2,a_weight_sep ! Weighted 
             ! size of partitions and separator
        TYPE (mc70_control), INTENT (IN) :: control
        INTEGER, INTENT (IN) :: lwork ! length of work array: must be atleast
           ! 9a_n + sumweight + a_ne/2
        INTEGER, INTENT (OUT) :: work(lwork) ! work array
        INTEGER, INTENT(INOUT) :: info1
        INTEGER, INTENT(IN) :: ml_max_levels ! no. levels in the multilevel grid (default 10)

        TYPE (mc70_multigrid) :: grid ! the multilevel of graphs (matrices)
        TYPE (zd11_type), TARGET :: matrix
        INTEGER :: i, inv1, inv2, ins
        INTEGER :: mp
        INTEGER :: mglevel_cur ! current level
        INTEGER :: err, print_level ! printing
        LOGICAL :: lerr
        INTEGER :: st ! stat parameter

        info1 = 0
! Set up printing
        IF (control%print_level<0) print_level = 0
! The default is control%print_level = 0
        IF (control%print_level==0) print_level = 1
        IF (control%print_level==1) print_level = 2
        IF (control%print_level>1) print_level = 3
        mp = control%unit_diagnostics
        IF (mp<0) print_level = 0
! Set error controls
        lerr = control%unit_error >= 0 .AND. print_level > 0
        err = control%unit_error


! construct the multigrid at this level
        matrix%n = a_n
        matrix%ne = a_ne
        ALLOCATE (matrix%ptr(a_n+1),matrix%col(a_ne),matrix%val(a_ne),STAT=st)
        IF (st/=0) info1 = mc70_err_memory_alloc
        IF (info1<0) THEN
          IF (lerr) CALL mc70_print_message(info1,err,' multilevel_partition')
          RETURN
        END IF
        matrix%ptr(1:a_n) = a_ptr(1:a_n)
        matrix%ptr(a_n+1) = a_ne+1
        matrix%col(1:a_ne) = a_row(1:a_ne)
        matrix%val(1:a_ne) = 1.0

        grid%size = a_n
        grid%level = 1
        grid%graph => matrix
        NULLIFY (grid%p)

        ALLOCATE (grid%where(a_n),grid%row_wgt(a_n),STAT=st)
        IF (st/=0) info1 = mc70_err_memory_alloc
        IF (info1<0) THEN
          IF (lerr) CALL mc70_print_message(info1,err,' multilevel_partition')
          RETURN
        END IF

! Initialise row weights
        grid%row_wgt(1:a_n) = a_weight(1:a_n)

! initialise mglevel_cur to the maximum number of levels
! allowed for this bisection
        mglevel_cur = ml_max_levels
           
        CALL multilevel(grid,control,sumweight,mglevel_cur,mp,print_level,&
           lwork,work,info1)

        IF (info1/=0) THEN
          IF (lerr) CALL mc70_print_message(info1,err,' multilevel_partition')
          RETURN
        END IF

        inv1 = 1
        inv2 = grid%part_div(1) + 1
        ins = grid%part_div(1) + grid%part_div(2) + 1

        a_weight_1 = 0
        a_weight_2 = 0
        a_weight_sep = 0
        DO i = 1, a_n
          SELECT CASE (grid%where(i))
          CASE (mc70_part1_flag) 
            partition(inv1) = i
            inv1 = inv1 + 1
            a_weight_1 = a_weight_1 + a_weight(i)
          CASE (mc70_part2_flag) 
            partition(inv2) = i
            inv2 = inv2 + 1
            a_weight_2 = a_weight_2 + a_weight(i)
          CASE default
            partition(ins) = i
            ins = ins + 1
            a_weight_sep = a_weight_sep + a_weight(i)
          END SELECT
        END DO

        a_n1 = grid%part_div(1)
        a_n2 = grid%part_div(2)

! deallocate the finest level
        CALL multigrid_deallocate_first(a_n,a_n,grid,info1)
        IF (info1/=0) THEN
          IF (lerr) CALL mc70_print_message(info1,err,' multilevel_partition')
          RETURN
        END IF

        DEALLOCATE (matrix%ptr,matrix%col,STAT=st)
        IF (st/=0) info1 = mc70_err_memory_dealloc
        IF (info1<0) THEN
          IF (lerr) CALL mc70_print_message(info1,err,' multilevel_partition')
          RETURN
        END IF

      END SUBROUTINE multilevel_partition

! ********************************************************

! main subroutine for computing multilevel structure.
! Offers heavy-edge collapsing and maximal independent vertex
! set for coarsening. We will need to test out to see
! which is better.

      RECURSIVE SUBROUTINE multilevel(grid,control,sumweight,&
          mglevel_cur,mp,print_level,lwork,work,info)

        REAL (myreal_mc70), PARAMETER :: half = 0.5_myreal_mc70
        REAL (myreal_mc70), PARAMETER :: one = 1.0_myreal_mc70

! Arguments
        TYPE (mc70_multigrid), INTENT (INOUT), TARGET :: grid ! this level
!  of matrix (grid)
        TYPE (mc70_control), INTENT (IN) :: control
        INTEGER, INTENT (IN) :: sumweight ! sum of weights (unchanged between 
                ! coarse and fine grid
        INTEGER, INTENT(INOUT) :: mglevel_cur ! current grid level
        INTEGER, INTENT (IN) :: mp, print_level ! diagnostic printing
        INTEGER, INTENT (IN) :: lwork ! length of work array
        INTEGER, INTENT(OUT) :: work(lwork) ! work array
        INTEGER, INTENT(INOUT) :: info ! Error flag

! Local variables
        TYPE (mc70_multigrid), POINTER :: cgrid ! the coarse level grid
        INTEGER :: cnvtx ! number of vertices (rows) in the coarse matrix
        TYPE (zd11_type), POINTER :: p ! the coarse grid prolongator
         
        INTEGER, DIMENSION (:), POINTER :: fwhere ! partition on fine grid
        INTEGER, DIMENSION (:), POINTER :: cwhere ! partition on coarse grid
        TYPE (zd11_type), POINTER :: cgraph ! the coarse graph
        TYPE (zd11_type), POINTER :: graph ! the fine graph
        INTEGER, DIMENSION (:), POINTER :: row_wgt ! fine  
! graph vertex weights
        INTEGER, DIMENSION (:), POINTER :: crow_wgt ! coarse 
! graph vertex weights
        REAL (myreal_mc70) :: grid_rdc_fac_min ! min grid reduction factor
        REAL (myreal_mc70) :: grid_rdc_fac_max ! max grid reduction factor
        REAL (myreal_mc70) :: one1
        INTEGER :: ml_switch ! controls when to stop coarsening
        INTEGER :: partition_ptr,work_ptr,a_ne
        INTEGER :: st ! stat parameter
        INTEGER :: i, j,k,l,a_weight_1,a_weight_2,a_weight_sep,ref_method
        INTEGER, ALLOCATABLE :: work1(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        info = 0
        one1 = 1.0

        ml_switch = max(2,control%ml_switch)
        IF (print_level==3) CALL level_print(mp,'size of grid on level ', &
          grid%level,' is ',real(grid%size,myreal_mc70))

        grid_rdc_fac_min = max(0.01_myreal_mc70,control%ml_min_reduction)
! max grid reduction factor must be at least half and at most one
        grid_rdc_fac_max = max(half,control%ml_max_reduction)
        grid_rdc_fac_max = min(one,grid_rdc_fac_max)

! Test to see if this is either the last level or
! if the matrix size too small
!!!write(*,*) 'grid',grid%level,mglevel_cur, grid%size
        IF (grid%level>=mglevel_cur .OR. grid%size<=ml_switch) THEN
          NULLIFY (grid%coarse)
          IF (print_level==3) CALL level_print(mp,'end of level ',grid%level)

! coarsest level in multilevel so compute separator
!!write(*,*) 'start 3520'

         ! CALL coarse_separator(grid%graph%n,grid%graph,grid%row_wgt, &
         !   grid%where,grid%part_div,control,info)
          a_ne = grid%graph%ptr(grid%graph%n+1)-1

          CALL mc70_coarse_partition(grid%graph%n,a_ne,grid%graph%ptr,&
            grid%graph%col,grid%row_wgt,sumweight,&
            grid%part_div(1),grid%part_div(2),grid%where,lwork,work,control,info)

        
          RETURN

        END IF

! Coarsest level not yet reached so carry on coarsening
        CALL coarsen_hec(mp,print_level,grid,info)
        IF (info<0) RETURN

        cgrid => grid%coarse
        cnvtx = cgrid%size
! allocate coarse grid quantities
        ALLOCATE (cgrid%where(cnvtx),cgrid%row_wgt(cnvtx),STAT=st)
        IF (st/=0) THEN
          info = -1
          RETURN
        END IF

! see if the grid reduction is achieved, if not, set the allowed
! maximum level to current level and partition this level
! deallocate the coarse grid quantities that haves been allocated so far
        IF (real(cgrid%size)/real(grid%size)>grid_rdc_fac_max .OR. &
            real(cgrid%size)/real(grid%size)<grid_rdc_fac_min .OR. cgrid%size<4) &
            THEN

          IF (print_level==3) THEN
        !  IF (.true.) THEN
            WRITE (mp,'(a,i10,a,f12.4,i4)') 'at level ', grid%level, &
              ' further coarsening gives reduction factor', &
              cgrid%size/real(grid%size)
            WRITE (mp,'(a,i10)') 'current size = ', grid%size
          END IF

! set current grid level and recurse
          mglevel_cur = grid%level

          CALL multilevel(grid,control,sumweight,mglevel_cur,mp, &
            print_level,lwork,work,info)

          IF (info<0) RETURN

          CALL multigrid_deallocate_last(cgrid,info)
          IF (info<0) RETURN

          RETURN

        END IF

! restriction ================

! form the coarse grid graph and matrix
! cmatrix = P^T*matrix = R*matrix
        p => cgrid%p
        graph => grid%graph
        cgraph => cgrid%graph

! get the coarse matrix
!!write(*,*) 'start 3594'
        CALL galerkin_graph(graph,p,cgraph,info)
!!write(*,*) 'end 3594'
        IF (info<0) RETURN

        IF (cgrid%graph%ptr(cgrid%graph%n+1)-1 .ge. cgrid%graph%n*(cgrid%graph%n-1)) &
            THEN
          IF (print_level==3) THEN
            WRITE (mp,'(a,i10,a,f12.4)') 'at level ', grid%level, &
              ' further coarsening gives reduction factor', &
              cgrid%size/real(grid%size)          
            WRITE (mp,'(a,i10)') 'current size = ', grid%size
            WRITE (mp,'(a,i10)') 'coarse size = ', cgrid%size
          END IF

! set current grid level and recurse
          mglevel_cur = grid%level-1

!!write(*,*) 'start 3568'
          CALL multilevel(grid,control,sumweight,mglevel_cur,mp, &
            print_level,lwork,work,info)
!!write(*,*) 'end 3568'

          IF (info<0) RETURN

          CALL multigrid_deallocate_last(cgrid,info)

          RETURN

        END IF
! row weight cw = R*w
        row_wgt => grid%row_wgt
        crow_wgt => cgrid%row_wgt
        CALL mc65_matrix_multiply_vector(p,row_wgt,crow_wgt,info,trans=.TRUE.)
        CALL multilevel(cgrid,control,sumweight,mglevel_cur,mp,&
           print_level,lwork,work,info)

! prolongation ================

! injection of the order from coarse grid to the
! fine grid, since cwhere(i) is the index of the
! i-th vertex in the new ordering, the order
! of this vertex should be where(i)
! grid%where = P*order_on_coarse_grid
! here P is a special matrix with only one non-zero entry per row
        fwhere => grid%where
        cwhere => cgrid%where
        grid%part_div(1:2) = 0
        CALL mc65_matrix_multiply_vector(p,cwhere,fwhere,info)

        DO i = 1, size(fwhere)
            IF (fwhere(i)==mc70_part1_flag) THEN
              grid%part_div(1) = grid%part_div(1) + 1
            ELSE
              IF (fwhere(i)==mc70_part2_flag) THEN
                grid%part_div(2) = grid%part_div(2) + 1
              END IF
            END IF
        END DO
        a_weight_1 = 0
        a_weight_2 = 0
        a_weight_sep = 0

        ! Set partition
        partition_ptr = 0
        work_ptr = partition_ptr+grid%graph%n
        i = 1
        j = grid%part_div(1)+1
        k = grid%part_div(1)+grid%part_div(2)+1
        DO l = 1,grid%graph%n
          SELECT CASE (grid%where(l))
           CASE (mc70_part1_flag)
            work(partition_ptr+i) = l
            a_weight_1 = a_weight_1 + grid%row_wgt(l)
            i = i+1
           CASE (mc70_part2_flag)
            work(partition_ptr+j) = l
            a_weight_2 = a_weight_2 + grid%row_wgt(l)
            j = j+1
           CASE (mc70_sep_flag)
            work(partition_ptr+k) = l
            a_weight_sep = a_weight_sep + grid%row_wgt(l)
            k = k+1
          END SELECT
        END DO
         a_ne = grid%graph%ptr(grid%graph%n+1)-1
         !  write(*,*) 'a_n1p =',grid%part_div(1),';'
         !  write(*,*) 'a_n2p =',grid%part_div(2),';'
         !  write(*,*) 'a_np =',grid%graph%n,';'
         !  write(*,*) 'a_partp =[',work(partition_ptr+1:partition_ptr+grid%graph%n),'];'
        
       IF (a_weight_sep .GT. 0) THEN
        ! Do not refine if separable graph
        IF (lwork .LT. work_ptr+8*grid%graph%n+sumweight+a_ne/2) THEN
         ALLOCATE (work1(8*grid%graph%n+sumweight+a_ne/2),stat = st)
         IF (st/=0) THEN
          info = mc70_err_memory_alloc
          RETURN
         END IF
         IF (control%refinement .GE. 2) THEN
              IF (min(a_weight_1,a_weight_2) + a_weight_sep .LT. &
                  max(a_weight_1,a_weight_2) ) THEN
                 ref_method = 2
              ELSE
                 ref_method = 1
              END IF
           ELSE  
             IF (control%refinement .LE. 0) THEN
                 ref_method = 1
             ELSE
                 ref_method = 2
             END IF
          END IF


         IF (ref_method .EQ. 1) THEN
            CALL mc70_refine_trim(grid%graph%n,a_ne,grid%graph%ptr,&
             grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1),&
             grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep,&
             work(partition_ptr+1:partition_ptr+grid%graph%n),&
             work1(1:8*grid%graph%n+sumweight+a_ne/2),&
             control) 
            CALL mc70_refine(grid%graph%n,a_ne,grid%graph%ptr,&
             grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1),&
             grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep,&
             work(partition_ptr+1:partition_ptr+grid%graph%n),&
             work1(1:8*grid%graph%n+sumweight+a_ne/2),&
             control,ref_method) 
         ELSE
            CALL mc70_refine(grid%graph%n,a_ne,grid%graph%ptr,&
             grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1),&
             grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep,&
             work(partition_ptr+1:partition_ptr+grid%graph%n),&
             work1(1:8*grid%graph%n+sumweight+a_ne/2),&
             control,ref_method) 
            CALL mc70_refine_trim(grid%graph%n,a_ne,grid%graph%ptr,&
             grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1),&
             grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep,&
             work(partition_ptr+1:partition_ptr+grid%graph%n),&
             work1(1:8*grid%graph%n+sumweight+a_ne/2),&
             control) 
         END IF

         DEALLOCATE (work1,stat = st)
         IF (st/=0) THEN
          info = mc70_err_memory_dealloc
          RETURN
         END IF
        ELSE
        
         IF (control%refinement .GE. 2) THEN
              IF (min(a_weight_1,a_weight_2)+a_weight_sep .LT. &
                  max(a_weight_1,a_weight_2)) THEN
                 ref_method = 2
              ELSE
                 ref_method = 1
              END IF
           ELSE  
             IF (control%refinement .LT. 1) THEN
                 ref_method = 1
             ELSE
                 ref_method = 2
             END IF
          END IF
         IF (ref_method .LE. 1) THEN
            CALL mc70_refine_trim(grid%graph%n,a_ne,grid%graph%ptr,&
             grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1),&
             grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep,&
             work(partition_ptr+1:partition_ptr+grid%graph%n),&
             work(work_ptr+1:work_ptr+8*grid%graph%n+sumweight+a_ne/2),&
             control) 
            CALL mc70_refine(grid%graph%n,a_ne,grid%graph%ptr,&
             grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1),&
             grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep,&
             work(partition_ptr+1:partition_ptr+grid%graph%n),&
             work(work_ptr+1:work_ptr+8*grid%graph%n+sumweight+a_ne/2),&
             control,ref_method) 
         ELSE
            CALL mc70_refine(grid%graph%n,a_ne,grid%graph%ptr,&
             grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1),&
             grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep,&
             work(partition_ptr+1:partition_ptr+grid%graph%n),&
             work(work_ptr+1:work_ptr+8*grid%graph%n+sumweight+a_ne/2),&
             control,ref_method) 
            CALL mc70_refine_trim(grid%graph%n,a_ne,grid%graph%ptr,&
             grid%graph%col,grid%row_wgt,sumweight,grid%part_div(1),&
             grid%part_div(2),a_weight_1,a_weight_2,a_weight_sep,&
             work(partition_ptr+1:partition_ptr+grid%graph%n),&
             work(work_ptr+1:work_ptr+8*grid%graph%n+sumweight+a_ne/2),&
             control) 
         END IF
        END IF
       END IF

        !   write(*,*) 'a_n1r =',grid%part_div(1),';'
        !   write(*,*) 'a_n2r =',grid%part_div(2),';'
        !   write(*,*) 'a_nr =',grid%graph%n,';'
        !   write(*,*) 'a_partr =[',work(partition_ptr+1:partition_ptr+grid%graph%n),'];'

        DO i = 1,grid%part_div(1)
          j = work(partition_ptr+i)
          grid%where(j) = mc70_part1_flag
        END DO
        DO i = grid%part_div(1)+1,grid%part_div(1)+grid%part_div(2)
          j = work(partition_ptr+i)
          grid%where(j) = mc70_part2_flag
        END DO
        DO i = grid%part_div(1)+grid%part_div(2)+1,grid%graph%n
          j = work(partition_ptr+i)
          grid%where(j) = mc70_sep_flag
        END DO

        IF (info<0) RETURN

        IF (print_level==3) CALL level_print(mp,' after post smoothing ', &
          grid%level)

! deallocate the previous level
        CALL multigrid_deallocate(cgrid,info)

      END SUBROUTINE multilevel

!***************************************************************
! ---------------------------------------------------
! mc70_partition matrix
! ---------------------------------------------------
! Partition the matrix and if one (or more) of the generated submatrices is 
! small enough, apply halo amd

      SUBROUTINE mc70_coarse_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,&
         a_n2,where1,lwork,work,control,info)

        INTEGER, INTENT(IN) :: a_n ! dimension of subproblem ND is applied to
        INTEGER, INTENT(IN) :: a_ne ! no. nonzeros of subproblem
        INTEGER, INTENT(INOUT) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. This is then 
             ! used to hold positions for submatrices after partitioning
        INTEGER, INTENT(INOUT) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.This is then used to hold row indices for
             ! submatrices after partitioning
        INTEGER, INTENT(INOUT) :: a_weight(a_n) ! On input a_weight(i) contains 
             ! the weight of column i. This is then used to hold the weights for
             ! the submatrices after partitioning.
        INTEGER, INTENT(IN) :: sumweight ! Sum entries in a_weight. Unchanged.
        INTEGER, INTENT(OUT) :: a_n1, a_n2 ! size of the two submatrices
        INTEGER, INTENT(OUT) :: where1(a_n) ! Computed partition
        INTEGER, INTENT(IN) :: lwork
        INTEGER, INTENT(OUT) :: work(lwork)
        TYPE (mc70_control), INTENT(IN) :: control 
        INTEGER, INTENT(INOUT) :: info
      !  REAL (myreal_mc70), OPTIONAL, INTENT(OUT) :: real_work(a_n)

! ---------------------------------------------
! Local variables
        INTEGER :: unit_diagnostics ! unit on which to print diagnostics
        LOGICAL :: printi, printd, use_multilevel
        INTEGER :: partition_ptr ! pointer into work array
        INTEGER :: work_ptr ! pointer into work array
        INTEGER :: partition_method
        INTEGER :: i,j,st,nstrt,nend
        INTEGER :: a_weight_1,a_weight_2,a_weight_sep, ref_method
        INTEGER, ALLOCATABLE :: work1(:)
        REAL(myreal_mc70) :: dummy

! ---------------------------------------------
! Printing levels
        unit_diagnostics = control%unit_diagnostics
        printi = (control%print_level==1 .AND. unit_diagnostics>=0)
        printd = (control%print_level>=2 .AND. unit_diagnostics>=0)
 
        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Start finding a partition'
        END IF

        ! Find the partition
        IF (control%partition_method .LE. 1) THEN
          partition_method = 1
        ELSE
          partition_method = 2
        END IF
         nstrt = -1
         nend = -1

        partition_ptr = 0 ! length a_n
        work_ptr = partition_ptr + a_n ! max length needed 11*a_n+a_ne
        IF (lwork .LT. work_ptr+9*a_n+sumweight+a_ne/2) THEN
         ALLOCATE (work1(9*a_n+sumweight+a_ne/2),stat = st)
         IF (st/=0) THEN
          info = mc70_err_memory_alloc
          RETURN
         END IF

   
         SELECT CASE (partition_method)
         CASE (1)
           ! Ashcraft method
           use_multilevel = .false.
           CALL mc70_ashcraft(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,2,nstrt,&
               nend,a_n1,a_n2,a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work1(1:9*a_n+sumweight+a_ne/2),control,info,dummy,&
               use_multilevel)
           !  write(*,*) 'dd',a_n1,a_n2,a_n
         CASE (2)
           ! Level set method
           use_multilevel = .false.
           CALL mc70_level_set(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,2,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n), &
               work1(1:4*a_n),control,info,dummy,use_multilevel)
         END SELECT

         IF (a_n1 .NE. 0 .AND. a_n2 .NE. 0 .AND. a_n .GT. 3) THEN
           IF (a_n1+a_n2 .LT. a_n) THEN
          ! Refine the partition
           IF (control%refinement .GE. 2) THEN
              IF (min(a_weight_1,a_weight_2) + a_weight_sep .LT. &
               max(a_weight_1,a_weight_2)) THEN
                 ref_method = 2
              ELSE
                 ref_method = 1
              END IF
           ELSE
              IF (control%refinement .LT. 1) THEN
                 ref_method = 1
              ELSE
                 ref_method = 2
              END IF
           END IF
           IF (ref_method .EQ. 1) THEN
            CALL mc70_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work1(1:8*a_n+sumweight+a_ne/2),control) 
             CALL mc70_refine(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work1(1:8*a_n+sumweight+a_ne/2),control,ref_method) 
           ELSE
           CALL mc70_refine(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work1(1:8*a_n+sumweight+a_ne/2),control,ref_method) 

            CALL mc70_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work1(1:8*a_n+sumweight+a_ne/2),control) 
           END IF


           END IF  
         ELSE
          GOTO 10
         END IF

         DO i = 1,a_n1
          j = work(partition_ptr+i)
          where1(j) = mc70_part1_flag
         END DO
         DO i = a_n1+1,a_n1+a_n2
          j = work(partition_ptr+i)
          where1(j) = mc70_part2_flag
         END DO
         DO i = a_n1+a_n2+1,a_n
          j = work(partition_ptr+i)
          where1(j) = mc70_sep_flag
         END DO
         DEALLOCATE (work1,stat = st)
         IF (st/=0) THEN
          info = mc70_err_memory_alloc
          RETURN
         END IF

        ELSE
         SELECT CASE (partition_method)
         CASE (1)
           ! Ashcraft method
           use_multilevel = .false.
           CALL mc70_ashcraft(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,2,&
               nstrt,nend,a_n1,&
               a_n2,a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work(work_ptr+1:work_ptr+9*a_n+sumweight+a_ne/2),control,info,&
               dummy,use_multilevel)
           !  write(*,*) 'dd',a_n1,a_n2,a_n
         CASE (2)
           ! Level set method
           CALL mc70_level_set(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,2,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n), &
               work(work_ptr+1:work_ptr+4*a_n),control,info,dummy,use_multilevel)
         END SELECT

         IF (a_n1 .NE. 0 .AND. a_n2 .NE. 0 .AND. a_n .GT. 3) THEN
           IF (a_n1+a_n2 .LT. a_n) THEN

           IF (control%refinement .GE. 2) THEN
              IF (min(a_weight_1,a_weight_2) + a_weight_sep .LT. &
               max(a_weight_1,a_weight_2)) THEN
                 ref_method = 2
              ELSE
                 ref_method = 1
              END IF
           ELSE
              IF (control%refinement .LT. 1) THEN
                 ref_method = 1
              ELSE
                 ref_method = 2
              END IF
           END IF
           IF (ref_method .EQ. 1) THEN
            CALL mc70_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control) 
             CALL mc70_refine(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control,ref_method) 
           ELSE
           CALL mc70_refine(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control,ref_method) 
            CALL mc70_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,a_n2,&
               a_weight_1,a_weight_2,a_weight_sep,&
               work(partition_ptr+1:partition_ptr+a_n),&
               work(work_ptr+1:work_ptr+8*a_n+sumweight+a_ne/2),control) 

           END IF


           END IF  
         ELSE
          GOTO 10
         END IF

         DO i = 1,a_n1
          j = work(partition_ptr+i)
          where1(j) = mc70_part1_flag
         END DO
         DO i = a_n1+1,a_n1+a_n2
          j = work(partition_ptr+i)
          where1(j) = mc70_part2_flag
         END DO
         DO i = a_n1+a_n2+1,a_n
          j = work(partition_ptr+i)
          where1(j) = mc70_sep_flag
         END DO
        END IF
 
        IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'Partition found'
        END IF
        GOTO 20

10      IF (printi .or. printd) THEN
          WRITE (unit_diagnostics,'(a)') ' '
          WRITE (unit_diagnostics,'(a)') 'No partition found'
        END IF

20        info = 0
        IF (printi .or. printd) THEN
          CALL mc70_print_message(info,unit_diagnostics, &
            'mc70_partition')
        END IF
        RETURN
      
      END SUBROUTINE mc70_coarse_partition


! *****************************************************************
      SUBROUTINE multigrid_deallocate(grid,info)
! deallocate a grid (at given level between last and first)
        TYPE (mc70_multigrid), POINTER :: grid
        INTEGER :: st, info, info65

        CALL mc65_matrix_destruct(grid%graph,info65)
        IF (info65/=0) THEN
          info = mc70_err_memory_dealloc
          RETURN
        END IF

        CALL mc65_matrix_destruct(grid%p,info65)
        IF (info65/=0) THEN
          info = mc70_err_memory_dealloc
          RETURN
        END IF

        DEALLOCATE (grid%graph,grid%p,grid%where,grid%row_wgt,STAT=st)
        IF (st/=0) THEN
          info = mc70_err_memory_dealloc
          RETURN
        END IF

        DEALLOCATE (grid,STAT=st)
        IF (st/=0) THEN
          info = mc70_err_memory_dealloc
          RETURN
        END IF

      END SUBROUTINE multigrid_deallocate

! *****************************************************************
      SUBROUTINE multigrid_deallocate_last(grid,info)

! deallocate a grid (at the last level). In this case the matrix grid%graph
! has not been formed yet
        TYPE (mc70_multigrid), POINTER :: grid
        INTEGER, INTENT (INOUT) :: info


        INTEGER :: ierr

        CALL mc65_matrix_destruct(grid%p,ierr)
        IF (ierr/=0) THEN
          info = mc70_err_memory_dealloc
          RETURN
        END IF

        DEALLOCATE (grid%graph,grid%p,grid%where,grid%row_wgt,STAT=ierr)
        IF (ierr/=0) THEN
          info = mc70_err_memory_dealloc
          RETURN
        END IF

        DEALLOCATE (grid,STAT=ierr)
        IF (ierr/=0) THEN
          info = mc70_err_memory_dealloc
          RETURN
        END IF

      END SUBROUTINE multigrid_deallocate_last
! *****************************************************************
      SUBROUTINE multigrid_deallocate_first(nvtx,n,grid,info)
! deallocate a grid (at the first level). In this case the matrix grid%p
! does not exist
        TYPE (mc70_multigrid) :: grid
        INTEGER, INTENT (INOUT) :: info
        INTEGER :: n, nvtx
        INTEGER :: ierr

        IF (nvtx/=n) THEN
          CALL mc65_matrix_destruct(grid%graph,ierr)
          IF (ierr/=0) THEN
            info = mc70_err_memory_dealloc
            RETURN
          END IF
        END IF

! in subroutine front, grid%graph is not allocated but is pointed to the
! finest level graph. so no need to deallocate
        NULLIFY (grid%graph)
        DEALLOCATE (grid%where,grid%row_wgt,STAT=ierr)
        IF (ierr/=0) info = mc70_err_memory_dealloc

      END SUBROUTINE multigrid_deallocate_first

!***************************************************************
      SUBROUTINE coarsen_hec(mp,print_level,grid,info)
! coarsen the grid using heavy-edge collapsing and set up the
! coarse grid equation, the prolongator and restrictor

        INTEGER, INTENT (IN) :: mp, print_level
        INTEGER, INTENT (INOUT) :: info
        TYPE (mc70_multigrid), INTENT (INOUT), TARGET :: grid
        TYPE (mc70_multigrid), POINTER :: cgrid

        ALLOCATE (cgrid)

        cgrid%fine => grid
        grid%coarse => cgrid

        IF (print_level==3) CALL level_print(mp,'before coarsening', &
          grid%level)

! find the prolongator
        CALL prolng_heavy_edge(grid,info)

        IF (print_level==3) CALL level_print(mp,'after coarsening ', &
          grid%level)

        cgrid%level = grid%level + 1

      END SUBROUTINE coarsen_hec

!********************************************
      SUBROUTINE prolng_heavy_edge(grid,info)

!    calculate the prolongator for heavy-edge collapsing:
!    match the vertices of the heaviest edges

        INTEGER, INTENT (INOUT) :: info
! input fine grid
        TYPE (mc70_multigrid), INTENT (INOUT) :: grid

! coarse grid based on the fine grid
        TYPE (mc70_multigrid), POINTER :: cgrid

! the fine grid row connectivity graph
        TYPE (zd11_type), POINTER :: graph

! the coarse grid prolongator
        TYPE (zd11_type), POINTER :: p

! the number of fine and coarse grid vertices
        INTEGER :: nvtx, cnvtx

! working variables
        INTEGER :: v, u, j, i, k, st
        INTEGER :: nz

        INTEGER, POINTER, DIMENSION (:) :: the_row
        REAL (myreal_mc70), POINTER, DIMENSION (:) :: the_row_val

! whether a vertex is matched already
        INTEGER, PARAMETER :: unmatched = -1

! matching status of each vertex
        INTEGER, POINTER, DIMENSION (:) :: match

! maximum weight and index of edges connected to the current vertex
        REAL (myreal_mc70) :: maxwgt
        INTEGER :: maxind

! order of which vertex is visited for matching
        INTEGER, ALLOCATABLE, DIMENSION (:) :: order
        INTEGER, ALLOCATABLE, DIMENSION (:) :: nrow

! allocate the prolongation matrix pointers
        cgrid => grid%coarse
        graph => grid%graph
        ALLOCATE (cgrid%p)

! allocate the graph and matrix pointer and the mincut pointer
! so that everything is defined
        ALLOCATE (cgrid%graph)

        p => cgrid%p
        nvtx = graph%n

! prolongator start here ================================

! initialise the matching status and randomly permute the vertex order
        ALLOCATE (match(nvtx),order(nvtx),nrow(grid%size),STAT=st)
        IF (st/=0) THEN
          info = mc70_err_memory_alloc
          RETURN
        END IF
        match = unmatched

        DO i = 1, nvtx
          order(i) = i
        END DO

! loop over each vertex and match along the heaviest edge
        cnvtx = 0
        nz = 0
        nrow = 0
        DO i = 1, nvtx
          v = order(i)
! If already matched, next vertex please
          IF (match(v)/=unmatched) CYCLE
! access the col. indices of row v
          CALL mc65_matrix_getrow(graph,v,the_row)
! access the entry values of row v
          CALL mc65_matrix_getrowval(graph,v,the_row_val)
          maxwgt = -huge(0.0_myreal_mc70)
! in the case no match is found then match itself
          maxind = v
! Loop over entries in row v
          DO j = 1, size(the_row)
! u is col index of en entry in row v (so u is neighbor of v)
            u = the_row(j)
! heavy edge matching
! if u is unmatched and value of the entry in col. u is greater
! than maxwgt, select u as the matching.
            IF (match(u)==unmatched .AND. maxwgt<abs(the_row_val(j))) THEN
              maxwgt = abs(the_row_val(j))
              maxind = u
            END IF
          END DO
! the neighbor with heaviest weight
          match(v) = maxind
! mark maxind as having been matched
          match(maxind) = v
! increase number of vertices in coarse graph by 1
          cnvtx = cnvtx + 1
! construct the prolongation matrix: find vertex v and maxind is linked
! with the coarse grid vertex cnvtx
          nz = nz + 1
          nrow(v) = nrow(v) + 1
          IF (maxind/=v) THEN
            nz = nz + 1
            nrow(maxind) = nrow(maxind) + 1
          END IF
        END DO

! storage allocation for col. indices and values of prolongation
! matrix P (order nvtx * cnvtx)
        CALL mc65_matrix_construct(p,nvtx,nz,info,n=cnvtx,type='general')
        p%val = 0.0_myreal_mc70

! Set column pointers for P
        p%ptr(1) = 1
        DO i = 1, nvtx
          p%ptr(i+1) = p%ptr(i) + nrow(i)
        END DO

! the row counter: nothing in matrix P filled yet
        nrow = 0

! Now fill in entries of matrix P
        match = unmatched
! loop over each vertex and match along the heaviest edge
        cnvtx = 0
        DO i = 1, nvtx
          v = order(i)
! if already matched, next vertex please
          IF (match(v)/=unmatched) CYCLE
          CALL mc65_matrix_getrow(graph,v,the_row)
          CALL mc65_matrix_getrowval(graph,v,the_row_val)
          maxwgt = -huge(0.0_myreal_mc70)
! in the case no match is found then match itself
          maxind = v
          DO j = 1, size(the_row)
            u = the_row(j)
! heavy edge matching
            IF (match(u)==unmatched .AND. maxwgt<abs(the_row_val(j))) THEN
              maxwgt = abs(the_row_val(j))
              maxind = u
            END IF
          END DO
! the neighbor with heaviest weight
          match(v) = maxind
          match(maxind) = v
          cnvtx = cnvtx + 1
! construct the prolongation matrix: vertex v and maxind are linked
! with the coarse grid vertex cnvtx
! k points to start of row v and now(v) holds number of
! entries in row v that have been filled (so k + nrow(v) is
! first free entry in row v)
          k = p%ptr(v)
          p%col(k+nrow(v)) = cnvtx
          p%val(k+nrow(v)) = 1.0_myreal_mc70
          nrow(v) = nrow(v) + 1

          IF (maxind/=v) THEN
            k = p%ptr(maxind)
            p%col(k+nrow(maxind)) = cnvtx
            p%val(k+nrow(maxind)) = 1.0_myreal_mc70
            nrow(maxind) = nrow(maxind) + 1
          END IF

        END DO

! deallocate the match pointer for vertices
        DEALLOCATE (match,STAT=st)
        IF (st/=0) THEN
          info = mc70_err_memory_dealloc
          RETURN
        END IF

        DEALLOCATE (order,nrow,STAT=st)
        IF (st/=0) THEN
          info = mc70_err_memory_dealloc
          RETURN
        END IF

! prolongator end
! size of coarse grid
        cgrid%size = cnvtx

! allocate coarse grid quantities
        ALLOCATE (cgrid%where(cnvtx),cgrid%row_wgt(cnvtx),STAT=st)
        IF (st/=0) info = mc70_err_memory_alloc

      END SUBROUTINE prolng_heavy_edge

!*******************************************************************
      SUBROUTINE level_print(mp,title1,level,title2,res)

        CHARACTER (len=*), INTENT (IN) :: title1
        INTEGER, INTENT (IN) :: mp, level
        REAL (myreal_mc70), OPTIONAL, INTENT (IN) :: res
        CHARACTER (len=*), OPTIONAL, INTENT (IN) :: title2
        INTEGER :: char_len1, char_len2

        char_len1 = len_trim(title1)

        IF (present(res) .AND. present(title2)) THEN
          char_len2 = len_trim(title2)
          WRITE (mp,'(a,i4,a,g14.3)') title1, level, title2, res
        ELSE IF (present(res)) THEN
          WRITE (mp,'(a,i4,g14.3)') title1, level,  res
        ELSE
          WRITE (mp,'(a,i4)') title1, level
        END IF

      END SUBROUTINE level_print

!*************************************************
      SUBROUTINE galerkin_graph(matrix,p,cmatrix,info)

! Given matrix on fine grid and a prolongation operator p,
! find the coarse matrix R*A*P

! matrix: fine grid matrix
        TYPE (zd11_type), INTENT (IN) :: matrix
! p: prolongation operator
        TYPE (zd11_type), INTENT (IN) :: p
! cmatrix: coarse grid matrix
        TYPE (zd11_type), INTENT (OUT) :: cmatrix
! r: restriction operator
        TYPE (zd11_type) :: r

! nvtx,cnvtx: size of fine and coarse grid
        INTEGER :: nvtx, cnvtx
        INTEGER :: nz

        INTEGER, INTENT (INOUT) :: info
        INTEGER :: info65

!!write(*,*) 'start 4637'
        CALL mc65_matrix_transpose(p,r,info65)
!!write(*,*) 'end 4637'
        IF (info65<0) THEN
          info = info65
          RETURN
        END IF
        nvtx = matrix%n
        cnvtx = p%n

! get the size of the coarse matrix first
!!write(*,*) 'start 4648'
        CALL galerkin_graph_rap_size(nvtx,cnvtx,nz,p%ptr(nvtx+1)-1,p%col, &
          p%ptr,matrix%ptr(nvtx+1)-1,matrix%col,matrix%ptr,r%ptr(cnvtx+1)-1, &
          r%col,r%ptr,info)
!!write(*,*) 'end 4648'
        IF (info<0) RETURN

!!write(*,*) 'start 4654'
        CALL mc65_matrix_construct(cmatrix,cnvtx,nz,info65)
!!write(*,*) 'end 4654'
        IF (info65<0) THEN
          info = info65
          RETURN
        END IF

!!write(*,*) 'start 4662'
        CALL galerkin_graph_rap(nvtx,cnvtx,p%ptr(nvtx+1)-1,p%val,p%col,p%ptr, &
          matrix%ptr(nvtx+1)-1,matrix%val,matrix%col,matrix%ptr, &
          r%ptr(cnvtx+1)-1,r%val,r%col,r%ptr,nz,cmatrix%val,cmatrix%col, &
          cmatrix%ptr,info)
!!write(*,*) 'end 4662'
        IF (info<0) RETURN

        CALL mc65_matrix_destruct(r,info65)
        IF (info65<0) THEN
          info = info65
          RETURN
        END IF

      END SUBROUTINE galerkin_graph

!*************************************************

      SUBROUTINE galerkin_graph_rap_size(nvtx,cnvtx,nz,nzp,pcol,pptr,nzaa, &
          acol,aptr,nzr,rcol,rptr,info)
! get the number of nonzeros in R*A*P
! nvtx: size of aa matrix
! cnvtx: size of ca matrix
        INTEGER, INTENT (IN) :: nvtx, cnvtx
! nz: number of nonzeros in R*A*P
        INTEGER, INTENT (OUT) :: nz

! P: matrix
        INTEGER, INTENT (IN) :: nzp
        INTEGER, INTENT (IN), DIMENSION (nzp) :: pcol
        INTEGER, INTENT (IN), DIMENSION (nvtx+1) :: pptr
! aa: matrix
        INTEGER, INTENT (IN) :: nzaa
        INTEGER, INTENT (IN), DIMENSION (nzaa) :: acol
        INTEGER, INTENT (IN), DIMENSION (nvtx+1) :: aptr
! R: matrix
        INTEGER, INTENT (IN) :: nzr
        INTEGER, INTENT (IN), DIMENSION (nzr) :: rcol
        INTEGER, INTENT (IN), DIMENSION (cnvtx+1) :: rptr

! mask: masking array to see if an entry has been seen before
        INTEGER, DIMENSION (:), ALLOCATABLE :: mask
! i,j,k: loop index
        INTEGER :: i, j, k
! nz: number of nonzeros so far in ca
        INTEGER :: nz1
! various neighbors
        INTEGER :: neigh, neighneigh

! col: column index of a row of r*matrix
        INTEGER, DIMENSION (:), ALLOCATABLE :: col

        INTEGER, INTENT (INOUT) :: info
        INTEGER st

        ALLOCATE (mask(nvtx),col(nvtx),STAT=st)
        IF (st/=0) THEN
          info = mc70_err_memory_alloc
          RETURN
        END IF
        mask = 0
        nz = 0
! loop over coarse grid points
        DO i = 1, cnvtx
! first form row i of (r*matrix)
          nz1 = 0
! for each vertex D that restricts to C (including itself).
          DO j = rptr(i), rptr(i+1) - 1
            neigh = rcol(j)
! find D's neighbor
            DO k = aptr(neigh), aptr(neigh+1) - 1
              neighneigh = acol(k)
              IF (mask(neighneigh)/=i) THEN
                nz1 = nz1 + 1
                col(nz1) = neighneigh
                mask(neighneigh) = i
              END IF
            END DO
          END DO
! form row i of (r*matrix)*p
          DO j = 1, nz1
            neigh = col(j)
            DO k = pptr(neigh), pptr(neigh+1) - 1
              neighneigh = pcol(k)
              IF (mask(neighneigh)/=-i .AND. neighneigh/=i) THEN
                nz = nz + 1
                mask(neighneigh) = -i
              END IF
            END DO
          END DO
        END DO
        DEALLOCATE (mask,col,STAT=st)
        IF (st/=0) info = mc70_err_memory_dealloc

      END SUBROUTINE galerkin_graph_rap_size
! ******************************************************
      SUBROUTINE galerkin_graph_rap(nvtx,cnvtx,nzp,pa,pcol,pptr,nzaa,aa,acol, &
          aptr,nzr,ra,rcol,rptr,nzca,ca,ccol,cptr,info)
! multiply R*A*P to get CA
! nvtx: size of aa matrix
! cnvtx: size of ca matrix
        INTEGER, INTENT (IN) :: nvtx, cnvtx

! p: matrix
        INTEGER, INTENT (IN) :: nzp
        REAL (myreal_mc70), INTENT (IN), DIMENSION (nzp) :: pa
        INTEGER, INTENT (IN), DIMENSION (nzp) :: pcol
        INTEGER, INTENT (IN), DIMENSION (nvtx+1) :: pptr
! aa: matrix
        INTEGER, INTENT (IN) :: nzaa
        REAL (myreal_mc70), INTENT (IN), DIMENSION (:) :: aa
        INTEGER, INTENT (IN), DIMENSION (nzaa) :: acol
        INTEGER, INTENT (IN), DIMENSION (:) :: aptr
! r: matrix
        INTEGER, INTENT (IN) :: nzr
        REAL (myreal_mc70), INTENT (IN), DIMENSION (nzr) :: ra
        INTEGER, INTENT (IN), DIMENSION (nzr) :: rcol
        INTEGER, INTENT (IN), DIMENSION (cnvtx+1) :: rptr
! ca: matrix
        INTEGER, INTENT (IN) :: nzca
        REAL (myreal_mc70), INTENT (INOUT), DIMENSION (nzca) :: ca
        INTEGER, INTENT (INOUT), DIMENSION (nzca) :: ccol
        INTEGER, INTENT (INOUT), DIMENSION (cnvtx+1) :: cptr


! mask: masking array to see if an entry has been seen before
        INTEGER, DIMENSION (:), ALLOCATABLE :: mask
! i,j,k,l: loop index
        INTEGER :: i, j, k
! nz: number of nonzeros so far in ca
        INTEGER :: nz, nzz, nz1
! various neighbors
        INTEGER :: neigh, neighneigh
! r_ij: (i,j) element of r
        REAL (myreal_mc70) :: r_ij
! col: column index of a row of r*matrix
! a: values of a row of r*matrix
        INTEGER, DIMENSION (:), ALLOCATABLE :: col
        REAL (myreal_mc70), DIMENSION (:), ALLOCATABLE :: a

        INTEGER, INTENT (INOUT) :: info
        INTEGER st

        ALLOCATE (col(nvtx),a(nvtx),mask(nvtx),STAT=st)
        IF (st/=0) THEN
          info = mc70_err_memory_alloc
          RETURN
        END IF
! now get the entries of the coarse matrix
        cptr(1) = 1
        mask = 0
        nz = 0
! loop over every coarse grid point
        DO i = 1, cnvtx
! first form row i of (r*matrix)
          nz1 = 0
! foreach each vertex D that restricts to C (including itself).
          DO j = rptr(i), rptr(i+1) - 1
            neigh = rcol(j)
            r_ij = ra(j)
! find D's neighbor
            DO k = aptr(neigh), aptr(neigh+1) - 1
              neighneigh = acol(k)
              nzz = mask(neighneigh)
              IF (nzz==0) THEN
                nz1 = nz1 + 1
                col(nz1) = neighneigh
                a(nz1) = r_ij*aa(k)
                mask(neighneigh) = nz1
              ELSE
                a(nzz) = a(nzz) + r_ij*aa(k)
              END IF
            END DO
          END DO
          mask(col(1:nz1)) = 0

! form row i of (r*matrix)*p
          DO j = 1, nz1
            neigh = col(j)
            r_ij = a(j)
            DO k = pptr(neigh), pptr(neigh+1) - 1
              neighneigh = pcol(k)
              IF (neighneigh==i) CYCLE
              nzz = mask(neighneigh)
              IF (nzz==0) THEN
                nz = nz + 1
                mask(neighneigh) = nz
                ca(nz) = r_ij*pa(k)
                ccol(nz) = neighneigh
              ELSE
                ca(nzz) = ca(nzz) + r_ij*pa(k)
              END IF
            END DO
          END DO

          mask(ccol(cptr(i):nz)) = 0
          cptr(i+1) = nz + 1
        END DO

        DEALLOCATE (col,a,mask,STAT=st)
        IF (st/=0) info = mc70_err_memory_dealloc

      END SUBROUTINE galerkin_graph_rap

! ---------------------------------------------------
! mc70_refine_trim
! ---------------------------------------------------
! Given a partition, trim the partition to make it minimal
      SUBROUTINE mc70_refine_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,&
               a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,control)
        INTEGER, INTENT(IN) :: a_n ! order of matrix
        INTEGER, INTENT(IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT(IN) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. 
        INTEGER, INTENT(IN) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.
        INTEGER, INTENT(IN) :: a_weight(a_n) ! On input a_weight(i) contains 
             ! the weight of column i 
        INTEGER, INTENT(IN) :: sumweight ! Sum of weights in a_weight
        INTEGER, INTENT(INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT(INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT(INOUT) :: a_weight_1,a_weight_2,a_weight_sep ! Weighted 
             ! size of partitions and separator
        INTEGER, INTENT(INOUT) :: partition(a_n) !First a_n1 entries contain
             ! list of (local) indices in partition 1; next a_n2 entries  
             ! contain list of (local) entries in partition 2; entries in 
             ! separator are listed at the end. This is updated to the new 
             ! partition
        INTEGER, INTENT(OUT) :: work(8*a_n+sumweight+a_ne/2) ! Work array
        TYPE(mc70_control), INTENT(IN) :: control

! ---------------------------------------------
! Local variables
        INTEGER :: work_part, work_next, work_prev, a_n1_orig, a_n2_orig
        INTEGER :: head1, head2, tail1, tail2
        INTEGER, PARAMETER :: sep1=-1
        INTEGER, PARAMETER :: sep2=-2
        INTEGER, PARAMETER :: sep3=-3
        INTEGER :: i,j,k,l,m,p,q,w1,w2
        LOGICAL :: next1, next2, imbal
        REAL(myreal_mc70) :: t1, t2
        REAL(myreal_mc70) :: ratio
      !  write(*,*) 'a_n1,a_n2',a_n1,a_n2
      !  write(*,*) 'partition'
      !  write(*,*) partition

        ratio = control%ratio
        IF (ratio.GT.real(sumweight-2)) THEN
           imbal = .FALSE.
        ELSE
           imbal = .TRUE.
        END IF
 
        ! Set work(work_part+1:work_part+a_n) to hold flags to indicate what 
        ! part of the partition the nodes are in
        work_part = 0
        k = 0
        l=0
        m=0
        DO i = 1, a_n1
            j = partition(i)
            work(work_part+j) = mc70_part1_flag
            k = k+1
        END DO
        DO i = a_n1+1, a_n1+a_n2
            j = partition(i)
            work(work_part+j) = mc70_part2_flag
            l = l+1
        END DO
        DO i = a_n1+a_n2+1, a_n
            j = partition(i)
            work(work_part+j) = mc70_sep_flag
            m = m+1
        END DO
        a_n1_orig = a_n1
        a_n2_orig = a_n2
        a_weight_sep = sumweight - a_weight_1 - a_weight_2

        ! Create two lists
        work_next = work_part + a_n
        work_prev = work_next + a_n
        head1 = 0
        head2 = 0
        tail1 = 0
        tail2 = 0
        DO i = a_n1+a_n2+1, a_n
            next1 = .false.
            next2 = .false.
            j = partition(i)
            IF (j .LT. a_n) THEN
             k = a_ptr(j+1)-1
            ELSE
             k = a_ne
            END IF
            DO l = a_ptr(j),k
              m = a_row(l)
              IF (work(work_part+m) .EQ. mc70_part1_flag) THEN
                next1 = .true.
              ELSE IF (work(work_part+m) .EQ. mc70_part2_flag) THEN
                next2 = .true.
              END IF
            END DO
         !   write(*,*) j, next1, next2
            IF (next1 .and. .not. next2) THEN
              ! Add to list 1
            !  write(*,*) 'Add',j,'list1'
              IF (head1.eq.0) THEN
                head1 = j
                work(work_next + j) = 0
                work(work_prev + j) = 0
                tail1 = j
              ELSE
                work(work_next + tail1) = j
                work(work_prev + j) = tail1
                work(work_next + j) = 0
                tail1 = j
              END IF
              work(work_part+j) = sep1
            ELSE IF (next2 .and. .not. next1) THEN
              ! Add to list 2
             ! write(*,*) 'Add',j,'list2'
              IF (head2.eq.0) THEN
                head2 = j
                work(work_next + j) = 0
                work(work_prev + j) = 0
                tail2 = j
              ELSE
                work(work_next + tail2) = j
                work(work_prev + j) = tail2
                work(work_next + j) = 0
                tail2 = j
              END IF
              work(work_part+j) = sep2
            ELSE IF (next1 .AND. next2) THEN
              work(work_part+j) = sep3
            !  write(*,*) 'Vertex',j,'next both'

            ELSE
              continue
            END IF
        END DO

        DO WHILE (head1.gt.0 .or. head2.gt.0)
           IF (head1 .gt. 0 .and. head2 .gt. 0) THEN
             w1 = a_weight(head1)
             w2 = a_weight(head2)
             CALL cost_function(a_weight_1+w1,a_weight_2,a_weight_sep-w1,&
              sumweight,ratio,imbal,t1)
             CALL cost_function(a_weight_1,a_weight_2+w2,a_weight_sep-w2,&
              sumweight,ratio,imbal,t2)

             IF (t1 .lt. t2) THEN
               GOTO 200
             ELSE
               GOTO 210
             END IF

           ELSE IF (head1 .gt.0) THEN
             GOTO 200
           ELSE 
             GOTO 210
           END IF


             ! move entry from separator to partition1
200          i = head1
            ! write(*,*) 'move ',i,'into partition 1'
             work(work_part+i) = mc70_part1_flag
             head1 = work(work_next+i)
             work(work_next+i) = 0
             a_n1 = a_n1+1
             a_weight_1 = a_weight_1 + a_weight(i)
             a_weight_sep = a_weight_sep - a_weight(i)
             ! update list
             IF (i .LT. a_n) THEN
              k = a_ptr(i+1)-1
             ELSE
              k = a_ne
             END IF
             DO l = a_ptr(i),k
              j = a_row(l)
              m = work(work_part+j)
              SELECT CASE (m)
               CASE (mc70_sep_flag)
                  ! Add to list 1
              !  write(*,*) 'add',j, 'to list 1'
                work(work_next + tail1) = j
                work(work_prev + j) = tail1
                work(work_next + j) = 0
                tail1 = j
                work(work_part + j) = sep1
                IF (head1 .EQ. 0) THEN
                    head1 = j
                END IF

               CASE (sep2)
                   ! Remove from list 2
               ! write(*,*) 'remove',j, 'from list 2'
                p = work(work_prev + j)
                q = work(work_next + j)

                IF (j .NE. head2 .AND. j .NE. tail2) THEN
                 work(work_prev + q) = p
                 work(work_next + p) = q
                 work(work_prev + j) = 0
                 work(work_next + j) = 0
                ELSE IF (j .NE. head2 .AND. j .EQ. tail2) THEN
                 work(work_next + p) = 0
                 work(work_prev + j) = 0
                 tail2 = p
                ELSE IF (j .NE. tail2 .AND. j .EQ. head2) THEN
                 work(work_prev + q) = p
                 work(work_next + j) = 0
                 head2 = q
                ELSE
                 head2 = 0
                 tail2 = 0

                END IF
                work(work_part + j) = sep3
                   
              END SELECT
            END DO
            GOTO 220

             ! move entry from separator to partition 2
210          i = head2
           !  write(*,*) 'move ',i,'into partition 2'
             work(work_part+i) = mc70_part2_flag
             head2 = work(work_next+i)
             work(work_next+i) = 0
            a_n2 = a_n2 + 1
            a_weight_2 = a_weight_2 + a_weight(i)
            a_weight_sep = a_weight_sep - a_weight(i)
             ! update list
             IF (i .LT. a_n) THEN
              k = a_ptr(i+1)-1
             ELSE
              k = a_ne
             END IF
             DO l = a_ptr(i),k
              j = a_row(l)
              m = work(work_part+j)
            !  write(*,*) j,m
              SELECT CASE (m)
               CASE (mc70_sep_flag)
                  ! Add to list 2
              !  write(*,*) 'add',j, 'to list 2'
                work(work_next + tail2) = j
                work(work_prev + j) = tail2
                work(work_next + j) = 0
                tail2 = j
                work(work_part + j) = sep2
                IF (head2 .EQ. 0) THEN
                    head2 = j
                END IF

               CASE (sep1)
                   ! Remove from list 1
              !  write(*,*) 'remove',j, 'from list 1'
                p = work(work_prev + j)
                q = work(work_next + j)
                
                IF (j .NE. head1 .AND. j .NE. tail1) THEN
                 work(work_prev + q) = p
                 work(work_next + p) = q
                 work(work_prev + j) = 0
                 work(work_next + j) = 0
                ELSE IF (j .NE. head1 .AND. j .EQ. tail1) THEN
                 work(work_next + p) = 0
                 work(work_prev + j) = 0
                 tail1 = p
                ELSE IF (j .NE. tail1 .AND. j .EQ. head1) THEN
                 work(work_prev + q) = p
                 work(work_next + j) = 0
                 head1 = q
                ELSE
                 head1 = 0
                 tail1 = 0
                END IF
                work(work_part + j) = sep3
              END SELECT
            END DO

220        CONTINUE

        END DO

      ! Check for any entries in separator that are still inside boundary and 
      ! move into a partition
      work(work_next+a_n1_orig + a_n2_orig +1:work_next+a_n) = 0
      
      DO i = a_n1_orig + a_n2_orig +1,a_n
        j = partition(i)
        IF (work(work_part+j) .EQ. mc70_sep_flag) THEN
          ! j is not on the boundary
   !       write(*,*) j, 'not boundary'
          IF (a_weight_1 .LT. a_weight_2) THEN
            ! Move j into partition 1
            work(work_part+j) = mc70_part1_flag
            a_n1 = a_n1 + 1
            a_weight_1 = a_weight_1 + a_weight(j)
            a_weight_sep = a_weight_sep - a_weight(j)
            
            head1 = j
            tail1 = j
            DO WHILE (head1 .GT. 0)
             q = head1
             IF (q .LT. a_n) THEN
              k = a_ptr(q+1)-1
             ELSE
              k = a_ne
             END IF
             DO l = a_ptr(q),k
               p = a_row(l)
               IF (work(work_part+p) .EQ. mc70_sep_flag) THEN
                 work(work_part+p) = mc70_part1_flag
                 a_n1 = a_n1 + 1
                 a_weight_1 = a_weight_1 + a_weight(q)
                 a_weight_sep = a_weight_sep - a_weight(q)
                 work(work_next+tail1) = p
                 tail1 = p
               END IF
             END DO
             IF (head1 .EQ. tail1) THEN
               head1 = 0
               tail1 = 0
             ELSE
               head1 = work(work_next+q)
               work(work_next+q) = 0
             END IF

            END DO

          ELSE
            ! Move j into partition 2
            work(work_part+j) = mc70_part2_flag
            a_n2 = a_n2 + 1
            a_weight_2 = a_weight_2 + a_weight(j)
            a_weight_sep = a_weight_sep - a_weight(j)
            head2 = j
            tail2 = j

            DO WHILE (head2 .GT. 0)
             q = head2
             IF (q .LT. a_n) THEN
              k = a_ptr(q+1)-1
             ELSE
              k = a_ne
             END IF
             DO l = a_ptr(q),k
               p = a_row(l)
               IF (work(work_part+p) .EQ. mc70_sep_flag) THEN
                 work(work_part+p) = mc70_part2_flag
                 a_n2 = a_n2 + 1
                 a_weight_2 = a_weight_2 + a_weight(q)
                 a_weight_sep = a_weight_sep - a_weight(q)
                 work(work_next+tail2) = p
                 tail2 = p
               END IF
             END DO
             IF (head2 .EQ. tail2) THEN
               head2 = 0
               tail2 = 0
             ELSE
               head2 = work(work_next+q)
               work(work_next+q) = 0
             END IF
            END DO
          END IF
        END IF
      END DO
      
      ! Reset partition matrix
      i = 1
      j = a_n1+1
      k = a_n1+a_n2+1
      DO l = 1,a_n
        m = work(work_part+l)
        SELECT CASE (m)
        CASE (mc70_part1_flag)
          partition(i) = l
          i = i+1
        CASE (mc70_part2_flag)
          partition(j) = l
          j = j+1
        CASE default
          partition(k) = l
          work(work_part+l) = mc70_sep_flag
          k = k+1
        END SELECT
      END DO
      a_weight_sep = sumweight - a_weight_1 - a_weight_2
     ! call check_partition1(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,work(work_part+1:work_part+a_n),a_weight_1,a_weight_2,a_weight)

      END SUBROUTINE mc70_refine_trim

! ---------------------------------------------------
! mc70_refine_block_trim
! ---------------------------------------------------
! Given a partition, trim the partition using blocks to make it minimal
      SUBROUTINE mc70_refine_block_trim(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,&
               a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,control)
        INTEGER, INTENT(IN) :: a_n ! order of matrix
        INTEGER, INTENT(IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT(IN) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. 
        INTEGER, INTENT(IN) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.
        INTEGER, INTENT(IN) :: a_weight(a_n) ! On input a_weight(i) contains 
             ! the weight of column i 
        INTEGER, INTENT(IN) :: sumweight ! Sum of weights in a_weight
        INTEGER, INTENT(INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT(INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT(INOUT) :: a_weight_1,a_weight_2,a_weight_sep ! Weighted 
             ! size of partitions and separator
        INTEGER, INTENT(INOUT) :: partition(a_n) !First a_n1 entries contain
             ! list of (local) indices in partition 1; next a_n2 entries  
             ! contain list of (local) entries in partition 2; entries in 
             ! separator are listed at the end. This is updated to the new 
             ! partition
        INTEGER, INTENT(OUT) :: work(8*a_n+sumweight+a_ne/2) ! Work array
        TYPE(mc70_control), INTENT(IN) :: control

! ---------------------------------------------
! Local variables
        INTEGER :: work_part, work_next, work_prev, a_n1_orig, a_n2_orig
        INTEGER :: head1, head2, tail1, tail2
        INTEGER, PARAMETER :: sep1=-1
        INTEGER, PARAMETER :: sep2=-2
        INTEGER, PARAMETER :: sep3=-3
        INTEGER :: i,j,k,l,m,p,q,w1,w2,l1,l2,l1new,l2new,ll
        LOGICAL :: next1, next2, imbal
        REAL(myreal_mc70) :: t1, t2
        REAL(myreal_mc70) :: ratio
      !  write(*,*) 'a_n1,a_n2',a_n1,a_n2
      !  write(*,*) 'partition'
      !  write(*,*) partition
        ratio = control%ratio
        IF (ratio.GT.real(sumweight-2)) THEN
           imbal = .FALSE.
        ELSE
           imbal = .TRUE.
        END IF
 
        ! Set work(work_part+1:work_part+a_n) to hold flags to indicate what 
        ! part of the partition the nodes are in
        work_part = 0
        k = 0
        l=0
        m=0
        DO i = 1, a_n1
            j = partition(i)
            work(work_part+j) = mc70_part1_flag
            k = k+1
        END DO
        DO i = a_n1+1, a_n1+a_n2
            j = partition(i)
            work(work_part+j) = mc70_part2_flag
            l = l+1
        END DO
        DO i = a_n1+a_n2+1, a_n
            j = partition(i)
            work(work_part+j) = mc70_sep_flag
            m = m+1
        END DO
        a_n1_orig = a_n1
        a_n2_orig = a_n2
        a_weight_sep = sumweight - a_weight_1 - a_weight_2

        ! Create two lists
        work_next = work_part + a_n
        work_prev = work_next + a_n
        head1 = 0
        head2 = 0
        tail1 = 0
        tail2 = 0
        l1 = 0
        l2 = 0
        DO i = a_n1+a_n2+1, a_n
            next1 = .false.
            next2 = .false.
            j = partition(i)
            IF (j .LT. a_n) THEN
             k = a_ptr(j+1)-1
            ELSE
             k = a_ne
            END IF
            DO l = a_ptr(j),k
              m = a_row(l)
              IF (work(work_part+m) .EQ. mc70_part1_flag) THEN
                next1 = .true.
              ELSE IF (work(work_part+m) .EQ. mc70_part2_flag) THEN
                next2 = .true.
              END IF
            END DO
         !   write(*,*) j, next1, next2
            IF (next1 .and. .not. next2) THEN
              ! Add to list 1
            !  write(*,*) 'Add',j,'list1'
              IF (head1.eq.0) THEN
                head1 = j
                work(work_next + j) = 0
                work(work_prev + j) = 0
                tail1 = j
              ELSE
                work(work_next + tail1) = j
                work(work_prev + j) = tail1
                work(work_next + j) = 0
                tail1 = j
              END IF
              work(work_part+j) = sep1
              l1 = l1 + 1
            ELSE IF (next2 .and. .not. next1) THEN
              ! Add to list 2
             ! write(*,*) 'Add',j,'list2'
              IF (head2.eq.0) THEN
                head2 = j
                work(work_next + j) = 0
                work(work_prev + j) = 0
                tail2 = j
              ELSE
                work(work_next + tail2) = j
                work(work_prev + j) = tail2
                work(work_next + j) = 0
                tail2 = j
              END IF
              work(work_part+j) = sep2
              l2 = l2 + 1
            ELSE IF (next1 .AND. next2) THEN
              work(work_part+j) = sep3
            !  write(*,*) 'Vertex',j,'next both'

            ELSE
              continue
            END IF
        END DO

        DO WHILE (head1.gt.0 .or. head2.gt.0)
           IF (head1 .gt. 0 .and. head2 .gt. 0) THEN
             j = head1
             w1 = 0
             DO WHILE (j .gt. 0)
               w1 = w1 + a_weight(j)
               j = work(work_next+j)
             END DO
             j = head2
             w2 = 0
             DO WHILE (j .gt. 0)
               w2 = w2 + a_weight(j)
               j = work(work_next+j)
             END DO
             CALL cost_function(a_weight_1+w1,a_weight_2,a_weight_sep-w1,&
              sumweight,ratio,imbal,t1)
             CALL cost_function(a_weight_1,a_weight_2+w2,a_weight_sep-w2,&
              sumweight,ratio,imbal,t2)

             IF (t1 .lt. t2) THEN
               GOTO 200
             ELSE
               GOTO 210
             END IF

           ELSE IF (head1 .gt.0) THEN
             GOTO 200
           ELSE 
             GOTO 210
           END IF


             ! move entries from separator to partition1
200         CONTINUE
            l1new = 0
            DO ll = 1,l1 
             i = head1
            ! write(*,*) 'move ',i,'into partition 1'
             work(work_part+i) = mc70_part1_flag
             head1 = work(work_next+i)
             work(work_next+i) = 0
             a_n1 = a_n1+1
             a_weight_1 = a_weight_1 + a_weight(i)
             a_weight_sep = a_weight_sep - a_weight(i)
             ! update list
             IF (i .LT. a_n) THEN
              k = a_ptr(i+1)-1
             ELSE
              k = a_ne
             END IF
             DO l = a_ptr(i),k
              j = a_row(l)
              m = work(work_part+j)
              SELECT CASE (m)
               CASE (mc70_sep_flag)
                  ! Add to list 1
              !  write(*,*) 'add',j, 'to list 1'
                work(work_next + tail1) = j
                work(work_prev + j) = tail1
                work(work_next + j) = 0
                tail1 = j
                work(work_part + j) = sep1
                IF (head1 .EQ. 0) THEN
                    head1 = j
                END IF
                l1new = l1new+1

               CASE (sep2)
                   ! Remove from list 2
               ! write(*,*) 'remove',j, 'from list 2'
                p = work(work_prev + j)
                q = work(work_next + j)

                IF (j .NE. head2 .AND. j .NE. tail2) THEN
                 work(work_prev + q) = p
                 work(work_next + p) = q
                 work(work_prev + j) = 0
                 work(work_next + j) = 0
                ELSE IF (j .NE. head2 .AND. j .EQ. tail2) THEN
                 work(work_next + p) = 0
                 work(work_prev + j) = 0
                 tail2 = p
                ELSE IF (j .NE. tail2 .AND. j .EQ. head2) THEN
                 work(work_prev + q) = p
                 work(work_next + j) = 0
                 head2 = q
                ELSE
                 head2 = 0
                 tail2 = 0

                END IF
                work(work_part + j) = sep3
                l2 = l2 -1
                   
              END SELECT
            END DO
           END DO
           l1 = l1new
           GOTO 220

             ! move entry from separator to partition 2
210        CONTINUE
           l2new = 0
           DO ll = 1,l2
            i = head2
           !  write(*,*) 'move ',i,'into partition 2'
             work(work_part+i) = mc70_part2_flag
             head2 = work(work_next+i)
             work(work_next+i) = 0
            a_n2 = a_n2 + 1
            a_weight_2 = a_weight_2 + a_weight(i)
            a_weight_sep = a_weight_sep - a_weight(i)
             ! update list
             IF (i .LT. a_n) THEN
              k = a_ptr(i+1)-1
             ELSE
              k = a_ne
             END IF
             DO l = a_ptr(i),k
              j = a_row(l)
              m = work(work_part+j)
            !  write(*,*) j,m
              SELECT CASE (m)
               CASE (mc70_sep_flag)
                  ! Add to list 2
              !  write(*,*) 'add',j, 'to list 2'
                work(work_next + tail2) = j
                work(work_prev + j) = tail2
                work(work_next + j) = 0
                tail2 = j
                work(work_part + j) = sep2
                IF (head2 .EQ. 0) THEN
                    head2 = j
                END IF
                l2new = l2new+1

               CASE (sep1)
                   ! Remove from list 1
              !  write(*,*) 'remove',j, 'from list 1'
                p = work(work_prev + j)
                q = work(work_next + j)
                
                IF (j .NE. head1 .AND. j .NE. tail1) THEN
                 work(work_prev + q) = p
                 work(work_next + p) = q
                 work(work_prev + j) = 0
                 work(work_next + j) = 0
                ELSE IF (j .NE. head1 .AND. j .EQ. tail1) THEN
                 work(work_next + p) = 0
                 work(work_prev + j) = 0
                 tail1 = p
                ELSE IF (j .NE. tail1 .AND. j .EQ. head1) THEN
                 work(work_prev + q) = p
                 work(work_next + j) = 0
                 head1 = q
                ELSE
                 head1 = 0
                 tail1 = 0
                END IF
                work(work_part + j) = sep3
                l1 = l1-1
              END SELECT
            END DO
           END DO
           l2 = l2new
220        CONTINUE

        END DO

      ! Check for any entries in separator that are still inside boundary and 
      ! move into a partition
      work(work_next+a_n1_orig + a_n2_orig +1:work_next+a_n) = 0
      
      DO i = a_n1_orig + a_n2_orig +1,a_n
        j = partition(i)
        IF (work(work_part+j) .EQ. mc70_sep_flag) THEN
          ! j is not on the boundary
   !       write(*,*) j, 'not boundary'
          IF (a_weight_1 .LT. a_weight_2) THEN
            ! Move j into partition 1
            work(work_part+j) = mc70_part1_flag
            a_n1 = a_n1 + 1
            a_weight_1 = a_weight_1 + a_weight(j)
            a_weight_sep = a_weight_sep - a_weight(j)
            
            head1 = j
            tail1 = j
            DO WHILE (head1 .GT. 0)
             q = head1
             IF (q .LT. a_n) THEN
              k = a_ptr(q+1)-1
             ELSE
              k = a_ne
             END IF
             DO l = a_ptr(q),k
               p = a_row(l)
               IF (work(work_part+p) .EQ. mc70_sep_flag) THEN
                 work(work_part+p) = mc70_part1_flag
                 a_n1 = a_n1 + 1
                 a_weight_1 = a_weight_1 + a_weight(q)
                 a_weight_sep = a_weight_sep - a_weight(q)
                 work(work_next+tail1) = p
                 tail1 = p
               END IF
             END DO
             IF (head1 .EQ. tail1) THEN
               head1 = 0
               tail1 = 0
             ELSE
               head1 = work(work_next+q)
               work(work_next+q) = 0
             END IF

            END DO

          ELSE
            ! Move j into partition 2
            work(work_part+j) = mc70_part2_flag
            a_n2 = a_n2 + 1
            a_weight_2 = a_weight_2 + a_weight(j)
            a_weight_sep = a_weight_sep - a_weight(j)
            head2 = j
            tail2 = j

            DO WHILE (head2 .GT. 0)
             q = head2
             IF (q .LT. a_n) THEN
              k = a_ptr(q+1)-1
             ELSE
              k = a_ne
             END IF
             DO l = a_ptr(q),k
               p = a_row(l)
               IF (work(work_part+p) .EQ. mc70_sep_flag) THEN
                 work(work_part+p) = mc70_part2_flag
                 a_n2 = a_n2 + 1
                 a_weight_2 = a_weight_2 + a_weight(q)
                 a_weight_sep = a_weight_sep - a_weight(q)
                 work(work_next+tail2) = p
                 tail2 = p
               END IF
             END DO
             IF (head2 .EQ. tail2) THEN
               head2 = 0
               tail2 = 0
             ELSE
               head2 = work(work_next+q)
               work(work_next+q) = 0
             END IF
            END DO
          END IF
        END IF
      END DO
      
      ! Reset partition matrix
      i = 1
      j = a_n1+1
      k = a_n1+a_n2+1
      DO l = 1,a_n
        m = work(work_part+l)
        SELECT CASE (m)
        CASE (mc70_part1_flag)
          partition(i) = l
          i = i+1
        CASE (mc70_part2_flag)
          partition(j) = l
          j = j+1
        CASE default
          partition(k) = l
          work(work_part+l) = mc70_sep_flag
          k = k+1
        END SELECT
      END DO
      a_weight_sep = sumweight - a_weight_1 - a_weight_2
     ! call check_partition1(a_n,a_ne,a_ptr,a_row,a_n1,a_n2,work(work_part+1:work_part+a_n),a_weight_1,a_weight_2,a_weight)

      END SUBROUTINE mc70_refine_block_trim

      SUBROUTINE expand_partition(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,&
               a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,control)
        INTEGER, INTENT(IN) :: a_n ! order of matrix
        INTEGER, INTENT(IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT(IN) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. 
        INTEGER, INTENT(IN) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.
        INTEGER, INTENT(IN) :: a_weight(a_n) ! On input a_weight(i) contains 
             ! the weight of column i 
        INTEGER, INTENT(IN) :: sumweight ! Sum of weights in a_weight
        INTEGER, INTENT(INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT(INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT(INOUT) :: a_weight_1,a_weight_2,a_weight_sep ! Weighted 
             ! size of partitions and separator
        INTEGER, INTENT(INOUT) :: partition(a_n) !First a_n1 entries contain
             ! list of (local) indices in partition 1; next a_n2 entries  
             ! contain list of (local) entries in partition 2; entries in 
             ! separator are listed at the end. This is updated to the new 
             ! partition
        INTEGER, INTENT(OUT) :: work(a_n) ! Work array
        TYPE(mc70_control), INTENT(IN) :: control

        ! Local variables
        INTEGER :: i,j,k,l,m,w,side
        INTEGER :: work_part
        ! Set work(work_part+1:work_part+a_n) to hold flags to indicate what 
        ! part of the partition the nodes are in
        work_part = 0
        DO i = 1, a_n1
            j = partition(i)
            work(work_part+j) = mc70_part1_flag
        END DO
        DO i = a_n1+1, a_n1+a_n2
            j = partition(i)
            work(work_part+j) = mc70_part2_flag
        END DO
        DO i = a_n1+a_n2+1, a_n
            j = partition(i)
            work(work_part+j) = mc70_sep_flag
        END DO

        IF (a_n1 .LT. a_n2) THEN
          side = mc70_part2_flag
        ELSE IF (a_n1 .GT. a_n2) THEN
          side = mc70_part1_flag
        ELSE
          side = mc70_sep_flag
        END IF

        DO i = a_n1+a_n2+1, a_n
            j = partition(i)
            ! search neighbours of j and add to separator
            IF (j.EQ.a_n) THEN
               k = a_ne
            ELSE
               k = a_ptr(j+1)-j
            END IF
            DO l=a_ptr(j),k
               m = a_row(l)
               IF (work(work_part+m) .EQ. mc70_part1_flag .AND. a_n1.GT.1) THEN
                 IF (side .EQ. mc70_part1_flag .OR. side .EQ. mc70_sep_flag) THEN
                  work(work_part+m) = mc70_sep_flag
                  a_n1 = a_n1 - 1
                  w = a_weight(m)
                  a_weight_1 = a_weight_1 - w
                  a_weight_sep = a_weight_sep + w
                 END IF
               ELSE IF (work(work_part+m) .EQ. mc70_part2_flag .AND. a_n2.GT.1) THEN
                 IF (side .EQ. mc70_part2_flag .OR. side .EQ. mc70_sep_flag) THEN
                  work(work_part+m) = mc70_sep_flag
                  a_n2 = a_n2 - 1
                  w = a_weight(m)
                  a_weight_2 = a_weight_2 - w
                  a_weight_sep = a_weight_sep + w
                 END IF
               END IF
            END DO
        END DO
        j =1
        k = j + a_n1
        l = k + a_n2
        DO i=1,a_n
          m = work(work_part+i)
          SELECT CASE (m)
          CASE (mc70_part1_flag)
             partition(j) =i
             j = j+1
          CASE (mc70_part2_flag) 
             partition(k) =i
             k = k+1
          CASE default
             partition(l) =i
             l = l+1
          END SELECT
        END DO
        
      END SUBROUTINE expand_partition

      SUBROUTINE expand_partition_simple(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,&
               a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,control)
        INTEGER, INTENT(IN) :: a_n ! order of matrix
        INTEGER, INTENT(IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT(IN) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. 
        INTEGER, INTENT(IN) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.
        INTEGER, INTENT(IN) :: a_weight(a_n) ! On input a_weight(i) contains 
             ! the weight of column i 
        INTEGER, INTENT(IN) :: sumweight ! Sum of weights in a_weight
        INTEGER, INTENT(INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT(INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT(INOUT) :: a_weight_1,a_weight_2,a_weight_sep ! Weighted 
             ! size of partitions and separator
        INTEGER, INTENT(INOUT) :: partition(a_n) !First a_n1 entries contain
             ! list of (local) indices in partition 1; next a_n2 entries  
             ! contain list of (local) entries in partition 2; entries in 
             ! separator are listed at the end. This is updated to the new 
             ! partition
        INTEGER, INTENT(OUT) :: work(a_n) ! Work array
        TYPE(mc70_control), INTENT(IN) :: control

        ! Local variables
        INTEGER :: i,j,k,l,m,w,side,t,s
        INTEGER :: work_part
        ! Set work(work_part+1:work_part+a_n) to hold flags to indicate what 
        ! part of the partition the nodes are in
        work_part = 0
        DO i = 1, a_n1
            j = partition(i)
            work(work_part+j) = mc70_part1_flag
        END DO
        DO i = a_n1+1, a_n1+a_n2
            j = partition(i)
            work(work_part+j) = mc70_part2_flag
        END DO
        DO i = a_n1+a_n2+1, a_n
            j = partition(i)
            work(work_part+j) = mc70_sep_flag
        END DO

    !    s = 0
        DO i = a_n1+a_n2+1, a_n
           j = partition(i)
           IF (j.EQ. a_n) THEN
              t = a_ne
           ELSE
              t = a_ptr(j+1)-1
           END IF
           DO k = a_ptr(j),t
              l = a_row(k)
              IF ((work(work_part+l).NE.mc70_sep_flag) ) THEN
               work(work_part+l) = mc70_sep_flag
                  ! Add l to list
           !       s=s+a_weight(l)
              END IF
           END DO
        END DO
        
        ! Update a_n1, a_n2, a_weight1, a_weight_2,a_weight_sep
        a_n1 = 0
        a_n2 = 0
        a_weight_1 = 0
        a_weight_2 = 0
        a_weight_sep = 0

        DO i = 1,a_n
          j = work(work_part+i)
          SELECT CASE (j)
          CASE (mc70_part1_flag)
            a_n1 = a_n1 + 1
            a_weight_1 = a_weight_1 + a_weight(i)
          CASE (mc70_part2_flag)
            a_n2 = a_n2 + 1
            a_weight_2 = a_weight_2 + a_weight(i)
          CASE (mc70_sep_flag)
            a_weight_sep = a_weight_sep + a_weight(i)
          END SELECT

        END DO

        j =1
        k = j + a_n1
        l = k + a_n2
        DO i=1,a_n
          m = work(work_part+i)
          SELECT CASE (m)
          CASE (mc70_part1_flag)
             partition(j) =i
             j = j+1
          CASE (mc70_part2_flag) 
             partition(k) =i
             k = k+1
          CASE default
             partition(l) =i
             l = l+1
          END SELECT
        END DO
        
      END SUBROUTINE expand_partition_simple

      SUBROUTINE expand_partition_kinks(a_n,a_ne,a_ptr,a_row,a_weight,&
               sumweight,radius,upper,ratio,a_n1,&
               a_n2,a_weight_1,a_weight_2,a_weight_sep,partition,work,control)
        INTEGER, INTENT(IN) :: a_n ! order of matrix
        INTEGER, INTENT(IN) :: a_ne ! number of entries in matrix
        INTEGER, INTENT(IN) :: a_ptr(a_n) ! On input a_ptr(i) contains 
             ! position in a_row that entries for column i start. 
        INTEGER, INTENT(IN) :: a_row(a_ne) ! On input a_row contains row 
             ! indices of the non-zero rows. Diagonal entries have been removed
             ! and the matrix expanded.
        INTEGER, INTENT(IN) :: a_weight(a_n) ! On input a_weight(i) contains 
             ! the weight of column i 
        INTEGER, INTENT(IN) :: sumweight ! Sum of weights in a_weight
        INTEGER, INTENT(IN) :: radius ! Check nodes at most distance radius 
                                ! from current node. Ball = these nodes
        REAL (myreal_mc70) :: ratio ! Current node in partition i.
             ! If |Ball \cap P_j|/|Ball \cap P_i| > ratio, current node will be 
             !    moved into expanded matrix
        REAL (myreal_mc70) :: upper ! Weight of expanded separator must not 
             ! exceed upper*(weight of initial separator)
        INTEGER, INTENT(INOUT) :: a_n1 ! Size of partition 1
        INTEGER, INTENT(INOUT) :: a_n2 ! Size of partition 2
        INTEGER, INTENT(INOUT) :: a_weight_1,a_weight_2,a_weight_sep ! Weighted 
             ! size of partitions and separator
        INTEGER, INTENT(INOUT) :: partition(a_n) !First a_n1 entries contain
             ! list of (local) indices in partition 1; next a_n2 entries  
             ! contain list of (local) entries in partition 2; entries in 
             ! separator are listed at the end. This is updated to the new 
             ! partition
        INTEGER, INTENT(OUT) :: work(5*a_n) ! Work array
        TYPE(mc70_control), INTENT(IN) :: control

        ! Local variables
        INTEGER :: i,j,jj,head,tail,headb,tailb,k,l,t,s
        INTEGER :: no_part1,no_part2,no_sep,w_sep_orig
        INTEGER :: work_part,work_next,work_mask,work_dist,work_nextb
        LOGICAL :: move

        ! Divide up workspace
        work_part = 0
        work_next = work_part + a_n
        work_mask = work_next + a_n
        work_dist = work_mask + a_n
        work_nextb = work_dist + a_n

        ! Set work(work_part+1:work_part+a_n) to hold flags to indicate what 
        ! part of the partition the nodes are in
        work_part = 0
        DO i = 1, a_n1
            j = partition(i)
            work(work_part+j) = mc70_part1_flag
        END DO
        DO i = a_n1+1, a_n1+a_n2
            j = partition(i)
            work(work_part+j) = mc70_part2_flag
        END DO
        DO i = a_n1+a_n2+1, a_n
            j = partition(i)
            work(work_part+j) = mc70_sep_flag
        END DO
        w_sep_orig = a_weight_sep

        ! Work through separator adding adjacent entries to list
        head = 0
        tail = 0
        work(work_next+1:work_next+a_n) = 0
        s=0
        DO i = a_n1+a_n2+1, a_n
           j = partition(i)
           IF (j.EQ. a_n) THEN
              t = a_ne
           ELSE
              t = a_ptr(j+1)-1
           END IF
           DO k = a_ptr(j),t
              l = a_row(k)
              IF ((work(work_part+l).NE.mc70_sep_flag) .AND. &
                  (work(work_next+l).EQ.0)) THEN
                  ! Add l to list
                  s=s+1
                  IF (tail.EQ.0) THEN
                     head = l
                     tail = l
                  ELSE
                     work(work_next+tail) = l
                     tail = l
                  END IF
              END IF
           END DO
        END DO
        work(work_nextb+1:work_nextb+a_n) = 0
        work(work_dist+1:work_dist+a_n) = -1
        work(work_mask+1:work_mask+a_n) = 0
        DO WHILE (head.NE.0)
          i = head
          IF ((work(work_part+i).EQ.mc70_part1_flag) .AND. a_n1.EQ.1) THEN
            GOTO 100
          END IF
          IF ((work(work_part+i).EQ.mc70_part2_flag) .AND. a_n2.EQ.1) THEN
            GOTO 100
          END IF
          IF ( real(a_weight_sep + a_weight(i)) .GT. upper*real(w_sep_orig) ) THEN
            GOTO 100
          END IF

          work(work_dist+i) = 0
          headb = 0
          tailb = 0
          ! create list with all entries in Ball
          IF (i.EQ.a_n) THEN
             jj = a_ne
          ELSE
             jj = a_ptr(i+1)-1
          END IF

          DO j = a_ptr(i),jj
             l = a_row(j)
             ! Add l to Ball list
             IF (tail.EQ.0) THEN
                 headb = l
                 tailb = l
             ELSE
                 work(work_nextb+tailb) = l
                 tailb = l
             END IF
             work(work_dist+l) = 1
          END DO

          k = headb
          DO WHILE (k.NE.0 .AND. work(work_dist+k).LT.radius)
            IF (k.EQ.a_n) THEN
              jj = a_ne
            ELSE
              jj = a_ptr(k+1)-1
            END IF
            DO j = a_ptr(i),jj
             l = a_row(j)
             IF (work(work_dist+l).EQ.-1) THEN
               ! Add l to Ball list
               IF (tail.EQ.0) THEN
                 headb = l
                 tailb = l
               ELSE
                 work(work_nextb+tailb) = l
                 tailb = l
               END IF
               work(work_dist+l) = work(work_dist+k)+1
             END IF
            END DO
            k = work(work_nextb+k)
          END DO


          ! count entries in ball that are in partition 1 and partition 2
          ! empty list as proceed
          no_part1=0
          no_part2=0
          no_sep = 0
          DO WHILE (headb.NE.0)
            k = headb
            IF (work(work_part+k).GE.mc70_part2_flag) THEN
              no_part2=no_part2+a_weight(k)             
            ELSE
              IF (work(work_part+k).LE.mc70_part1_flag) THEN
                no_part1=no_part1+a_weight(k)  
              ELSE
                no_sep = no_sep + a_weight(k)  
              END IF
            END IF
            headb = work(work_nextb+k)
            work(work_dist+k) = -1
            work(work_nextb+k) = 0
            k = headb
          END DO

          move = .FALSE.
          IF (work(work_part+i).EQ.mc70_part2_flag) THEN
             !If |Ball \cap (P_1 \cup S)|/|Ball \cap P_2| > ratio
             IF (no_part2.EQ.0) THEN
               move = .TRUE.
               a_n2 = a_n2 - 1
               a_weight_2 = a_weight_2 - a_weight(i)
               a_weight_sep = a_weight_sep + a_weight(i)
             ELSE
                IF (real(no_part1+no_sep)/real(no_part2).GT.ratio) THEN
                 move = .TRUE.
                 a_n2 = a_n2 - 1
               a_weight_2 = a_weight_2 - a_weight(i)
               a_weight_sep = a_weight_sep + a_weight(i)
                END IF
             END IF
          ELSE
             !If |Ball \cap P_2|/|Ball \cap P_1| > ratio
             IF (no_part1.EQ.0) THEN
               move = .TRUE.
               a_n1 = a_n1 - 1
               a_weight_1 = a_weight_1 - a_weight(i)
               a_weight_sep = a_weight_sep + a_weight(i)
             ELSE
                IF (real(no_part2+no_sep)/real(no_part1).GT.ratio) THEN
                 move = .TRUE.
                 a_n1 = a_n1 - 1
               a_weight_1 = a_weight_1 - a_weight(i)
               a_weight_sep = a_weight_sep + a_weight(i)
                END IF
             END IF           
          END IF    
          IF (move) THEN
            k = work(work_part+i)
            work(work_part+i) = (k - 1)**2
            j = i
            IF (j.EQ. a_n) THEN
              t = a_ne
            ELSE
              t = a_ptr(j+1)-1
            END IF
            DO k = a_ptr(j),t
              l = a_row(k)
              IF ((work(work_part+l).NE.mc70_sep_flag) .AND. &
                  (work(work_next+l).EQ.0).AND. work(work_mask+l).EQ.0) THEN
                  ! Add l to list (note: list is non-empty)
                  work(work_next+tail) = l
                  tail = l
              END IF
            END DO
          END IF       

          ! remove i from list
100       head = work(work_next+i)
          work(work_mask+i) = 1
          work(work_next+i) = 0
          IF (head.EQ.0) THEN
            tail = 0
          END IF
          work(work_dist+i) = -1
        END DO


        j =1
        k = j + a_n1
        l = k + a_n2 
        t=0
        DO i=1,a_n
          jj = work(work_part+i)
          IF (jj .EQ. mc70_part1_flag) THEN
             partition(j) =i
             j = j+1
          ELSE  IF (jj .EQ. mc70_part2_flag) THEN
             partition(k) =i
             k = k+1
          ELSE
             partition(l) =i
             l = l+1
             t = t + a_weight(i)
          END IF
        END DO
        
      END SUBROUTINE expand_partition_kinks


      SUBROUTINE cost_function(a_weight_1,a_weight_2,a_weight_sep,sumweight,&
          ratio,imbal,tau)
       
        INTEGER, INTENT(IN) :: a_weight_1,a_weight_2,a_weight_sep ! Weighted 
             ! size of partitions and separator
        INTEGER, INTENT(IN) :: sumweight
        REAL(myreal_mc70), INTENT(IN) :: ratio
        LOGICAL, INTENT(IN) :: imbal ! Use penalty function?
        REAL(myreal_mc70), INTENT(OUT) :: tau
        REAL(myreal_mc70) :: beta
        beta = 0.000375

       IF (.true.) THEN
        IF (imbal .AND. real(max(a_weight_1,a_weight_2))/ &
                   real(min(a_weight_1,a_weight_2)) .GE. ratio ) THEN
           tau = real(sumweight-2) + real(max(a_weight_1,a_weight_2))/ &
                   real(min(a_weight_1,a_weight_2))
        ELSE
           tau = ((real(a_weight_sep)**1.0)/real(a_weight_1))/real(a_weight_2)
        END IF
      ELSE
        IF (imbal .AND. real(max(a_weight_1,a_weight_2))/ &
                   real(min(a_weight_1,a_weight_2)) .GE. ratio ) THEN
           tau = real(sumweight)*(1.0_myreal_mc70 + beta*real(sumweight)) 
        ELSE
           tau = real(a_weight_sep)*(1.0_myreal_mc70 + &
                        beta*real(abs(a_weight_1 - a_weight_2))  )
        END IF

      END IF
      !  IF (tau<0.0) THEN
      !    write(*,*) a_weight_1,a_weight_2,a_weight_sep
      !  END IF

      END SUBROUTINE cost_function

END MODULE hsl_mc70_double
