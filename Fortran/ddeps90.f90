! COPYRIGHT (c) 2011 Science and Technology Facilities Council
! Original date 24 January 2010. Version 1.0.0
!
! Written by: Sue Thorne
!
! Version 1.1.0
! See ChangeLog for version history
!

    MODULE hsl_mc79_integer

      IMPLICIT NONE
      PRIVATE

! ---------------------------------------------------
! Precision
! ---------------------------------------------------
      INTEGER, PARAMETER :: myint_mc79 = kind(1)

! ---------------------------------------------------
! Error flags
! ---------------------------------------------------
      INTEGER (myint_mc79), PARAMETER :: mc79_err_memory_alloc = -1, & ! memory
! alloc
! error
        mc79_err_memory_dealloc = -2, & ! memory dealloc error
        mc79_err_m = -3, & ! m<1
        mc79_err_n = -4 ! n<1


! ---------------------------------------------------
! Derived type definitions
! ---------------------------------------------------

      TYPE, PUBLIC :: mc79_control
        INTEGER :: lp = 6 ! stream number for error messages
        INTEGER :: wp = 6 ! stream number for warning messages
        INTEGER :: mp = 6 ! stream number for diagnostic messages
        INTEGER :: print_level = 0 ! amount of informational output required
      END TYPE mc79_control

      TYPE, PUBLIC :: mc79_info
        INTEGER :: flag = 0 ! error/warning flag
        INTEGER :: hz_comps = -1 ! number of horizontal comps after fine perm
        INTEGER :: vt_comps = -1 ! number of vertical comps after fine perm
        INTEGER :: sq_comps = -1 ! number of comps in square after fine perm
        INTEGER :: m1 = 0 ! number of rows in R1
        INTEGER :: m2 = 0 ! number of rows in R2
        INTEGER :: m3 = 0 ! number of rows in R3
        INTEGER :: mbar = 0 ! number of rows in Rbar
        INTEGER :: n1 = 0 ! number of columns in C1
        INTEGER :: n2 = 0 ! number of columns in C2
        INTEGER :: n3 = 0 ! number of columns in C3
        INTEGER :: nbar = 0 ! number of columns in Cbar
        INTEGER :: stat = 0 ! holds Fortran stat parameter

      END TYPE mc79_info

      PUBLIC mc79_matching, mc79_coarse, mc79_fine

    CONTAINS

! ---------------------------------------------------------------

      SUBROUTINE mc79_matching(m,n,ptr,row,rowmatch,colmatch,control,info)
! subroutine mc79_matching(m,n,ptr,row,rowmatch,colmatch,control,info)
!  constructs a matching for a matrix A. The matrix is stored using CSC format.

! m: is an integer with intent(in). It holds the number of rows in A
        INTEGER (myint_mc79), INTENT (IN) :: m

! n: is an integer with intent(in). It holds the number of columns in A
        INTEGER (myint_mc79), INTENT (IN) :: n

! ptr: holds pointer for A
        INTEGER (myint_mc79), INTENT (IN) :: ptr(n+1)

! row: holds indices of A
        INTEGER (myint_mc79), INTENT (IN) :: row(:)

! rowmatch: is an integer array with intent(out) of size m. It holds
! the row matching
        INTEGER (myint_mc79), INTENT (OUT) :: rowmatch(m)


! colmatch: is an integer array with intent(out) of size m. It holds
! the column matching
        INTEGER (myint_mc79), INTENT (OUT) :: colmatch(n)

! control: is of derived type mc79_control with intent(in). Controls
! action
        TYPE (mc79_control), INTENT (IN) :: control

! info: is of derived type mc79_info with intent(out).
! info%flag
! = 0 if successful
! = mc79_ERR_MEMORY_ALLOC if memory allocation failed
! = mc79_ERR_MEMORY_DEALLOC if memory deallocation failed
! = mc79_ERR_M if m<1
! = mc79_ERR_N if n<1
        TYPE (mc79_info), INTENT (OUT) :: info


! ---------------------------------------------
! Local allocatables
! ---------------------------------------------
        INTEGER, ALLOCATABLE, DIMENSION (:) :: ptr_temp, row_temp

! ---------------------------------------------
! Local variables
! ---------------------------------------------

        LOGICAL :: printe ! errors to be printed?
        LOGICAL :: printi ! basic diagnostic to be printed?
        LOGICAL :: printd ! additional diagnostic to be printed?
        INTEGER :: lp ! stream for printing errors
        INTEGER :: mp ! stream for printing diagnostic info
        INTEGER :: ne ! number of nonzero entries in A
        INTEGER :: i, j, l

        info%flag = 0
! ---------------------------------------------
! Set stream numbers
! ---------------------------------------------
        lp = control%lp
        mp = control%mp

! ---------------------------------------------
! Printing levels
! ---------------------------------------------
        printe = (control%print_level>=0 .AND. lp>=0)
        printi = (control%print_level>=1 .AND. mp>=0)
        printd = (control%print_level>=2 .AND. mp>=0)

        IF (printi) THEN
          WRITE (mp,'(a)') ' '
          WRITE (mp,'(a15)') 'MC79_matching:'
        END IF

        info%flag = 0

! Check that m, n and ne have valid entries
        IF (m<1) THEN
          info%flag = mc79_err_m
          GO TO 300
        END IF

        IF (n<1) THEN
          info%flag = mc79_err_n
          GO TO 300
        END IF

        IF (printi) THEN
          WRITE (mp,'(a4,i15)') 'm = ', m
          WRITE (mp,'(a4,i15)') 'n = ', n
        END IF

! Build temp = A^T if necessary
        ne = ptr(n+1) - 1
        IF (m<n) THEN

          IF (printd) THEN
            WRITE (mp,'(a35)') 'm smaller n: build transpose of A'
          END IF

          ALLOCATE (row_temp(ne),ptr_temp(m+1),STAT=info%stat)
          IF (info%stat>0) GO TO 100

          ptr_temp(1:m) = 0
          ptr_temp(m+1) = ne + 1

          DO i = 1, ne
            j = row(i)
            ptr_temp(j) = ptr_temp(j) + 1
          END DO

          ptr_temp(1) = ptr_temp(1) + 1
          DO i = 2, m
            ptr_temp(i) = ptr_temp(i) + ptr_temp(i-1)
          END DO

          DO i = 1, n
            DO j = ptr(i), ptr(i+1) - 1
              l = row(j)
              ptr_temp(l) = ptr_temp(l) - 1
              row_temp(ptr_temp(l)) = i
            END DO
          END DO

        END IF

! Find the maximum matching
        IF (m>=n) THEN
          CALL find_max_match(m,n,ne,ptr,row,rowmatch,colmatch,info%flag, &
            info%stat)
        ELSE
          CALL find_max_match(n,m,ne,ptr_temp,row_temp,colmatch,rowmatch, &
            info%flag,info%stat)
        END IF
        IF (info%flag<0) GO TO 300

! Run through and count number of unmatched rows and columns
        info%mbar = 0
        DO i = 1, m
          IF (rowmatch(i)==0) info%mbar = info%mbar + 1
        END DO

        info%nbar = 0
        DO i = 1, n
          IF (colmatch(i)==0) info%nbar = info%nbar + 1
        END DO

        IF (m<n) THEN
          DEALLOCATE (row_temp,ptr_temp,STAT=info%stat)
          IF (info%stat>0) GO TO 200
        END IF

        IF (printd) THEN
          CALL mc79_print_message(info%flag,lp,context='mc79_matching')
! ---------------------------------------------
! Print out perm
! ---------------------------------------------
          WRITE (mp,'(a10)') 'rowmatch = '
          WRITE (mp,'(5i15)') (rowmatch(i),i=1,m)
          WRITE (mp,'(a10)') 'colmatch = '
          WRITE (mp,'(5i15)') (colmatch(i),i=1,n)
        ELSE IF (printi) THEN
! ---------------------------------------------
! Print out first few entries of perm
! ---------------------------------------------
          WRITE (mp,'(a10)') 'rowmatch(1:min(5,m)) = '
          WRITE (mp,'(5i15)') (rowmatch(i),i=1,min(5,m))
          WRITE (mp,'(a10)') 'colmatch(1:min(5,n)) = '
          WRITE (mp,'(5i15)') (colmatch(i),i=1,min(5,n))
        END IF

        RETURN

100     info%flag = mc79_err_memory_alloc
        IF (printe) CALL mc79_print_message(info%flag,lp, &
          context='mc79_matching')

200     info%flag = mc79_err_memory_dealloc
        IF (printe) CALL mc79_print_message(info%flag,lp, &
          context='mc79_matching')
        RETURN

300     IF (printe) CALL mc79_print_message(info%flag,lp, &
          context='mc79_matching')
        RETURN


      END SUBROUTINE mc79_matching

! ---------------------------------------------------------------

      SUBROUTINE mc79_coarse(m,n,ptr,row,rowperm,colperm,control,info)
! subroutine mc79_matching(m,n,ne,lptr,ptr,idx,csr,rowperm,colperm,control,info)
! constructs a matching for a matrix A. If (csr), then A is stored in CSR 
! format; otherwise, A is stored in CSC format.

! m: is an integer with intent(in). It holds the number of rows in A
        INTEGER (myint_mc79), INTENT (IN) :: m

! n: is an integer with intent(in). It holds the number of columns in B
        INTEGER (myint_mc79), INTENT (IN) :: n

! ptr: holds row pointer for B
        INTEGER (myint_mc79), INTENT (IN) :: ptr(n+1)

! row: holds indices of B
        INTEGER (myint_mc79), INTENT (IN) :: row(:)


! rowperm: is an integer array with intent(out) of size m. It holds
! the row permutation. If rowperm(i) = j, row i in the original matrix 
! becomes row j in the perumted matrix 
        INTEGER (myint_mc79), INTENT (OUT) :: rowperm(m)


! colperm: is an integer array with intent(out) of size n. It holds
! the column permutation. If colperm(i) = j, column i in the original 
! matrix becomes column j in the perumted matrix 
        INTEGER (myint_mc79), INTENT (OUT) :: colperm(n)

! control: is of derived type mc79_control with intent(in). Controls
! action
        TYPE (mc79_control), INTENT (IN) :: control

! info: is of derived type mc79_info with intent(out).
! info%flag
! = 0 if successful
! = mc79_ERR_MEMORY_ALLOC if memory allocation failed
! = mc79_ERR_MEMORY_DEALLOC if memory deallocation failed
! = mc79_ERR_M if m<1
! = mc79_ERR_N if n<1
        TYPE (mc79_info), INTENT (OUT) :: info


! ---------------------------------------------
! Local allocatables
! ---------------------------------------------
        INTEGER, ALLOCATABLE, DIMENSION (:) :: ptr_temp, row_temp, rowmatch, &
          colmatch, col_markers, row_markers

! ---------------------------------------------
! Local variables
! ---------------------------------------------

        LOGICAL :: printe ! errors to be printed?
        LOGICAL :: printi ! info to be printed?
        LOGICAL :: printd ! extended info to be printed?
        INTEGER :: lp ! stream for printing errors
        INTEGER :: mp ! stream for printing info
        INTEGER :: i, nhrows, nhcols, nvrows, nvcols, vindex, sqindex, hindex, &
          j, l, ne, in1, in2, in3, inbar

! ---------------------------------------------
! Set stream numbers
! ---------------------------------------------
        lp = control%lp
        mp = control%mp

! ---------------------------------------------
! Printing levels
! ---------------------------------------------
        printe = (control%print_level>=0 .AND. lp>=0)
        printi = (control%print_level>=1 .AND. mp>=0)
        printd = (control%print_level>=2 .AND. mp>=0)

        info%flag = 0
        IF (printi) THEN
          WRITE (mp,'(a)') ' '
          WRITE (mp,'(a15)') 'MC79_coarse:'
        END IF

! Check that m, n and ne have valid entries
        IF (m<1) THEN
          info%flag = mc79_err_m
          GO TO 300
        END IF

        IF (n<1) THEN
          info%flag = mc79_err_n
          GO TO 300
        END IF

        IF (printi) THEN
          WRITE (mp,'(a4,i15)') 'm = ', m
          WRITE (mp,'(a4,i15)') 'n = ', n
        END IF

! Check whether there are any nonzero entries
        IF (ptr(n+1)==1) THEN
          info%m1 = 0
          info%n1 = n
          info%m2 = 0
          info%n2 = 0
          info%m3 = m
          info%n3 = 0
          info%mbar = info%m3 - info%n3
          info%nbar = info%n1 - info%m1
          rowperm = (/ (i,i=1,m) /)
          colperm = (/ (i,i=1,n) /)
          RETURN
        END IF

        vindex = -1
        sqindex = -2
        hindex = -3

! Build temp = A^T 
        ne = ptr(n+1) - 1

        ALLOCATE (row_temp(ne),ptr_temp(m+1),STAT=info%stat)
        IF (info%stat>0) GO TO 100
        ptr_temp(1:m) = 0
        ptr_temp(m+1) = ne + 1

        DO i = 1, ne
          j = row(i)
          ptr_temp(j) = ptr_temp(j) + 1
        END DO

        ptr_temp(1) = ptr_temp(1) + 1
        DO i = 2, m
          ptr_temp(i) = ptr_temp(i) + ptr_temp(i-1)
        END DO

        DO i = 1, n
          DO j = ptr(i), ptr(i+1) - 1
            l = row(j)
            ptr_temp(l) = ptr_temp(l) - 1
            row_temp(ptr_temp(l)) = i
          END DO
        END DO

        ALLOCATE (rowmatch(m),colmatch(n),STAT=info%stat)
        IF (info%stat>0) GO TO 100

! Find the maximum matching

        IF (m>=n) THEN
          CALL find_max_match(m,n,ne,ptr,row,rowmatch,colmatch,info%flag, &
            info%stat)
        ELSE
          CALL find_max_match(n,m,ne,ptr_temp,row_temp,colmatch,rowmatch, &
            info%flag,info%stat)
        END IF
        IF (info%flag<0) GO TO 300

        ALLOCATE (row_markers(m),col_markers(n),STAT=info%stat)
        IF (info%stat>0) GO TO 100

        row_markers(:) = sqindex
        col_markers(:) = sqindex

        CALL find_rect(n,m,ne,ptr_temp,row_temp,vindex,sqindex,rowmatch, &
          colmatch,row_markers,col_markers,nvcols,nvrows,info%flag,info%stat)
        IF (info%flag<0) GO TO 300
        CALL find_rect(m,n,ne,ptr,row,hindex,sqindex,colmatch,rowmatch, &
          col_markers,row_markers,nhrows,nhcols,info%flag,info%stat)
        IF (info%flag<0) GO TO 300

        info%m1 = nhrows
        info%n1 = nhcols
        info%m2 = m - nhrows - nvrows
        info%n2 = info%m2
        info%m3 = nvrows
        info%n3 = nvcols
        info%mbar = info%m3 - info%n3
        info%nbar = info%n1 - info%m1

! Run through rows and count how many not matched. At same time, place
!  in rowperm

        in1 = 1
        in2 = in1 + info%m1
        in3 = in2 + info%m2
        inbar = m

        info%mbar = 0
        DO i = 1, m
          IF (rowmatch(i)==0) THEN
            rowperm(inbar) = i
            inbar = inbar - 1
          ELSE IF (row_markers(i)==hindex) THEN
            rowperm(in1) = i
            in1 = in1 + 1

          ELSE IF (row_markers(i)==sqindex) THEN
            rowperm(in2) = i
            in2 = in2 + 1
          ELSE
            rowperm(in3) = i
            in3 = in3 + 1
          END IF
        END DO

! Run through columns and place in colperm
        inbar = 1
        in1 = inbar + info%nbar
        in2 = inbar + info%n1
        in3 = in2 + info%n2

        DO i = 1, n
          IF (colmatch(i)==0) THEN
            colperm(inbar) = i
            inbar = inbar + 1
          ELSE IF (col_markers(i)==hindex) THEN
            colperm(in1) = i
            in1 = in1 + 1

          ELSE IF (col_markers(i)==sqindex) THEN
            colperm(in2) = i
            in2 = in2 + 1
          ELSE
            colperm(in3) = i
            in3 = in3 + 1
          END IF

        END DO
        info%hz_comps = -1
        info%vt_comps = -1
        info%sq_comps = -1

        DEALLOCATE (rowmatch,colmatch,row_markers,col_markers,ptr_temp, &
          row_temp,STAT=info%stat)
        IF (info%stat>0) GO TO 200
        IF (printd) THEN
          CALL mc79_print_message(info%flag,lp,context='mc79_coarse')
! ---------------------------------------------
! Print out perm
! ---------------------------------------------
          WRITE (mp,'(a10)') 'info%m1 = '
          WRITE (mp,'(i15)') info%m1
          WRITE (mp,'(a10)') 'info%m2 = '
          WRITE (mp,'(i15)') info%m2
          WRITE (mp,'(a10)') 'info%m3 = '
          WRITE (mp,'(i15)') info%m3
          WRITE (mp,'(a10)') 'info%mbar = '
          WRITE (mp,'(i15)') info%mbar
          WRITE (mp,'(a10)') 'info%n1 = '
          WRITE (mp,'(i15)') info%n1
          WRITE (mp,'(a10)') 'info%n2 = '
          WRITE (mp,'(i15)') info%n2
          WRITE (mp,'(a10)') 'info%n3 = '
          WRITE (mp,'(i15)') info%n3
          WRITE (mp,'(a10)') 'info%nbar = '
          WRITE (mp,'(i15)') info%nbar
          WRITE (mp,'(a10)') 'rowperm = '
          WRITE (mp,'(5i15)') (rowperm(i),i=1,m)
          WRITE (mp,'(a10)') 'colperm = '
          WRITE (mp,'(5i15)') (colperm(i),i=1,n)
        ELSE IF (printi) THEN
! ---------------------------------------------
! Print out first few entries of perm
! ---------------------------------------------
          WRITE (mp,'(a10)') 'info%m1 = '
          WRITE (mp,'(i15)') info%m1
          WRITE (mp,'(a10)') 'info%m2 = '
          WRITE (mp,'(i15)') info%m2
          WRITE (mp,'(a10)') 'info%m3 = '
          WRITE (mp,'(i15)') info%m3
          WRITE (mp,'(a10)') 'info%mbar = '
          WRITE (mp,'(i15)') info%mbar
          WRITE (mp,'(a10)') 'info%n1 = '
          WRITE (mp,'(i15)') info%n1
          WRITE (mp,'(a10)') 'info%n2 = '
          WRITE (mp,'(i15)') info%n2
          WRITE (mp,'(a10)') 'info%n3 = '
          WRITE (mp,'(i15)') info%n3
          WRITE (mp,'(a10)') 'info%nbar = '
          WRITE (mp,'(i15)') info%nbar
          WRITE (mp,'(a10)') 'rowperm(1:min(5,m)) = '
          WRITE (mp,'(5i15)') (rowperm(i),i=1,min(5,m))
          WRITE (mp,'(a10)') 'colperm(1:min(5,n)) = '
          WRITE (mp,'(5i15)') (colperm(i),i=1,min(5,n))
        END IF

        RETURN

100     info%flag = mc79_err_memory_alloc
        IF (printe) CALL mc79_print_message(info%flag,lp, &
          context='mc79_coarse')

200     info%flag = mc79_err_memory_dealloc
        IF (printe) CALL mc79_print_message(info%flag,lp, &
          context='mc79_coarse')
        RETURN

300     IF (printe) CALL mc79_print_message(info%flag,lp, &
          context='mc79_coarse')

      END SUBROUTINE mc79_coarse


! ---------------------------------------------------------------

      SUBROUTINE mc79_fine(m,n,ptr,row,rowperm,colperm,rowcomp,colcomp, &
          control,info)
! subroutine mc79_matching(m,n,ptr,row,rowperm,colperm,control,info)
! constructs a matching for a matrix A. A is stored in CSC format.

! m: is an integer with intent(in). It holds the number of rows in A
        INTEGER (myint_mc79), INTENT (IN) :: m

! n: is an integer with intent(in). It holds the number of columns in B
        INTEGER (myint_mc79), INTENT (IN) :: n

! ptr: holds row pointer for A
        INTEGER (myint_mc79), INTENT (IN) :: ptr(n+1)

! row: holds indices of A
        INTEGER (myint_mc79), INTENT (IN) :: row(:)

! rowperm: is an integer array with intent(out) of size m. It holds
! the row permutation. If rowperm(i) = j, the row j in the original matrix 
! becomes row i in the perumted matrix 
        INTEGER (myint_mc79), INTENT (OUT) :: rowperm(m)


! colperm: is an integer array with intent(out) of size n. It holds
! the column permutation. If colperm(i) = j, the column j in the original matrix
! becomes column i in the permuted matrix 
        INTEGER (myint_mc79), INTENT (OUT) :: colperm(n)


        INTEGER (myint_mc79), INTENT (OUT) :: rowcomp(m+2)

        INTEGER (myint_mc79), INTENT (OUT) :: colcomp(n+2)

! control: is of derived type mc79_control with intent(in). Controls
! action
        TYPE (mc79_control), INTENT (IN) :: control

! info: is of derived type mc79_info with intent(out).
! info%flag
! = 0 if successful
! = mc79_ERR_MEMORY_ALLOC if memory allocation failed
! = mc79_ERR_MEMORY_DEALLOC if memory deallocation failed
! = mc79_ERR_M if m<1
! = mc79_ERR_N if n<1
        TYPE (mc79_info), INTENT (OUT) :: info

! ---------------------------------------------
! Local allocatables
! ---------------------------------------------
        INTEGER, ALLOCATABLE, DIMENSION (:) :: ptr_temp, row_temp, rowmatch, &
          colmatch, col_markers, row_markers

! ---------------------------------------------
! Local variables
! ---------------------------------------------

        LOGICAL :: printe ! errors to be printed?
        LOGICAL :: printi ! info to be printed?
        LOGICAL :: printd ! extended info to be printed?
        INTEGER :: lp ! stream for printing errors
        INTEGER :: mp ! stream for printing info
        INTEGER :: i, l, ne, nhrows, nhcols, nvrows, nvcols, inbar, cmbase, &
          rnbase, cnbase, mn, vindex, sqindex, hindex, j

        info%flag = 0

! ---------------------------------------------
! Set stream numbers
! ---------------------------------------------
        lp = control%lp
        mp = control%mp

! ---------------------------------------------
! Printing levels
! ---------------------------------------------
        printe = (control%print_level>=0 .AND. lp>=0)
        printi = (control%print_level>=1 .AND. mp>=0)
        printd = (control%print_level>=2 .AND. mp>=0)

        info%flag = 0

        IF (printi) THEN
          WRITE (mp,'(a)') ' '
          WRITE (mp,'(a15)') 'MC79_fine:'
        END IF

! Check that m, n and ne have valid entries
        IF (m<1) THEN
          info%flag = mc79_err_m
          GO TO 300
        END IF

        IF (n<1) THEN
          info%flag = mc79_err_n
          GO TO 300
        END IF

        IF (printi) THEN
          WRITE (mp,'(a4,i15)') 'm = ', m
          WRITE (mp,'(a4,i15)') 'n = ', n
        END IF

        rowcomp(:) = 0
        colcomp(:) = 0
! Check whether there are any nonzero entries
        IF (ptr(n+1)==1) THEN
          info%m1 = 0
          info%n1 = n
          info%m2 = 0
          info%n2 = 0
          info%m3 = m
          info%n3 = 0
          info%mbar = info%m3 - info%n3
          info%nbar = info%n1 - info%m1
          info%hz_comps = 1
          rowcomp(1) = 1
          rowcomp(2) = 1
          colcomp(1) = 1
          colcomp(2) = 1 + info%n1
          info%sq_comps = 0
          info%vt_comps = 1
          rowcomp(info%hz_comps+info%sq_comps+info%vt_comps+1) &
            = rowcomp(info%hz_comps+info%sq_comps+1) + info%m3
          colcomp(info%hz_comps+info%sq_comps+info%vt_comps+1) &
            = colcomp(info%hz_comps+info%sq_comps+1) + info%n3
          rowperm = (/ (i,i=1,m) /)
          colperm = (/ (i,i=1,n) /)
          RETURN
        END IF

        vindex = -1
        sqindex = -2
        hindex = -3

! Build temp = A^T 
        ne = ptr(n+1) - 1

        ALLOCATE (row_temp(ne),ptr_temp(m+1),STAT=info%stat)
        IF (info%stat>0) GO TO 100
        ptr_temp(1:m) = 0
        ptr_temp(m+1) = ne + 1

        DO i = 1, ne
          j = row(i)
          ptr_temp(j) = ptr_temp(j) + 1
        END DO

        ptr_temp(1) = ptr_temp(1) + 1
        DO i = 2, m
          ptr_temp(i) = ptr_temp(i) + ptr_temp(i-1)
        END DO

        DO i = 1, n
          DO j = ptr(i), ptr(i+1) - 1
            l = row(j)
            ptr_temp(l) = ptr_temp(l) - 1
            row_temp(ptr_temp(l)) = i
          END DO
        END DO


        ALLOCATE (rowmatch(m),colmatch(n),STAT=info%stat)
        IF (info%stat>0) GO TO 100

! Find the maximum matching

        IF (m>=n) THEN
          CALL find_max_match(m,n,ne,ptr,row,rowmatch,colmatch,info%flag, &
            info%stat)
        ELSE
          CALL find_max_match(n,m,ne,ptr_temp,row_temp,colmatch,rowmatch, &
            info%flag,info%stat)
        END IF

        IF (info%flag<0) GO TO 300

        ALLOCATE (row_markers(m),col_markers(n),STAT=info%stat)
        IF (info%stat>0) GO TO 100

        row_markers(:) = sqindex
        col_markers(:) = sqindex

        CALL find_rect(m,n,ne,ptr,row,hindex,sqindex,colmatch,rowmatch, &
          col_markers,row_markers,nhrows,nhcols,info%flag,info%stat)

        IF (info%flag<0) GO TO 300

        CALL find_rect(n,m,ne,ptr_temp,row_temp,vindex,sqindex,rowmatch, &
          colmatch,row_markers,col_markers,nvcols,nvrows,info%flag,info%stat)

        IF (info%flag<0) GO TO 300

        info%m1 = nhrows
        info%n1 = nhcols
        info%m2 = m - nhrows - nvrows
        info%n2 = info%m2
        info%m3 = nvrows
        info%n3 = nvcols
        info%mbar = info%m3 - info%n3
        info%nbar = info%n1 - info%m1



! Find connected components in the horizontal submatrix
        mn = max(m,n)

        colperm(:) = 0
        rowperm(:) = 0
        IF (nhcols>0) THEN
          cmbase = 0
          rnbase = 0
          cnbase = 0

          IF (nhrows>0) THEN
            CALL find_conn_comps(n,m,ne,ptr,row,ptr_temp,row_temp,cmbase, &
              cnbase,rnbase,hindex,nhcols,nhrows,row_markers,col_markers, &
              rowcomp,colcomp,rowperm,colperm,info%hz_comps,info%flag, &
              info%stat)
            IF (info%flag<0) GO TO 300


          ELSE
            inbar = cnbase + 1

            info%nbar = 0
            DO i = 1, n
              IF (colmatch(i)==0) THEN
                info%nbar = info%nbar + 1
                colperm(inbar) = i
                inbar = inbar + 1
              END IF
            END DO

            info%hz_comps = 1
            rowcomp(1) = 1
            rowcomp(2) = 1
            colcomp(1) = 1
            colcomp(2) = 1 + nhcols

          END IF
! IF (nhrows .eq. 0) info%hz_comps = 0

        ELSE
          info%hz_comps = 0

        END IF

        IF (info%m2>0) THEN

          CALL find_conn_comp_sq(n,m,ne,ptr,row,nhrows,nhcols,info%m2,sqindex, &
            info%hz_comps,rowmatch,row_markers,rowcomp,colcomp,rowperm, &
            colperm,info%sq_comps,info%flag,info%stat)
          IF (info%flag<0) GO TO 300

        ELSE
          info%sq_comps = 0

        END IF
        IF (nvrows>0) THEN
          cmbase = info%hz_comps + info%sq_comps
          rnbase = nhrows + info%m2
          cnbase = nhcols + info%m2
          IF (nvcols>0) THEN

            CALL find_conn_comps(m,n,ne,ptr_temp,row_temp,ptr,row,cmbase, &
              rnbase,cnbase,vindex,nvrows,nvcols,col_markers,row_markers, &
              colcomp,rowcomp,colperm,rowperm,info%vt_comps,info%flag, &
              info%stat)

            IF (info%flag<0) GO TO 300

          ELSE
            inbar = rnbase + 1

            info%mbar = 0
            DO i = 1, m
              IF (rowmatch(i)==0) THEN
                info%mbar = info%mbar + 1
                rowperm(inbar) = i
                inbar = inbar + 1
              END IF
            END DO

            info%vt_comps = 1
            rowcomp(info%hz_comps+info%sq_comps+info%vt_comps+1) &
              = rowcomp(info%hz_comps+info%sq_comps+1) + nvrows
            colcomp(info%hz_comps+info%sq_comps+info%vt_comps+1) &
              = colcomp(info%hz_comps+info%sq_comps+1) + nvcols

          END IF

        ELSE
          info%vt_comps = 0

        END IF


        DEALLOCATE (rowmatch,colmatch,row_markers,col_markers,ptr_temp, &
          row_temp,STAT=info%stat)
        IF (info%stat>0) GO TO 200
        IF (printd) THEN
          CALL mc79_print_message(info%flag,lp,context='mc79_coarse')
! ---------------------------------------------
! Print out perm
! ---------------------------------------------
          WRITE (mp,'(a10)') 'info%m1 = '
          WRITE (mp,'(i15)') info%m1
          WRITE (mp,'(a10)') 'info%m2 = '
          WRITE (mp,'(i15)') info%m2
          WRITE (mp,'(a10)') 'info%m3 = '
          WRITE (mp,'(i15)') info%m3
          WRITE (mp,'(a10)') 'info%mbar = '
          WRITE (mp,'(i15)') info%mbar
          WRITE (mp,'(a10)') 'info%n1 = '
          WRITE (mp,'(i15)') info%n1
          WRITE (mp,'(a10)') 'info%n2 = '
          WRITE (mp,'(i15)') info%n2
          WRITE (mp,'(a10)') 'info%n3 = '
          WRITE (mp,'(i15)') info%n3
          WRITE (mp,'(a10)') 'info%nbar = '
          WRITE (mp,'(i15)') info%nbar
          l = info%hz_comps + info%sq_comps + info%vt_comps + 1
          WRITE (mp,'(a10)') 'rowcomp = '
          WRITE (mp,'(5i15)') (rowcomp(i),i=1,l)
          WRITE (mp,'(a10)') 'colcomp = '
          WRITE (mp,'(5i15)') (colcomp(i),i=1,l)
          WRITE (mp,'(a10)') 'rowperm = '
          WRITE (mp,'(5i15)') (rowperm(i),i=1,m)
          WRITE (mp,'(a10)') 'colperm = '
          WRITE (mp,'(5i15)') (colperm(i),i=1,n)
        ELSE IF (printi) THEN
! ---------------------------------------------
! Print out first few entries of perm
! ---------------------------------------------
          WRITE (mp,'(a10)') 'info%m1 = '
          WRITE (mp,'(i15)') info%m1
          WRITE (mp,'(a10)') 'info%m2 = '
          WRITE (mp,'(i15)') info%m2
          WRITE (mp,'(a10)') 'info%m3 = '
          WRITE (mp,'(i15)') info%m3
          WRITE (mp,'(a10)') 'info%mbar = '
          WRITE (mp,'(i15)') info%mbar
          WRITE (mp,'(a10)') 'info%n1 = '
          WRITE (mp,'(i15)') info%n1
          WRITE (mp,'(a10)') 'info%n2 = '
          WRITE (mp,'(i15)') info%n2
          WRITE (mp,'(a10)') 'info%n3 = '
          WRITE (mp,'(i15)') info%n3
          WRITE (mp,'(a10)') 'info%nbar = '
          WRITE (mp,'(i15)') info%nbar
          l = info%hz_comps + info%sq_comps + info%vt_comps + 1
          WRITE (mp,'(a10)') 'rowcomp = '
          WRITE (mp,'(5i15)') (rowcomp(i),i=1,min(5,l))
          WRITE (mp,'(a10)') 'colcomp = '
          WRITE (mp,'(5i15)') (colcomp(i),i=1,min(5,l))
          WRITE (mp,'(a10)') 'rowperm(1:min(5,m)) = '
          WRITE (mp,'(5i15)') (rowperm(i),i=1,min(5,m))
          WRITE (mp,'(a10)') 'colperm(1:min(5,n)) = '
          WRITE (mp,'(5i15)') (colperm(i),i=1,min(5,n))
        END IF

        RETURN

100     info%flag = mc79_err_memory_alloc
        IF (printe) CALL mc79_print_message(info%flag,lp,context='mc79_fine')

200     info%flag = mc79_err_memory_dealloc
        IF (printe) CALL mc79_print_message(info%flag,lp,context='mc79_fine')
        RETURN

300     IF (printe) CALL mc79_print_message(info%flag,lp,context='mc79_fine')


      END SUBROUTINE mc79_fine



! ---------------------------------------------------------------


      SUBROUTINE mc79_print_message(flag,unit,context)
! Prints out errors and warnings according to value of flag

! flag: is an integer scaler of intent(in). It is the information flag
! whose corresponding error message is printed
        INTEGER (myint_mc79), INTENT (IN) :: flag

! unit: is an integer scaler of intent(in). It is the unit number the
! error message should be printed on
        INTEGER (myint_mc79), INTENT (IN) :: unit

! context: is an optional assumed size character array of intent(in).
! It describes the context under which the error occured
        CHARACTER (len=*), OPTIONAL, INTENT (IN) :: context

        INTEGER (myint_mc79) :: length

        IF (unit<=0) RETURN

        IF (flag<0) THEN
          WRITE (unit,advance='yes',fmt='('' ERROR: '')')
        END IF

        IF (present(context)) THEN
          length = len_trim(context)
          WRITE (unit,advance='no',fmt='('' '', a,'': '')') context(1:length)
        END IF

        SELECT CASE (flag)
        CASE (0)
          WRITE (unit,'(A)') 'successful completion'

        CASE (mc79_err_memory_alloc)
          WRITE (unit,'(A)') 'memory allocation failure'

        CASE (mc79_err_memory_dealloc)
          WRITE (unit,'(A)') 'memory deallocation failure'

        CASE (mc79_err_m)
          WRITE (unit,'(A)') 'restriction m>=1 violated'

        CASE (mc79_err_n)
          WRITE (unit,'(A)') 'restriction n>=1 violated'
        END SELECT

      END SUBROUTINE mc79_print_message



!------------------------------------------------------------------------------
      SUBROUTINE find_max_match(m,n,ne,ptr,row,rowset,colset,flag,stat)

! find_max_match uses depth-first search to find an augmenting path from
!  each column node to get the maximum matching.

! m is a scalar INTEGER with INTENT(IN). It holds the number of rows in A

! n is a scalar INTEGER with INTENT(IN). It holds the number of columns in A

! ptr points to the start of the columns in A

! row contains the row indices

! D is a scalar of type zd11_type and INTENT(IN). It holds the transpose of 
!  the graph A for which we want to find a matching.

! rowset is an INTEGER array of length m and INTENT(OUT). It describes the 
!  matching. If rowset(i) = j > 0, then column j in A is matched to row i in 
!  A. If rowset(i) = 0, then row i is unmatched.

! colset is an INTEGER array of length n and INTENT(OUT). It describes the 
!  matching. If colset(j) = i > 0, then row i in A is matched to column j in 
!  A. If colset(j) = 0, then column j is unmatched.

! flag is a scalar INTEGER with INTENT(OUT). It is equal to 0 if no errors or 
!   warnings were encountered.

! stat is a scalar INTEGER with INTENT(OUT). It holds the stat argument from 
!   allocations and deallocations.

! -------------------------------------
!     Scalar Arguments
! -------------------------------------
        INTEGER, INTENT (IN) :: m, n, ne, ptr(n+1), row(ne)
        INTEGER, INTENT (OUT) :: flag, stat

! -------------------------------------
!     Array Arguments
! -------------------------------------
        INTEGER, INTENT (OUT) :: colset(n), rowset(m)

! -------------------------------------
!     Local allocatables
! -------------------------------------
        INTEGER, ALLOCATABLE, DIMENSION (:) :: mask
!      mask(i) contains the index of the root of the current depth-first 
!      search. Row i has been visited in current pass when equality holds.
        INTEGER, ALLOCATABLE, DIMENSION (:) :: next_chp
!      next_chp(j) is a pointer into row to the next row in A to be explored 
!      from column j in A for the cheap assignment. It is set to -1 when all 
!      rows have been considered for cheap assignment
        INTEGER, ALLOCATABLE, DIMENSION (:) :: prev_col
!      prev_col contains pointers toward the root of the depth-first search from
!      a column of A to a column of A. The pair (prevrw,prev_col) represent a
!      matched pair.
        INTEGER, ALLOCATABLE, DIMENSION (:) :: prevrw
!      prevrw contains pointers toward the root of the depth-first search from 
!      a column of A to a row of A. The pair (prevrw,prev_col) represent a
!      matched pair.
        INTEGER, ALLOCATABLE, DIMENSION (:) :: try_row
!      try_row (j) is a pointer into row to the next row to be explored from 
!      column j in the depth-first search.

! -------------------------------------
!     Local Scalars 
! -------------------------------------
        INTEGER j, last_row, next_row, nodec, next_col, pcol, prow, i, xrow

! Allocate all of the arrays
        flag = 0
        stat = 0

        ALLOCATE (mask(m),next_chp(n),prev_col(n),prevrw(n),try_row(n), &
          STAT=stat)
        IF (stat>0) THEN
          flag = mc79_err_memory_alloc
          RETURN
        END IF

        rowset(:) = 0
        colset(:) = 0
        mask(:) = 0
        next_chp(:) = 0
        prev_col(:) = 0
        prevrw(:) = 0
        try_row(:) = 0


        DO nodec = 1, n

!       Initialize node 'j' as the root of the path.
          j = nodec
          prevrw(j) = 0
          prev_col(j) = 0
          next_chp(j) = ptr(j)
10        CONTINUE

!       main loop begins here. Each time through, try to find a
!            cheap assignment from node j.
          next_row = next_chp(j)
          last_row = ptr(j+1) - 1
          IF (next_row>0) THEN
            DO xrow = next_row, last_row
              i = row(xrow)
              IF (rowset(i)==0) GO TO 50
            END DO
!         mark column when all adjacent rows have been considered for cheap
!           assignment.
            next_chp(j) = -1
          END IF

!       Each time through, take a step forward if possible, or backtrack. Quit 
!         when backtracking takes us back to the beginning of the search.
          try_row(j) = ptr(j)
          next_row = try_row(j)
          IF (last_row>=next_row) THEN
            DO xrow = next_row, last_row
              i = row(xrow)
              IF (mask(i)<nodec) GO TO 20
            END DO
            GO TO 30

!         row i is currently not visited, so move forward
20          try_row(j) = xrow + 1
            mask(i) = nodec
            next_col = rowset(i)
            IF (next_col<0 .OR. next_col==j) THEN
              GO TO 70
            ELSE IF (next_col>0) THEN
!           The forward step led to a matched row. Try to extend augmenting 
!             path from the column matched by this row.
              prev_col(next_col) = j
              prevrw(next_col) = i
              try_row(next_col) = ptr(next_col)
              j = next_col
              GO TO 10
            ELSE
              GO TO 40 ! Should never be reached. Added as safe guard
            END IF
          END IF
!       No forward step so backtrack. If we backtrack all the way, then the 
!        search is done

30        j = prev_col(j)
          IF (j>0) THEN
            GO TO 10
          ELSE
            CYCLE
          END IF
!       We have an unmatched row

40        CONTINUE

!       Update the matching by alternating the matching edge backward toward 
!         the root
50        rowset(i) = j
          prow = prevrw(j)
          pcol = prev_col(j)
          DO WHILE (pcol>0)
            rowset(prow) = pcol
            j = pcol
            prow = prevrw(pcol)
            pcol = prev_col(pcol)
          END DO
60        CONTINUE
        END DO

!     Compute the matching from the columns of A
        DO i = 1, m
          j = rowset(i)
          IF (j>0) colset(j) = i
        END DO


70      CONTINUE
        DEALLOCATE (mask,next_chp,prev_col,prevrw,try_row,STAT=stat)
        IF (stat>0) THEN
          flag = mc79_err_memory_dealloc
          RETURN
        END IF

        RETURN


      END SUBROUTINE find_max_match

! ---------------------------------------------

      SUBROUTINE find_rect(m,n,ne,ptr,col,marked,unmarked,colset,rowset, &
          col_markers,row_markers,nhrows,nhcols,flag,stat)

! Finds a rectangular portion of matrix A by depth-first search. Uses a 
!   depth-first serch to find all the rows and columns, which can be reached 
!   via alternating paths beginning from all the unmatched columns.

! m holds the number of rows in A

! n holds the number of columns in A

! ne holds the number of nonzero entries in A

! ptr holds points to the start of rows in A^T

! col holds the column indices of entries in A^T

! D is a scalar of type zd11_type and INTENT(IN). It holds the transpose of 
!  the graph A for which we want to find a matching.

! marked is an INTEGER scalar with INTENT(IN). It holds the value to store in 
!  marker vectors to indicate that row/column has been reached and is
!  therefore in the horizontal block

! unmarked is an INTEGER scalar with INTENT(IN). It holds the initial value in 
!  marker vectors and indicates that the row/column is free to be chosen

! colset is an INTEGER array of length n and INTENT(IN). It describes the 
!  matching. If colset(j) = i > 0, then row i in A is matched to column j in 
!  A. If colset(j) = 0, then column j is unmatched.

! rowset is an INTEGER array of length m and INTENT(IN). It describes the 
!  matching. If rowset(i) = j > 0, then column j in A is matched to row i in 
!  A. If rowset(i) = 0, then row i is unmatched.

! col_markers is an INTEGER array of length n and INTENT(OUT). It is used to 
!  mark whether a column is free for use

! row_markers is an INTEGER array of length m and INTENT(OUT). It is used to 
!  mark whether a row is free for use

! nhrows is an INTEGER scalar with INTENT(OUT). It holds the number of rows in
!  the horizontal block.

! nhcols is an INTEGER scalar with INTENT(OUT). It holds the number of columns 
!  in the horizontal block.

! flag is a scalar INTEGER with INTENT(OUT). It is equal to 0 if no errors or 
!   warnings were encountered.

! stat is a scalar INTEGER with INTENT(OUT). It holds the stat argument from 
!   allocations and deallocations.


! -------------------------------------
!     Scalar Arguments
! -------------------------------------

        INTEGER, INTENT (IN) :: marked, unmarked, m, n, ne
        INTEGER, INTENT (OUT) :: nhcols, nhrows, flag, stat

! -------------------------------------
!     Array Arguments
! -------------------------------------
        INTEGER, INTENT (IN) :: colset(n), rowset(m), col(ne), ptr(n+1)
        INTEGER, INTENT (INOUT) :: col_markers(n), row_markers(m)

! -------------------------------------
!     Local Scalars
! -------------------------------------
        INTEGER :: j, fromc, next_col, next_row, p, i, xrow

! -------------------------------------
!     Local Allocatables
! -------------------------------------
        INTEGER, ALLOCATABLE, DIMENSION (:) :: prev_col, try_row
!     try_row (j) is a pointer into  col to the next row to be explored from 
!      col 'j' in the search.
!     prev_col contains pointers toward the root of the search from
!                   column to column.

        flag = 0
        flag = 0
        stat = 0

        ALLOCATE (prev_col(n),try_row(n),STAT=stat)
        IF (stat>0) THEN
          flag = mc79_err_memory_alloc
          RETURN
        END IF
        prev_col(:) = 0
        try_row(:) = 0

        nhcols = 0
        nhrows = 0
        DO p = 1, n
!        Search for an unmatched column to start the alternating path.
          IF (colset(p)==0) THEN
            fromc = p
!         The path starts from unmatched column "fromc". Put fromc into 
!         horizontal set "hc" and indicate fromc is the root of the path.
            nhcols = nhcols + 1
            col_markers(fromc) = marked
            try_row(fromc) = ptr(fromc)
            prev_col(fromc) = 0
            j = fromc
10          CONTINUE
!          The main depth-first search loop begins here.
!          If possible, we take a step forward; otherwise, we backtrack.
!          We quit loop when we backtrack to the beginning of the search.

!          Find a forward step from column 'j' to an unmarked row.
            next_row = try_row(j)
            DO xrow = next_row, ptr(j+1) - 1
              IF (row_markers(col(xrow))==unmarked) GO TO 20
            END DO

!         No forward step available so backtrack.  If we backtrack all the way,
!         we have completed all searchs beginning at column 'p'.
            j = prev_col(j)
            IF (j/=0) THEN
              GO TO 10
            ELSE
              CYCLE
            END IF

!         Take a double forward step from 'j' to 'i' and THEN via matching 
!         edge from 'i' to column 'next_col'.
20          try_row(j) = xrow + 1
            i = col(xrow)
            row_markers(i) = marked
            nhrows = nhrows + 1
            next_col = rowset(i)
            nhcols = nhcols + 1
            col_markers(next_col) = marked
            prev_col(next_col) = j
            try_row(next_col) = ptr(next_col)
            j = next_col
            GO TO 10
          END IF
30        CONTINUE
        END DO

        DEALLOCATE (prev_col,try_row,STAT=stat)
        IF (stat>0) THEN
          flag = mc79_err_memory_dealloc
          RETURN
        END IF

      END SUBROUTINE find_rect

! --------------------------------------

      SUBROUTINE find_conn_comps(m,n,ne,a_ptr,a_col,d_ptr,d_col,cmbase,rnbase, &
          cnbase,vindex,nvrows,nvcols,col_markers,row_markers,colcomp,rowcomp, &
          col_perm,row_perm,num_comp,flag,stat)

!   Finds the connected components in the subgraph spanned by the rows and 
!    columns in the vertical block. 

! m holds the number of rows in A

! n holds the number of columns in A

! ne holds the number of nonzero entries in A

! a_ptr holds points to the start of rows in A

! a_col holds the column indices of entries in A

! d_ptr holds points to the start of rows in A^T

! d_col holds the column indices of entries in A^T

! cmbase is an INTEGER scalar of INTENT(IN). It holds the number of components 
!  found in previous fine analysis of the coarse partition

! rnbase is an INTEGER scalar of INTENT(IN). It holds the number of rows in 
!  earlier numbered partitions (0 for the horizontal block, nhrows+nsrows for
!  the vertical partition)

! cnbase is an INTEGER scalar of INTENT(IN). It holds the number of columns in 
!  earlier numbered partitions

! vindex is an INTEGER scalar of INTENT(IN). It is used to check whether the 
!  nodes belong in the vertical block

! nvrows is an INTEGER scalar of INTENT(IN). It holds the number of rows in 
!  the vertical block

! nvcols is an INTEGER scalar of INTENT(IN). It holds the number of columns in 
!  the vertical block

! num_comp is an INTEGER scalar of INTENT(OUT). It holds the number of connected
!  components

! col_markers is an INTEGER array with length  n and INTENT(INOUT). Initially,
!        col_markers(i) = vindex if i belongs to the vertical component;
!                  < 0 otherwise.
!  During execution,
!        col_markers(i) = j, if i belongs to the jth component.
!  After execution, original values restored

! row_markers is an INTEGER array with length  m and INTENT(INOUT). Initially,
!        row_markers(i) = vindex if i belongs to vr;
!                  < 0  otherwise.
!  During execution,
!        row_markers(i) = j, if i belongs to the jth component;
!                  < 0 otherwise.
!  After execution, original values restored

!  colcomp is an rank-one INTEGER array with length  n and INTENT(INOUT). It 
!    holds the address (in the new ordering) of the first column in each 
!    component.

!  rowcomp is an rank-one INTEGER array with length  m and INTENT(INOUT). It 
!    holds the address (in the new ordering) of the first row in each component.

!  col_perm is an rank-one INTEGER array with length  n and INTENT(INOUT). 
!    If col_perm(j) = k, then the j-th column in the permuted matrix 
!    corresponds to the k-th column in A.

!  row_perm is an rank-one INTEGER array with length  n and INTENT(INOUT). If 
!    col_perm(j) = k, then the j-th row in the permuted matrix corresponds to 
!    the k-th row in A.

! flag is a scalar INTEGER with INTENT(OUT). It is equal to 0 if no errors or 
!   warnings were encountered.

! stat is a scalar INTEGER with INTENT(OUT). It holds the stat argument from 
!   allocations and deallocations.


! -------------------------------------
!     Scalar arguments
! -------------------------------------
        INTEGER, INTENT (IN) :: cmbase, cnbase, nvcols, nvrows, rnbase, &
          vindex, m, n, ne, a_ptr(m+1), a_col(ne), d_ptr(n+1), d_col(ne)
        INTEGER, INTENT (OUT) :: flag, stat, num_comp

! -------------------------------------
!     Array Arguments
! -------------------------------------
        INTEGER, INTENT (INOUT) :: colcomp(n+2), rowcomp(m+2), col_perm(n), &
          col_markers(n), row_perm(m), row_markers(m)

! -------------------------------------
!     Local Allocatables
! -------------------------------------
        INTEGER, ALLOCATABLE, DIMENSION (:) :: ctab, rtab, nextcl, nextrw, &
          predcl, predrw
!     ctab and rtab hold a temporary copy of the address (in the new ordering)
!       of the first column/row in each component
!     working variables:

!     predrw and predcl hold the stack of paths
!        predrw(i) = j means that we have in the path an edge leaving from row 
!             j to column node i.
!        predcl(i) = j means that we have in the path an edge leaving from 
!             column node j to row node i.
!     nextcl(i) is index of first unsearched edge leaving from column node i.
!     nextrw(i) is index of first unsearched edge leaving from row node i.

! -------------------------------------
!     Local Scalars
! -------------------------------------
        INTEGER cn, j, compn, p, rn, i, xcol, xrow, mn
!     cn holds the number of the scanned columns
!     rn holds the number of the scanned rows
        flag = 0
        mn = max(n,m)

        ALLOCATE (ctab(mn),rtab(mn),nextcl(n),nextrw(m),predcl(m),predrw(n), &
          STAT=stat)
        IF (stat>0) THEN
          flag = mc79_err_memory_alloc
          RETURN
        END IF
        ctab(:) = 0
        rtab(:) = 0
        cn = 0
        rn = 0
        num_comp = 0

!     Note that the number of vertical rows > number of vertical columns.
!         Start each search for a connected component with an unmarked row in 
!         the vertical block.
        j = 0
        DO p = 1, m
          IF (row_markers(p)==vindex) THEN
            i = p
            IF (a_ptr(i)==a_ptr(i+1)) THEN
              IF (num_comp>0) THEN
                row_markers(i) = num_comp
                rn = rn + 1
              ELSE
                num_comp = num_comp + 1
                ctab(num_comp) = cnbase + 1 + cn
                rtab(num_comp) = rnbase + 1 + cn
                row_markers(i) = num_comp
                colcomp(cmbase+num_comp) = ctab(num_comp)
                rowcomp(cmbase+num_comp) = rtab(num_comp)
                rn = rn + 1
              END IF

            ELSE
!         Update the value of the current working component. Put 'i' into the 
!         new component as the root of path
              IF (i==1) THEN
                num_comp = num_comp + 1
                ctab(num_comp) = cnbase + 1 + cn
                rtab(num_comp) = rnbase + 1 + rn
                colcomp(cmbase+num_comp) = ctab(num_comp)
                rowcomp(cmbase+num_comp) = rtab(num_comp)
              ELSE
                IF (num_comp>1 .OR. a_ptr(i)/=a_ptr(i-1)) THEN
                  num_comp = num_comp + 1
                  ctab(num_comp) = cnbase + 1 + cn
                  rtab(num_comp) = rnbase + 1 + rn
                  colcomp(cmbase+num_comp) = ctab(num_comp)
                  rowcomp(cmbase+num_comp) = rtab(num_comp)


                END IF

              END IF
              row_markers(i) = num_comp
              rn = rn + 1
              nextrw(i) = a_ptr(i)
              predcl(i) = 0
10            CONTINUE
!         From row node to col node - try to find a forward step if possible
!               else backtrack
              DO xcol = nextrw(i), a_ptr(i+1) - 1
                j = a_col(xcol)
                IF (col_markers(j)==vindex) GO TO 20
              END DO
!         Go backward for one step  (back to col node)
              nextrw(i) = a_ptr(i+1)
              j = predcl(i)
              IF (j==0) THEN
                CYCLE
              ELSE
                GO TO 30
              END IF
!         Go forward one step : find a forward step from row 'i' to column
!           'j'.  put 'j' into the current component
20            nextrw(i) = xcol + 1
              col_markers(j) = num_comp
              cn = cn + 1
              nextcl(j) = d_ptr(j)
              predrw(j) = i

30            DO xrow = nextcl(j), d_ptr(j+1) - 1
                i = d_col(xrow)
                IF (row_markers(i)==vindex) GO TO 40
              END DO
!         Go back one step  (back to row node)
              nextcl(j) = d_ptr(j+1)
              i = predrw(j)
              GO TO 10
!         Go forward one step : find a forward step from column 'j' to
!            row 'i'.  put row into the current component
40            nextcl(j) = xrow + 1
              row_markers(i) = num_comp
              rn = rn + 1
              nextrw(i) = a_ptr(i)
              predcl(i) = j
              GO TO 10
            END IF
          END IF
50        CONTINUE
        END DO
!     Generate the column and row permutations (col_perm and row_perm) so that
!       each component is numbered consecutively
        colcomp(cmbase+1+num_comp) = cnbase + 1 + nvcols
        rowcomp(cmbase+1+num_comp) = rnbase + 1 + nvrows
        DO j = 1, n
          compn = col_markers(j)
          IF (compn>0) THEN
            col_perm(ctab(compn)) = j
            ctab(compn) = ctab(compn) + 1
            col_markers(j) = vindex
          END IF
        END DO
        DO i = 1, m
          compn = row_markers(i)
          IF (compn>0) THEN
            row_perm(rtab(compn)) = i
            rtab(compn) = rtab(compn) + 1
            row_markers(i) = vindex
          END IF
        END DO
        DEALLOCATE (ctab,rtab,nextcl,nextrw,predcl,predrw,STAT=stat)
        IF (stat>0) THEN
          flag = mc79_err_memory_dealloc
          RETURN
        END IF
      END SUBROUTINE find_conn_comps

      SUBROUTINE find_conn_comp_sq(m,n,ne,ptr,idx,nhcols,nhrows,nscols,sqindx, &
          hrzcmp,colset,col_markers,colcomp,rowcomp,colperm,rowperm,sqcmpn, &
          flag,stat)
!    find_conn_comp_sq finds the lower block triangular form of the square
!     submatrix in the general block triangular form. The square submatrix
!     consists entirely of matched rows and columns.  

!     If non-matching edges are directed from rows to columns, and matching
!     edges are shrunk into single vertices, the resulting directed graph has
!     strongly connected components.

!     find_conn_comp_sq uses Tarjan's algorithm to find the strongly connected
!     components by depth-first search. All the pairs have been visited
!     will be labeled in the order they are visited, and associated a
!     lowlink for each vertex, stored in the stack - lowlnk.

! m holds the number of rows in A

! n holds the number of columns in A

! ne holds the number of nonzero entries in A

! ptr holds points to the start of rows in A

! idx holds the column indices of entries in A

! nhrows is an INTEGER scalar of INTENT(IN). It holds the number of rows in 
!  the horizontal block

! nhcols is an INTEGER scalar of INTENT(IN). It holds the number of columns in 
!  the horizontal block

! nscols is an INTEGER scalar of INTENT(IN). It holds the number of columns and
!  rows in the square partition

! sqindx is an INTEGER scalar of INTENT(IN). It holds the value used to 
!  indentify rows and columns in the square partition

! hrzcmp is an INTEGER scalar of INTENT(IN). It holds the  number of 
!  components in vertical partition

! colset is an INTEGER array of length n and INTENT(IN). It describes the 
!  matching. If colset(j) = i > 0, then row i in A is matched to column j in 
!  A. If colset(j) = 0, then column j is unmatched.

! col_markers is an INTEGER array with length n and INTENT(INOUT). Initially,
!        col_markers(i) = vindex if i belongs to the vertical component;
!                  < 0 otherwise.
!  During execution,
!        col_markers(i) = j, if i belongs to the jth component.
!  After execution, original values restored

! sqcmpn is an INTEGER scalar of INTENT(OUT). It holds the number of 
!    components in the square partition

! colcomp is an INTEGER array with length n and INTENT(INOUT). It 
!    holds the address (in the new ordering) of the first column in each 
!    component.

! rowcomp is an INTEGER array with length m and INTENT(INOUT). It 
!    holds the address (in the new ordering) of the first row in each 
!    component.

!  col_perm is an rank-one INTEGER array with length n and INTENT(INOUT). 
!    If col_perm(j) = k, then the j-th column in the permuted matrix 
!    corresponds to the k-th column in A.

!  row_perm is an rank-one INTEGER array with length n and INTENT(INOUT). If 
!    col_perm(j) = k, then the j-th row in the permuted matrix corresponds to 
!    the k-th row in A.

! flag is a scalar INTEGER with INTENT(OUT). It is equal to 0 if no errors or 
!   warnings were encountered.

! stat is a scalar INTEGER with INTENT(OUT). It holds the stat argument from 
!   allocations and deallocations.


! -------------------------------------
!     Scalar arguments
! -------------------------------------
        INTEGER, INTENT (IN) :: hrzcmp, nhcols, nhrows, nscols, sqindx, m, n, &
          ne, ptr(m+1), idx(ne)
        INTEGER, INTENT (OUT) :: sqcmpn, stat, flag

! -------------------------------------
!     Array arguments
! -------------------------------------
        INTEGER, INTENT (INOUT) :: colcomp(n+2), colperm(n), col_markers(n), &
          colset(n), rowcomp(m+2), rowperm(m)


! -------------------------------------
!     Local scalars
! -------------------------------------
        INTEGER cmpbeg, col, compnt, counter, fcol, fnlpos, frow, j, pair, &
          passes, rootcl, scol, stackp, xcol

!     fnlpos holds the number of pairs whose positions in final ordering
!             have been found.
!     counter holds the number of pairs on the stack  (stack pointer)

! -------------------------------------
!     Local Allocatables
! -------------------------------------
        INTEGER, ALLOCATABLE, DIMENSION (:) :: cbegin, lowlnk, prev, trycol
!        trycol(i) contains a pointer to next unsearched column for row i
!        cbegin(i) contains the beginning of component i
!        lowlnk stores the lowlink for each pair. This is the smallest stack 
!          position of any pair to which a path from pair i has been found. It 
!          is set to n+1 when pair i is removed from the stack.
!        prev(i) holds the pair at the end of the path when pair i was
!                  placed on the stack

        flag = 0

        ALLOCATE (cbegin(nscols),lowlnk(n),prev(n),trycol(m),STAT=stat)
        IF (stat>0) THEN
          flag = mc79_err_memory_alloc
          RETURN
        END IF

        cbegin(:) = 0
        lowlnk(:) = 0
        prev(:) = 0
        trycol(:) = 0

        fnlpos = 0
        sqcmpn = 0
        DO col = 1, n
          IF (col_markers(col)==sqindx) col_markers(col) = 0
        END DO
        DO j = 1, m
          trycol(j) = ptr(j)
        END DO

!     Find a starting pair
        DO rootcl = 1, n
          IF (col_markers(rootcl)==0) THEN
!       Place the pair (rootcl, colset(rootcl)) at beginning of stack
            fcol = rootcl
            counter = 1
            lowlnk(fcol) = counter
            col_markers(fcol) = counter
            cbegin(nscols) = fcol

!         Either put a new pair on the stack or backtrack
            DO passes = 1, 2*nscols - 1
              frow = colset(fcol)

!           Check whether all edges leaving pair (frow,fcol)  been searched?
              IF (trycol(frow)>0) THEN
!             Check all edges leaving from row "frow" until we find a new 
!             column "scol" that has not yet been encountered or until all 
!             edges are exhausted.
                DO xcol = trycol(frow), ptr(frow+1) - 1
                  scol = idx(xcol)
                  IF (col_markers(scol)==0) THEN
                    GO TO 10
                  ELSE IF (col_markers(scol)>0) THEN
!                  If scol has been on stack already, update value of low(fcol)
!                  if necessary 
                    IF (lowlnk(scol)<lowlnk(fcol)) lowlnk(fcol) = lowlnk(scol)
                  END IF
                END DO
!             There are no more edges leaving frow
                trycol(frow) = -1
                GO TO 20
!             Add the new pair  (scol, colset(scol)) to the stack
10              trycol(frow) = xcol + 1
                prev(scol) = fcol
                fcol = scol
                counter = counter + 1
                lowlnk(fcol) = counter
                col_markers(fcol) = counter
                cbegin(nscols+1-counter) = fcol
                CYCLE
              END IF

20            IF (lowlnk(fcol)>=col_markers(fcol)) THEN

!              If frow is the root of a block, then we have found a component. 
!              Order the nodes in this  block by starting at the top of the 
!              stack and working down to the root of the block
                sqcmpn = sqcmpn + 1
                cmpbeg = fnlpos + 1
                DO stackp = nscols + 1 - counter, nscols
                  pair = cbegin(stackp)
                  fnlpos = fnlpos + 1
                  col_markers(pair) = fnlpos
                  counter = counter - 1
                  lowlnk(pair) = nscols + 1
                  IF (pair==fcol) EXIT
                END DO

30              cbegin(sqcmpn) = cmpbeg
!             If there are any pairs left on the stack, then backtrack. 
!             Otherwise, check if all the pairs have been ordered
                IF (counter==0) GO TO 50
              END IF

!           Backtrack to previous pair on path
              scol = fcol
              fcol = prev(fcol)
              IF (lowlnk(scol)<lowlnk(fcol)) lowlnk(fcol) = lowlnk(scol)
40            CONTINUE
            END DO
            CYCLE ! Should never be reached. Added as safe guard
50          IF (fnlpos>=nscols) GO TO 70
          END IF
        END DO

70      DO compnt = 1, sqcmpn
          colcomp(compnt+hrzcmp) = (cbegin(compnt)+nhcols)
          rowcomp(compnt+hrzcmp) = (cbegin(compnt)+nhcols) - (nhcols-nhrows)
        END DO
        colcomp(hrzcmp+sqcmpn+1) = nhcols + nscols + 1
        rowcomp(hrzcmp+sqcmpn+1) = nhrows + nscols + 1

        DO col = 1, n
          j = col_markers(col)
          IF (j>0) THEN
            colperm(nhcols+j) = col
            rowperm(nhrows+j) = colset(col)
            col_markers(col) = sqindx
          END IF
        END DO
        DEALLOCATE (cbegin,lowlnk,prev,trycol,STAT=stat)
        IF (stat>0) THEN
          flag = mc79_err_memory_dealloc
          RETURN
        END IF
      END SUBROUTINE find_conn_comp_sq


    END MODULE hsl_mc79_integer
! (c) STFC 2010-2011
! Originating author: Jonathan Hogg
!
! Given a pivot order, this package performs common tasks
! required in the analyse phase of a symmetric sparse direct solver.
! Either the entire analyse may be performed or individual tasks.
! The matrix may be hled in assembled form or in elemental form.
!
! Version 1.2.0
! See ChangeLog for version history

! To convert to long:
! s/_integer/_long_integer
! Set pkg_type to long
module hsl_mc78_integer
   implicit none

   private
   public :: mc78_control
   public :: mc78_analyse, mc78_supervars, mc78_compress_by_svar, mc78_etree, &
      mc78_elt_equiv_etree, mc78_postorder, mc78_col_counts, mc78_supernodes, &
      mc78_stats, mc78_row_lists, mc78_optimize_locality

   integer, parameter :: dp = kind(0d0) ! not package type
   integer, parameter :: long = selected_int_kind(18)

   integer, parameter :: minsz_ms = 16 ! minimum size to use merge sort

   integer, parameter :: pkg_type = kind(0) ! package type - integer or long

   type mc78_control
      integer :: heuristic = 1 ! 1=ma77 2=cholmod
      integer :: nrelax(3) = (/ 4, 16, 48 /) ! CHOLMOD-like
      real(dp) :: zrelax(3) = (/ 0.8, 0.1, 0.05 /) ! CHOLMOD-like
      integer :: nemin = 16  ! Node amalgamation parameter

      integer :: unit_error = 6
      integer :: unit_warning = 6
      logical :: ssa_abort = .false. ! If .true., then return with an error if
         ! an assembled matrix is detected as symbolically singular (we do
         ! not garuntee to detect _all_ symbolically singular matrices).
         ! If .false., then a warning is raised instead.

      logical :: svar = .false. ! If .true. then supervariables are used in
         ! the assembled case, otherwise they are not. Supervaraibles are
         ! always used in the elemental case.
      logical :: sort = .false. ! If .true. then entries within each supernode's
         ! row lists are sorted. Otherwise they might not be.
      logical :: lopt = .false. ! If .true. then variable ordering is optimized
         ! for cache locality. Otherwise it is not.
   end type mc78_control

   integer, parameter :: MC78_ERROR_ALLOC = -1 ! allocate error
   integer, parameter :: MC78_ERROR_SSA   = -2 ! symbolically singular assembled
   integer, parameter :: MC78_ERROR_ROW_SMALL = -3 ! supplied row array to short
   integer, parameter :: MC78_ERROR_UNKNOWN = -99 ! internal/unknown error

   ! Warning flags are treated as bit masks, add together if multiple occour
   integer, parameter :: MC78_WARNING_SSA = 1 ! symbolically singular assembled
   integer, parameter :: MC78_WARNING_BLK_SVAR = 2 ! svar and blk pivs requested

   interface mc78_analyse
      module procedure mc78_analyse_assembled_integer
      module procedure mc78_analyse_elemental_integer
   end interface mc78_analyse

   interface mc78_supervars
      module procedure mc78_supervars_integer
   end interface mc78_supervars

   interface mc78_compress_by_svar
      module procedure mc78_compress_by_svar_integer
   end interface mc78_compress_by_svar

   interface mc78_etree
      module procedure mc78_etree_integer
   end interface mc78_etree

   interface mc78_elt_equiv_etree
      module procedure mc78_elt_equiv_etree_integer
   end interface

   interface mc78_postorder
      ! Note: cannot distinguish postorder_std between integer and long versions
      module procedure mc78_postorder_std
      module procedure mc78_postorder_detect
   end interface mc78_postorder

   interface mc78_col_counts
      module procedure mc78_col_counts_integer
   end interface mc78_col_counts

   ! Note: cannot distinguish mc78_supernodes between integer and long versions
   ! Note: cannot distinguish mc78_stats between integer and long versions

   interface mc78_row_lists
      module procedure mc78_row_lists_nosvar_integer
      module procedure mc78_row_lists_svar_integer
   end interface mc78_row_lists

   ! Note: cannot distinguish mc78_optimize_locality between integer and long
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            Main analysis routines   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! For assembled matrix input, this subroutine performs a full analysis.
! This is essentially a wrapper around the rest of the package.
!
! Performance might be improved by:
! * Improving the sort algorithm used in find_row_idx
!
subroutine mc78_analyse_assembled_integer(n, ptr, row, perm, nnodes, sptr, &
      sparent, rptr, rlist, control, info, stat, nfact, nflops, piv_size)
   integer, intent(in) :: n ! Dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, intent(out) :: nnodes ! number of supernodes found
   integer, dimension(:), allocatable, intent(out) :: sptr ! supernode pointers
   integer, dimension(:), allocatable, intent(out) :: sparent ! assembly tree
   integer(long), dimension(:), allocatable, intent(out) :: rptr
      ! pointers to rlist
   integer, dimension(:), allocatable, intent(out) :: rlist ! row lists
   ! For details of control, info : see derived type descriptions
   type(mc78_control), intent(in) :: control
   integer, intent(out) :: info
   integer, optional, intent(out) :: stat
   integer(long), optional, intent(out) :: nfact
   integer(long), optional, intent(out) :: nflops
   integer, dimension(n), optional, intent(inout) :: piv_size ! If
      ! present, then matches matrix order and specifies block pivots. 
      ! piv_size(i) is the number of entries pivots in the block pivot
      ! containing column i.

   integer(pkg_type), dimension(:), allocatable :: bptr ! copy of matrix with
      ! added entries for block pivots - column pointers
   integer, dimension(:), allocatable :: brow ! copy of matrix with added
      ! entries for block pivots - row indices
   integer :: flag ! return status flag for call to compress_by_svar
   integer(pkg_type) :: sz
   integer, dimension(:), allocatable :: svara ! array for supervariables
   integer :: i
   integer, dimension(:), allocatable :: invp ! inverse permutation of perm
   integer, dimension(:), allocatable :: sinvp
   integer :: j
   integer :: k
   integer :: nsvar ! number of supervariables
   integer, dimension(:), allocatable :: sperm
   integer(pkg_type), dimension(:), allocatable :: ptr2
   integer :: realn ! number of variables with an actual entry present
   integer, dimension(:), allocatable :: row2
   integer :: st ! stat argument in allocate calls
   logical :: svar_r

   integer :: svar_type ! 0=none, 1=col, 2=compressed form

   integer, dimension(:), allocatable :: scc

   ! initialise
   info = 0

   svar_r = control%svar

   ! Ensure allocatable output arguments are deallocated
   deallocate(sptr, stat=st)
   deallocate(sparent, stat=st)
   deallocate(rptr, stat=st)
   deallocate(rlist, stat=st)

   ! Initialise inverse permutation and check for duplicates
   allocate(invp(n), stat=st)
   if(st.ne.0) goto 490
   do i = 1, n
      j = perm(i)
      invp(j) = i
   end do

   if(present(piv_size)) then
      call convert_to_blk_piv(n, invp, piv_size)
      allocate(bptr(n+1), brow(ptr(n+1)-1+2*n), stat=st)
      if(st.ne.0) goto 490
      call mc78_block_prep(n, ptr, row, bptr, brow, perm, invp, piv_size, &
         st)
      if(st.ne.0) goto 490
      if(svar_r) then
         ! Supervariables don't interact well with block pivots, so don't do it
         svar_r = .false.
         info = info + MC78_WARNING_BLK_SVAR
      endif
   endif

   ! Determine supervariables (if required)
   if(svar_r) then
      allocate(svara(n), stat=st)
      if(st.ne.0) goto 490
      realn = n
      call mc78_supervars(realn, ptr, row, perm, invp, nsvar, svara, st)
      if(st.ne.0) goto 490
      if(n.ne.realn) then
         if(control%ssa_abort) then
            if(control%unit_error.gt.0) write(control%unit_error, "(2a)") &
               "HSL_MC78: Error, matrix is symbolically singular and ", &
               "control%ssa_abort=.true.."
            info = MC78_ERROR_SSA
            return
         else
            if(control%unit_warning.gt.0) write(control%unit_warning, "(a)") &
               "HSL_MC78: Warning, matrix is symbolically singular."
            info = info + MC78_WARNING_SSA
         endif
      endif
      if(3*nsvar.lt.2*n) then
         svar_type = 2 ! Use compressed form
      else
         svar_type = 0 ! Do not use supervariables
      endif

      select case(svar_type)
      case(0) ! do not use supervariables
         ! release resources
         deallocate(svara, stat=st)
      case(2) ! Compressed form
         ! It is worth using the compressed form
         ! Determine upper bound on size of data for compressed array
         sz = 0
         k = 1
         do i = 1, nsvar
            j = invp(k)
            sz = sz + ptr(j+1) - ptr(j)
            k = k + svara(i)
         end do
         allocate(ptr2(nsvar+1), row2(sz), sperm(nsvar), sinvp(nsvar), stat=st)
         if(st.ne.0) goto 490
         call mc78_compress_by_svar(n, ptr, row, invp, nsvar, svara, &
            ptr2, sz, row2, flag, st)
         select case(flag)
         case(0) ! Everything OK
            ! Do nothing
         case(-1) ! Allocate failure
            goto 490
         case default ! Should never happen
            info = MC78_ERROR_UNKNOWN
            return
         end select
         ! Compressed matrix is in pivot order already
         do i = 1, nsvar
            sperm(i) = i
            sinvp(i) = i
         end do
      end select
   else
      svar_type = 0
      realn = n ! Assume full rank
   endif

   select case(svar_type)
   case(0)
      if(present(piv_size)) then
         call mc78_inner_analyse(n, realn, bptr, brow, perm, invp, nnodes, &
            sptr, sparent, scc, rptr, rlist, control, info, st, &
            block_pivots=piv_size)
      else
         call mc78_inner_analyse(n, realn, ptr, row, perm, invp, nnodes, &
            sptr, sparent, scc, rptr, rlist, control, info, st)
      endif
   case(2)
      call mc78_inner_analyse(nsvar, i, ptr2, row2, sperm, sinvp, nnodes, &
         sptr, sparent, scc, rptr, rlist, control, info, st, wt=svara, &
         block_pivots=piv_size)
      if(st.ne.0) goto 490
      if(info.lt.0) return
      if(i.ne.nsvar) then
         ! Note: This code should NEVER execute
         if(control%unit_error.gt.0) &
            write(control%unit_error, "(a,2(a,i8))") "MC78_ANALYSE Internal ", &
               "Error: supervariable matrix is rank deficient: i = ", i, &
               "nsvar = ", nsvar
         info = MC78_ERROR_UNKNOWN
         return
      endif
      call svar_unmap(n, nsvar, svara, perm, invp, nnodes, sinvp, sptr, st)
   end select
   if(st.ne.0) goto 490
   if(info.lt.0) return

   ! Calculate info%num_factor and info%num_flops
   call mc78_stats(nnodes, sptr, scc, nfact=nfact, nflops=nflops)

   if(control%lopt) then
      ! Reorder elimination variables for better cache locality
      if(control%sort) then
         call mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, &
            sparent, rptr, rlist, st, sort=.true.)
      else
         call mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, &
            sparent, rptr, rlist, st)
      endif
      if(st.ne.0) goto 490
   else if(control%sort) then
      call dbl_tr_sort(n, nnodes, rptr, rlist, st)
      if(st.ne.0) goto 490
   endif

   if(present(piv_size)) &
      call convert_from_blk_piv(n, perm, piv_size)

   !print *, "n, realn = ", n, realn
   !print *, "ptr = ", ptr
   !do i = 1, n
   !   print *, "row(", i, ") = ", row(ptr(i):ptr(i+1)-1)
   !end do
   !print *, "nnodes = ", nnodes
   !print *, "sptr = ", sptr(1:nnodes+1)
   !print *, "sparent = ", sparent(1:nnodes)
   !print *, "rptr = ", rptr(1:nnodes+1)
   !print *, "rlist = ", rlist(1:rptr(nnodes+1)-1)
   !print *, "invp = ", invp
   !print *, "perm = ", perm
   !print *, "piv_size = ", piv_size

   return

   !!!!!!!!!!!!!!!!!!
   ! Error handlers !
   !!!!!!!!!!!!!!!!!!

   490 continue
   info = MC78_ERROR_ALLOC
   if(present(stat)) stat = st
   return
end subroutine mc78_analyse_assembled_integer

!
! Inner core for assembled analyse routine, used to make calls in compressed
! (supervariable) case and standard case uniform
!
subroutine mc78_inner_analyse(n, realn, ptr, row, perm, invp, nnodes, sptr, &
      sparent, scc, rptr, rlist, control, info, st, wt, block_pivots)
   integer, intent(in) :: n ! Dimension of system
   integer, intent(out) :: realn ! Symbolic dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, intent(out) :: nnodes ! number of supernodes found
   integer, dimension(:), allocatable, intent(out) :: sptr ! supernode pointers
   integer, dimension(:), allocatable, intent(out) :: sparent ! assembly tree
   integer, dimension(:), allocatable, intent(out) :: scc ! supernodal col cnt
   integer(long), dimension(:), allocatable, intent(out) :: rptr
      ! pointers to rlist
   integer, dimension(:), allocatable, intent(out) :: rlist ! row lists
   ! For details of control, info : see derived type descriptions
   type(mc78_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(out) :: st
   integer, dimension(n), optional, intent(in) :: wt ! Weights of columns
      ! (i.e. size of each supervariable they represent)
   integer, dimension(n), optional, intent(inout) :: block_pivots ! If
      ! present, then matches pivot order and specifies block pivots. 
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot

   integer, dimension(:), allocatable :: cc ! number of entries in each column
   integer :: ntot ! total number of variables
   integer, dimension(:), allocatable :: parent ! parent of each node in etree
   integer, dimension(:), allocatable :: tperm ! temporary permutation vector

   ! Build elimination tree
   allocate(parent(n), stat=st)
   if(st.ne.0) return
   call mc78_etree(n, ptr, row, perm, invp, parent, st)
   if(st.ne.0) return

   ! Postorder tree (modifies perm!)
   call mc78_postorder(n, realn, ptr, perm, invp, parent, st, block_pivots)
   if(st.ne.0) return

   if(n.ne.realn) then
      if(control%ssa_abort) then
         if(control%unit_error.gt.0) write(control%unit_error, "(2a)") &
            "HSL_MC78: Error, matrix is symbolically singular and ", &
            "control%ssa_abort=.true.."
         info = MC78_ERROR_SSA
         return
      else
         if(control%unit_warning.gt.0) write(control%unit_warning, "(a)") &
            "HSL_MC78: Warning, matrix is symbolically singular."
         info = info + MC78_WARNING_SSA
      endif
   endif

   ! Determine column counts
   allocate(cc(n+1), stat=st)
   if(st.ne.0) return
   call mc78_col_counts(n, ptr, row, perm, invp, parent, cc, st, wt=wt)
   if(st.ne.0) return

   ! Identify supernodes
   allocate(tperm(n), sptr(n+1), sparent(n), scc(n), stat=st)
   if(st.ne.0) return
   call mc78_supernodes(n, realn, parent, cc, tperm, nnodes, sptr, sparent, &
      scc, invp, control, info, st, wt=wt, block_pivots=block_pivots)
   if(info.lt.0) return

   ! Apply permutation to obtain final elimination order
   call apply_perm(n, tperm, perm, invp, cc, block_pivots=block_pivots)

   ! Determine column patterns - keep%nodes(:)%index
   allocate(rptr(nnodes+1), rlist(sum(scc(1:nnodes))), stat=st)
   if(st.ne.0) return
   if(present(wt)) then
      ntot = sum(wt)
      call mc78_row_lists(n, wt, ntot, ptr, row, perm, invp, nnodes, sptr, &
         sparent, scc, rptr, rlist, control, info, st)
   else
      call mc78_row_lists(n, ptr, row, perm, invp, nnodes, sptr, &
         sparent, scc, rptr, rlist, control, info, st)
   endif
   if(st.ne.0) return
end subroutine mc78_inner_analyse

!
! This subroutine performs full analyse when A is in elemental form.
! This is essentially a wrapper around the rest of the package.
!
subroutine mc78_analyse_elemental_integer(n, nelt, starts, vars, perm, &
      eparent, nnodes, sptr, sparent, rptr, rlist, control, info, stat, &
      nfact, nflops, piv_size)
   integer, intent(in) :: n ! Maximum integer used to index an element
   integer, intent(in) :: nelt ! Number of elements
   integer(pkg_type), dimension(nelt+1), intent(in) :: starts ! Element pointers
   integer, dimension(starts(nelt+1)-1), intent(in) :: vars ! Variables
      !assoicated with each element. Element i has vars(starts(i):starts(i+1)-1)
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(nelt), intent(out) :: eparent ! On exit, eparent(i) holds
      ! node of assembly that element i is a child of.
   integer, intent(out) :: nnodes ! number of supernodes found
   integer, dimension(:), allocatable, intent(out) :: sptr ! supernode pointers
   integer, dimension(:), allocatable, intent(out) :: sparent ! assembly tree
   integer(long), dimension(:), allocatable, intent(out) :: rptr
      ! pointers to rlist
   integer, dimension(:), allocatable, intent(out) :: rlist ! row lists
   ! For details of control, info : see derived type descriptions
   type(mc78_control), intent(in) :: control
   integer, intent(out) :: info
   integer, optional, intent(out) :: stat
   integer(long), optional, intent(out) :: nfact ! If present, then on exit
      ! contains the number of entries in L
   integer(long), optional, intent(out) :: nflops ! If present, then on exit
      ! contains the number of floating point operations in factorize.
   integer, dimension(n), optional, intent(inout) :: piv_size ! If
      ! present, then matches matrix order and specifies block pivots. 
      ! piv_size(i) is the number of entries pivots in the block pivot
      ! containing column i.

   integer, dimension(:), allocatable :: cc ! number of entries in each column
   integer :: i
   integer, dimension(:), allocatable :: invp ! inverse permutation of perm
   integer :: j
   integer :: nsvar ! number of supervariables
   integer, dimension(:), allocatable :: parent ! parent of each node in etree
   integer, dimension(:), allocatable :: perm2 ! temporary permutation vector
   integer(pkg_type), dimension(:), allocatable :: ptr ! column pointers for
      ! equivilent matrix
   integer :: realn ! Set to actual number of variables present
   integer, dimension(:), allocatable :: row ! row indices for equivilent matrix
   integer :: st ! stat argument in allocate calls
   integer, dimension(:), allocatable :: svara ! array for supervariables
   integer, dimension(:), allocatable :: sinvp ! inverse permutation of svars
   integer, dimension(:), allocatable :: sperm ! permutation vector of svars
   integer(long) :: sz ! temporary var for size of arrays at allocation

   integer, dimension(:), allocatable :: scc

   ! initialise
   info = 0

   ! Ensure allocatable output arguments are deallocated
   deallocate(sptr, stat=st)
   deallocate(sparent, stat=st)
   deallocate(rptr, stat=st)
   deallocate(rlist, stat=st)

   ! Initialise inverse permutation and check for duplicates
   allocate(invp(n), stat=st)
   if(st.ne.0) goto 490
   do i = 1, n
      j = perm(i)
      invp(j) = i
   end do

   if(present(piv_size)) &
      call convert_to_blk_piv(n, invp, piv_size)

   ! Determine supernodes, build equivilant lwr matrix and find elimination tree
   sz = starts(nelt+1)-1
   if(present(piv_size)) sz = sz + n
   allocate(ptr(n+2), row(sz), svara(n+1), parent(n), stat=st)
   if(st.ne.0) goto 490
   realn = n
   call mc78_elt_equiv_etree(realn, nelt, starts, vars, perm, invp, nsvar, &
      svara, ptr, row, eparent, parent, st, block_pivots=piv_size)
   if(st.ne.0) goto 490

   ! Set up permutations of supervariables (initially the idenity)
   allocate(sperm(nsvar), sinvp(nsvar), stat=st)
   if(st.ne.0) goto 490
   sperm(1:nsvar) = (/ (i, i=1,nsvar) /)
   sinvp(1:nsvar) = (/ (i, i=1,nsvar) /)

   ! Postorder tree (modifies perm!)
   call mc78_postorder(nsvar, sperm, sinvp, parent, st, &
      block_pivots=piv_size)
   if(st.ne.0) goto 490

   ! Determine column counts
   allocate(cc(nsvar+1), stat=st)
   if(st.ne.0) goto 490
   call mc78_col_counts(nsvar, ptr, row, sperm, sinvp, parent, cc, st, wt=svara)
   if(st.ne.0) goto 490

   ! Identify supernodes
   allocate(perm2(nsvar), sptr(nsvar+1), sparent(nsvar), scc(nsvar), stat=st)
   if(st.ne.0) goto 490
   call mc78_supernodes(nsvar, nsvar, parent, cc, perm2, nnodes, sptr, &
      sparent, scc, sinvp, control, info, st, wt=svara, block_pivots=piv_size)
   if(info.eq.MC78_ERROR_ALLOC) goto 490
   if(info.lt.0) return

   ! Apply permutation to obtain final elimination order
   call apply_perm(nsvar, perm2, sperm, sinvp, cc, block_pivots=piv_size)

   ! Determine column patterns - keep%nodes(:)%index
   allocate(rptr(nnodes+1), rlist(sum(scc(1:nnodes))), stat=st)
   if(st.ne.0) goto 490
   call mc78_row_lists(nsvar, svara, n, ptr, row, sperm, sinvp, &
      nnodes, sptr, sparent, scc, rptr, rlist, control, info, st)
   if(st.ne.0) goto 490

   ! Unmap from supervariables to real variables
   call svar_unmap(n, nsvar, svara, perm, invp, nnodes, sinvp, &
      sptr, st)
   if(st.ne.0) goto 490

   ! Calculate info%num_factor and info%num_flops
   call mc78_stats(nnodes, sptr, scc, nfact=nfact, nflops=nflops)

   if(control%lopt) then
      ! Reorder elimination variables for better cache locality
      if(control%sort) then
         call mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, &
            sparent, rptr, rlist, st, sort=.true.)
      else
         call mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, &
            sparent, rptr, rlist, st)
      endif
      if(st.ne.0) goto 490
   else if(control%sort) then
      call dbl_tr_sort(n, nnodes, rptr, rlist, st)
      if(st.ne.0) goto 490
   endif

   ! Adjust eparent with final permuation. Use svara to contain a mapping
   ! from original variables to supernodes
   do i = 1, nnodes
      do j = sptr(i), sptr(i+1)-1
         svara(invp(j)) = i
      end do
   end do
   svara(n+1) = n+1
   do i = 1, nelt
      eparent(i) = svara(eparent(i))
   end do

   if(present(piv_size)) &
      call convert_from_blk_piv(n, perm, piv_size)

   !print *, "n, realn, nelt = ", n, realn, nelt
   !print *, "starts = ", starts
   !print *, "vars = ", vars
   !print *, "nnodes = ", nnodes
   !print *, "sptr = ", sptr(1:nnodes+1)
   !print *, "sparent = ", sparent(1:nnodes)
   !print *, "eparent = ", eparent(:)
   !print *, "rptr = ", rptr(1:nnodes+1)
   !print *, "rlist = ", rlist(1:rptr(nnodes+1)-1)
   !print *, "invp = ", invp
   !print *, "perm = ", perm

   return

   !!!!!!!!!!!!!!!!!!
   ! Error handlers !
   !!!!!!!!!!!!!!!!!!

   490 continue
   info = MC78_ERROR_ALLOC
   if(present(stat)) stat = st
   return
end subroutine mc78_analyse_elemental_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Supervariable routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine find supervariables of A using the algorithm of [1].
!
! [1] I.S. Duff and J.K. Reid. "Exploiting zeros on the diagonal in the direct
!     solution of indefinite sparse symmetric linear systems".
!     ACM TOMS 22(2), pp227-257. 1996.
!
subroutine mc78_supervars_integer(n, ptr, row, perm, invp, nsvar, svar, st)
   integer, intent(inout) :: n ! Dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, intent(out) :: nsvar ! number of supervariables
   integer, dimension(n), intent(out) :: svar ! number of vars in each svar
   integer, intent(out) :: st

   logical :: full_rank ! flags if supervariable 1 has ever become empty.
      ! If it has not, then the varaibles in s.v. 1 are those that never
      ! occur
   integer :: i
   integer(long) :: ii
   integer :: j
   integer :: idx ! current index
   integer :: next_sv ! head of free sv linked list
   integer :: nsv ! new supervariable to move j to
   integer :: piv ! current pivot
   integer :: col ! current column of A
   integer :: sv ! current supervariable
   integer :: svc ! temporary holding supervariable count
   integer, dimension(:), allocatable :: sv_new  ! Maps each supervariable to
      ! a new supervariable with which it is associated.
   integer, dimension(:), allocatable :: sv_seen ! Flags whether svariables have
      ! been seen in the current column. sv_seen(j) is set to col when svar j
      ! has been encountered.
   integer, dimension(:), allocatable :: sv_count ! number of variables in sv.

   allocate(sv_new(n+1), sv_seen(n+1), sv_count(n+1), stat=st)
   if(st.ne.0) return

   svar(:) = 1
   sv_count(1) = n
   sv_seen(1) = 0

   ! Setup linked list of free super variables
   next_sv = 2
   do i = 2, n
      sv_seen(i) = i+1
   end do
   sv_seen(n+1) = -1

   ! Determine supervariables using modified Duff and Reid algorithm
   full_rank = .false.
   do col = 1, n
      if(ptr(col+1).ne.ptr(col)) then
         ! If column is not empty, add implicit diagonal entry
         j = col
         sv = svar(j)
         if(sv_count(sv).eq.1) then ! Are we only (remaining) var in sv
            full_rank = full_rank .or. (sv.eq.1)
            ! MUST BE the first time that sv has been seen for this
            ! column, so just leave j in sv, and go to next variable.
            ! (Also there can be no other vars in this block pivot)
         else
            ! There is at least one other variable remaining in sv
            ! MUST BE first occurence of sv in the current row/column,
            ! so define a new supervariable and associate it with sv.
            sv_seen(sv) = col
            sv_new(sv) = next_sv
            nsv = next_sv
            next_sv = sv_seen(next_sv)
            sv_new(nsv) = nsv ! avoids problems with duplicates
            sv_seen(nsv) = col
            ! Now move j from sv to nsv
            nsv = sv_new(sv)
            svar(j) = nsv
            sv_count(sv) = sv_count(sv) - 1
            sv_count(nsv) = 1
            ! This sv cannot be empty as initial sv_count was > 1
         endif
      endif
      do ii = ptr(col), ptr(col+1)-1
         j = row(ii)
         sv = svar(j)
         if(sv_count(sv).eq.1) then ! Are we only (remaining) var in sv
            full_rank = full_rank .or. (sv.eq.1)
            ! If so, and this is first time that sv has been seen for this
            ! column, then we can just leave j in sv, and go to next variable.
            if(sv_seen(sv).lt.col) cycle
            ! Otherwise, we have already defined a new supervariable associated
            ! with sv. Move j to this variable, then retire (now empty) sv.
            nsv = sv_new(sv)
            if(sv.eq.nsv) cycle
            svar(j) = nsv
            sv_count(nsv) = sv_count(nsv) + 1
            ! Old sv is now empty, add it to top of free stack
            sv_seen(sv) = next_sv
            next_sv = sv
         else
            ! There is at least one other variable remaining in sv
            if(sv_seen(sv).lt.col) then
               ! this is the first occurence of sv in the current row/column,
               ! so define a new supervariable and associate it with sv.
               sv_seen(sv) = col
               sv_new(sv) = next_sv
               sv_new(next_sv) = next_sv ! avoids problems with duplicates
               next_sv = sv_seen(next_sv)
               sv_count(sv_new(sv)) = 0
               sv_seen(sv_new(sv)) = col
            endif
            ! Now move j from sv to nsv
            nsv = sv_new(sv)
            svar(j) = nsv
            sv_count(sv) = sv_count(sv) - 1
            sv_count(nsv) = sv_count(nsv) + 1
            ! This sv cannot be empty as sv_count was > 1
         endif
      end do
   end do

   ! Note: block pivots do not mix well with supervariables as any significant
   ! number (unless aligned to s.v.s) will demolish any gain from using them.
   ! Converting vlock pivots to s.v.s results in potentially large amount of
   ! unneeded fillin to left of block pivot.
   !! If block pivots are being used, we force all pivots of a block pivot
   !! to be in either the same supervariable, or in supervariables of size 1
   !if(present(block_pivots)) then
   !   piv = 1
   !   do while(piv.le.n)
   !      ! Check if we need to split pivots
   !      split = .false.
   !      sv = svar(piv)
   !      do i = piv+1, piv+block_pivots(piv)
   !         j = invp(i)
   !         if(svar(j).ne.sv) then
   !            split = .true.
   !            exit
   !         endif
   !      end do
   !      ! Do split if required
   !      if(split) then
   !         j = invp(i)
   !         do i = piv, piv+block_pivots(piv)
   !            sv = svar(j)
   !            if(sv_count(sv).eq.1) cycle ! Already a singleton
   !            ! Otherwise create a new sv and move j to it
   !            nsv = next_sv
   !            next_sv = sv_seen(next_sv)
   !            svar(j) = nsv
   !            sv_count(nsv) = 1
   !            sv_count(sv) = sv_count(sv) - 1
   !         end do
   !      endif
   !      piv = piv + block_pivots(piv) + 1
   !   end do
   !endif

   ! Now modify pivot order such that all variables in each supervariable are
   ! consecutive. Do so by iterating over pivots in elimination order. If a
   ! pivot has not already been listed, then order that pivot followed by
   ! any other pivots in that supervariable.

   ! We will build a new inverse permutation in invp, and then find perm
   ! afterwards. First copy invp to perm:
   perm(:) = invp(:)
   ! Next we iterate over the pivots that have not been ordered already
   ! Note: as we begin, all entries of sv_seen are less than or equal to n+1
   ! hence we can use <=n+1 or >n+1 as a flag to indicate if a variable has been
   ! ordered.
   idx = 1
   nsvar = 0
   do piv = 1, n
      if(sv_seen(piv).gt.n+1) cycle ! already ordered
      ! Record information for supervariable
      sv = svar(perm(piv))
      if(.not.full_rank .and. sv.eq.1) cycle ! Don't touch unused vars
      nsvar = nsvar + 1
      svc = sv_count(sv)
      sv_new(nsvar) = svc ! store # vars in s.v. to copy to svar later
      j = piv
      ! Find all variables that are members of sv and order them.
      do while(svc.gt.0)
         do j = j, n
            if(svar(perm(j)).eq.sv) exit
         end do
         sv_seen(j) = n+2 ! flag as ordered
         invp(idx) = perm(j)
         idx = idx + 1
         svc = svc - 1
         j = j + 1
      end do
   end do
   ! Push unused variables to end - these are those vars still in s.v. 1
   if(.not.full_rank) then
      svc = sv_count(1)
      ! Find all variables that are members of sv and order them.
      j = 1
      do while(svc.gt.0)
         do j = j, n
            if(svar(perm(j)).eq.1) exit
         end do
         invp(idx) = perm(j)
         idx = idx + 1
         svc = svc - 1
         j = j + 1
      end do
      n = n - sv_count(1)
   end if
   ! Recover perm as inverse of invp
   do piv = 1, n
      perm(invp(piv)) = piv
   end do
   ! sv_new has been used to store number of variables in each svar, copy into
   ! svar where it is returned.
   svar(1:nsvar) = sv_new(1:nsvar)
end subroutine mc78_supervars_integer

!
! This subroutine takes a set of supervariables and compresses the supplied
! matrix using them.
!
! As we would need a full scan of the matrix to calculate the correct size of
! row2, we instead allow the user to make a guess at a good size and return
! an error if this turns out to be incorrect. An upper bound on the required
! size may be obtained by summing the number of entries in the first column of
! each supervariable.
!
! Error returns:
!   MC78_ERROR_ALLOC      Failed to allocate memory
!   MC78_ERROR_ROW_SMALL  row2 too small
subroutine mc78_compress_by_svar_integer(n, ptr, row, invp, nsvar, svar, ptr2, &
      lrow2, row2, info, st)
   integer, intent(in) :: n ! Dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(in) :: invp ! inverse of perm
   integer, intent(in) :: nsvar
   integer, dimension(nsvar), intent(in) :: svar ! super variables of A
   integer(pkg_type), dimension(nsvar+1), intent(out) :: ptr2
   integer(pkg_type), intent(in) :: lrow2
   integer, dimension(lrow2), intent(out) :: row2
   integer, intent(out) :: info
   integer, intent(out) :: st

   integer :: piv, svc, sv, col
   integer(pkg_type) :: j, idx
   integer, dimension(:), allocatable :: flag, sv_map

   info = 0 ! by default completed succefully

   allocate(flag(nsvar), sv_map(n), stat=st)
   if(st.ne.0) then
      info = MC78_ERROR_ALLOC
      return
   endif
   flag(:) = 0

   ! Setup sv_map
   piv = 1
   do svc = 1, nsvar
      do piv = piv, piv + svar(svc) - 1
         sv_map( invp(piv) ) = svc
      end do
   end do

   piv = 1
   idx = 1
   do svc = 1, nsvar
      col = invp(piv)
      ptr2(svc) = idx
      do j = ptr(col), ptr(col+1)-1
         sv = sv_map(row(j))
         if(flag(sv).eq.piv) cycle ! Already dealt with this supervariable
         if(idx.gt.lrow2) then
            ! oops, row2 is too small
            info = MC78_ERROR_ROW_SMALL
            return
         endif
         ! Add row entry for this sv
         row2(idx) = sv
         flag(sv) = piv
         idx = idx + 1
      end do
      piv = piv + svar(svc)
   end do
   ptr2(svc) = idx
end subroutine mc78_compress_by_svar_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Elimination tree routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine determines the elimination tree of a PAP^T where A is a
! sparse symmetric matrix stored in compressed sparse column form with
! entries both above and below the diagonal present in the argument matrix.
! P is a permutation stored in order such that order(i) gives the pivot
! position of column i. i.e. order(3) = 5 means that the fifth pivot is
! A_33.
!
! The elimination tree is returned in the array parent. parent(i) gives the
! parent in the elimination tree of pivot i.
!
! The algorithm used is that of Liu [1].
!
! [1] Liu, J. W. 1986. A compact row storage scheme for Cholesky factors using
!     elimination trees. ACM TOMS 12, 2, 127--148.
!
subroutine mc78_etree_integer(n, ptr, row, perm, invp, parent, st)
   integer, intent(in) :: n ! dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! column pointers of A
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! row indices of A
   integer, dimension(n), intent(in) :: perm ! perm(i) is the pivot position
      ! of column i
   integer, dimension(n), intent(in) :: invp ! inverse of perm
   integer, dimension(n), intent(out) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls

   integer(pkg_type) :: i ! next index into row
   integer :: j ! current entry in row
   integer :: k ! current ancestor
   integer :: l ! next ancestor
   integer :: piv ! current pivot
   integer :: rowidx ! current column of A = invp(piv)
   integer, dimension(:), allocatable :: vforest ! virtual forest, used for
      ! path compression (shortcuts to top of each tree)

   ! Allocate virtual forest and initialise it
   allocate(vforest(n), stat=st)
   if(st.ne.0) return
   vforest(:) = n+1

   ! Loop over rows of A in pivot order
   piv = 1
   do while(piv.le.n)
      !print *, "row ", piv
      rowidx = invp(piv)
      ! Loop over entries in row in lower triangle of PAP^T
      do i = ptr(rowidx), ptr(rowidx+1)-1
         j = perm(row(i))
         if(j.ge.piv) cycle ! not in lower triangle
         !print *, "  entry ", j
         k = j
         do while(vforest(k).lt.piv)
            l = vforest(k)
            vforest(k) = piv
            k = l
         end do
         ! Check if we have already done this pivot
         if(vforest(k).eq.piv) cycle 
         parent(k) = piv
         vforest(k) = piv
      end do
      parent(piv) = n + 1 ! set to be a root if not overwritten
      piv = piv + 1 ! move on to next pivot
   end do
end subroutine mc78_etree_integer

!
! This subroutine identifies supervariables of A using a modified variant of
! the algorithm of Duff and Reid [1]. A lower triangular equivilant matrix
! is returned that is expressed in terms of these supervariables. The grouping
! of variables into supervaribles is returned through a modified pivot order
! and an array specifying the number of variables in each supervariable in
! elimination order. Finally the vector eparent is also returned. This contains
! the variable (in natural numbering) that corresponds to the least pivot in
! each supervariable.
!
! [1] I.S. Duff and J.K. Reid. "Exploiting zeros on the diagonal in the direct
!     solution of indefinite sparse symmetric linear systems".
!     ACM TOMS 22(2), pp227-257. 1996.
!
! Note: If block pivots are present they have priority over supervariables - 
! members of same block pivot must remain in same supervariable. This is
! enforced by moving them all at once.
subroutine mc78_elt_equiv_etree_integer(n, nelt, starts, vars, perm, invp, &
      nsvar, svar, ptr, row, eparent, parent, st, block_pivots)
   integer, intent(inout) :: n ! dimension of system
   integer, intent(in) :: nelt ! number of elements
   integer(pkg_type), dimension(nelt+1), intent(in) :: starts ! variable
      ! pointers of elements
   integer, dimension(starts(nelt+1)-1), intent(in) :: vars ! variables of
      ! elements
   integer, dimension(n), intent(inout) :: perm ! perm(i) is the pivot position
      ! of column i
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, intent(out) :: nsvar ! number of supervariables found
   integer, dimension(n), intent(out) :: svar ! size of each supervariable
   integer(pkg_type), dimension(n+1), intent(out) :: ptr ! column pointers
      ! for equivilant lower triangular form
   integer, dimension(:), intent(out) :: row ! row indices
      ! for equivilant lower triangular form
   integer, dimension(nelt), intent(out) :: eparent ! parent nodes of each
      ! element - i.e. the least pivot in each element
   integer, dimension(n), intent(out) :: parent ! parent(i) is parent of node
      ! i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls
   integer, dimension(n), optional, intent(inout) :: block_pivots ! Used for
      ! block pivots, see description in analyse phase.

   integer :: csv ! column supervariable in loop
   integer :: elt ! current element
   logical :: full_rank ! flags if supervariable 1 has ever become empty.
      ! If it has not, then the varaibles in s.v. 1 are those that never
      ! occour
   integer :: i
   integer(pkg_type) :: ii
   integer(pkg_type) :: idx ! Next insert position
   integer :: j
   integer :: k
   integer :: minpiv ! minimum pivot of current element
   integer, dimension(:), allocatable :: mp_head, mp_next ! mp_head and mp_next
      ! store a linked list of elements for which a given variable is the
      ! minimum pivot.
   integer :: next_sv ! Top of stack of free supervariables (stored as a linked
      ! list in unused part of sv_seen)
   integer :: orign ! original system dimension
   integer :: nsv ! temporary variable storing new supervariable to move var to
   integer :: piv ! current pivot
   integer :: sv ! current supervariable
   integer :: svc ! temporary variable storing supervariable count remaining
   integer, dimension(:), allocatable :: sv_count ! sv_count(s) is the number
      ! of variables in supervariable s.
   integer, dimension(:), allocatable :: sv_map ! sv_map(v) is the current
      ! supervariable to which variable v belongs.
   integer, dimension(:), allocatable :: sv_new ! sv_map(s) is new
      ! supervariable for variables currently in supervariable s.
   integer, dimension(:), allocatable :: sv_seen ! sv_seen(s) is used to flag
      ! if supervariable s has been found in the current element before.
      ! In addition the part corresponding to unused supervariables is used
      ! to store a stack (as a linked list) of empty supervariables.
   integer(pkg_type), dimension(:), allocatable :: uprptr ! column pointers for
      ! upper triangular equivilant form.
   integer, dimension(:), allocatable :: uprrow ! row indices for upper
      ! triangular equivilant form.
   logical :: used ! flag if a variable has been used

   ! Initialise supervariable representation
   allocate(sv_new(max(nelt+1,n+1)), sv_seen(max(nelt,n)+1), &
      sv_map(max(nelt,n)+1), sv_count(max(nelt,n)+1), stat=st)
   if(st.ne.0) return
   sv_map(:) = 1 ! All vars are intially in supervariable 1
   sv_count(1) = n ! ... which thus has all variables
   sv_seen(1) = 0 ! Flag supervariable 1 as unseen on first iteration
   orign = n

   if(present(block_pivots)) then
      ! Do not mix supervariables and block pivots

      ! Still need to determine minimum pivots and rank
      sv_seen(:) = 0
      do elt = 1, nelt
         minpiv = n+1
         do ii = starts(elt), starts(elt+1)-1
            j = vars(ii)
            ! Mark variable as used
            sv_seen(j) = 1
            ! Determine minimum pivot
            minpiv = min(minpiv, perm(j))
         end do
         ! Store the minimum pivot as original variable (avoids messy remapping)
         if(minpiv.le.n) eparent(elt) = invp(minpiv)
      end do

      ! Build invp that pushes unsued vars to the end. Be careful of unused vars
      ! that are in fact part of block pivots, split them out.
      perm(:) = invp(:)
      piv = 1
      j = 1
      ! Handle variables that are actually used
      do while(piv.le.n)
         used = .true.
         do i = piv, n
            used = used .and. (sv_seen(perm(i)).eq.1)
            if(block_pivots(i).ge.2) exit ! end of block pivot
         end do

         if(used) then
            ! Block pivot is entirely composed of used variables
            do piv = piv, i
               invp(j) = perm(piv)
               sv_seen(perm(piv)) = 2
               j = j + 1
            end do
         else
            ! Block pivot has some unused variables in it
            k = 0
            do piv = piv, i
               if(sv_seen(perm(piv)).eq.1) then
                  invp(j) = perm(piv)
                  sv_seen(perm(piv)) = 2
                  j = j + 1
                  if(k.eq.0) then
                     ! This is the new start of the block pivot
                     select case(block_pivots(piv))
                     case(0) ! was in the middle. now a start
                        block_pivots(piv) = 1
                     case(2) ! was the end. now a 1x1
                        block_pivots(piv) = 3
                     end select
                  endif
                  k = piv
               endif
            end do
            if(k.ne.0) then
               ! The was at least one used variable in the block pivot
               select case(block_pivots(k))
               case(0) ! was the middle. now an end
                  block_pivots(k) = 2
               case(1) ! was the start. now a 1x1
                  block_pivots(k) = 3
               end select
            endif
         endif
         piv = i + 1
      end do
      ! Handle unused variables; build supervariables
      nsvar = 0
      do piv = 1, n
         i = perm(piv)
         if(sv_seen(i).eq.0) then
            invp(j) = i
            j = j + 1
            block_pivots(piv) = 3 ! Force to 1x1
         else
            nsvar = nsvar + 1
            svar(nsvar) = 1
            sv_new(i) = nsvar
         endif
      end do
      ! Map block_pivots in original variable order into sv_map
      do i = 1, n
         sv_map(perm(i)) = block_pivots(i)
      end do
      ! Map sv_map in new pivot order back into block_pivots
      do i = 1, n
         block_pivots(i) = sv_map(invp(i))
      end do
      ! Reestablish perm
      do i = 1, n
         perm(invp(i)) = i
      end do
   else
      ! Setup linked list of free supervariables - we utilise the unused part
      ! of sv_seen for this
      next_sv = 2
      do i = 2, n
         sv_seen(i) = i+1
      end do
      sv_seen(n+1) = -1
      
      ! Determine supervariables. At the same time find the least pivot
      ! associated with each element.
      nsvar = 1
      full_rank = .false.
      do elt = 1, nelt
         minpiv = n+1
         do ii = starts(elt), starts(elt+1)-1
            j = vars(ii)
            ! Determine minimum pivot
            minpiv = min(minpiv, perm(j))
            sv = sv_map(j)
            if(sv_count(sv).eq.1) then ! Are we only (remaining) var in sv
               full_rank = full_rank .or. (sv.eq.1)
               ! If so, and this is first time that sv has been seen for this
               ! element, then we can just leave j in sv, and go to next
               ! variable.
               if(sv_seen(sv).lt.elt) cycle
               ! Otherwise, we have already defined a new supervariable
               ! associated with sv. Move j to this variable, then retire (now
               ! empty) sv.
               ! Note: as only var in sv, cannot have fellows in block pivot
               nsv = sv_new(sv)
               if(sv.eq.nsv) cycle  ! don't delete a variable because of a
                                    ! duplicate
               sv_map(j) = nsv
               sv_count(nsv) = sv_count(nsv) + 1
               ! Old sv is now empty, add it to top of free stack
               sv_seen(sv) = next_sv
               next_sv = sv
               nsvar = nsvar - 1
            else
               ! There is at least one other variable remaining in sv
               if(sv_seen(sv).lt.elt) then
                  ! this is the first occurence of sv in the current element,
                  ! so define a new supervariable and associate it with sv.
                  sv_seen(sv) = elt
                  sv_new(sv) = next_sv
                  sv_new(next_sv) = next_sv  ! ensure we are tolerant of
                                             ! duplicates
                  next_sv = sv_seen(next_sv)
                  sv_count(sv_new(sv)) = 0
                  sv_seen(sv_new(sv)) = elt
                  nsvar = nsvar + 1
               endif
               ! Now move j from sv to nsv
               nsv = sv_new(sv)
               sv_map(j) = nsv
               sv_count(sv) = sv_count(sv) - 1
               sv_count(nsv) = sv_count(nsv) + 1
               ! We know sv can't be empty, so it doesn't need adding to free
               ! stack
            endif
         end do
         ! Store the minimum pivot as original variable (avoids messy remapping)
         if(minpiv.le.n) then
            eparent(elt) = invp(minpiv)
         else
            eparent(elt) = n+1
         endif
      end do

      ! Now modify pivot order such that all variables in each supervariable are
      ! consecutive. Do so by iterating over pivots in elimination order. If a
      ! pivot has not already been listed, then order that pivot followed by
      ! any other pivots in that supervariable.

      ! We will build a new inverse permutation in invp, and then find perm
      ! afterwards. First copy invp to perm:
      perm(:) = invp(:)
      ! Next we iterate over the pivots that have not been ordered already
      ! Note: as we begin, all entries of sv_seen are less than or equal to n+1
      ! hence we can use <=n+1 or >n+1 as a flag to indicate if a variable has
      ! been ordered.
      idx = 1
      nsvar = 0
      do piv = 1, n
         if(sv_seen(piv).gt.n+1) cycle ! already ordered
         ! Record information for supervariable
         sv = sv_map(perm(piv))
         if(.not.full_rank .and. sv.eq.1) cycle ! Don't touch unused vars
         svc = sv_count(sv)
         nsvar = nsvar + 1
         svar(nsvar) = svc
         ! Find all variables that are members of sv and order them.
         j = piv
         do while(svc.gt.0)
            do j = j, n
               if(sv_map(perm(j)).eq.sv) exit
            end do
            sv_seen(j) = n+2 ! flag as ordered
            sv_new(perm(j)) = nsvar ! new mapping to sv
            invp(idx) = perm(j)
            idx = idx + 1
            svc = svc - 1
            j = j + 1
         end do
      end do
      sv_new(n+1) = nsvar+1
      ! Push unused variables to end - these are those vars still in s.v. 1
      if(.not.full_rank) then
         svc = sv_count(1)
         ! Find all variables that are members of sv and order them.
         j = 1
         do while(svc.gt.0)
            do j = j, n
               if(sv_map(perm(j)).eq.1) exit
            end do
            invp(idx) = perm(j)
            idx = idx + 1
            svc = svc - 1
            j = j + 1
         end do
         n = n - sv_count(1)
      end if
      ! Recover perm as inverse of invp
      do piv = 1, n
         perm(invp(piv)) = piv
      end do
   endif

   ! build linked lists by supervariable
   allocate(mp_head(nsvar+1), mp_next(nelt), stat=st)
   if(st.ne.0) return
   mp_head(:) = -1
   do elt = 1, nelt
      if(eparent(elt).gt.orign) cycle
      minpiv = sv_new(eparent(elt))
      if(present(block_pivots) .and. minpiv.ne.1) then
         do while(block_pivots(minpiv-1).lt.2)
            minpiv = minpiv - 1
            if(minpiv.eq.1) exit
         end do
      endif
      ! Store element in linked list for minpiv
      mp_next(elt) = mp_head(minpiv)
      mp_head(minpiv) = elt
   end do

   ! Iterate over columns in pivot order, storing the lower triangular
   ! equivilant matrix as we go. At the same time, build the column counts for
   ! the upper triangle in uprptr, but offset by 2 (ie uprptr(i+2) for col i).
   ! Observe that all the pivots associated with the supervariable
   ! to which minpiv belongs _must_ appear in each element that minpiv does.
   ! Note: This only generates the lower triangular part of the matrix!
   allocate(uprptr(nsvar+2), stat=st)
   if(st.ne.0) return
   uprptr(:) = 0
   sv_seen(:) = 0
   idx = 1
   ptr(:) = -1
   do csv = 1, nsvar
      elt = mp_head(csv)
      ptr(csv) = idx
      sv_seen(csv) = csv ! Mark diagonal as seen, as it is implicit.
      do while(elt.ne.-1)
         do ii = starts(elt), starts(elt+1)-1
            sv = sv_new(vars(ii))
            ! Skip this sv if it is already included (or is implicit)
            if(sv_seen(sv).ge.csv) cycle
            sv_seen(sv) = csv ! Mark as seen
            ! If we can't skip it, then add entry (sv,csv) to lwr matrix
            row(idx) = sv
            idx = idx + 1
            ! Add count in upper triangle for (csv,sv)
            uprptr(sv+2) = uprptr(sv+2) + 1
         end do
         ! Move on to next element for which this is the minimum pivot
         elt = mp_next(elt)
      end do
      if(present(block_pivots) .and. csv.ne.nsvar) then
         ! Add entry (csv+1,csv) to ensure elimination tree correct
         if(block_pivots(csv).lt.2 .and. sv_seen(csv+1).ne.csv) then
            sv_seen(csv+1) = csv
            row(idx) = csv+1
            idx = idx + 1
            ! Add count in upper triangle for (csv, csv+1)
            uprptr(csv+1+2) = uprptr(csv+1+2) + 1
         endif
      endif
   end do
   ptr(nsvar+1) = idx

   ! Build upper form - work out column start for col i in uprptr(i+1)
   uprptr(1:2) = 1
   do i = 1, nsvar
      uprptr(i+2) = uprptr(i+1) + uprptr(i+2)
   end do

   ! Now iterate over lwr form, droppping entries into upr form
   allocate(uprrow(uprptr(nsvar+2)), stat=st)
   if(st.ne.0) return
   do csv = 1, nsvar
      do ii = ptr(csv), ptr(csv+1) - 1
         sv = row(ii)
         uprrow(uprptr(sv+1)) = csv
         uprptr(sv+1) = uprptr(sv+1) + 1
      end do
   end do

   ! Now determine supervariable elimination tree
   call etree_no_perm(nsvar, uprptr, uprrow, parent, st)
   if(st.ne.0) return
end subroutine mc78_elt_equiv_etree_integer

!
! Specialised version of mc78_etree that assumes elimination order is identity
!
subroutine etree_no_perm(n, ptr, row, parent, st)
   integer, intent(in) :: n ! dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! column pointers of A
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! row indices of A
   integer, dimension(n), intent(out) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls

   integer(pkg_type) :: ii ! next index into row
   integer :: k ! current ancestor
   integer :: l ! next ancestor
   integer :: piv ! current pivot
   integer, dimension(:), allocatable :: vforest ! virtual forest, used for
      ! path compression (shortcuts to top of each tree)

   ! Allocate virtual forest and initialise it
   allocate(vforest(n), stat=st)
   if(st.ne.0) return
   vforest(:) = n+1

   ! Loop over rows of A in pivot order
   do piv = 1, n
      ! Loop over entries in row in lower triangle of PAP^T
      do ii = ptr(piv), ptr(piv+1)-1
         k = row(ii)
         do while(vforest(k).lt.piv)
            l = vforest(k)
            vforest(k) = piv
            k = l
         end do
         if(vforest(k).eq.piv) cycle ! Already done from here, don't overwrite
         parent(k) = piv
         vforest(k) = piv
      end do
      parent(piv) = n + 1 ! set to be a root if not overwritten
   end do
end subroutine etree_no_perm

!
! This subroutine will postorder the elimination tree. That is to say it will
! reorder the nodes of the tree such that they are in depth-first search order.
!
! This is done by performing a depth-first search to identify mapping from the
! original pivot order to the new one. This map is then applied to order, invp
! and parent to enact the relabelling.
!
subroutine mc78_postorder_std(n, perm, invp, parent, st, block_pivots)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: perm ! perm(i) is the pivot
      ! position of column i
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, dimension(n), intent(inout) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls
   integer, dimension(n), optional, intent(inout) :: block_pivots ! If
      ! present, then matches pivot order and specifies block pivots. 
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot

   integer, dimension(:), allocatable :: chead ! chead(i) is first child of i
   integer, dimension(:), allocatable :: cnext ! cnext(i) is next child of i
   integer :: i
   integer :: id
   integer :: j
   integer, dimension(:), allocatable :: map ! mapping from original pivot
      ! order to new one
   integer :: node
   integer :: shead ! pointer to top of stack
   integer, dimension(:), allocatable :: stack ! stack for depth first search

   !
   ! Build linked lists of children for each node
   !
   allocate(chead(n+1), cnext(n+1), stat=st)
   if(st.ne.0) return
   chead(:) = -1 ! no parent if necessary
   do i = n, 1, -1 ! do in reverse order so they come off in original order
      j = parent(i)
      cnext(i) = chead(j)
      chead(j) = i
   end do

   !
   ! Perform depth first search to build map
   !
   allocate(map(n+1), stack(n), stat=st)
   if(st.ne.0) return
   ! Place virtual root on top of stack
   shead = 1
   stack(shead) = n+1
   id = n + 1 ! next node id
   do while(shead.ne.0)
      ! Get node from top of stack
      node = stack(shead)
      shead = shead - 1

      ! Number it
      map(node) = id
      id = id - 1

      ! Place all its children on the stack such that the last child is
      ! at the top of the stack and first child closest to the bottom
      i = chead(node)
      do while(i.ne.-1)
         shead = shead + 1
         stack(shead) = i
         i = cnext(i)
      end do
   end do

   !
   ! Apply map to perm, invp and parent (and block_pivots if present)
   !

   ! invp is straight forward, use stack as a temporary
   stack(1:n) = invp(1:n)
   do i = 1, n
      j = map(i)
      invp(j) = stack(i)
   end do

   ! perm can be easily done as the inverse of invp
   do i = 1, n
      perm(invp(i)) = i
   end do

   ! parent is done in two stages. The first copies it to stack and permutes
   ! parent(i), but not the locations. i.e. if 1 is a parent of 3, and
   ! map(1)=2 and map(3)=4, then the first stage sets stack(1) = 4.
   ! The second stage then permutes the entries of map back into parent
   do i = 1, n
      stack(i) = map(parent(i))
   end do
   do i = 1, n
      parent(map(i)) = stack(i)
   end do

   ! permute block_pivots if required
   if(present(block_pivots)) then
      stack(1:n) = block_pivots(1:n)
      do i = 1, n
         block_pivots(map(i)) = stack(i)
      end do
   endif
end subroutine mc78_postorder_std

!
! This subroutine will postorder the elimination tree. That is to say it will
! reorder the nodes of the tree such that they are in depth-first search order.
!
! This is done by performing a depth-first search to identify mapping from the
! original pivot order to the new one. This map is then applied to order, invp
! and parent to enact the relabelling.
!
subroutine mc78_postorder_detect(n, realn, ptr, perm, invp, parent, st, &
      block_pivots)
   integer, intent(in) :: n
   integer, intent(out) :: realn
   integer(pkg_type), dimension(n+1), intent(in) :: ptr
   integer, dimension(n), intent(inout) :: perm ! perm(i) is the pivot
      ! position of column i
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, dimension(n), intent(inout) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls
   integer, dimension(n), optional, intent(inout) :: block_pivots ! If
      ! present, then matches pivot order and specifies block pivots. 
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot

   integer, dimension(:), allocatable :: chead ! chead(i) is first child of i
   integer, dimension(:), allocatable :: cnext ! cnext(i) is next child of i
   integer :: i
   integer :: id
   integer :: j
   integer, dimension(:), allocatable :: map ! mapping from original pivot
      ! order to new one
   integer :: node
   integer :: shead ! pointer to top of stack
   integer, dimension(:), allocatable :: stack ! stack for depth first search

   realn = n

   !
   ! Build linked lists of children for each node
   !
   allocate(chead(n+1), cnext(n+1), stat=st)
   if(st.ne.0) return
   chead(:) = -1 ! no parent if necessary
   do i = n, 1, -1 ! do in reverse order so they come off in original order
      j = parent(i)
      cnext(i) = chead(j)
      chead(j) = i
   end do

   !
   ! Perform depth first search to build map
   !
   allocate(map(n+1), stack(n), stat=st)
   if(st.ne.0) return
   ! Place virtual root on top of stack
   shead = 1
   stack(shead) = n+1
   id = n + 1 ! next node id
   do while(shead.ne.0)
      ! Get node from top of stack
      node = stack(shead)
      shead = shead - 1

      ! Number it
      map(node) = id
      id = id - 1

      ! Place all its children on the stack such that the last child is
      ! at the top of the stack and first child closest to the bottom
      if(node.eq.n+1) then
         ! Virtual root node, detect children with no entries at same time
         ! placing those that are empty at the top of the stack
         ! First do those which are proper roots
         i = chead(node)
         do while(i.ne.-1)
            if(ptr(invp(i)+1)-ptr(invp(i)).eq.0) then
               i = cnext(i)
               cycle
            endif
            shead = shead + 1
            stack(shead) = i
            i = cnext(i)
         end do
         ! Second do those which are null roots
         i = chead(node)
         do while(i.ne.-1)
            if(ptr(invp(i)+1)-ptr(invp(i)).ne.0) then
               i = cnext(i)
               cycle
            endif
            realn = realn - 1
            shead = shead + 1
            stack(shead) = i
            i = cnext(i)
         end do
      else ! A normal node
         i = chead(node)
         do while(i.ne.-1)
            shead = shead + 1
            stack(shead) = i
            i = cnext(i)
         end do
      endif
   end do

   !
   ! Apply map to perm, invp and parent (and block_pivots if present)
   !

   ! invp is straight forward, use stack as a temporary
   stack(1:n) = invp(1:n)
   do i = 1, n
      j = map(i)
      invp(j) = stack(i)
   end do

   ! perm can be easily done as the inverse of invp
   do i = 1, n
      perm(invp(i)) = i
   end do

   ! parent is done in two stages. The first copies it to stack and permutes
   ! parent(i), but not the locations. i.e. if 1 is a parent of 3, and
   ! map(1)=2 and map(3)=4, then the first stage sets stack(1) = 4.
   ! The second stage then permutes the entries of map back into parent
   do i = 1, n
      stack(i) = map(parent(i))
   end do
   do i = 1, n
      parent(map(i)) = stack(i)
   end do

   ! permute block_pivots if required
   if(present(block_pivots)) then
      stack(1:n) = block_pivots(1:n)
      do i = 1, n
         block_pivots(map(i)) = stack(i)
      end do
   endif
end subroutine mc78_postorder_detect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Column count routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine determines column counts given the elimination tree and
! pattern of the matrix PAP^T.
!
! The algorithm is a specialisation of that given by Gilbert, Ng and Peyton [1],
! to only determine column counts. It is also described in Section 4.4 "Row
! counts" of [2].
!
! The essential technique is to determine the net number of entries introduced
! at a node (the "weight" in [1]). This is composed over the following terms:
!  wt[i] = [ - #children of node i
!            - #common indices between children
!            + #additional "new" row indices from column of A ]
!
! The clever part of this algorithm is how to determine the number of common
! indices between the children. This is accomplished by storing the last column
! at which an index was encountered, and a partial elimination tree. This
! partial elimination tree consists of all nodes processed so far, plus their
! parents. As we have a postorder on the tree, the current top of the tree
! containing node i is the least common ancestor of node i and the current node.
! We then observe that the first time an index will be double counted is at the
! least common ancestor of the current node and the last node where it was
! encountered.
!
! [1] Gilbert, Ng, Peyton, "An efficient algorithm to compute row and column
!     counts for sparse Cholesky factorization", SIMAX 15(4) 1994.
!
! [2] Tim Davis's book "Direct Methods for Sparse Linear Systems", SIAM 2006.
!
subroutine mc78_col_counts_integer(n, ptr, row, perm, invp, parent, cc, st, wt)
   integer, intent(in) :: n ! dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! column pointers of A
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! row indices of A
   integer, dimension(n), intent(in) :: perm ! perm(i) is the pivot
      ! position of column i
   integer, dimension(n), intent(in) :: invp ! inverse of perm
   integer, dimension(n), intent(in) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, dimension(n+1), intent(out) :: cc ! On exit, cc(i) is the
      ! number of entries in the lower triangular part of L (includes diagonal)
      ! for the column containing pivot i. For most of the routine however, it
      ! is used as a work space to track the net number of entries appearing
      ! for the first time at node i of the elimination tree (this may be
      ! negative).
   integer, intent(out) :: st ! stat parmeter for allocate calls
   integer, dimension(:), optional, intent(in) :: wt ! weights (eg number of
      ! variables in each supervariable)
   
   integer :: col ! column of matrix associated with piv
   integer, dimension(:), allocatable :: first ! first descendants
   integer :: i
   integer(pkg_type) :: ii
   integer :: totalwt
   integer, dimension(:), allocatable :: last_nbr ! previous neighbour
   integer, dimension(:), allocatable :: last_p ! previous p?
   integer :: par ! parent node of piv
   integer :: piv ! current pivot
   integer :: pp ! last pivot where u was encountered
   integer :: lca ! least common ancestor of piv and pp
   integer :: u ! current entry in column col
   integer :: uwt ! weight of u
   integer, dimension(:), allocatable :: vforest ! virtual forest

   !
   ! Determine first descendants, and set cc = 1 for leaves and cc = 0 for
   ! non-leaves.
   !
   allocate(first(n+1), stat=st)
   if(st.ne.0) return
   do i = 1, n+1
      first(i) = i
   end do
   if(present(wt)) then
      totalwt = 0 ! Find sum of weights so we can determine non-physical value
      do i = 1, n
         par = parent(i)
         first(par) = min(first(i), first(par)) ! first descendant
         if(first(i).eq.i) then ! is it a leaf or not?
            cc(i) = wt(invp(i))
         else
            cc(i) = 0
         endif
         totalwt = totalwt + wt(invp(i))
      end do
      cc(n+1) = totalwt + 1 ! Set to non-physical value
   else
      do i = 1, n
         par = parent(i)
         first(par) = min(first(i), first(par)) ! first descendant
         if(first(i).eq.i) then ! is it a leaf or not?
            cc(i) = 1
         else
            cc(i) = 0
         endif
      end do
      cc(n+1) = n + 1 ! Set to non-physical value
   endif

   !
   ! We store the partial elimination trees in a virtual forest. It is
   ! initialised such that each node is in its own tree to begin with.
   !
   allocate(vforest(n+1), stat=st)
   if(st.ne.0) return
   vforest(:) = 0

   !
   ! Initialise previous pivot and neightbour arrays to indicate no previous
   ! pivot or neightbour.
   !
   allocate(last_p(n+1), last_nbr(n+1), stat=st)
   if(st.ne.0) return
   last_p(:) = 0
   last_nbr(:) = 0

   !
   ! Determine cc(i), the number of net new entries to pass up tree from
   ! node i.
   !
   do piv = 1, n
      ! Loop over entries in column below the diagonal
      col = invp(piv)
      do ii = ptr(col), ptr(col+1)-1
         u = perm(row(ii))
         if(u.le.piv) cycle ! not in lower triangular part

         ! Check if entry has been seen by a descendant of this pivot, if
         ! so we skip the tests that would first add one to the current
         ! pivot's weight before then subtracting it again.
         if(first(piv).gt.last_nbr(u)) then
            ! Count new entry in current column
            uwt = 1
            if(present(wt)) uwt = wt(invp(u))
            cc(piv) = cc(piv) + uwt

            ! Determine least common ancestor of piv and the node at which
            ! u was last encountred
            pp = last_p(u)
            if(pp.ne.0) then
               ! u has been seen before, find top of partial elimination
               ! tree for node pp
               lca = FIND(vforest, pp)
               ! prevent double counting of u at node lca
               cc(lca) = cc(lca) - uwt
            endif

            ! Update last as u has now been seen at piv.
            last_p(u) = piv
         endif

         ! Record last neighbour of u so we can determine if it has been
         ! seen in this subtree before
         last_nbr(u) = piv
      end do
      ! Pass uneliminated variables up to parent
      par = parent(piv)
      if(present(wt)) then
         cc(par) = cc(par) + cc(piv) - wt(invp(piv))
      else
         cc(par) = cc(par) + cc(piv) - 1
      endif

      ! place the parent of piv into the same partial elimination tree as piv
      vforest(piv) = par ! operation "UNION" from [1]
   end do
end subroutine mc78_col_counts_integer

! Return top most element of tree containing u.
! Implements path compression to speed up subsequent searches.
integer function FIND(vforest, u)
   integer, dimension(:), intent(inout) :: vforest
   integer, intent(in) :: u

   integer :: current, prev

   prev = -1
   current = u
   do while(vforest(current).ne.0)
      prev = current
      current = vforest(current)
      if(vforest(current).ne.0) vforest(prev) = vforest(current)
   end do

   FIND = current
end function FIND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Supernode amalgamation routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine identifies (relaxed) supernodes from the elimination tree
! and column counts.
!
! A node, u, and its parent, v, are merged if:
! (a) No new fill-in is introduced i.e. cc(v) = cc(u)-1
! (b) The number of columns in both u and v is less than nemin
!
! Note: assembly tree must be POSTORDERED on output
subroutine mc78_supernodes(n, realn, parent, cc, sperm, nnodes, sptr, sparent, &
      scc, invp, control, info, st, wt, block_pivots)
   integer, intent(in) :: n
   integer, intent(in) :: realn
   integer, dimension(n), intent(in) :: parent ! parent(i) is the
      ! parent of supernode i in the elimination/assembly tree. 
   integer, dimension(n), intent(in) :: cc ! cc(i) is the column count
      ! of supernode i, including elements eliminated at supernode i.
   integer, dimension(n), intent(out) :: sperm ! on exit contains a permutation
      ! from pivot order to a new pivot order with contigous supernodes
   integer, intent(out) :: nnodes ! number of supernodes
   integer, dimension(n+1), intent(out) :: sptr
   integer, dimension(n), intent(out) :: sparent
   integer, dimension(n), intent(out) :: scc
   integer, dimension(n), intent(in) :: invp
   type(mc78_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(out) :: st ! stat paremter from allocate calls
   integer, dimension(n), optional, intent(in) :: wt ! weights (number of vars
      ! in each initial node)
   integer, dimension(n), optional, intent(in) :: block_pivots ! If
      ! present, then matches pivot order and specifies block pivots. 
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot

   integer :: i, j, k
   integer :: flag
   integer, dimension(:), allocatable :: height ! used to track height of tree
   logical, dimension(:), allocatable :: mark ! flag array for nodes to finalise
   integer, dimension(:), allocatable :: map ! map vertex idx -> supernode idx
   integer, dimension(:), allocatable :: nelim ! number of eliminated variables
   integer, dimension(:), allocatable :: nvert ! number of elimd supervariables
   integer :: node
   integer, dimension(:), allocatable :: npar ! temporary array of snode pars
   integer :: par ! parent of current node
   integer :: shead ! current head of stack
   integer, dimension(:), allocatable :: stack ! used to navigate tree
   integer :: v
   integer, dimension(:), allocatable :: vhead ! heads of vertex linked lists
   integer, dimension(:), allocatable :: vnext ! next element in linked lists
   integer(long), dimension(:), allocatable :: ezero ! number of explicit zeros
   integer, dimension(:), allocatable :: chead ! chead(i) is first child of i
   integer, dimension(:), allocatable :: cnext ! cnext(i) is next child of i
   integer, dimension(:), allocatable :: child
   integer :: nchild
   integer :: start ! First pivot in block pivot
   integer :: totalwt ! sum of weights

   !
   ! Initialise supernode representation
   !
   allocate(nelim(n+1), nvert(n+1), vhead(n+1), vnext(n+1), stack(n), &
      height(n+1), mark(n), stat=st)
   if(st.ne.0) goto 490
   vnext(:) = -1
   vhead(:) = -1
   height(:) = 1

   ! Initialise number of variables in each node
   if(present(wt)) then
      totalwt = 0
      do i = 1, n
         nelim(i) = wt(invp(i))
         totalwt = totalwt + wt(invp(i))
      end do
   else ! All nodes initially contain a single variable
      nelim(:) = 1
      totalwt = n
   endif
   nvert(:) = 1

   allocate(map(n+1), npar(n+1), ezero(n+1), stat=st)
   if(st.ne.0) goto 490

   ezero(:) = 0 ! Initially no explicit zeros
   ezero(n+1) = huge(ezero) ! ensure root is not merged

   ! Ensure virtual root never gets amlgamated
   nelim(n+1) = totalwt+1 + control%nemin

   !
   ! Build child linked lists for nodes; merge block pivots if needed
   !
   allocate(chead(n+1), cnext(n+1), child(n), stat=st)
   if(st.ne.0) goto 490
   if(present(block_pivots)) then
      chead(:) = -1 ! no parent if necessary
      do i = realn, 1, -1 ! do in reverse order so come off in original order
         if(block_pivots(i).lt.2) cycle
         j = parent(i)
         if(j.ne.n+1) then
            do while(block_pivots(j).lt.2)
               j = parent(j)
            end do
         end if
         cnext(i) = chead(j)
         chead(j) = i
      end do
   else
      chead(:) = -1 ! no parent if necessary
      do i = realn, 1, -1 ! do in reverse order so come off in original order
         j = parent(i)
         cnext(i) = chead(j)
         chead(j) = i
      end do
   endif

   !
   ! Merge supernodes.
   !
   flag = 0
   v = 1
   nnodes = 0
   start=n+2
   do par = 1, n+1
      if(present(block_pivots) .and. par.lt.n+1) then
         if(block_pivots(par).lt.2) then
            if(start.ge.n+1) start = par
            cycle
         endif
         ! Merge pivots start to par, but don't add to vertex list (yet)

         do node = start, par-1
            ! Add together eliminated variables
            nelim(par) = nelim(par) + nelim(node)
            nvert(par) = nvert(par) + nvert(node)

            ! nodes have same height
            height(par) = max(height(par), height(node))
         end do
      endif

      nchild = 0
      node = chead(par)
      do while(node.ne.-1)
         nchild = nchild + 1
         child(nchild) = node
         node = cnext(node)
      end do
      call sort_by_val(nchild, child, cc, st)
      if(st.ne.0) goto 490

      do j = 1, nchild
         node = child(j)
         if(do_merge(node, par, nelim, cc, ezero, control, invp, flag, wt)) then
            ! Merge contents of node into par. Delete node.
            call merge_nodes(node, par, nelim, nvert, vhead, vnext, height, &
               ezero, cc)
            mark(node) = .false.
         else
            mark(node) = .true.
         endif
      end do

      if(present(block_pivots) .and. par.lt.n+1) then
         if(block_pivots(par).ge.2) then
            ! Add vertices start to par-1 into par
            do node = start, par-1
               vnext(node) = vhead(par)
               vhead(par) = node
               mark(node) = .false.
            end do
            start = n+2
         endif
      endif
   end do

   if(flag.ne.0) then
      if(control%unit_error.gt.0) write(control%unit_error, "(a)") &
         "MC78 Internal Error: Unrecognised amalgamation heuristic."
      info = MC78_ERROR_UNKNOWN
      return
   endif

   do node = 1, realn
      if(.not.mark(node)) cycle
      ! Node not merged, now a complete supernode

      ! Record start of supernode
      nnodes = nnodes + 1
      sptr(nnodes) = v
      npar(nnodes) = parent(node)
      if(present(wt)) then
         scc(nnodes) = cc(node) + nelim(node) - wt(invp(node))
      else
         scc(nnodes) = cc(node) + nelim(node) - 1
      endif

      ! Record height in tree of parent vertices
      height(parent(node)) = max(height(parent(node)), height(node) + 1)

      ! Determine last vertex of node so we can number backwards
      v = v + nvert(node)
      k = v

      ! Loop over member vertices of node and number them
      shead = 1
      stack(shead) = node
      do while(shead.gt.0)
         i = stack(shead)
         shead = shead - 1

         ! Order current vertex
         k = k - 1
         sperm(i) = k
         map(i) = nnodes

         ! Stack successor, if any
         if(vnext(i).ne.-1) then
            shead = shead + 1
            stack(shead) = vnext(i)
         endif

         ! Descend into tree rooted at i
         if(vhead(i).ne.-1) then
            shead = shead + 1
            stack(shead) = vhead(i)
         endif
      end do
   end do
   sptr(nnodes+1) = v ! Record end of final supernode
   map(n+1) = nnodes + 1 ! virtual root vertex maps to virtual root sn
   npar(nnodes+1) = n + 1

   ! Handle permutation of empty columns
   do i = realn+1, n
      sperm(i) = i
   end do

   ! Allocate arrays for return and copy data into them correctly
   do node = 1, nnodes
      par = npar(node) ! parent /vertex/ of supernode
      par = map(par)   ! parent /node/   of supernode
      sparent(node) = par ! store parent
   end do

   return

   490 continue
   info = MC78_ERROR_ALLOC
   return
end subroutine mc78_supernodes

!
! Sort n items labelled by idx into decreasing order of val(idx(i))
!
recursive subroutine sort_by_val(n, idx, val, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: idx
   integer, dimension(:), intent(in) :: val
   integer, intent(out) :: st

   integer :: ice_idx, ice_val, ik_idx, ik_val
   integer :: klo,kor,k,kdummy

   st = 0

   if(n.ge.minsz_ms) then
      call sort_by_val_ms(n, idx, val, st)
   else
      klo = 2
      kor = n
      do kdummy = klo,n
         ! items kor, kor+1, .... ,nchild are in order
         ice_idx = idx(kor-1)
         ice_val = val(ice_idx)
         do k = kor,n
            ik_idx = idx(k)
            ik_val = val(ik_idx)
            if (ice_val >= ik_val) exit
            idx(k-1) = ik_idx
         end do
         idx(k-1) = ice_idx
         kor = kor - 1
      end do
   endif
end subroutine sort_by_val

! Sort n items labelled by idx into decreasing order of val(idx(i))
!
! Merge sort version, dramatically improves performance for nodes with large
! numbers of children
! (Passes to simple sort for small numbers of entries)
recursive subroutine sort_by_val_ms(n, idx, val, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: idx
   integer, dimension(:), intent(in) :: val
   integer, intent(out) :: st

   integer :: i, j, jj, jj2, k, kk, kk2
   integer :: mid
   integer, dimension(:), allocatable :: work

   if(n.le.1) return
   if(n.lt.minsz_ms) then
      call sort_by_val(n, idx, val, st)
      return
   endif
   mid = (n-1)/2 + 1

   ! Recurse to order half lists
   call sort_by_val_ms(mid, idx(1:mid), val, st)
   if(st.ne.0) return
   call sort_by_val_ms(n - mid, idx(mid+1:n), val, st)
   if(st.ne.0) return

   ! Merge two half lists
   ! (Take a copy of the first half list so we don't overwrite it)
   allocate(work(mid), stat=st)
   if(st.ne.0) return
   work(:) = idx(1:mid)
   j = 1
   k = mid+1
   jj = work(j)
   jj2 = val(jj)
   kk = idx(k)
   kk2 = val(kk)
   do i = 1, n
      if(jj2.ge.kk2) then
         idx(i) = jj
         j = j + 1
         if(j.gt.mid) exit
         jj = work(j)
         jj2 = val(jj)
      else
         idx(i) = kk
         k = k + 1
         if(k.gt.n) exit
         kk = idx(k)
         kk2 = val(kk)
      endif
   end do
   if(j.le.mid) idx(i+1:n) = work(j:mid)
end subroutine sort_by_val_ms

!
! Return .true. if we should merge node and par, .false. if we should not
!
logical function do_merge(node, par, nelim, cc, ezero, control, invp, info, wt)
   integer, intent(in) :: node ! node to merge and delete
   integer, intent(in) :: par ! parent to merge into
   integer, dimension(:), intent(in) :: nelim
   integer, dimension(:), intent(in) :: cc
   integer(long), dimension(:), intent(in) :: ezero
   type(mc78_control), intent(in) :: control
   integer, dimension(:), intent(in) :: invp
   integer, intent(out) :: info
   integer, dimension(:), optional, intent(in) :: wt

   real(dp) :: z, ne

   info = 0

   if(ezero(par).eq.huge(ezero)) then
      do_merge = .false.
      return
   endif

   select case(control%heuristic)
   case(1)
      !
      ! HSL_MA77 style nemin
      !
      if(present(wt)) then
         do_merge = (cc(par).eq.cc(node)-wt(invp(node)) .and. &
            nelim(par).eq.wt(invp(par))) .or. &
            (nelim(par).lt.control%nemin .and. nelim(node).lt.control%nemin)
      else
         do_merge = (cc(par).eq.cc(node)-1 .and. nelim(par).eq.1) .or. &
            (nelim(par).lt.control%nemin .and. nelim(node).lt.control%nemin)
      endif
   case(2)
      !
      ! CHOLMOD style nrelax/zrelax
      !
      ! FIXME: currently assumes nodes are square, not trapezoidal

      ! calculate number of non-zeros in new node
      z = ezero(par) + ezero(node) + &
         (cc(par)-1+nelim(par) - cc(node)+1) * nelim(par)
      ! find this as a fraction of total non-zeros in new node
      ne = nelim(par) + nelim(node)
      z = z / ( (cc(par)-1+ne)*ne )

      do_merge = (ne .le. control%nrelax(1)) .or. &
         (cc(par).eq.cc(node)-1 .and. nelim(par).eq.1) .or. &
         (ne .le. control%nrelax(2) .and. z .lt. control%zrelax(1)) .or. &
         (ne .le. control%nrelax(3) .and. z .lt. control%zrelax(2)) .or. &
         (z .lt. control%zrelax(3))
   case default
      ! Note: This bit of code should NEVER execute
      do_merge = .false.
      info = MC78_ERROR_UNKNOWN
   end select
end function do_merge

!
! This subroutine merges node with its parent, deleting node in the process.
!
subroutine merge_nodes(node, par, nelim, nvert, vhead, vnext, height, ezero, cc)
   integer, intent(in) :: node ! node to merge and delete
   integer, intent(in) :: par ! parent to merge into
   integer, dimension(:), intent(inout) :: nelim
   integer, dimension(:), intent(inout) :: nvert
   integer, dimension(:), intent(inout) :: vhead
   integer, dimension(:), intent(inout) :: vnext
   integer, dimension(:), intent(inout) :: height
   integer(long), dimension(:), intent(inout) :: ezero
   integer, dimension(:), intent(in) :: cc

   ! Add node to list of children merged into par
   vnext(node) = vhead(par)
   vhead(par) = node

   ! Work out number of explicit zeros in new node
   ! FIXME: probably wrong now with weights and block pivots
   ezero(par) = ezero(par) + ezero(node) + &
      (cc(par)-1+nelim(par) - cc(node) + 1_long) * nelim(par)

   ! Add together eliminated variables
   nelim(par) = nelim(par) + nelim(node)
   nvert(par) = nvert(par) + nvert(node)

   ! nodes have same height
   height(par) = max(height(par), height(node))
end subroutine merge_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Statistics routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine merely calculates interesting statistics
!
subroutine mc78_stats(nnodes, sptr, scc, nfact, nflops)
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: scc
   integer(long), optional, intent(out) :: nfact
   integer(long), optional, intent(out) :: nflops

   integer :: j
   integer :: m ! number of entries in retangular part of ndoe
   integer :: nelim ! width of node
   integer :: node ! current node of assembly tree
   integer(long) :: r_nfact, r_nflops

   if(.not.present(nfact) .and. .not.present(nflops)) return ! nothing to do

   r_nfact = 0
   r_nflops = 0
   do node = 1, nnodes
      nelim = sptr(node+1) - sptr(node)
      m = scc(node) - nelim

      ! number of entries
      r_nfact = r_nfact + (nelim * (nelim+1)) / 2 ! triangular block
      r_nfact = r_nfact + nelim * m ! below triangular block

      ! flops
      do j = 1, nelim
         r_nflops = r_nflops + (m+j)**2
      end do
   end do

   if(present(nfact)) nfact = r_nfact
   if(present(nflops)) nflops = r_nflops

   !print *, "n = ", n
   !print *, "nnodes = ", nnodes
   !print *, "nfact = ", nfact
   !print *, "sum cc=", sum(cc(1:n))
   !print *, "nflops = ", nflops
end subroutine mc78_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Row list routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine determines the row indices for each supernode
!
subroutine mc78_row_lists_nosvar_integer(n, ptr, row, perm, invp, nnodes, &
      sptr, sparent, scc, rptr, rlist, control, info, st)
   integer, intent(in) :: n
   integer(pkg_type), dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n), intent(in) :: invp
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer, dimension(nnodes), intent(in) :: scc
   integer(long), dimension(nnodes+1), intent(out) :: rptr
   integer, dimension(sum(scc(1:nnodes))), intent(out) :: rlist
   type(mc78_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(out) :: st

   integer :: child ! current child of node
   integer :: col ! current column of matrix corresponding to piv
   integer(long) :: i
   integer(long) :: idx ! current insert position into nodes(node)%index
   integer :: j
   integer :: node ! current node of assembly tree
   integer :: piv ! current pivot position
   integer, dimension(:), allocatable :: seen ! tag last time index was seen
   integer, dimension(:), allocatable :: chead ! head of child linked lists
   integer, dimension(:), allocatable :: cnext ! pointer to next child

   ! Allocate and initialise memory
   allocate(seen(n), chead(nnodes+1), cnext(nnodes+1), stat=st)
   if(st.ne.0) then
      info = MC78_ERROR_ALLOC
      return
   endif
   seen(:) = 0
   chead(:) = -1

   ! Build child linked lists (backwards so pop off in good order)
   do node = nnodes, 1, -1
      i = sparent(node)
      cnext(node) = chead(i)
      chead(i) = node
   end do

   ! Loop over nodes from bottom up building row lists.
   rptr(1) = 1
   do node = 1, nnodes

      ! Allocate space for row indices
      rptr(node+1) = rptr(node) + scc(node)
      idx = rptr(node) ! insert position

      ! Add entries eliminated at this node
      do piv = sptr(node), sptr(node+1)-1
         seen(piv) = node
         rlist(idx) = piv
         idx = idx + 1
      end do

      ! Find indices inherited from children
      child = chead(node)
      do while (child.ne.-1)
         do i = rptr(child), rptr(child+1)-1
            j = rlist(i)
            if(j.lt.sptr(node)) cycle ! eliminated
            if(seen(j).eq.node) cycle ! already seen
            seen(j) = node
            rlist(idx) = j
            idx = idx + 1
         end do
         child = cnext(child)
      end do

      ! Find new indices from A
      do piv = sptr(node), sptr(node+1)-1
         col = invp(piv)
         do i = ptr(col), ptr(col+1)-1
            j = perm(row(i))
            if(j.lt.piv) cycle ! in upper triangle
            if(seen(j).eq.node) cycle ! already seen in this snode
            ! Otherwise, this is a new entry
            seen(j) = node
            rlist(idx) = j
            idx = idx + 1
         end do
      end do

      ! Note: following error check won't work with block pivots
      !if(idx .ne. rptr(node+1)) then
      !   ! Note: This bit of code should NEVER execute
      !  if(control%unit_error.gt.0) write(control%unit_error, "(3(a,i8))") &
      !      "MC78 Internal Error: node ", node, ": found ", idx-rptr(node), &
      !      " entries, but expected to find ", rptr(node+1)-rptr(node)
      !   info = MC78_ERROR_UNKNOWN
      !   !print *, rlist(1:idx-1)
      !   return
      !endif
   end do
end subroutine mc78_row_lists_nosvar_integer

subroutine mc78_row_lists_svar_integer(nsvar, svar, n, ptr, row, perm, invp, &
      nnodes, sptr, sparent, scc, rptr, rlist, control, info, st)
   integer, intent(in) :: nsvar
   integer, dimension(nsvar), intent(in) :: svar
   integer, intent(in) :: n
   integer(pkg_type), dimension(nsvar+1), intent(in) :: ptr
   integer, dimension(ptr(nsvar+1)-1), intent(in) :: row
   integer, dimension(nsvar), intent(in) :: perm
   integer, dimension(nsvar), intent(in) :: invp
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer, dimension(nnodes), intent(in) :: scc
   integer(long), dimension(nnodes+1), intent(out) :: rptr
   integer, dimension(sum(scc(1:nnodes))), intent(out) :: rlist
   type(mc78_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(out) :: st

   integer :: child ! current child of node
   integer :: col ! current column of matrix corresponding to piv
   integer(long) :: i
   integer(long) :: idx ! current insert position into nodes(node)%index
   integer :: j
   integer :: k
   integer :: node ! current node of assembly tree
   integer :: piv ! current pivot position
   integer, dimension(:), allocatable :: seen ! tag last time index was seen
   integer, dimension(:), allocatable :: chead ! head of child linked lists
   integer, dimension(:), allocatable :: cnext ! pointer to next child
   integer, dimension(:), allocatable :: svptr ! pointers for row list starts

   ! Allocate and initialise memory
   allocate(seen(n), chead(nnodes+1), cnext(nnodes+1), svptr(nsvar+1), stat=st)
   if(st.ne.0) return
   seen(:) = 0
   chead(:) = -1

   ! Build svptr array
   svptr(1) = 1
   do i = 1, nsvar
      svptr(i+1) = svptr(i) + svar(invp(i))
   end do

   ! Build child linked lists (backwards so pop off in good order)
   do node = nnodes, 1, -1
      i = sparent(node)
      cnext(node) = chead(i)
      chead(i) = node
   end do

   ! Loop over nodes from bottom up building row lists.
   rptr(1) = 1
   do node = 1, nnodes

      ! Allocate space for row indices
      rptr(node+1) = rptr(node) + scc(node)
      idx = rptr(node) ! insert position

      ! Add entries eliminated at this node
      do i = sptr(node), sptr(node+1)-1
         do piv = svptr(i), svptr(i+1)-1
            seen(piv) = nnodes+1
            rlist(idx) = piv
            idx = idx + 1
         end do
      end do

      ! Find indices inherited from children
      child = chead(node)
      do while (child.ne.-1)
         do i = rptr(child), rptr(child+1)-1
            j = rlist(i)
            if(seen(j).ge.node) cycle ! already seen (or eliminated already)
            seen(j) = node
            rlist(idx) = j
            idx = idx + 1
         end do
         child = cnext(child)
      end do

      ! Find new indices from A
      do piv = sptr(node), sptr(node+1)-1
         col = invp(piv)
         do i = ptr(col), ptr(col+1)-1
            j = perm(row(i))
            if(seen(svptr(j)).ge.node) cycle ! already seen (or eliminated)
            ! Otherwise, this is a new entry
            ! Iterate over variables in supervariable
            do k = svptr(j), svptr(j+1)-1
               seen(k) = node
               rlist(idx) = k
               idx = idx + 1
            end do
         end do
      end do

      if(idx .ne. rptr(node+1)) then
         ! Note: This bit of code should NEVER execute
         if(control%unit_error.gt.0) write(control%unit_error, "(3(a,i8))") &
            "MC78 Internal Error: node ", node, ": found ", idx-rptr(node), &
            " entries, but expected to find ", rptr(node+1)-rptr(node)
         info = MC78_ERROR_UNKNOWN
         !print *, rlist(1:idx-1)
         return
      endif
   end do
end subroutine mc78_row_lists_svar_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Optimize cache locality routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! The following subroutine reorders the elimination order within each node in
! such a way that the order of the variables in
! the "primary" child has ordering agreement with the node. Consider
! the following tree:
!          A
!         / \
!        B   C
! 
! and assume B has more partially summed variables than C. Then
! B is the primary child and the ordering of the
! corresponding fully summed variables in the parent A matches the ordering
! of the partially summed variables in B (but not in C). Any additional
! partially summed variables present in C but not in B are then ordered in A
! such that they match C.
!
! This is done by two passes of the tree.
! The first builds a map from variables to the nodes at which
! they are eliminated, and orders the children of each node such
! that the first has the largest number of partially summed variables.
! The second pass uses a depth first search of the now ordered tree. It loops
! over non-fully summed variables and when it first encounters each it will
! place it as the next variable at its elimination node.
!
subroutine mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, sparent, &
      rptr, rlist, st, sort)

   integer, intent(in) :: n ! dimension of system
   integer, intent(in) :: realn ! symbolic dimension of system
   integer, dimension(n), intent(inout) :: perm ! on exit, will have been
      ! reordered for better cache locality
   integer, dimension(n), intent(inout) :: invp ! inverse of perm. on exit
      ! will have been changed to match new perm
   integer, intent(in) :: nnodes ! number of supernodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), intent(inout) :: rlist
   integer, intent(out) :: st ! stat parameter
   logical, optional, intent(in) :: sort

   integer :: i
   integer(long) :: ii ! loop index
   integer :: id ! current insert position into list
   integer :: k
   integer :: j ! temporary variable
   integer, allocatable :: list(:) ! list of nodes in a weighted depth-first
      ! order such that children are visitied in decreasing number of
      ! partially summed variables (ie child with most p.s.v. visited first)
   integer, allocatable :: map(:) ! maps variables to nodes where
      ! they are eliminated
   integer :: node ! current node
   integer, allocatable :: ord(:) ! tracks number of variables ordered at
      ! each node
   integer, allocatable :: perm2(:) ! permutation to apply to perm
   integer :: pnode ! parent node
   integer :: shead ! current top of stack
   integer, allocatable :: stack(:) ! used for depth first walk of tree
   integer :: start ! first entry on stack of child from current node
   integer, allocatable :: chead(:) ! heads of child linked lists
   integer, allocatable :: cnext(:) ! tails of child linked lists

   ! Allocate arrays for depth first search of tree
   allocate (map(n), list(nnodes+1), stack(nnodes), chead(nnodes+1), &
      cnext(nnodes), stat=st)
   if (st /= 0) return

   !
   ! Build elimination map
   !
   do node = 1, nnodes
      do ii = sptr(node), sptr(node+1)-1
         map(ii) = node
      end do
   end do

   !
   ! Build child linked lists
   !
   chead(:) = -1 ! no child if necessary
   do i = nnodes, 1, -1 ! do in reverse order so they come off in original order
      j = sparent(i)
      cnext(i) = chead(j)
      chead(j) = i
   end do

   !
   ! Perform depth first search of tree, such that children of a node are
   ! visited in order of the number of partially summer variables, largest
   ! first.
   !
   shead = 1
   stack(shead) = nnodes + 1
   id = nnodes+1
   do while(shead.ne.0)
      ! Get node from top of stack
      node = stack(shead)
      shead = shead - 1

      ! Number it
      list(id) = node
      id = id - 1

      ! Place all its children on the stack
      start = shead + 1
      i = chead(node)
      do while(i.ne.-1)
         shead = shead + 1
         stack(shead) = i
         i = cnext(i)
      end do
      ! Order children just placed on stack such that child with least partially
      ! summed variables is at the top
      call order_children(shead-start+1, stack(start:shead), nnodes, sptr, &
         rptr, st)
      if(st.ne.0) return
   end do

   !
   ! Next loop over children reordering partially summed variables.
   !
   allocate(ord(nnodes),perm2(n),stat=st)
   if (st.ne.0) return

   do node = 1, nnodes
      ord(node) = sptr(node)
   end do

   do k = 1, nnodes
      node = list(k)

      ! Order variables first encountered at this node
      do ii = rptr(node), rptr(node+1)-1
         j = rlist(ii)
         pnode = map(j)
         if(pnode .ne. -1) then ! check if we have ordered j already
            ! order at parent
            perm2(j) = ord(pnode)
            ord(pnode) = ord(pnode) + 1
            map(j) = -1 ! mark as ordered
         endif
         rlist(ii) = perm2(j)
      end do
   end do

   do i = realn+1, n
      perm2(i) = i
   end do

   !
   ! Apply permutation to perm and invp
   !
   ! Use perm as a temporary variable to permute invp.
   perm(1:n) = invp(1:n)
   do i = 1, n
      j = perm2(i)
      invp(j) = perm(i)
   end do

   ! Recover invp as inverse of perm
   do i = 1, n
      perm(invp(i)) = i
   end do

   if(present(sort)) then
      if(sort) then
         call dbl_tr_sort(n, nnodes, rptr, rlist, st)
         if(st.ne.0) return
      endif
   endif
end subroutine mc78_optimize_locality

! Orders nodes stored in child(1:n) such that the number of partially summed
! variables at each node is decreasing (ie one with least is in posn n)
!
! Simple sort version, good for nodes with small numbers of children
! (Passes to mergesort for large numbers of entries)
recursive subroutine order_children(n, child, nnodes, sptr, rptr, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: child
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, intent(out) :: st

   integer :: nelim, m
   integer :: k, kdummy, klo, kor
   integer :: ice_idx, ice_psum, ik_idx, ik_psum

   st = 0

   if(n.ge.minsz_ms) then
      call order_children_ms(n, child, nnodes, sptr, rptr, st)
   else
      klo = 2
      kor = n
      do kdummy = klo, n
         ! items kor, kor+1, .... ,n are in order
         ice_idx = child(kor-1)
         nelim = sptr(ice_idx+1) - sptr(ice_idx)
         m = int(rptr(ice_idx+1) - rptr(ice_idx))
         ice_psum = m - nelim
         do k = kor, n
            ik_idx = child(k)
            nelim = sptr(ik_idx+1) - sptr(ik_idx)
            m = int(rptr(ik_idx+1) - rptr(ik_idx))
            ik_psum = m - nelim
            if (ice_psum .ge. ik_psum) exit
            child(k-1) = ik_idx
         end do
         child(k-1) = ice_idx
         kor = kor - 1
      end do
   endif
end subroutine order_children

! Orders nodes stored in child(1:n) such that the number of partially summed
! variables at each node is decreasing (ie one with least is in posn n)
!
! Merge sort version, dramatically improves performance for nodes with large
! numbers of children
! (Passes to simple sort for small numbers of entries)
recursive subroutine order_children_ms(n, child, nnodes, sptr, rptr, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: child
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, intent(out) :: st

   integer :: i, j, jj, jj2, k, kk, kk2, m, nelim
   integer :: mid
   integer, dimension(:), allocatable :: work

   if(n.le.1) return
   if(n.lt.minsz_ms) then
      call order_children(n, child, nnodes, sptr, rptr, st)
      return
   endif
   mid = (n-1)/2 + 1

   ! Recurse to order half lists
   call order_children_ms(mid, child(1:mid), nnodes, sptr, rptr, st)
   if(st.ne.0) return
   call order_children_ms(n - mid, child(mid+1:n), nnodes, sptr, rptr, st)
   if(st.ne.0) return

   ! Merge two half lists
   ! (Take a copy of the first half list so we don't overwrite it)
   allocate(work(mid), stat=st)
   if(st.ne.0) return
   work(:) = child(1:mid)
   j = 1
   k = mid+1
   jj = work(j)
   nelim = sptr(jj+1) - sptr(jj)
   m = int(rptr(jj+1) - rptr(jj))
   jj2 = m - nelim
   kk = child(k)
   nelim = sptr(kk+1) - sptr(kk)
   m = int(rptr(kk+1) - rptr(kk))
   kk2 = m - nelim
   do i = 1, n
      if(jj2.ge.kk2) then
         child(i) = jj
         j = j + 1
         if(j.gt.mid) exit
         jj = work(j)
         nelim = sptr(jj+1) - sptr(jj)
         m = int(rptr(jj+1) - rptr(jj))
         jj2 = m - nelim
      else
         child(i) = kk
         k = k + 1
         if(k.gt.n) exit
         kk = child(k)
         nelim = sptr(kk+1) - sptr(kk)
         m = int(rptr(kk+1) - rptr(kk))
         kk2 = m - nelim
      endif
   end do
   if(j.le.mid) child(i+1:n) = work(j:mid)
end subroutine order_children_ms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assorted auxilary routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Converts piv_size to block_pivots:
! piv_size(i) is size of block pivot containing column i of A
! block_pivots(j) is a flag for pivot j of L and has one of the following values
!  0 - pivot i is in the middle of a block pivot
!  1 - pivot i is the first pivot of a block pivot
!  2 - pivot i is the last pivot of a block pivot
!  3 - pivot i is a 1x1 pivot
!
subroutine convert_to_blk_piv(n, invp, block)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: invp
   integer, dimension(n), intent(inout) :: block

   integer, dimension(:), allocatable :: blk2
   integer :: i, cnt

   allocate(blk2(n))

   ! Take a copy of block so we can change it
   blk2(1:n) = block(1:n)

   ! Iterate over pivots in elimination order, recording starts and ends of blks
   cnt = blk2(invp(1))-1 ! Initialise for first pivot
   block(1) = 1 ! First pivot is start of a block
   do i = 2, n
      block(i) = 0
      if(cnt.eq.0) then
         ! this is first pivot of a block, previous is last pivot of a block
         cnt = blk2(invp(i))
         block(i-1) = block(i-1) + 2
         block(i) = block(i) + 1
      endif
      cnt = cnt - 1
   end do
   block(n) = block(n) + 2 ! end of matrix must end a block pivot

end subroutine convert_to_blk_piv

!
! Converts block_pivots back to piv_size:
! block_pivots(j) is a flag for pivot j of L and has one of the following values
!  0 - pivot i is in the middle of a block pivot
!  1 - pivot i is the first pivot of a block pivot
!  2 - pivot i is the last pivot of a block pivot
!  3 - pivot i is a 1x1 pivot
! piv_size(i) is size of block pivot containing column i of A
!
subroutine convert_from_blk_piv(n, perm, block)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n), intent(inout) :: block

   integer, dimension(:), allocatable :: blk2
   integer :: i, sa, cnt

   allocate(blk2(n))

   ! convert first/last notation to block size notation
   cnt = -1; sa = -1 ! these values should never actually be used
   do i = 1, n
      select case(block(i))
      case (0) ! middle pivot of a block
         cnt  = cnt + 1
      case (1) ! first pivot of a block
         sa = i
         cnt = 1
      case (2) ! end pivot of a block
         cnt = cnt + 1
         block(sa:i) = cnt
      case (3) ! only pivot of a block
         block(i) = 1
      end select
   end do

   ! Permute back to original matrix order
   blk2(1:n) = block(1:n)
   do i = 1, n
      block(i) = blk2(perm(i))
   end do
end subroutine convert_from_blk_piv

!
! This subroutine copies a matrix pattern and adds subdiagonal entries as
! needed to force block pivots to have a parent-child relation in the
! elimination tree.
!
subroutine mc78_block_prep(n, ptr, row, bptr, brow, perm, invp, block_pivots, &
      st)
   integer, intent(in) :: n ! Dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer(pkg_type), dimension(n+1), intent(out) :: bptr ! Column pointers
   integer, dimension(:), intent(out) :: brow ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(n), intent(inout) :: invp ! inverse permutation of perm
   integer, dimension(n), intent(inout) :: block_pivots ! Matches pivot order
      ! and specifies block pivots.
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot
   integer, intent(out) :: st

   integer :: col
   integer :: i, j, k
   integer(pkg_type) :: ii
   integer :: idx
   integer :: piv
   integer, dimension(:), allocatable :: seen

   allocate(seen(n), stat=st)
   if(st.ne.0) return

   ! First pass through ptr and ensure that all block pivots contain no
   ! empty columns
   perm(:) = invp(:)
   piv = 1
   j = 1
   ! Handle variables that are actually used
   do while(piv.le.n)
      do i = piv, n
         if(block_pivots(i).ge.2) exit ! end of block pivot
      end do

      k = 0
      do piv = piv, i
         if(ptr(perm(piv)).eq.ptr(perm(piv)+1)) cycle
         invp(j) = perm(piv)
         j = j + 1
         if(k.eq.0) then
            ! This is the new start of the block pivot
            select case(block_pivots(piv))
            case(0) ! was in the middle. now a start
               block_pivots(piv) = 1
            case(2) ! was the end. now a 1x1
               block_pivots(piv) = 3
            end select
         endif
         k = piv
      end do
      if(k.ne.0) then
         ! The was at least one used variable in the block pivot
         select case(block_pivots(k))
         case(0) ! was the middle. now an end
            block_pivots(k) = 2
         case(1) ! was the start. now a 1x1
            block_pivots(k) = 3
         end select
      endif
      piv = i + 1
   end do
   ! Handle unused variables
   do piv = 1, n
      i = perm(piv)
      if(ptr(i).eq.ptr(i+1)) then
         invp(j) = i
         j = j + 1
         block_pivots(piv) = 3 ! Force to 1x1
      endif
   end do
   ! Map block_pivots in original variable order into sv_map
   do i = 1, n
      seen(perm(i)) = block_pivots(i)
   end do
   ! Map sv_map in new pivot order back into block_pivots
   do i = 1, n
      block_pivots(i) = seen(invp(i))
   end do
   ! Reestablish perm
   do i = 1, n
      perm(invp(i)) = i
   end do

   ! Now iterate over cleaned up block pivot sequence
   seen(:) = 0
   piv = 1
   idx = 1
   do col = 1, n
      piv = perm(col)
      bptr(col) = idx
      if(block_pivots(piv).eq.3) then
         ! 1x1 pivot, just copy the column
         idx = idx + int(ptr(col+1) - ptr(col))
         brow(bptr(col):idx-1) = row(ptr(col):ptr(col+1)-1)
      else
         ! copy the column, but add an entry on subdiagonal(s)
         do ii = ptr(col), ptr(col+1)-1
            j = row(ii)
            seen(j) = col
            brow(idx) = j
            idx = idx + 1
         end do
         if(block_pivots(piv).ne.1) then
            ! Not the first column, add an entry above the diagonal
            j = invp(piv-1)
            if(seen(j).lt.col) then
               brow(idx) = j
               idx = idx + 1
            endif
         endif
         if(block_pivots(piv).ne.2) then
            ! Not the last column, add an entry below the diagonal
            j = invp(piv+1)
            if(seen(j).lt.col) then
               brow(idx) = j
               idx = idx + 1
            endif
         endif
      endif
   end do
   bptr(n+1) = idx
end subroutine mc78_block_prep

!
! This subroutine will take information concerning a compressed matrix and a
! supervariable map, and will decompress the information so it relates to
! the original matrix
!
subroutine svar_unmap(n, nsvar, svar, perm, invp, nnodes, sinvp, &
      snptr, st)
   integer, intent(in) :: n
   integer, intent(in) :: nsvar
   integer, dimension(nsvar), intent(in) :: svar
   integer, dimension(n), intent(out) :: perm
   integer, dimension(n), intent(inout) :: invp
   integer, intent(in) :: nnodes
   integer, dimension(nsvar), intent(in) :: sinvp
   integer, dimension(nnodes+1), intent(inout) :: snptr
   integer, intent(out) :: st

   integer, dimension(:), allocatable :: svptr
   integer :: i, j, k
   integer :: j1, j2
   integer :: idx

   ! Set up svptr
   allocate(svptr(nsvar+1), stat=st)
   if(st.ne.0) return
   svptr(1) = 1
   do i = 1, nsvar
      svptr(i+1) = svptr(i) + svar(i)
   end do

   ! Take a copy of invp in perm to ease remapping
   perm(:) = invp(:)

   ! Remap invp
   idx = 1
   do i = 1, nsvar
      j = sinvp(i)
      do k = svptr(j), svptr(j+1)-1
         invp(idx) = perm(k)
         idx = idx + 1
      end do
   end do

   ! Expand supernode pointer
   j1 = snptr(1)
   do i = 1, nnodes
      j2 = snptr(i+1)
      snptr(i+1) = snptr(i)
      do j = j1, j2-1
         snptr(i+1) = snptr(i+1) + svar(sinvp(j))
      end do
      j1 = j2
   end do

   ! Finally, recover perm as inverse of invp
   do i = 1, n
      perm(invp(i)) = i
   end do
end subroutine svar_unmap


!
! This subroutine performs a double transpose sort on the row indices of sn
!
subroutine dbl_tr_sort(n, nnodes, rptr, rlist, st)
   integer, intent(in) :: n
   integer, intent(in) :: nnodes
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), intent(inout) :: rlist
   integer, intent(out) :: st

   integer :: node
   integer :: i, j
   integer(long) :: ii, jj
   integer(long), dimension(:), allocatable :: ptr
   integer(long), dimension(:), allocatable :: nptr
   integer, dimension(:), allocatable :: col

   allocate(ptr(n+2), stat=st)
   if(st.ne.0) return
   ptr(:) = 0

   ! Count number of entries in each row. ptr(i+2) = #entries in row i
   do node = 1, nnodes
      do ii = rptr(node), rptr(node+1)-1
         j = rlist(ii) ! row entry
         ptr(j+2) = ptr(j+2) + 1
      end do
   end do

   ! Determine row starts. ptr(i+1) = start of row i
   ptr(1:2) = 1
   do i = 1, n
      ptr(i+2) = ptr(i+1) + ptr(i+2)
   end do

   jj = ptr(n+2)-1 ! total number of entries
   allocate(col(jj), stat=st)
   if(st.ne.0) return

   ! Now fill in col array
   do node = 1, nnodes
      do ii = rptr(node), rptr(node+1)-1
         j = rlist(ii) ! row entry
         col( ptr(j+1) ) = node
         ptr(j+1) = ptr(j+1) + 1
      end do
   end do

   ! Finally transpose back into nodes
   allocate(nptr(nnodes))
   nptr(:) = rptr(1:nnodes)
   do i = 1, n
      do jj = ptr(i), ptr(i+1)-1
         node = col(jj)
         rlist(nptr(node)) = i
         nptr(node) = nptr(node) + 1
      end do
   end do
end subroutine dbl_tr_sort

!
! This subroutine applies the permutation perm to order, invp and cc
!
subroutine apply_perm(n, perm, order, invp, cc, block_pivots)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n), intent(inout) :: order
   integer, dimension(n), intent(inout) :: invp
   integer, dimension(n), intent(inout) :: cc
   integer, dimension(n), optional, intent(inout) :: block_pivots

   integer :: i
   integer :: j

   ! Use order as a temporary variable to permute cc. Don't care about cc(n+1)
   order(1:n) = cc(1:n)
   do i = 1, n
      j = perm(i)
      cc(j) = order(i)
   end do

   ! Use order as a temporary variable to permute invp.
   order(1:n) = invp(1:n)
   do i = 1, n
      j = perm(i)
      invp(j) = order(i)
   end do

   ! Use order as a temporary variable to permute block_pivots if present
   if(present(block_pivots)) then
      order(1:n) = block_pivots(1:n)
      do i = 1, n
         j = perm(i)
         block_pivots(j) = order(i)
      end do
   endif

   ! Recover order as inverse of invp
   do i = 1, n
      order(invp(i)) = i
   end do
end subroutine apply_perm

end module hsl_mc78_integer
! COPYRIGHT (c) 2006 Council for the Central Laboratory
!                    of the Research Councils
! This package may be copied and used in any application, provided no
! changes are made to these or any other lines.
! Original date 21 February 2006. Version 1.0.0.
! 6 March 2007 Version 1.1.0. Argument stat made non-optional

MODULE HSL_ZD11_double

!  ==========================
!  Sparse matrix derived type
!  ==========================

  TYPE, PUBLIC :: ZD11_type
    INTEGER :: m, n, ne
    CHARACTER, ALLOCATABLE, DIMENSION(:) :: id, type
    INTEGER, ALLOCATABLE, DIMENSION(:) :: row, col, ptr
    REAL ( KIND( 1.0D+0 ) ), ALLOCATABLE, DIMENSION(:) :: val
  END TYPE

CONTAINS

   SUBROUTINE ZD11_put(array,string,stat)
     CHARACTER, allocatable :: array(:)
     CHARACTER(*), intent(in) ::  string
     INTEGER, intent(OUT) ::  stat

     INTEGER :: i,l

     l = len_trim(string)
     if (allocated(array)) then
        deallocate(array,stat=stat)
        if (stat/=0) return
     end if
     allocate(array(l),stat=stat)
     if (stat/=0) return
     do i = 1, l
       array(i) = string(i:i)
     end do

   END SUBROUTINE ZD11_put

   FUNCTION ZD11_get(array)
     CHARACTER, intent(in):: array(:)
     CHARACTER(size(array)) ::  ZD11_get
! Give the value of array to string.

     integer :: i
     do i = 1, size(array)
        ZD11_get(i:i) = array(i)
     end do

   END FUNCTION ZD11_get

END MODULE HSL_ZD11_double


! COPYRIGHT (c) 2001 Council for the Central Laboratory
!                    of the Research Councils
!
! Version 2.2.0
! See ChangeLog for version history.
!
 module HSL_MC65_double

  use hsl_zd11_double
  implicit none
  private

  INTEGER, PARAMETER :: myreal = KIND( 1.0D0 )
  INTEGER, PARAMETER :: myint = KIND( 1)

  integer (kind = myint), public, parameter :: &
       MC65_ERR_MEMORY_ALLOC = -1, &  ! memory alloc failure
       MC65_ERR_MEMORY_DEALLOC = -2, & ! memory deallocate failure
       ! dimension mismatch in matrix_sum
       MC65_ERR_SUM_DIM_MISMATCH = -3,&
       ! dimension mismatch in matrix_multiply
       MC65_ERR_MATMUL_DIM_MISMATCH = -4, &
       ! dimension mismatch in matrix_multiply_graph
       MC65_ERR_MATMULG_DIM_MISMATCH = -5, &
       ! y = Ax with A pattern only and x real.
       MC65_ERR_MATVEC_NOVALUE = -6, &
       ! no vacant unit has been found MC65_MATRIX_WRITE
       MC65_ERR_NO_VACANT_UNIT = -7,&
       ! if in MC65_MATRIX_READ, the file <file_name> does not exist.
       MC65_ERR_READ_FILE_MISS = -8,&
       ! if in MC65_MATRIX_READ, opening of
       ! the file <file_name> returns iostat /= 0
       MC65_ERR_READ_OPEN = -9, &
       ! in MC65_matrix_read, string length of matrix type
       ! exceeds maximum allowed (2000)
       MC65_ERR_READ_MAXLEN = -10, &
       ! in MC65_matrix_read, the particular FORM
       ! is not supported
       MC65_ERR_READ_WRONGFORM = -11, &
       ! in MC65_MATRIX_CONDENSE, try to condense
       ! a matrix that is not of type "pattern"
       MC65_ERR_CONDENSE_NOPAT = -12,&
       ! if PTR(1:M+1) is not monotonically increasing
       ! MATRIX_CONSTRUCT
       MC65_ERR_RANGE_PTR = -13

  integer (kind = myint), public, parameter :: &
       !  COL is not within the range of [1,n] MATRIX_CONSTRUCT,
       ! such entries are excluded in the construct of the matrix
       MC65_WARN_RANGE_COL = 1, &
       ! IRN is not within the  range of {\tt [1,M]}. In
       ! matrix_construct using coordinate formatted input data.
       ! related entries excluded
       MC65_WARN_RANGE_IRN = 2, &
       ! JCN is not within the  range of {\tt [1,N]}. In
       ! matrix_construct using coordinate formatted input data.
       ! related entries excluded
       MC65_WARN_RANGE_JCN = 3,&
       ! IN MC65_MATRIX_DIAGONAL_FIRST, some diagonal elements
       ! are missing!
       MC65_WARN_MV_DIAG_MISSING = 4,&
       ! if duplicate entries were
       !  found when cleaning the matrix
       MC65_WARN_DUP_ENTRY = 5,&
       ! if both out-of-range and duplicate entries are found. In
       ! matrix construct
       MC65_WARN_ENTRIES = 6,&
       ! if both out-of-range IRN and JCN entries found, such entries
       ! are excluded in the construct of the matrix
       MC65_WARN_RANGE_BOTH = 7
  public &
       MC65_print_message, &
       MC65_matrix_construct, &
       MC65_matrix_destruct, &
       MC65_matrix_reallocate, &
       MC65_matrix_transpose, &
       MC65_matrix_copy, &
       MC65_matrix_clean, &
       MC65_matrix_sort, &
       MC65_matrix_sum, &
       MC65_matrix_symmetrize, &
       MC65_matrix_getrow, &
       MC65_matrix_getrowval, &
       MC65_matrix_is_symmetric, &
       MC65_matrix_is_pattern,&
       MC65_matrix_is_different, &
       MC65_matrix_multiply, &
       MC65_matrix_multiply_graph, &
       MC65_matrix_multiply_vector, &
       MC65_matrix_to_coo, &
       MC65_matrix_remove_diagonal,&
       MC65_matrix_diagonal_first, &
       MC65_matrix_write, &
       MC65_matrix_read,&
       MC65_matrix_condense

!       MC65_matrix_fill, &
!       MC65_matrix_component,&
!       MC65_matrix_crop_unsym, &
!       MC65_matrix_crop



  interface MC65_print_message
     module procedure csr_print_message
  end interface

  interface MC65_matrix_construct
     module procedure csr_matrix_construct
     module procedure csr_to_csr_matrix
     module procedure coo_to_csr_format
  end interface
  interface MC65_matrix_is_symmetric
     module procedure csr_matrix_is_symmetric
  end interface
  interface MC65_matrix_is_pattern
     module procedure csr_matrix_is_pattern
  end interface
  interface MC65_matrix_sort
     module procedure csr_matrix_sort
  end interface
  interface MC65_matrix_clean
     module procedure csr_matrix_clean
  end interface
  interface MC65_matrix_multiply
     module procedure csr_matrix_multiply
  end interface
  interface MC65_matrix_multiply_graph
     module procedure csr_matrix_multiply_graph
  end interface
  interface MC65_matrix_to_coo
     module procedure csr_matrix_to_coo
  end interface
  interface MC65_matrix_write
     module procedure csr_matrix_write
  end interface
  interface MC65_matrix_read
     module procedure csr_matrix_read
  end interface
  interface MC65_matrix_reallocate
     module procedure csr_matrix_reallocate
  end interface

  interface MC65_matrix_destruct
     module procedure csr_matrix_destruct
  end interface
  interface MC65_matrix_transpose
     module procedure csr_matrix_transpose
  end interface
  interface MC65_matrix_symmetrize
     module procedure csr_matrix_symmetrize
  end interface
  interface MC65_matrix_getrow
     module procedure csr_matrix_getrow
  end interface
  interface MC65_matrix_getrowval
     module procedure csr_matrix_getrowval
  end interface
  interface MC65_matrix_sum
     module procedure csr_matrix_sum
  end interface
  interface MC65_matrix_copy
     module procedure csr_matrix_copy
  end interface
  interface MC65_matrix_multiply_vector
     module procedure csr_matrix_multiply_rvector
     module procedure csr_matrix_multiply_ivector
  end interface
  interface MC65_matrix_condense
     module procedure csr_matrix_condense
  end interface
  interface MC65_matrix_diagonal_first
     module procedure csr_matrix_diagonal_first
  end interface
  interface MC65_matrix_remove_diagonal
     module procedure csr_matrix_remove_diagonal
  end interface
!!$  interface MC65_matrix_fill
!!$     module procedure csr_matrix_fill
!!$  end interface
!!$  interface MC65_matrix_component
!!$     module procedure csr_matrix_component
!!$  end interface
!!$  interface MC65_matrix_crop
!!$     module procedure csr_matrix_crop
!!$  end interface
!!$  interface MC65_matrix_crop_unsym
!!$     module procedure csr_matrix_crop_unsym
!!$  end interface
  interface MC65_matrix_is_different
     module procedure csr_matrix_diff
  end interface
  interface expand
     module procedure expand1
     module procedure iexpand1
  end interface


contains

  subroutine csr_print_message(info,stream,context)
    ! printing error message
    ! info: is an integer scaler of INTENT (IN).
    ! It is the information flag
    !       whose corresponding error message is to be printed.
    ! stream: is an OPTIONAL integer scaler of INTENT (IN).
    !         It is the unit number
    !         the user wish to print the error message to.
    !         If this number
    !         is negative, printing is supressed. If not supplied,
    !         unit 6 is the default unit to print the error message.
    ! context: is an OPTIONAL assumed size CHARACTER array
    !          of INTENT (IN).
    !          It describes the context under which the error occurred.
    integer (kind = myint), intent (in) :: info
    integer (kind = myint), intent (in), optional :: stream
    character (len = *), optional, intent (in) :: context
    integer (kind = myint) :: unit,length



    if (present(stream)) then
       unit = stream
    else
       unit = 6
    end if
    if (unit <= 0) return

    if (info > 0) then
       write(unit,advance = "yes", fmt = "(' WARNING: ')")
    else if (info < 0) then
       write(unit,advance = "yes", fmt = "(' ERROR: ')")
    end if

    if (present(context)) then
       length = len_trim(context)
       write(unit,advance = "no", fmt = "(a,' : ')") context(1:length)
    end if

    select case (info)
    case (0)
       write(unit,"(a)") "successful completion"
    case (MC65_ERR_MEMORY_ALLOC)
       write(unit,"(A)") "memory allocation failure failure"
    case(MC65_ERR_MEMORY_DEALLOC)
       write(unit,"(A)") "memory deallocate failure"
    case(MC65_ERR_SUM_DIM_MISMATCH)
       write(unit,"(A)") &
            "dimension mismatch of matrices in MC65_matrix_sum"
    case(MC65_ERR_MATMUL_DIM_MISMATCH)
       write(unit,"(A)") "dimension mismatch in matrix_multiply"
    case(MC65_ERR_MATMULG_DIM_MISMATCH)
       write(unit,"(A)") "dimension mismatch in matrix_multiply_graph"
    case(MC65_ERR_MATVEC_NOVALUE)
       write(unit,"(A)") "try to compute y = Ax or y = A^Tx "
       write(unit,"(A)") "with A pattern only and x real"
    case(MC65_ERR_RANGE_PTR)
       write(unit,"(A)") "PTR(1:M+1) is not monotonically &
            &increasing  in MATRIX_CONSTRUCT"
    case(MC65_ERR_NO_VACANT_UNIT)
       write(unit,"(A)") "no vacant I/O unit has been &
            &found MC65_MATRIX_WRITE"
    case(MC65_ERR_READ_FILE_MISS)
       write(unit,"(A)") "in MC65_MATRIX_READ, &
            &the file to be read does not exist"
    case(MC65_ERR_READ_OPEN)
       write(unit,"(A)") "error opening file in MC65_MATRIX_READ"
    case(MC65_ERR_READ_MAXLEN)
       write(unit,"(A)") "in MC65_matrix_read, &
            &string length of matrix type"
       write(unit,"(A)") "exceeds maximum allowed (2000)"
    case(MC65_ERR_READ_WRONGFORM)
       write(unit,"(A)") "in MC65_matrix_read, &
            &the supplied input format is not supported"
    case(MC65_ERR_CONDENSE_NOPAT)
       write(unit,"(A)") "in MC65_MATRIX_CONDENSE, try to condense"
       write(unit,"(A)") "a matrix that is not of type 'pattern'"
       ! warnings ===========
    case(MC65_WARN_MV_DIAG_MISSING)
       write(unit,"(A)") "IN MC65_MATRIX_DIAGONAL_FIRST, &
            &some diagonal elements"
       write(unit,"(A)") "are missing"
    case(MC65_WARN_RANGE_COL)
       write(unit,"(A)") "Some column indices in COL array are not "
       write(unit,"(A)") "within the range of [1,N] &
            &in MATRIX_CONSTRUCT"
       write(unit,"(A)") "such entries are excluded &
            &in the constructed matrix"
    case(MC65_WARN_RANGE_IRN)
       write(unit,"(A)") "IRN is not within the &
            &range of [1,M] in MATRIX_CONSTRUCT"
       write(unit,"(A)") "using coordinate formatted input data. &
            &Such entries are excluded"
    case(MC65_WARN_RANGE_JCN)
       write(unit,"(A)") "JCN is not within the &
            &range of [1,N] in matrix_construct"
       write(unit,"(A)") " using coordinate formatted &
            &input data. Such entries are excluded"
    case (MC65_WARN_DUP_ENTRY)
       write(unit,"(A)") "duplicate entries were found&
            &and had been merged (summed)"
    case (MC65_WARN_ENTRIES)
       write(unit,"(A)") "both duplicate entries and out-of-range entries &
            &were found &
            &and had been treated"
    case (MC65_WARN_RANGE_BOTH)
       write(unit,"(A)") "both out-of-range IRN and JCN entries found. &
            &Such entries are excluded in the constructed matrix"
    case default
       write(*,"(a,i10,a)") "you have supplied an info flag of ",&
            info," this is not a recognized info flag"
    end select
  end subroutine csr_print_message


  subroutine csr_matrix_construct(matrix,m,nz,info,n,type,stat)
    ! subroutine csr_matrix_construct(matrix,m,nz,info[,n,type])
    ! 1) initialise the matrix row/column dimensions and
    ! matrix type, allocate row pointers
    ! by default n = m and type = "general"
    ! 2) allocate a space of nz for the column indices matrix%col
    ! and when matrix%type /= "pattern",
    ! also allocate space for entry values matrix%val.

    ! matrix: of the derived type ZD11_type, INTENT (inOUT),
    !         the sparse matrix in compressed sparse row format
    type (zd11_type), intent (inout) :: matrix

    ! m: is an integer scaler of INTENT (IN), hold the row dimension
    integer (kind = myint), intent (in) :: m

    ! nz: integer scaler of INTENT (IN).
    !     It contains the storage needs to be allocated for the
    !     entries of matrix.
    !     including the storage needs to be allocated for
    !     the array that holds
    !     the column indices; and when type /= "pattern",
    !     also the space allocated to hold the values of the matrix.
    !     $nz$ must be greater than or equal to the number of
    !     entries in the matrix
    !     if $nz < 0$, a storage  of 0 entry is allocated.
    integer (kind = myint), intent (in) :: nz

    ! info: is an integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory allocation failed
    integer (kind = myint), intent (out) :: info

    ! n: is an optional integer scaler of INTENT (IN). If supplied
    !    it contains the column dimension of the matrix
    integer (kind = myint), optional, intent (in) :: n

    ! type: is an optional character array of unspecified length.
    !       INTENT (IN).
    !       If supplied gives the type of the matrix.
    !       If not supplied,
    !       the type of the matrix will be set to "general"
    character (len=*), optional, intent (in) :: type

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! ==================== locals =========================
    ! ierr: error tag for deallocation
    integer (kind = myint) :: ierr
    ! nz1: space allocated
    integer (kind = myint) :: nz1

    if (present(stat)) stat = 0
    info = 0
    ierr = 0

    matrix%m = m
    if (present(n)) then
       matrix%n = n
    else
       matrix%n = m
    end if

    ! CHECK THAT storage is not already existing and if so
    ! see if it is already large enough
    if (allocated(matrix%ptr)) then
       if (size(matrix%ptr) < m+1) then
          deallocate(matrix%ptr,stat = ierr)
          if (present(stat)) stat = ierr
          if (ierr /= 0) then
             info = MC65_ERR_MEMORY_DEALLOC
             return
          end if
          allocate(matrix%ptr(m+1),stat = ierr)
          if (present(stat)) stat = ierr
       end if
    else
       allocate(matrix%ptr(m+1),stat = ierr)
       if (present(stat)) stat = ierr
    end if
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

    ! get rid of the storage of type anyway without
    ! checking to see if
    ! it has large enough storage.
    if (allocated(matrix%type)) then
       deallocate(matrix%type, stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          info = MC65_ERR_MEMORY_DEALLOC
          return
       end if
    end if

    ! by default the matrix has entry values.
    if (present(type)) then
       call zd11_put(matrix%type,type,stat = ierr)
    else
       call zd11_put(matrix%type,"general",stat = ierr)
    end if
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

    ! allocate space for the column indices and possibly the values
    nz1 = nz
    if (nz1 < 0) then
       nz1 = 0
    end if

    ! check to see if col is already allocated and
    ! if so if the storage
    ! is large enough
    if (allocated(matrix%col)) then
       if (size(matrix%col) < nz1) then
          deallocate(matrix%col,stat = ierr)
          if (present(stat)) stat = ierr
          if (ierr /= 0) then
             info = MC65_ERR_MEMORY_DEALLOC
             return
          end if
          allocate(matrix%col(nz1),stat = ierr)
          if (present(stat)) stat = ierr
       end if
    else
       allocate(matrix%col(nz1),stat = ierr)
       if (present(stat)) stat = ierr
    end if
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

    if (ZD11_get(matrix%type) == "pattern")  then
       if (allocated(matrix%val)) then
          deallocate(matrix%val,stat = ierr)
          if (present(stat)) stat = ierr
          if (ierr /= 0) then
             info = MC65_ERR_MEMORY_DEALLOC
          end if
       end if
       return
    end if

    ! check if matrix%val is allocated already, if so
    ! cehck if the size is large enough
    if (allocated(matrix%val)) then
       if (size(matrix%val) < nz1) then
          deallocate(matrix%val,stat = ierr)
          if (present(stat)) stat = ierr
          if (ierr /= 0) then
             info = MC65_ERR_MEMORY_DEALLOC
             return
          end if
          allocate(matrix%val(nz1),stat = ierr)
          if (present(stat)) stat = ierr
       end if
    else
       allocate(matrix%val(nz1),stat = ierr)
       if (present(stat)) stat = ierr
    end if
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

  end subroutine csr_matrix_construct

  subroutine csr_matrix_destruct(matrix,info,stat)
    ! subroutine csr_matrix_destruct(matrix,info):
    !
    ! destruct the matrix object by deallocating all
    ! space occupied by
    ! matrix. including matrix%ptr, matrix%col and matrix%type.
    ! when matrix%type /= "pattern", also
    ! deallocate matrix%val

    ! matrix: is of the derived type ZD11_type,
    !         with INTENT (INOUT). It
    !         the sparse matrix object to be destroyed.
    type (zd11_type), intent (inout) :: matrix

    ! info: is an integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
    integer (kind = myint), intent (out) :: info

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! ===================== local variables =============
    ! ierr: error tag for deallocation
    integer (kind = myint) :: ierr

    info = 0
    if (present(stat)) stat = 0

    deallocate(matrix%col,matrix%ptr, stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if


    if (ZD11_get(matrix%type) == "pattern")  then
       deallocate(matrix%type, stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          info = MC65_ERR_MEMORY_DEALLOC
       end if
       return
    end if

    deallocate(matrix%type, matrix%val, stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if

  end subroutine csr_matrix_destruct

  subroutine csr_matrix_reallocate(matrix,nz,info,stat)
    ! subroutine csr_matrix_reallocate(matrix,nz,info)
    !
    ! reallocate a space of nz for the column index array
    ! and when matrix%type /= "pattern", for the value array
    ! of the sparse matrix in the compact sparse row format.
    ! if nz < 0, a space of 0 will be allocated.
    ! The original context of the array(s),
    ! will be copied to the beginning of the reallocated
    ! array(s). When nz is smaller than the size of
    ! the said array(s), only part of the original
    ! array(s) up to the nz-th element is copied.

    ! matrix: is of the derived type ZD11_type with INTENT (INOUT),
    !         this is the sparse matrix to be reallocated.
    type (zd11_type), intent (inout) :: matrix

    ! nz:  is an integer scaler of INTENT (IN). It holds the
    !     space to be reallocated for the array of
    !     the column indices; and when the matrix
    !     is of type "pattern",
    !     also holds the space to be reallocated
    !     for the array of the entry
    !     values of the matrix.
    !     If nz < 0, a space of 0 will be allocated.

    integer (kind = myint), intent (in) :: nz

    ! info: is an integer scaler of INTENT (OUT).
    !       = 0  if the subroutine returns successfully;
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed.
    integer (kind = myint), intent (out) :: info

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! ================= local variables ================
    ! nz1: actual space allocated
    integer (kind = myint) :: nz1, ierr

    info = 0
    if (present(stat)) stat = 0

    if (nz < 0) then
       nz1 = 0
    else
       nz1 = nz
    end if

    if (size(matrix%col) /= nz1) then
       call expand(matrix%col,nz1,info,ierr)
       if (present(stat)) stat = ierr
       if (info < 0) then
          return
       end if
    end if

    if (ZD11_get(matrix%type) == "pattern")  return

    if (size(matrix%val) /= nz1) then
       call expand(matrix%val,nz1,info,ierr)
       if (present(stat)) stat = ierr
    end if


  end subroutine csr_matrix_reallocate

  subroutine csr_matrix_transpose(MATRIX1,MATRIX2,info,merge,pattern,stat)
    ! subroutine csr_matrix_transpose(MATRIX1,MATRIX2[,merge])
    !
    ! find the transpose of matrix MATRIX1 and
    ! put it in matrix MATRIX2. MATRIX2 = MATRIX1^T
    ! If merge is present and merge = .true.,
    ! repeated entries in MATRIX1 will be summed
    ! in MATRIX2
    ! If pattern is present and pattern = .true.,
    ! the new matrix will be of type "pattern" and
    ! will have no values

    ! MATRIX1: of the derived type ZD11_type, INTENT (INOUT),
    !    the matrix to be transposed
    type (zd11_type), intent (in) :: MATRIX1
    ! MATRIX2: of the derived type ZD11_type, INTENT (INOUT),
    !    MATRIX2 = MATRIX1^T. If merge = .true.,
    !    repeated entries of MATRIX1 are
    !    summed in MATRIX2.
    type (zd11_type), intent (inout) :: MATRIX2

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
    integer (kind = myint), intent (out) :: info

    ! merge: optional logical of INTENT (IN).
    !        if merge = .true, repeated entries of
    !        MATRIX1 will be merged in MATRIX2
    !        by summing the entry values. By default,
    !        no merge is performed.
    logical, intent (in), optional :: merge

    ! pattern: optional logical of INTENT (IN).
    !          if pattern = .true., MATRIX2 will have only
    !          the pattern of MATRIX1^T
    logical, intent (in), optional :: pattern

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! ================== local variables ===============
    !
    ! rowsz: number of entries in each row of MATRIX2
    integer (kind = myint), allocatable, dimension (:) :: rowsz


    ! mask: a masking array to see if a particular entry in a row
    !    has already appeared before.
    integer (kind = myint), dimension (:), allocatable :: mask

    integer (kind = myint) :: n,m,ierr,nz

    ! to_merge: whether repeated entries should be
    ! merged by summing the entry values
    logical :: to_merge

    ! pattern_only: whether the transpose matrix will be
    !    pattern only.
    logical :: pattern_only

    info = 0
    if (present(stat)) stat = 0

    to_merge = .false.
    pattern_only = .false.
    if (present(pattern)) then
       pattern_only = pattern
    end if
    if (ZD11_get(MATRIX1%type) == "pattern") pattern_only = .true.

    if (present(merge)) then
       to_merge = merge
    end if

    m = MATRIX1%n
    n = MATRIX1%m

    ! number of entries in each row of MATRIX2
    allocate(rowsz(m),stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if
    rowsz = 0

    if (to_merge) then
       allocate(mask(m),stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          info = MC65_ERR_MEMORY_ALLOC
          return
       end if
       mask = 0
       CALL csr_matrix_transpose_rowsz(n,MATRIX1%ptr,MATRIX1%col,rowsz,mask)
    else
       CALL csr_matrix_transpose_rowsz(n,MATRIX1%ptr,MATRIX1%col,rowsz)
    end if

    nz = sum(rowsz)

    if (pattern_only) then
       call csr_matrix_construct(MATRIX2,m,nz,info,n = n,&
            type = "pattern",stat=ierr)
    else
       call csr_matrix_construct(MATRIX2,m,nz,info,n = n,&
            type = ZD11_get(MATRIX1%type),stat=ierr)
    end if
    if (present(stat)) stat = ierr

    if (info < 0 ) return

    if (pattern_only) then
       if (to_merge) then
          call csr_matrix_transpose_pattern(m,n,MATRIX1%ptr,MATRIX1%col, &
                                  MATRIX2%ptr,MATRIX2%col,rowsz,mask)
       else
          call csr_matrix_transpose_pattern(m,n,MATRIX1%ptr,MATRIX1%col, &
                                  MATRIX2%ptr,MATRIX2%col,rowsz)
       end if
    else
       if (to_merge) then
          call csr_matrix_transpose_values(m,n,MATRIX1%ptr,MATRIX1%col, &
                   MATRIX1%val,MATRIX2%ptr,MATRIX2%col,MATRIX2%val,rowsz,mask)
       else
          call csr_matrix_transpose_values(m,n,MATRIX1%ptr,MATRIX1%col, &
                   MATRIX1%val,MATRIX2%ptr,MATRIX2%col,MATRIX2%val,rowsz)
       end if
    end if

    if (to_merge) then
       deallocate(mask,stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          info = MC65_ERR_MEMORY_DEALLOC
          return
       end if
    end if

    deallocate(rowsz,stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
    end if

  end subroutine csr_matrix_transpose


  subroutine csr_matrix_transpose_rowsz(n,ia,ja,rowsz,mask)
    ! get the number of elements in each row
    ! this is an internal private routine
    ! called by csr_matrix_transpose

    ! n: row dimension of A
    ! ia: row pointer of A
    ! ja: column indices of A
    ! rowsz: the number of elements in each row
    ! mask: zero on entry. On exit mask <= 0. If mask(k) = -i,
    !       column index
    !       k has been seen in row i of A.
    !       this argument is optional, and only present
    !       if repeated entries in A is to be merged when forming A^T
    integer (kind = myint), intent (in) :: n,ia(*),ja(*)
    integer (kind = myint), intent (inout) :: rowsz(*)
    integer (kind = myint), intent (inout), optional :: mask(*)
    integer (kind = myint) :: i,j,k

    if (present(mask)) then
       do i=1,n
          do j = ia(i),ia(i+1)-1
             k = ja(j)
             if (mask(k) /= -i) then
                mask(k) = -i
                rowsz(k) = rowsz(k) + 1
             end if
          end do
       end do
    else
       do i=1,n
          do j = ia(i),ia(i+1)-1
             k = ja(j)
             rowsz(k) = rowsz(k) + 1
          end do
       end do
    end if
  end subroutine csr_matrix_transpose_rowsz



  subroutine csr_matrix_transpose_values(m,n,ia,ja,aa,ib,jb,&
       bb,rowsz,mask)
    ! transpose the matrix A=(ia,ja,aa) that is not pattern only
    ! to give B=(ib,jb,bb). If mask is present,
    ! repeated entries in A will be merged
    ! this is an internal private routine
    ! called by csr_matrix_transpose

    ! n: row dimension of A
    ! m: column dimension of A
    ! ia: row pointer of A
    ! ja: column indices of A
    ! aa: values of A
    ! ib: row pointer of A
    ! jb: column indices of A
    ! bb: values of B
    ! rowsz: the number of elements in each row
    ! mask: zero on entry. On exit mask <= 0. If mask(k) = -i,
    !       column index
    !       k has been seen in row i of A.
    !       this argument is optional, and only present
    !       if repeated entries in A is to be merged when forming A^T
    integer (kind = myint), intent (in) :: m,n,ia(*),ja(*)
    integer (kind = myint), intent (out) :: ib(*),jb(*)
    real (kind = myreal), intent(in) :: aa(*)
    real (kind = myreal), intent(out) :: bb(*)
    integer (kind = myint), intent (inout) :: rowsz(*)
    integer (kind = myint), intent (inout), optional :: mask(*)

    integer (kind = myint) :: i,j,k,l1,l2

    ib(1) = 1
    do i = 1, m
       ib(i+1) = ib(i) + rowsz(i)
    end do
    rowsz(1:m) = ib(1:m)

    if (present(mask)) then
       do i=1,n
          l1 = ia(i); l2 = ia(i+1)-1
          do j = l1,l2
             k = ja(j)
             if (mask(k) <= 0) then
                ! record the position of (k,i) in the
                ! transpose matrix
                mask(k) = rowsz(k)
                jb(rowsz(k)) = i
                bb(rowsz(k)) = aa(j)
                rowsz(k) = rowsz(k) + 1
             else
                bb(mask(k)) = bb(mask(k))+aa(j)
             end if
          end do
          ! Note: Can't use following array statement instead of the next loop
          ! as there may be repeated entries in a column.
          ! mask(ja(l1:l2)) = 0
          do j = l1, l2
             mask(ja(j)) = 0
          end do
       end do

    else
       do i=1,n
          do j = ia(i),ia(i+1)-1
             k = ja(j)
             jb(rowsz(k)) = i
             bb(rowsz(k)) = aa(j)
             rowsz(k) = rowsz(k) + 1
          end do
       end do
    end if

  end subroutine csr_matrix_transpose_values




  subroutine csr_matrix_transpose_pattern(m,n,ia,ja,ib,jb,rowsz,mask)
    ! transpose the matrix A=(ia,ja) that is  pattern only
    ! to give B=(ib,jb). If mask is present,
    ! repeated entries in A will be merged
    ! this is an internal private routine
    ! called by csr_matrix_transpose

    ! n: row dimension of A
    ! m: column dimension of A
    ! ia: row pointer of A
    ! ja: column indices of A
    ! aa: values of A
    ! ib: row pointer of A
    ! jb: column indices of A
    ! bb: values of B
    ! rowsz: the number of elements in each row
    ! mask: zero on entry. On exit mask <= 0. If mask(k) = -i,
    !       column index
    !       k has been seen in row i of A.
    !       this argument is optional, and only present
    !       if repeated entries in A is to be merged when forming A^T
    integer (kind = myint), intent (in) :: m,n,ia(*),ja(*)
    integer (kind = myint), intent (out) :: ib(*),jb(*)
    integer (kind = myint), intent (inout) :: rowsz(*)
    integer (kind = myint), intent (inout), optional :: mask(*)

    integer (kind = myint) :: i,j,k,l1,l2
    ib(1) = 1
    do i = 1, m
       ib(i+1) = ib(i) + rowsz(i)
    end do
    rowsz(1:m) = ib(1:m)

    if (present(mask)) then
       do i=1,n
          l1 = ia(i); l2 = ia(i+1)-1
          do j = l1,l2
             k = ja(j)
             if (mask(k) <= 0) then
                ! record the position of (k,i) in the
                ! transpose matrix
                mask(k) = rowsz(k)
                jb(rowsz(k)) = i
                rowsz(k) = rowsz(k) + 1
             end if
          end do
          ! Note: Can't use following array statement instead of the next loop
          ! as there may be repeated entries in a column.
          ! mask(ja(l1:l2)) = 0
          do j = l1, l2
             mask(ja(j)) = 0
          end do
       end do

    else
       do i=1,n
          do j = ia(i),ia(i+1)-1
             k = ja(j)
             jb(rowsz(k)) = i
             rowsz(k) = rowsz(k) + 1
          end do
       end do
    end if

  end subroutine csr_matrix_transpose_pattern


  subroutine csr_matrix_copy(MATRIX1,MATRIX2,info,stat)
    ! subroutine csr_matrix_copy(MATRIX1,MATRIX2,info)
    !
    ! Subroutine {\tt MC65\_MATRIX\_COPY} creates a new sparse matrix
    ! {\tt B} which is a copy of the existing sparse matrix
    ! {\tt MATRIX1}.
    ! Storage is allocated for {\tt MATRIX2} and this storage should
    ! be deallocated when {\tt MATRIX2} is no longer used by calling
    ! subroutine {\tt MC65_MATRIX_DESTRUCT}.


    ! MATRIX1: of the derived type ZD11_type, INTENT (IN),
    !    the matrix to be copied from
    type (zd11_type), intent (in) :: MATRIX1

    ! MATRIX2: of the derived type ZD11_type, INTENT (IN),
    !    the matrix to be copied to
    type (zd11_type), intent (inout) :: MATRIX2

    ! info: integer scaler of INTENT (INOUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
    integer (kind = myint), intent (out) :: info

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! n: column size
    ! m: row size
    ! nz: number of entries in a.
    integer (kind = myint) :: nz,n,m

    if (present(stat)) stat = 0

    nz = MATRIX1%ptr(MATRIX1%m+1)-1
    n = MATRIX1%n
    m = MATRIX1%m
    if (present(stat)) then
       call csr_matrix_construct(MATRIX2,m,nz,info, n = n,&
            type = zd11_get(MATRIX1%type),stat=stat)
    else
       call csr_matrix_construct(MATRIX2,m,nz,info, n = n,&
            type = zd11_get(MATRIX1%type))
    end if
    if (info < 0) return
    MATRIX2%ptr(1:m+1) = MATRIX1%ptr(1:m+1)
    MATRIX2%col(1:nz) = MATRIX1%col(1:nz)
    if (ZD11_get(MATRIX1%type) == "pattern")  return
    MATRIX2%val(1:nz) = MATRIX1%val(1:nz)
  end subroutine csr_matrix_copy



  subroutine csr_matrix_clean(matrix,info,realloc,stat)
    ! subroutine csr_matrix_clean(matrix,info[,realloc])
    !
    ! cleaning the matrix by removing redundant entries
    ! for pattern only matrix, or
    ! by summing it with existing ones for non-pattern only matrix

    ! matrix:  is of the derived type {\tt ZD11\_type} with
    !     INTENT (INOUT),
    !     it is the  the sparse matrix to be cleaned
    type (zd11_type), intent (inout) :: matrix

    ! info: is a integer scaler of INTENT (OUT).
    ! {\tt INFO = 0} if the subroutine completes successfully;
    ! {\tt INFO = MC65\_ERR\_MEMORY\_ALLOC} if memory allocation
    !     failed; and
    ! {\tt INFO = MC65\_ERR\_MEMORY\_DEALLOC} if memory
    !     deallocation failed;
    ! INFO = MC65_WARN_DUP_ENTRY if duplicate entries were
    !  found when cleaning the matrix
    ! Error message can be printed by calling subroutine
    !    {\tt MC65\_ERROR\_MESSAGE}
    integer (kind = myint), intent (out) :: info

    ! realloc:  is an optional real scaler of INTENT (IN).
    !     It is used to control the reallocation of storage.
    ! \begin{itemize}
    !  \item{} if {\tt REALLOC < 0}, no reallocation.
    !  \item{} if {\tt REALLOC == 0}, reallocation is carried out
    !  \item{} if {\tt REALLOC > 0}, reallocation iff
    !       memory saving is greater than
    !      {\tt REALLOC*100}\%. For example,
    !      if {\tt REALLOC = 0.5}, reallocation will be carried out
    !      if saving in storage is greater than 50\%
    !  \end{itemize}
    ! If {\tt REALLOC} is not present, no reallocation
    !    is carried out.
    real (kind = myreal), INTENT(IN), optional :: realloc

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! m: row dimension of matrix
    ! n: column dimension of the matrix
    ! nz: number of entries after the matrix is cleaned
    integer (kind = myint) :: n,m,nz

    ! ia: row pointer of the original matrix
    ! ja: column entries
    ! ib:  row pointer of the cleaned matrix
    integer (kind = myint), allocatable, dimension(:) :: ib

    integer (kind = myint) :: ierr


    ! mask: mask array which tells if an entry has been seen before
    integer (kind = myint), allocatable, dimension (:) :: mask

    ! nz_old: space occupied before cleaning
    ! nentry_old: number of entries in the matrix before cleaning
    integer (kind = myint) :: nz_old, nentry_old

    info = 0
    if (present(stat)) stat = 0
    nz_old = size(matrix%col)
    m = matrix%m
    n = matrix%n
    nentry_old = matrix%ptr(m+1)-1

    allocate(mask(n), ib(m+1), stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if
    mask = 0

    if (ZD11_get(matrix%type) /= "pattern" ) then
       call csr_matrix_clean_private(m,nz,ib,mask,matrix%ptr,matrix%col, &
                                     a = matrix%val)
    else
       call csr_matrix_clean_private(m,nz,ib,mask,matrix%ptr,matrix%col)
    end if

    if (present(realloc)) then
       if (realloc == 0.0_myreal.or.(realloc > 0.and.&
            &nz_old > (1.0_myreal + realloc)*nz)) then
          call csr_matrix_reallocate(matrix,nz,info,stat=ierr)
          if (present(stat)) stat = ierr
          if (info < 0) return
       end if
    end if

    deallocate(mask,ib,stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
    end if

    if (nentry_old /= matrix%ptr(m+1)-1) then
       INFO = MC65_WARN_DUP_ENTRY
    END if

  end subroutine csr_matrix_clean


  subroutine csr_matrix_clean_private(m,nz,ib,mask,ia,ja,a)
    ! clean the matrix a = (ia,ja,a) by summing repeated entries
    ! or clean the matrix a = (ia,ja) by removing repeated entries
    ! this is a private routine called by csr_matrix_clean.
    !
    ! m: row dimension
    ! nz: number of nonzeros in the cleaned matrix
    ! ia: row pointer
    ! ja: column indices
    ! a: values. Optional
    ! ib: new row pointer for cleaned matrix
    ! mask: naming array set to zero on entry
    integer (kind = myint), intent (in) :: m
    integer (kind = myint), intent(inout) :: ia(*),ja(*),mask(*)
    integer (kind = myint), intent(out) :: nz,ib(*)
    real (kind = myreal), intent (inout), optional :: a(*)

    integer (kind = myint) :: i,j,jj

    ib(1) = 1
    nz = 1
    if (present(a)) then
       do i = 1, m
          do j = ia(i),ia(i+1)-1
             jj = ja(j)
             if (mask(jj) == 0) then
                mask(jj) = nz
                ja(nz) = jj
                a(nz) = a(j)
                nz = nz + 1
             else
                a(mask(jj)) = a(mask(jj)) + a(j)
             end if
          end do
          ib(i+1) = nz
          do j = ib(i),nz-1
             mask(ja(j)) = 0
          end do
       end do
    else
       do i = 1, m
          do j = ia(i),ia(i+1)-1
             jj = ja(j)
             if (mask(jj) == 0) then
                mask(jj) = nz
                ja(nz) = jj
                nz = nz + 1
             end if
          end do
          ib(i+1) = nz
          do j = ib(i),nz-1
             mask(ja(j)) = 0
          end do
       end do
    end if
    nz = nz - 1
    ia(1:m+1) = ib(1:m+1)
  end subroutine csr_matrix_clean_private


  subroutine csr_matrix_sort(matrix,info,stat)
    ! subroutine csr_matrix_sort(matrix,info)
    !
    ! sort all rows of a matrix into increasing column indices
    ! This is achieved by doing transpose twice.
    ! Repeated entries will be merged.

    ! matrix: is of the derived type ZD11_type, INTENT (inOUT),
    !         the sparse matrix in compressed sparse row format
    type (zd11_type), intent (inout) :: matrix

    ! info: integer scaler of INTENT (OUT).
    !      INFO = 0 if the subroutine completes successfully;
    !      INFO = MC65_ERR_MEMORY_ALLOC if memory
    !      allocation failed; and
    !      INFO = MC65_ERR_MEMORY_DEALLOC if memory
    !      deallocation failed.
    integer (kind = myint), intent (out) :: info

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! temporary copy of transpose of the matrix
    type (zd11_type) :: matrix_trans

    if (present(stat)) stat = 0

    if ( present(stat) ) then
       call csr_matrix_transpose(matrix,matrix_trans,info,merge=.true.,&
            stat=stat)
    else
       call csr_matrix_transpose(matrix,matrix_trans,info,merge=.true.)
    end if
    if (info < 0) return

    if ( present(stat) )then
       call csr_matrix_destruct(matrix,info,stat=stat)
    else
       call csr_matrix_destruct(matrix,info)
    end if
    if (info < 0) return

    if ( present(stat) ) then
       call csr_matrix_transpose(matrix_trans,matrix,info,stat=stat)
    else
       call csr_matrix_transpose(matrix_trans,matrix,info)
    end if
    if (info < 0) return

    if ( present(stat) ) then
       call csr_matrix_destruct(matrix_trans,info,stat=stat)
    else
       call csr_matrix_destruct(matrix_trans,info)
    end if

  end subroutine csr_matrix_sort

  subroutine csr_matrix_sum(matrix1,matrix2,result_matrix,&
       info,scaling,graph,stat)
    ! subroutine csr_matrix_sum(matrix1,matrix2,result_matrix,
    !       info[,scaling,graph])
    ! Subroutine {\tt MC65\_MATRIX\_SUM} adds two matrix
    !   {\tt MATRIX1} and {\tt MATRIX2} together,
    !    and return the resulting matrix in {\tt RESULT\_MATRIX}.
    !   If the {\tt OPTIONAL} argument {\tt GRAPH}
    !   is present and is {\tt .TRUE.},
    !   the {\tt RESULT\_MATRIX} is the sum of {\tt MATRIX1}
    !   and {\tt MATRIX2}, but with diagonal removed and
    !   all other entries set to 1.0
    !   (1.0D0 for {\tt HSL\_MC65\_DOUBLE}).
    !    Otherwise if the {\tt OPTIONAL} argument
    !   {\tt SCALING} is present,
    !    the scaled summation {\tt Matrix1+SCALING*MATRISX2}
    !    is calculated.

    ! matrix1: of the derived type ZD11_type, INTENT (IN),
    !         a sparse matrix in compressed sparse row format
    ! matrix2: of the derived type ZD11_type, INTENT (IN),
    !         a sparse matrix in compressed sparse row format
    type (zd11_type), intent (in) :: matrix1,matrix2

    ! \itt{result_matrix} is of the derived type
    !   {\tt ZD11\_type} with INTENT (inOUT),
    !   result_matrix is one of the following
    ! \begin{itemize}
    !     \item{} If {\tt GRAPH} is present and is
    !       {\tt .TRUE.}, {\tt RESULT\_MATRIX}
    !       is a matrix of type ``general'', and is set to the sum
    !       of {\tt MATRIX1} and {\tt MATRIX2},
    !       with the diagonal removed
    !       and with all entry values set to 1.0
    !     \item{} Otherwise {\tt RESULT\_MATRIX =
    !       MATRIX1 + MATRIX2}. The type of {\tt RESULT\_MATRIX}
    !       is set to ``pattern'' either {\tt MATRIX1}
    !       or {\tt MATRIX2} is of this type, or else it is
    !       set to the type of {\tt MATRIX1} if
    !       both {\tt MATRIX1} and {\tt MATRIX2}
    !       have the same type, otherwise
    !       it is set as ``general''.
    !     \item{} If the {\tt OPTIONAL} argument
    !       {\tt SCALING} is present,
    !       {\tt RESULT\_MATRIX = MATRIX1 + SCALING*MATRIX2}
    ! \end{itemize}
    type (zd11_type), intent (inout) :: result_matrix

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
    !       = MC65_ERR_SUM_DIM_MISMATCH if row (column)
    !         dimension of matrix1 and matrix2 does not match.
    integer (kind = myint), intent (out) :: info

    ! scaling:  matrix2 is to be scaled by *scaling is
    real (kind = myreal), intent (in), optional :: scaling

    ! graph: optional logical scaler of intent (in).
    !        When present and is .true.
    !        the diagonals of the summation is removed, and
    !        entry values set to one if the result_matrix is
    !        not pattern only.
    logical, intent(in), optional :: graph

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! m: row dimension of matrix1
    ! n: column dimension of matrix1
    ! nz: number of entries in the summation matrix1+matrix2
    integer (kind = myint) :: n,m,nz

    integer (kind = myint) :: ierr

    ! mask: mask array to indicate if a particular column index
    !       has been seen
    integer (kind = myint), allocatable, dimension (:) :: mask
    logical :: is_graph

    info = 0
    if (present(stat)) stat = 0

    is_graph = .false.
    if (present(graph)) then
       is_graph = graph
    end if

    m = matrix1%m
    n = matrix1%n
    if (m /= matrix2%m .or. n /= matrix2%n) then
       info = MC65_ERR_SUM_DIM_MISMATCH
       return
    end if

    allocate (mask(n), stat = ierr)
    if ( present(stat) ) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

    mask = 0

    ! work out the number of entries of the summation.
    if (is_graph) then
       call csr_matrix_sum_getnz(m,matrix1%ptr,matrix1%col,matrix2%ptr, &
                       matrix2%col,mask,nz,graph = is_graph)
    else
       call csr_matrix_sum_getnz(m,matrix1%ptr,matrix1%col,matrix2%ptr, &
                       matrix2%col,mask,nz)
    end if

    ! construct the new matrix
    ! if one of the matrix is pattern only, then
    ! the result matrix has to be pattern only. If
    ! a graph is to be formed, the matrix is always of type "general"
    ! and filled with values 1.0
    if (is_graph) then
       call csr_matrix_construct(result_matrix,m,nz,info,&
            n = n,type = "general",stat=ierr)
    else if (ZD11_get(matrix1%type) == "pattern".or.&
         &ZD11_get(matrix2%type) == "pattern") then
       call csr_matrix_construct(result_matrix,m,nz,info,&
            n = n,type = "pattern",stat=ierr)
    else if (ZD11_get(matrix1%type) == ZD11_get(matrix2%type)) then
       call csr_matrix_construct(result_matrix,m,nz,info,n = n,&
            type = ZD11_get(matrix1%type),stat=ierr)
    else
       call csr_matrix_construct(result_matrix,m,nz,info,&
            n = n,type = "general",stat=ierr)
    end if
    if (present(stat)) stat=ierr
    if (info < 0) then
       return
    end if

    if (ZD11_get(result_matrix%type) /= "pattern" ) then
       IF (PRESENT(SCALING).and.(.not.is_graph)) THEN
          call csr_matrix_sum_values(m,nz,matrix1%ptr,matrix1%col,matrix1%val, &
                   matrix2%ptr,matrix2%col,matrix2%val,result_matrix%ptr, &
                   result_matrix%col,result_matrix%val,mask,scaling = scaling)
       ELSE if (is_graph) then
          call csr_matrix_sum_graph(m,nz,matrix1%ptr,matrix1%col, &
                   matrix2%ptr,matrix2%col,result_matrix%ptr, &
                   result_matrix%col,result_matrix%val,mask)
       else
          call csr_matrix_sum_values(m,nz,matrix1%ptr,matrix1%col,matrix1%val, &
                   matrix2%ptr,matrix2%col,matrix2%val,result_matrix%ptr, &
                   result_matrix%col,result_matrix%val,mask)
       END IF
    else
       call csr_matrix_sum_pattern(m,nz,matrix1%ptr,matrix1%col, &
                   matrix2%ptr,matrix2%col,result_matrix%ptr, &
                   result_matrix%col,mask)
    end if

    deallocate(mask,stat = ierr)
    if ( present(stat) ) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
    end if

  end subroutine csr_matrix_sum


  subroutine csr_matrix_sum_getnz(m,ia,ja,ib,jb,mask,nz,graph)
    ! get the number of nonzeros in the sum of two matrix
    ! (ia,ja) and (ib,jb), and return that number in nz
    ! this is a private routine called by csr_matrix_sum
    !
    ! m: row dimension of A
    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    ! mask: masking array to show whether a particular entry
    !       in the summed matrix has already been see. Initialised
    !       to zero on entry, on exit <= 0.
    ! nz: number of nonzeros in the summed matrix
    ! GRAPH: optional. If present and is true, diagonals
    !        in the summation matrix is not counted
    integer (kind = myint), intent (in) :: m,ia(*),ib(*),ja(*),jb(*)
    integer (kind = myint), intent (inout) :: mask(*)
    integer (kind = myint), intent (out) :: nz
    logical, intent (in),optional :: graph

    logical :: is_graph

    integer (kind = myint) :: i,j,jj

    is_graph = .false.
    if (present(graph)) then
       is_graph = graph
    end if

    if (is_graph) then
       nz = 0
       do i = 1, m
          do j = ia(i), ia(i+1)-1
             jj = ja(j)
             if (jj == i) cycle
             if (mask(jj) /= -i) then
                mask(jj) = -i
                nz = nz + 1
             end if
          end do
          do j = ib(i), ib(i+1)-1
             jj = jb(j)
             if (jj == i) cycle
             if (mask(jj) /= -i) then
                mask(jj) = -i
                nz = nz + 1
             end if
          end do
       end do
    else

       nz = 0
       do i = 1, m
          do j = ia(i), ia(i+1)-1
             jj = ja(j)
             if (mask(jj) /= -i) then
                mask(jj) = -i
                nz = nz + 1
             end if
          end do
          do j = ib(i), ib(i+1)-1
             jj = jb(j)
             if (mask(jj) /= -i) then
                mask(jj) = -i
                nz = nz + 1
             end if
          end do
       end do
    end if
  end subroutine csr_matrix_sum_getnz


  subroutine csr_matrix_sum_values(m,nz,ia,ja,a,ib,jb,b,ic,jc,c,&
       mask,scaling)
    ! sum the two matrix A=(ia,ja,a) and B=(ib,jb,b) into C=(ic,jc,c)
    ! if scaling is present, B will be scaled by it before summing
    !
    !
    ! m: row dimension of A
    ! nz: number of nonzeros in the summation C
    ! ia: row pointer of A
    ! ja: column indices of A
    ! a: values of A
    ! ib: row pointer of B
    ! jb: column indices of B
    ! b: values of B
    ! ic: row pointer of C
    ! jc: column indices of C
    ! c: values of C
    ! mask: masking array to show whether a particular entry
    !       in the summed matrix has already been see. Initialised
    !       to <= 0 on entry
    integer (kind = myint), intent (in) :: m,ia(*),ib(*),ja(*),jb(*)
    real (kind = myreal), intent (in) :: a(*),b(*)
    integer (kind = myint), intent (out) :: ic(*),jc(*)
    real (kind = myreal), intent (out) :: c(*)
    integer (kind = myint), intent (inout) :: mask(*)
    real (kind = myreal), intent (in), optional :: scaling
    integer (kind = myint), intent (out) :: nz


    integer (kind = myint) :: i,j,jj


    ic(1) = 1
    nz = 1
    if (.not.present(scaling)) then
       do i = 1, m
          do j = ia(i), ia(i+1)-1
             jj = ja(j)
             if (mask(jj) <= 0) then
                mask(jj) = nz
                jc(nz) = jj
                c(nz) = a(j)
                nz = nz + 1
             else
                c(mask(jj)) = c(mask(jj)) +  a(j)
             end if
          end do
          do j = ib(i),ib(i+1)-1
             jj = jb(j)
             if (mask(jj) <= 0) then
                mask(jj) = nz
                jc(nz) = jj
                c(nz) = b(j)
                nz = nz + 1
             else
                c(mask(jj)) = c(mask(jj)) +  b(j)
             end if
          end do
          ic(i+1) = nz
          do j = ic(i),nz-1
             mask(jc(j)) = 0
          end do
       end do
    else
       do i = 1, m
          do j = ia(i), ia(i+1)-1
             jj = ja(j)
             if (mask(jj) <= 0) then
                mask(jj) = nz
                jc(nz) = jj
                c(nz) = a(j)
                nz = nz + 1
             else
                c(mask(jj)) = c(mask(jj)) +  a(j)
             end if
          end do
          do j = ib(i),ib(i+1)-1
             jj = jb(j)
             if (mask(jj) <= 0) then
                mask(jj) = nz
                jc(nz) = jj
                c(nz) = scaling*b(j)
                nz = nz + 1
             else
                c(mask(jj)) = c(mask(jj)) +  scaling*b(j)
             end if
          end do
          ic(i+1) = nz
          do j = ic(i),nz-1
             mask(jc(j)) = 0
          end do
       end do

    end if
    nz = nz-1
  end subroutine csr_matrix_sum_values



  subroutine csr_matrix_sum_graph(m,nz,ia,ja,ib,jb,ic,jc,c,mask)
    ! sum the two matrix A=(ia,ja,a) and B=(ib,jb,b) into C=(ic,jc,c)
    ! diagonals of C is ignored and
    ! all entries of C set to one
    !
    ! m: row dimension of A
    ! nz: number of nonzeros in the summation C
    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    ! ic: row pointer of C
    ! jc: column indices of C
    ! c: values of C
    ! mask: masking array to show whether a particular entry
    !       in the summed matrix has already been see. Initiualised
    !       to <= 0 on entry
    integer (kind = myint), intent (in) :: m,ia(*),ib(*),ja(*),jb(*)
    integer (kind = myint), intent (out) :: ic(*),jc(*)
    real (kind = myreal), intent (out) :: c(*)
    integer (kind = myint), intent (inout) :: mask(*)
    integer (kind = myint), intent (out) :: nz


    integer (kind = myint) :: i,j,jj


    ic(1) = 1
    nz = 1
    do i = 1, m
       do j = ia(i), ia(i+1)-1
          jj = ja(j)
          if (jj == i) cycle
          if (mask(jj) <= 0) then
             mask(jj) = nz
             jc(nz) = jj
             c(nz) =  1.0_myreal
             nz = nz + 1
          else
             c(mask(jj)) = 1.0_myreal
          end if
       end do
       do j = ib(i),ib(i+1)-1
          jj = jb(j)
          if (jj == i) cycle
          if (mask(jj) <= 0) then
             mask(jj) = nz
             jc(nz) = jj
             c(nz) = 1.0_myreal
             nz = nz + 1
          else
             c(mask(jj)) = 1.0_myreal
          end if
       end do
       ic(i+1) = nz
       do j = ic(i),nz-1
          mask(jc(j)) = 0
       end do
    end do
    nz = nz-1
  end subroutine csr_matrix_sum_graph


  subroutine csr_matrix_sum_pattern(m,nz,ia,ja,ib,jb,ic,jc,mask)
    ! sum the two matrix A=(ia,ja) and B=(ib,jb, into C=(ic,jc),
    !    patterns only.
    ! m: row dimension of A
    ! nz: number of nonzeros in the summation C
    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    ! ic: row pointer of C
    ! jc: column indices of C
    ! mask: masking array to show whether a particular entry
    !       in the summed matrix has already been see. Initiualised
    !       to <= 0 on entry
    ! graph: optional logical, when present and is .true.,
    !        the diagonals of the C will not be generated.
    integer (kind = myint), intent (in) :: m,ia(*),ib(*),ja(*),jb(*)
    integer (kind = myint), intent (out) :: ic(*),jc(*)
    integer (kind = myint), intent (inout) :: mask(*)
    integer (kind = myint), intent (out) :: nz

    integer (kind = myint) :: i,j,jj

!!$    is_graph = .false.
!!$    if (present(graph)) then
!!$       is_graph = graph
!!$    end if


    ic(1) = 1
    nz = 1

!!$    if (is_graph) then
!!$       do i = 1, m
!!$          do j = ia(i), ia(i+1)-1
!!$             jj = ja(j)
!!$             if (jj == i) cycle
!!$             if (mask(jj) <= 0) then
!!$                mask(jj) = nz
!!$                jc(nz) = jj
!!$                nz = nz + 1
!!$             end if
!!$          end do
!!$          do j = ib(i),ib(i+1)-1
!!$             jj = jb(j)
!!$             if (jj == i) cycle
!!$             if (mask(jj) <= 0) then
!!$                mask(jj) = nz
!!$                jc(nz) = jj
!!$                nz = nz + 1
!!$             end if
!!$          end do
!!$          ic(i+1) = nz
!!$          do j = ic(i),nz-1
!!$             mask(jc(j)) = 0
!!$          end do
!!$       end do
!!$    else
       do i = 1, m
          do j = ia(i), ia(i+1)-1
             jj = ja(j)
             if (mask(jj) <= 0) then
                mask(jj) = nz
                jc(nz) = jj
                nz = nz + 1
             end if
          end do
          do j = ib(i),ib(i+1)-1
             jj = jb(j)
             if (mask(jj) <= 0) then
                mask(jj) = nz
                jc(nz) = jj
                nz = nz + 1
             end if
          end do
          ic(i+1) = nz
          do j = ic(i),nz-1
             mask(jc(j)) = 0
          end do
       end do
!    end if
    nz = nz-1
  end subroutine csr_matrix_sum_pattern

  subroutine csr_matrix_symmetrize(matrix,info,graph,stat)
    ! The subroutine {\tt MC65\_MATRIX\_SYMMETRIZE} make the
    !    matrix symmetric by summing it with its transpose.
    !    If the optional argument {\tt GRAPH}
    !    is present and is set to {\tt .TRUE.}, the diagonal of the
    !    summation is not included, and whenever entry values are
    !    available, they are set to {\tt 1.0}(or {\tt 1.0D0}for
    !    {\tt HSL\_MC65\_double}.

    ! matrix: of the derived type ZD11_type, INTENT (INOUT),
    !         the sparse matrix in compressed sparse row format
    !         to be made symmetric by matrix+matrix^T
    type (zd11_type), intent (inout) :: matrix

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
    !       = MC65_ERR_SUM_DIM_MISMATCH if trying to symmetrize
    !          a nonsquare matrix.
    integer (kind = myint), intent (out) :: info

    ! graph: optional logical scaler. when present and is .true.,
    !        the diagonals of matrix+matrix^T will be removed,
    !        and entry values set to 1.0_myreal matrix is
    !        not pattern only.
    logical, intent (in), optional :: graph

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! matrix1: transposed of matrix
    ! matrix2: A+A^T
    type (zd11_type) :: matrix1, matrix2

    ! is_graph: whether we are generating a graph by the summation
    logical :: is_graph

    if ( present(stat) ) stat = 0

    is_graph = .false.
    if (present(graph)) is_graph = graph
    if (present(stat)) then
       call csr_matrix_transpose(matrix,matrix1,info,stat=stat)
    else
       call csr_matrix_transpose(matrix,matrix1,info)
    end if
    if (info < 0) return

    if (is_graph) then
       if (present(stat)) then
          call csr_matrix_sum(matrix,matrix1,matrix2,info,graph=.true.,&
               &stat=stat)
       else
          call csr_matrix_sum(matrix,matrix1,matrix2,info,graph=.true.)
       end if
    else
       if (present(stat)) then
          call csr_matrix_sum(matrix,matrix1,matrix2,info,stat=stat)
       else
          call csr_matrix_sum(matrix,matrix1,matrix2,info)
       end if
    end if

    if (info < 0) return

    if (present(stat)) then
       call csr_matrix_copy(matrix2,matrix,info,stat=stat)
    else
       call csr_matrix_copy(matrix2,matrix,info)
    end if

    if (info < 0) return

    if (present(stat)) then
       call csr_matrix_destruct(matrix2,info,stat=stat)
    else
       call csr_matrix_destruct(matrix2,info)
    end if

    if (info < 0) return

    if (present(stat)) then
       call csr_matrix_destruct(matrix1,info,stat=stat)
    else
       call csr_matrix_destruct(matrix1,info)
    end if

    if (info < 0) return

  end subroutine csr_matrix_symmetrize


  subroutine csr_matrix_getrow(a,i,col)

    ! get the column indices of the i-th row of matrix a

    ! a: of the derived type ZD11_type, INTENT (IN), TARGET
    !    the sparse matrix in compressed sparse row format
    type (zd11_type), intent (in), target :: a

    ! i: is an integer of intent (in), it is the
    ! the index of a row of A to be extracted.
    integer (kind = myint), intent (in) :: i

    ! col: is the array of column indices for the row
    integer (kind = myint), pointer, dimension (:) :: col

    col => a%col(a%ptr(i):a%ptr(i+1)-1)

  end subroutine csr_matrix_getrow


  subroutine csr_matrix_getrowval(a,i,val)

    ! get the values of the i-th row of matrix a

    ! a: of the derived type ZD11_type, INTENT (IN), TARGET
    !    the sparse matrix in compressed sparse row format
    type (zd11_type), intent (in), target :: a

    ! i: the row index of A which the column values are
    !     to be extracted.
    integer (kind = myint), intent (in) :: i

    ! val: the real (or double) array of
    !     values of the i-th row of the matrix
    real (kind = myreal), pointer, dimension (:) :: val

    val => a%val(a%ptr(i):a%ptr(i+1)-1)

  end subroutine csr_matrix_getrowval


  subroutine csr_matrix_is_symmetric(matrix,symmetric,info,&
       pattern,tol,stat)
    ! function csr_matrix_is_symmetric(matrix,symmetric,&
    !    info[,pattern,tol,stat])
    !
    ! The subroutine {\tt MC65\_MATRIX\_IS\_SYMMETRIC}
    !   check whether a matrix is symmetric or not.
    !   {\tt symmetric = .TRUE.} if the matrix is symmetric,
    !   and {\tt symmetric= .FALSE.} if not. When the matrix
    !   is of type ``pattern'', only structural symmetry is
    !   checked, otherwise the matrix is said to be symmetric
    !   if it is both structural symmetric
    !   and that for each entry {\tt MATRIX(I,J)}, the inequality
    !  {\tt MATRIX(I,J) - MATRIX(J,I) <= E} holds.
    !  Here {\tt E} is a tolerance set by default to
    !  {\tt E = 0.0} (or {\tt E = 0.0D0) for HSL\_MC65\_DOUBLE}.

    ! Two optional arguments are supplied. If the optional argument
    !   {\tt PATTERN} is present and is set to {\tt .TRUE.},
    !   only structural symmetry is checked. If the optional
    !   argument {\tt TOL} is present, the internal tolerance
    !   {\tt E} is set to {\tt TOL}.



    ! matrix:  is of the derived type ZD11_type with INTENT (IN).
    !     It is the sparse matrix whose symmetry is to be checked.
    type (zd11_type), intent (in) :: matrix

    ! symmetric: is logical of intent (out).
    !             It is set to {\tt .TRUE.} only in one of
    !            the following cases
    !            1) When {\tt PATTERN} is not present, or
    !               is present but {\tt PATTERN /= .TRUE.}
    !               1a) the matrix is not of type "pattern" and
    !                   is symmetric both structurally
    !                   and numerically.
    !               1b) the matrix is of type "pattern" and
    !                   is structurally symmetric.
    !            2) when {\tt PATTERN} is present and
    !               {\tt PATTERN = .TRUE.}, and
    !               the matrix is structurally symmetric
    logical, intent (out) :: symmetric


    ! info: integer scaler of INTENT (OUT).
    !      {\tt INFO = 0} if the subroutine returns successfully;
    !      {\tt INFO = MC65_ERR_MEMORY_ALLOC} if memory
    !         allocation failed; and
    !      {\tt INFO = MC65_ERR_MEMORY_DEALLOC} if memory
    !        deallocation failed
    integer (kind = myint), intent (out) :: info

    ! pattern:  is an optional logical scaler of INTENT (IN).
    !         when it is present and is {\tt .TRUE.},
    !         only structural symmetry is checked.
    logical, optional, intent (in) :: pattern

    ! tol: an optional real scaler of INTENT (IN).
    !      tol is the tolerance used to decide if two
    !      entries are the same
    real (kind = myreal), optional, intent (in) :: tol

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! matrix_trans: the transpose of matrix
    type (zd11_type) :: matrix_trans

    ! m: row dimension of matrix
    ! n: column dimension of matrix
    ! nz: number of nonzero entries in the matrix
    integer (kind = myint) :: n,m


    integer (kind = myint), dimension (:), allocatable :: mask
    integer (kind = myint) :: ierr
    real (kind = myreal), dimension (:), allocatable :: val

    ! structural: = .true. if only structural symmetry is checked
    logical :: structural

    ! e: error allowed for judging symmetry on matrix that
    !    is not pattern only.
    real (kind = myreal) :: e

    info = 0
    if (present(stat)) stat = 0

    e = real(0,myreal)
    if (present(tol)) then
       e = tol
    end if

    structural = .false.
    if (present(pattern)) then
       structural = pattern
    end if
    if (zd11_get(matrix%type) == "pattern") structural = .true.

    symmetric = .true.
    m = matrix%m
    n = matrix%n
    if (n /= m) then
       symmetric = .false.
       return
    end if

    ! transpose matrix
    if (structural) then
       if (present(stat)) then
          call csr_matrix_transpose(matrix,matrix_trans,info,&
               merge = .true.,pattern = .true., stat = stat)
       else
          call csr_matrix_transpose(matrix,matrix_trans,info,&
               merge = .true.,pattern = .true.)
       end if
    else
       if (present(stat)) then
          call csr_matrix_transpose(matrix,matrix_trans,info,&
               merge = .true., stat=stat)
       else
          call csr_matrix_transpose(matrix,matrix_trans,info,&
               merge = .true.)
       end if
    end if
    if (info < 0) return

    allocate(mask(m),stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
    end if

    if (structural) then
       symmetric = csr_matrix_is_same_pattern(m,n,&
            ia = matrix%ptr,ja = matrix%col,&
            ib = matrix_trans%ptr,jb = matrix_trans%col,mask = mask)
    else
       allocate(val(m),stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          info = MC65_ERR_MEMORY_ALLOC
       end if
       symmetric = csr_matrix_is_same_values(m,n,&
            ia = matrix%ptr,ja = matrix%col,a = matrix%val, &
            ib = matrix_trans%ptr,jb = matrix_trans%col,&
            b = matrix_trans%val, &
            mask = mask,val = val,tol = e)
       ! check only structural symmetry as well as values
       deallocate(val,stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          info = MC65_ERR_MEMORY_DEALLOC
          return
       end if
    end if

    deallocate(mask,stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if
    call csr_matrix_destruct(matrix_trans,info)

  end subroutine csr_matrix_is_symmetric


  function csr_matrix_is_same_pattern(m,n,ia,ja,ib,jb,mask) &
       result(thesame)
    ! check if two matrix, A=(ia,ja) and B = (ib,jb),
    ! is structurally the same. Here B is assumed to have
    ! no repeated entries. A and B is assumed to have the same row
    ! size.
    ! This is a private function called by csr_matrix_is_symmetric
    !
    ! m: row size of A
    ! n: column size of A
    ! ia: row pointer array of A
    ! ja: column indices of A
    ! ib: row pointer array of B
    ! jb: column indices of B
    ! mask: a mask which indicate if a certain column entry k
    !       has been seen in this row i (if mask(k) = i). size >= n
    ! thesame: A and B is the same pattern
    integer (kind = myint), intent (in) :: m,n,ia(*),ja(*),ib(*),jb(*)
    integer (kind = myint), intent (inout) :: mask(*)
    logical :: thesame
    !
    ! na: number of distinctive entries in a row of A
    integer (kind = myint) :: na
    integer (kind = myint) :: i,j,k,l1,l2

    thesame = .true.

    mask(1:n) = 0
    do i = 1, m
       na = 0
       do j = ia(i),ia(i+1)-1
          k = ja(j)
          if (mask(k) /= i) then
             mask(k) = i
             na = na + 1
          end if
       end do
       l1 = ib(i); l2 = ib(i+1) - 1
       ! make sure the number of distinctive entries
       ! of the A and B is the same. Since B is assumed to have
       ! no repeated entries, its number of distinctive entries
       ! is simply l2-l1+1
       if (na /= l2 - l1 + 1) then
          thesame = .false.
          return
       end if
       do j = l1,l2
          if (mask(jb(j)) /= i) then
             thesame = .false.
             return
          end if
       end do
    end do
  end function csr_matrix_is_same_pattern





  function csr_matrix_is_same_values(m,n,ia,ja,a,ib,jb,b,mask,&
       val,tol) result(thesame)
    ! check if two matrix, A=(ia,ja,a) and B = (ib,jb,b),
    ! is the same structurally as well as in values.
    ! Here B is assumed to have
    ! no repeated entries. A and B is assumed to have the same row
    ! size as well as column size.
    ! This is a private function called by csr_matrix_is_symmetric
    !
    ! m: row size of A
    ! n: column size of A.
    ! ia: row pointer array of A
    ! ja: column indices of A
    ! a: entry values of A
    ! ib: row pointer array of B
    ! jb: column indices of B
    ! b: entry values of B
    ! mask: a mask which indicate if a certain column entry k
    !       has been seen in this row i (if mask(k) = i). Initialised
    !       to zero. size >= n
    ! val: the entry values in a row of A. If A has repeated entries,
    !      val will contain the sum of the repeated entry values.
    ! thesame: A and B is the same pattern
    ! tol: error allowed for judging symmetry on matrix that
    !     is not pettern only.
    integer (kind = myint), intent (in) :: m,n,ia(*),ja(*),ib(*),jb(*)
    integer (kind = myint), intent (inout) :: mask(*)
    real (kind = myreal), intent (in) :: a(*),b(*)
    real (kind = myreal), intent (inout) :: val(*)
    logical :: thesame
    real (kind = myreal), intent (in) :: tol
    !
    ! na: number of distinctive entries in a row of A
    integer (kind = myint) :: na
    integer (kind = myint) :: i,j,k,l1,l2

    thesame = .true.
    mask(1:n) = 0
    do i = 1, m
       na = 0
       do j = ia(i),ia(i+1)-1
          k = ja(j)
          if (mask(k) /= i) then
             mask(k) = i
             na = na + 1
             val(k) = a(j)
          else
             val(k) = val(k) + a(j)
          end if
       end do
       l1 = ib(i); l2 = ib(i+1) - 1
       ! make sure the number of distinctive entries
       ! of the A and B is the same. Since B
       ! is assumed to have no repeated entries,
       ! its number of distinctive entries
       ! is simply l2-l1+1
       if (na /= l2 - l1 + 1) then
          thesame = .false.
          return
       end if
       do j = l1,l2
          k = jb(j)
          if (mask(k) /= i) then
             thesame = .false.
             return
          else
             if (abs(val(k)-b(j)) > tol) then
                thesame = .false.
                return
             end if
          end if
       end do
    end do
  end function csr_matrix_is_same_values


  subroutine csr_matrix_diff(MATRIX1,MATRIX2,diff,info,pattern,stat)
    !   subroutine csr_matrix_diff(MATRIX1,MATRIX2,diff,&
    !      info[,pattern,stat])

    ! return the difference between two matrix MATRIX1 and MATRIX2.
    !    Here MATRIX1 and MATRIX2 are assumed to have not
    !    repeated entries.
    !
    ! 1) if the dimension of the matrix or the number of entries
    !    is different, huge(0.0_myreal) is returned.
    ! 2) else if the two matrices have different patterns,
    !    huge(0.0_myreal) is returned
    ! 3) else if both matrices are not pattern only, f1 is returned
    ! here: f1=sum(abs((a-b)%val))
    ! 4) else if either are pattern only, 0.0_myreal is returned
    !
    ! If the optional argument pattern is present and is .true.,
    !   both matrices are taken as pattern only.

    ! MATRIX1,MATRIX2: of the derived type ZD11_type, INTENT (IN),
    !      the sparse matrix in compressed sparse row format
    !      Two matrices to be compared.
    type (zd11_type), intent (in) :: MATRIX1,MATRIX2

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
    integer (kind = myint), intent(out) :: info

    ! pattern: optional logical scaler of INTENT (IN). If present and
    !          is .true., only the difference in pattern is returned.
    logical, optional, INTENT (IN) :: pattern


    ! diff:  is a real (double precision real in HSL_MC65_double)
    !         scaler returned on successful completion of the
    !         function. It is difference in the matrices
    !         {\tt MATRIX1} and {\tt MATRIX2} defined as follows
    !   \begin{itemize}
    !   \item{} If the row or column dimension of {\tt MATRIX1}
    !           and {\tt MATRIX2} does not match,
    !           {\tt HUGE\{0.0\}} ({\tt HUGE\{0.0D0\}} for
    !           {\tt HSL_MC65_double}) is returned;
    !   \item{} Else if {\tt MATRIX1} and {\tt MATRIX2} has different
    !            pattern, {\tt HUGE\{0.0\}} ({\tt HUGE\{0.0D0\}} for
    !           {\tt HSL_MC65_double}) is returned;
    !   \item{} Else if either {\tt MATRIX1} or {\tt MATRIX2}
    !            is of type {\tt ``pattern''},
    !            0.0 (0.0D0 for {\tt HSL_MC65_double}) is returned.
    !   \item{} Else if neither {\tt MATRIX1} nor {\tt MATRIX2}
    !            is of the type
    !            {\tt ``pattern''}, the sum of the absolute value
    !             of the difference in entry value,
    !            {\tt SUM(ABS((MATRIX1\%VAL - MATRIX2\%VAL))},
    !            is returned.
    !   \end{itemize}

    real (kind = myreal), intent (out) :: diff

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat


    ! ma,ma: row dimension of a and b
    ! na,nb: column dimension of a and b
    ! nza,nzb: number of entries of a and b and c
    integer (kind = myint) :: ma,mb,na,nb,nza,nzb


    ! pattern_only: only find the difference in patterns
    logical :: pattern_only,pattern_the_same

    integer (kind = myint), allocatable :: mask(:)

    integer (kind = myint) :: ierr


    info = 0
    if (present(stat)) stat = 0

    pattern_only = .false.
    if (present(pattern)) then
       pattern_only = pattern
    end if
    if (ZD11_get(MATRIX1%type) == "pattern".or.&
         ZD11_get(MATRIX2%type) == "pattern") &
         pattern_only = .true.


    ma = MATRIX1%m
    mb = MATRIX2%m
    na = MATRIX1%n
    nb = MATRIX2%n
    nza = MATRIX1%ptr(ma+1)-1
    nzb = MATRIX2%ptr(mb+1)-1
    if (ma /= mb .or. na /= nb .or. nza /= nzb) then
       diff = huge(0.0_myreal)
       return
    end if

    allocate(mask(na),stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if
    if (pattern_only) then
       pattern_the_same = csr_matrix_is_same_pattern(m = ma,n = na,&
            ia = MATRIX1%ptr,ja = MATRIX1%col,&
            ib = MATRIX2%ptr,jb = MATRIX2%col,mask = mask)
       if (pattern_the_same) then
          diff = 0.0_myreal
       else
          diff = huge(0.0_myreal)
       end if
    else
       diff = csr_matrix_diff_values(m = ma,n = na,&
            ia = MATRIX1%ptr,ja = MATRIX1%col,a=MATRIX1%val,&
            ib = MATRIX2%ptr,jb = MATRIX2%col,b=MATRIX2%val,&
            mask = mask)
    end if


    deallocate(mask, stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if

  end subroutine csr_matrix_diff


  function csr_matrix_diff_values(m,n,ia,ja,a,ib,jb,b,mask) &
       result(diff)
    ! find the difference in values of two matrices,
    !   A=(ia,ja,a) and B = (ib,jb,b),
    ! If the two are of different pattern,
    !   the difference is huge(0.0_myreal).
    ! Here A and B are assumed to have
    ! no repeated entries, they are assumed to have the same row
    ! size as well as column size.
    ! This is a private function called by csr_matrix_diff
    !
    ! m: row size of A
    ! n: column size of A.
    ! ia: row pointer array of A
    ! ja: column indices of A
    ! a: entry values of A
    ! ib: row pointer array of B
    ! jb: column indices of B
    ! b: entry values of B
    ! mask: a mask which indicate if a certain column entry k
    !       has been seen in this row i (if mask(k) = i). Initialised
    !       to zero. size >= n
    ! diff: difference in values between A and B.
    !       If A and B is of different pattern,
    !       diff - huge(0.0_myreal),
    !       else diff = sum(abs(A-B))
    real (kind = myreal) :: diff


    integer (kind = myint), intent (in) :: m,n,ia(*),ja(*),ib(*),jb(*)
    integer (kind = myint), intent (inout) :: mask(*)
    real (kind = myreal), intent (in) :: a(*),b(*)

    !
    integer (kind = myint) :: i,j,k,la1,la2,lb1,lb2,jj

    mask(1:n) = 0
    diff = 0
    do i = 1, m
       la1 = ia(i);la2 = ia(i+1)-1
       do j = la1,la2
          k = ja(j)
          if (mask(k) == 0) then
             mask(k) = j
          end if
       end do

       ! of the A and B is the same. Since neither
       ! is assumed to have repeated entries,
       ! its number of distinctive entries
       ! is simply ia(i+1)-ia(i)
       lb1 = ib(i);lb2 = ib(i+1)-1
       if (la2-la1 /= lb2-lb1) then
          diff = huge(0.0_myreal)
          return
       end if
       do j = lb1,lb2
          k = jb(j)
          jj = mask(k)
          if (jj == 0) then
             diff = huge(0.0_myreal)
             return
          else
             diff = diff + abs(b(j) - a(jj))
          end if
       end do
       do j = la1,la2
          mask(ja(j)) = 0
       end do
    end do
  end function csr_matrix_diff_values





  subroutine csr_matrix_multiply(matrix1,matrix2,result_matrix,info,stat)
    ! subroutine csr_matrix_multiply(matrix1,matrix2,&
    !     result_matrix,info)
    !
    ! multiply two matrices matrix1 and matrix2,
    ! and put the resulting matrix in result_matrix
    ! The storage is allocated inside this routine for
    ! the sparse matrix object {tt RESULT\_MATRIX},
    ! and when result_matrix is no longer needed,
    ! this memory should be deallocated using
    ! {\tt MATRIX\_DESTRUCT}


    ! matrix1: is of the derived type ZD11_type, INTENT (IN),
    !         the sparse matrix in compressed sparse row format
    ! matrix2: is of the derived type ZD11_type, INTENT (IN),
    !         the sparse matrix in compressed sparse row format
    type (zd11_type), intent (in) :: matrix1,matrix2

    ! result_matrix: of the derived type ZD11_type, INTENT (INOUT)
    !                result_matrix=matrix1*matrix2
    type (zd11_type), intent (inout) :: result_matrix

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed.
    !       = MC65_ERR_MATMUL_DIM_MISMATCH dimension mismatch
    integer (kind = myint), intent (out) :: info


    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! ================ local variables =================
    ! m: row size of the matrix1
    ! nz: number of nonzeros in the product matrix
    ! n: column side of matrix2
    integer (kind = myint) :: n,nz,m

    ! mask: an array with a size the same as the
    !       column size of matrix2,
    !       used to indicating whether an entry in the product
    !       has been seen or not.
    integer (kind = myint), allocatable, dimension(:) :: mask

    integer (kind = myint) :: ierr
    if (present(stat)) stat = 0

    info = 0

    m = matrix1%m
    n = matrix2%n

    ! check if the dimension match
    if (matrix1%n /= matrix2%m) then
       info = MC65_ERR_MATMUL_DIM_MISMATCH
       return
    end if

    allocate(mask(n), stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

    ! now find the size
    call matmul_size(m,n,matrix1%ptr,matrix1%col,matrix2%ptr,matrix2%col,&
         nz,mask)

    ! construct the new matrix
    ! if one of the matrix is pattern only, then
    ! the result matrix has to be pettern only
    if (ZD11_get(matrix1%type) == "pattern".or.&
         ZD11_get(matrix2%type) == "pattern") then
       call csr_matrix_construct(result_matrix,m,nz,info,n=n, &
            type = "pattern",stat=ierr)
       if (present(stat)) stat = ierr
       if (info < 0) return
    else if (ZD11_get(matrix1%type) == ZD11_get(matrix2%type)) then
       call csr_matrix_construct(result_matrix,m,nz,info,n=n, &
            type = ZD11_get(matrix1%type),stat=ierr)
       if (present(stat)) stat = ierr
       if (info < 0) return
    else
       call csr_matrix_construct(result_matrix,m,nz,info,n=n, &
            type = "general",stat=ierr)
       if (present(stat)) stat = ierr
       if (info < 0) return
    end if

    if (ZD11_get(result_matrix%type) == "pattern")  then
       call matmul_noval(m,n,matrix1%ptr,matrix1%col,matrix2%ptr,matrix2%col, &
                       result_matrix%ptr,result_matrix%col,mask)
    else
       call matmul_normal(m,n,matrix1%ptr,matrix1%col,matrix1%val, &
                       matrix2%ptr,matrix2%col,matrix2%val,result_matrix%ptr, &
                       result_matrix%col,result_matrix%val,mask)
    end if

    deallocate(mask, stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if

  end subroutine csr_matrix_multiply


  subroutine matmul_size(m,n,ia1,ja1,ia2,ja2,nz,mask)
    ! work out the number of nonzeros in the
    ! product of two matrix A=(ia1,ja1) and B=(ia2,ja2),

    ! m: row size of A
    ! n: column size of B
    integer (kind = myint), intent (in) :: n,m
    ! ia1: row pointer of A
    ! ja1: column indices of A
    ! ia2: row pointer of B
    ! ja2: column indices of B
    integer (kind = myint), intent (in), dimension (*) :: &
         ia1,ja1,ia2,ja2

    ! nz: number of nonzeros in the product matrix
    integer (kind = myint), intent (inout) :: nz

    ! mask: working array
    integer (kind = myint), intent (inout) :: mask(*)
    ! i,j,k,neigh,icol_add: working integers
    integer (kind = myint) :: i,j,k,neigh,icol_add

    ! initialise the mask array which is an array that has
    ! value i if column index already exist in row i of
    ! the new matrix
    mask(1:n) = 0
    nz = 0
    do i=1,m
       do j = ia1(i),ia1(i+1)-1
          neigh = ja1(j)
          do k = ia2(neigh),ia2(neigh+1)-1
             icol_add = ja2(k)
             if (mask(icol_add) /= i) then
                nz = nz + 1
                mask(icol_add) = i     ! add mask
             end if
          end do
       end do
    end do
  end subroutine matmul_size


  subroutine matmul_normal(m,n,ia,ja,a,ib,jb,b,ic,jc,c,mask)
    ! calculate the product of two matrix A=(ia,ja,a) and
    ! B = (ib,jb,b).

    ! m: row size of A
    ! n: column size of B
    integer (kind = myint), intent (in) :: n,m

    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    integer (kind = myint), intent (in) :: ia(*),ib(*),ja(*),jb(*)
    ! ic: row pointer of C
    ! jc: column indices of C
    integer (kind = myint), intent (out) :: ic(*),jc(*)
    ! a: entry values of A
    ! b: entry values of B
    ! c: entry values of C
    real (kind = myreal), intent (in) :: a(*),b(*)
    real (kind = myreal), intent (out) :: c(*)
    ! working array
    integer (kind = myint), intent (inout) :: mask(*)
    ! working variables
    integer (kind = myint) :: nz,i,j,k,icol,icol_add,neigh
    real (kind = myreal) :: aij

    ! initialise the mask array which is an array that
    ! has non-zero value if the
    ! column index already exist, in which case
    ! the value is the index of that column
    mask(1:n) = 0
    nz = 0
    ic(1) = 1
    do i=1,m
       do j = ia(i),ia(i+1)-1
          aij = a(j)
          neigh = ja(j)
          do k = ib(neigh),ib(neigh+1)-1
             icol_add = jb(k)
             icol = mask(icol_add)
             if (icol == 0) then
                nz = nz + 1
                jc(nz) = icol_add
                c(nz) =  aij*b(k)
                mask(icol_add) = nz     ! add mask
             else
                c(icol) = c(icol) + aij*b(k)
             end if
          end do
       end do
       do j = ic(i),nz
          ! done this row i, so set mask to zero again
          mask(jc(j)) = 0
       end do
       ic(i+1) = nz+1
    end do
  end subroutine matmul_normal




  subroutine matmul_noval(m,n,ia,ja,ib,jb,ic,jc,mask)
    ! calculate the product of two matrix A=(ia,ja) and
    ! B = (ib,jb). Sparse patterns only is calculated.

    ! m: row size of A
    ! n: column size of B
    integer (kind = myint), intent (in) :: n,m
    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    integer (kind = myint), intent (in) :: ia(*),ib(*),ja(*),jb(*)
    ! ic: row pointer of C
    ! jc: column indices of C
    integer (kind = myint), intent (out) :: ic(*),jc(*)

    ! working array
    integer (kind = myint), intent (inout) :: mask(*)
    ! working variables
    integer (kind = myint) :: nz,i,j,k,icol,icol_add,neigh

    ! initialise the mask array which is an array
    ! that has non-zero value if the
    ! column index already exist, in which case
    ! the value is the index of that column
    mask(1:n) = 0
    nz = 0
    ic(1) = 1
    do i=1,m
       do j = ia(i),ia(i+1)-1
          neigh = ja(j)
          do k = ib(neigh),ib(neigh+1)-1
             icol_add = jb(k)
             icol = mask(icol_add)
             if (icol == 0) then
                nz = nz + 1
                jc(nz) = icol_add
                mask(icol_add) = nz     ! add mask
             else
                cycle
             end if
          end do
       end do
       do j = ic(i),nz
          ! done this row i, so set mask to zero again
          mask(jc(j)) = 0
       end do
       ic(i+1) = nz+1
    end do
  end subroutine matmul_noval

  subroutine csr_matrix_multiply_graph(matrix1,matrix2,&
       result_matrix,info,col_wgt,stat)
    ! subroutine csr_matrix_multiply_graph(matrix1,matrix2,&
    !    result_matrix,info[,col_wgt])
    !
    ! Given two matrix, find
    ! result_matrix = (\bar matrix1)*(\bar matrix2) where
    ! (\bar A) is the matrix A with all non-zero entries set to 1.
    ! The diagonal entries of the result_matrix will be discarded.
    ! the routine is used for finding the row connectivity graph
    ! A, which is (\bar A)*(\bar A)^T, with diagonal elements removed

    ! matrix1: of the derived type ZD11_type, INTENT (IN),
    !         the sparse matrix in compressed sparse row format
    ! matrix2: of the derived type ZD11_type, INTENT (IN),
    !         the sparse matrix in compressed sparse row format
    type (zd11_type), intent (in) :: matrix1,matrix2

    ! result_matrix: of the derived type ZD11_type, INTENT (INOUT)
    !                result_matrix=matrix1*matrix2
    type (zd11_type), intent (inout) :: result_matrix

    ! col_wgt: optional real array of INTENT (IN).
    !          the super-column weight for the columns of
    !          matrix1. this is used when calculating A*A^T
    !          with A having super-columns The effect is the same as
    !          result_matrix =
    !          (\bar matrix1)*DIAG(col_wgt)*(\bar matrix2)
    real (kind = myreal), dimension (matrix1%n), intent (in), &
         optional :: col_wgt

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed.
    !       = MC65_ERR_MATMULG_DIM_MISMATCH dimension mismatch
    integer (kind = myint), intent (out), optional :: info

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! ================= local variables =================
    ! m: row size of matrix1
    ! n: column size of matrix2
    ! nz: number of nonzeros in the product matrix
    integer (kind = myint) :: n,nz,m

    ! mask: an array with a size the same as the column
    !       size of matrix2, used to indicating whether an
    !       entry in the product has been seen or not.
    integer (kind = myint), allocatable, dimension(:) :: mask

    integer (kind = myint) :: ierr

    info = 0
    if (present(stat)) stat = 0

    m = matrix1%m

    n = matrix2%n

    ! Check if dimensions match
    if ( matrix1%n /= matrix2%m) then
       info = MC65_ERR_MATMULG_DIM_MISMATCH
       return
    end if

    allocate(mask(n), stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
    end if

    ! now find the size
    call matmul_size_graph(m,n,matrix1%ptr,matrix1%col,matrix2%ptr,&
         matrix2%col,nz,mask)

    ! initialise the new matrix
    call csr_matrix_construct(result_matrix,m,nz,info,n = n,&
         type = "general",stat=ierr)
    if (present(stat)) stat = ierr
    if (info < 0) return

    if (present(col_wgt)) then
       call matmul_wgraph(m,n,matrix1%ptr,matrix1%col,matrix2%ptr, &
                  matrix2%col,result_matrix%ptr,result_matrix%col, &
                  result_matrix%val,mask,col_wgt)
    else
       call matmul_graph(m,n,matrix1%ptr,matrix1%col,matrix2%ptr, &
                  matrix2%col,result_matrix%ptr,result_matrix%col, &
                  result_matrix%val,mask)
    end if

    deallocate(mask, stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if

  end subroutine csr_matrix_multiply_graph


  subroutine matmul_size_graph(m,n,ia1,ja1,ia2,ja2,nz,mask)
    ! work out the number of nonzeros in the
    ! product of two matrix A=(ia1,ja1) and B=(ia2,ja2),
    ! diagonal element of A*B is not counted
    ! because this routine is used in finding row graphs.
    ! (a graph is represented as a symm. matrix with no diag. entries)
    ! m: row size of A
    ! n: column size of B
    integer (kind = myint), intent (in) :: n,m
    ! ia1: row pointer of A
    ! ja1: column indices of A
    ! ia2: row pointer of B
    ! ja2: column indices of B
    integer (kind = myint), intent (in), dimension (*) :: &
         ia1,ja1,ia2,ja2

    ! nz: number of nonzeros in the product matrix
    integer (kind = myint), intent (inout) :: nz

    ! mask: working array
    integer (kind = myint), intent (inout) :: mask(*)

    ! i,j,k,neigh,icol_add: working integers
    integer (kind = myint) :: i,j,k,neigh,icol_add

    ! initialise the mask array which is an array that has
    !  value i if column index already exist in row i
    !  of the new matrix
    mask(1:n) = 0
    nz = 0
    do i=1,m
       do j = ia1(i),ia1(i+1)-1
          neigh = ja1(j)
          do k = ia2(neigh),ia2(neigh+1)-1
             icol_add = ja2(k)
             if (icol_add /= i .and. mask(icol_add) /= i) then
                nz = nz + 1
                mask(icol_add) = i     ! add mask
             end if
          end do
       end do
    end do
  end subroutine matmul_size_graph

  subroutine matmul_graph(m,n,ia,ja,ib,jb,ic,jc,c,mask)
    ! given two matrix, A=(ia,ja) and B = (ib,jb), find
    ! C = (\bar A)*(\bar B) where
    ! (\bar A) is the matrix A with all non-zero entries set to 1.
    ! The diagonal entries of the c matrix will be discarded.
    ! the routine is mainly used for finding the row
    ! connectivity graph
    ! of a matrix A, which is (\bar A)*(\bar A)^T
    ! it is assumed that the size of arrays (ic,jc,c)
    ! used to hold C is large enough.

    ! m: row size of A
    ! n: column size of B
    integer (kind = myint), intent (in) :: n,m
    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    integer (kind = myint), intent (in), dimension (*) :: &
         ia,ib,ja,jb
    ! ic: row pointer of C
    ! jc: column indices of C
    integer (kind = myint), intent (out), dimension (*)  :: ic,jc
    ! c: entry values of C
    real (kind = myreal), intent (out), dimension (*)  :: c
    ! working array
    integer (kind = myint), intent (inout) :: mask(*)
    ! working variables
    integer (kind = myint) :: nz,i,j,k,icol,icol_add,neigh

    ! initialise the mask array which is an array
    ! that has non-zero value if the
    ! column index already exist, in which case
    ! the value is the index of that column
    mask(1:n) = 0
    nz = 0
    ic(1) = 1
    do i=1,m
       do j = ia(i),ia(i+1)-1
          neigh = ja(j)
          do k = ib(neigh),ib(neigh+1)-1
             icol_add = jb(k)
             if (icol_add /= i) then ! do not consider diagonal
                icol = mask(icol_add)
                if (icol == 0) then
                   nz = nz + 1
                   jc(nz) = icol_add
                   c(nz) =  1.0_myreal
                   mask(icol_add) = nz     ! add mask
                else
                   c(icol) = c(icol) + 1.0_myreal
                end if
             end if
          end do
       end do
       do j = ic(i),nz
          ! done this row i, so set mask to zero again
          mask(jc(j)) = 0
       end do
       ic(i+1) = nz+1
    end do
  end subroutine matmul_graph


  subroutine matmul_wgraph(m,n,ia,ja,ib,jb,ic,jc,c,mask,col_wgt)
    ! given two matrix, find
    ! C = (\bar A)*W*(\bar B) where
    ! (\bar A) is the matrix A, with all non-zero entries set to 1.
    ! W is a diagonal matrix consists of the column weights of A
    ! The diagonal entries of the c matrix will be discarded.
    ! the routine is mainly used for finding the
    ! row connectivity graph of
    ! a matrix A with column weights, which is (\bar A)*W*(\bar A)^T
    ! it is assumed that the size of arrays (ic,jc,c)
    ! used to hold C is large enough.

    ! m: row size of A
    ! n: column size of B
    integer (kind = myint), intent (in) :: n,m

    ! ia: row pointer of A
    ! ja: column indices of A
    ! ib: row pointer of B
    ! jb: column indices of B
    integer (kind = myint), intent (in), dimension (*) :: &
         ia,ib,ja,jb
    ! col_wgt: the column weight, or the diagonal matrix W
    real(kind = myreal), intent (in), dimension (*) :: col_wgt
    ! ic: row pointer of C
    ! jc: column indices of C
    integer (kind = myint), intent (out), dimension (*)  :: ic,jc
    ! c: entry values of C
    real (kind = myreal), intent (out), dimension (*)  :: c
    ! working array
    integer (kind = myint), intent (inout) :: mask(*)
    ! working variables
    integer (kind = myint) :: nz,i,j,k,icol,icol_add,neigh

    ! initialise the mask array which is an array
    ! that has non-zero value if the
    ! column index already exist, in which case
    ! the value is the index of that column
    mask(1:n) = 0
    nz = 0
    ic(1) = 1
    do i=1,m
       do j = ia(i),ia(i+1)-1
          neigh = ja(j)
          do k = ib(neigh),ib(neigh+1)-1
             icol_add = jb(k)
             if (icol_add /= i) then ! do not consider diagonal
                icol = mask(icol_add)
                if (icol == 0) then
                   nz = nz + 1
                   jc(nz) = icol_add
                   c(nz) =  col_wgt(neigh)
                   mask(icol_add) = nz     ! add mask
                else
                   c(icol) = c(icol) + col_wgt(neigh)
                end if
             end if
          end do
       end do
       do j = ic(i),nz
          ! done this row i, so set mask to zero again
          mask(jc(j)) = 0
       end do
       ic(i+1) = nz+1
    end do
  end subroutine matmul_wgraph



  subroutine csr_matrix_multiply_rvector(matrix,x,y,info,trans)
    ! subroutine csr_matrix_multiply_rvector(matrix,x,y,info)
    !
    ! y = matrix*x where x and y are real vectors. Dimension of y
    ! is checked and returned if it is smaller than the row dimension
    ! of x
    !
    ! if trans is present and == 1;
    ! y = matrix^T*x where x and y are real vectors.


    ! matrix: of the derived type ZD11_type, INTENT (IN),
    !         the sparse matrix in compressed sparse row format
    type (zd11_type), intent (in) :: matrix

    ! x: real array of intent (IN), a vector to be multiplied
    !    with the matrix
    real (kind = myreal), intent (in), dimension (*) :: x

    ! y: real array of intent (OUT), the result of matrix*x
    !    or matrix^T*x
    real (kind = myreal), intent (out), dimension (*) :: y

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MATVEC_NOVALUE, &  if A is of type
    !         pattern only and x real.
    integer (kind = myint), intent (out) :: info

    ! trans: optional logical scalar. If present and trans = true,
    ! the transpose of the matrix is multiplied with the vector
    logical, optional, intent (in) :: trans
    ! ========================local variables=======================
    integer (kind = myint) :: m,i,l1,l2,j,n,jc
    real (kind = myreal) :: xx
    ! transpose: whether it is a matrix transpose multiplying a vector
    logical :: transpose

    transpose = .false.

    if (present(trans)) then
       if (trans) transpose = .true.
    end if

    info = 0

    m = matrix%m; n = matrix%n

    if (ZD11_get(matrix%type) == "pattern") then
       info = MC65_ERR_MATVEC_NOVALUE
       return
    end if
    if (transpose) then
       y(1:n) = 0
       do i = 1, m
          xx = x(i)
          do j = matrix%ptr(i),matrix%ptr(i+1)-1
             jc = matrix%col(j)
             y(jc) = y(jc) + matrix%val(j)*xx
          end do
       end do
    else
       do i = 1, m
          l1 = matrix%ptr(i); l2 = matrix%ptr(i+1)-1
          y(i) = dot_product(matrix%val(l1:l2),x(matrix%col(l1:l2)))
       end do
    end if

  end subroutine csr_matrix_multiply_rvector


  subroutine csr_matrix_multiply_ivector(matrix,x,y,info,trans)
    ! subroutine csr_matrix_multiply_ivector(matrix,x,y,info)
    !
    ! y = matrix*x where x and y are integer vectors. Entries of
    ! matrix is assumed to be one. Dimension of y
    ! is checked and returned if it is smaller than the row dimension
    ! of x
    !
    ! if trans is present and == 1;
    ! y = matrix^T*x where x and y are integers and entries of
    ! the matrix are assumed to have the values one.

    ! matrix: of the derived type ZD11_type, INTENT (IN),
    !         the sparse matrix in compressed sparse row format
    type (zd11_type), intent (in) :: matrix

    ! x: integer array of intent (IN), a vector to be
    !    multiplied with the matrix
    integer (kind = myint), intent (in), dimension (*) :: x

    ! y: integer array of intent (OUT), the result of
    !    matrix*x or matrix^T*x
    integer (kind = myint), intent (out), dimension (*) :: y

    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    integer (kind = myint), intent (out) :: info

    ! trans: optional integer scalar. If present and trans = true,
    ! transpose of the matrix is multiplied with the vector
    logical, optional, intent (in) :: trans

    ! local ==========
    integer (kind = myint) :: m,n,i,l1,l2,j,xx,jc
    ! transpose: whether it is a matrix transpose
    !   multiplying a vector
    logical :: transpose


    info = 0

    transpose = .false.

    if (present(trans)) then
       if (trans) transpose = .true.
    end if

    m = matrix%m; n = matrix%n

    if (transpose) then
       y(1:n) = 0
       do i = 1, m
          xx = x(i)
          do j = matrix%ptr(i),matrix%ptr(i+1)-1
             jc = matrix%col(j)
             y(jc) = y(jc) + xx
          end do
       end do
    else
       do i = 1, m
          l1 = matrix%ptr(i); l2 = matrix%ptr(i+1)-1
          y(i) = sum(x(matrix%col(l1:l2)))
       end do
    end if
  end subroutine csr_matrix_multiply_ivector

  subroutine csr_to_csr_matrix(matrix,m,n,ptr,col,&
       info,val,type,checking,stat)
    !   subroutine csr_to_csr_matrix(matrix,m,n,ptr,col,
    !      info[,val,type,copying,checking])
    ! convert ptr, col and a arrays into the compressed sparse row
    ! matrix format. New storage is allocated for MATRIX and
    ! the content of ptr,col (and val) is copied into MATRIX.

    ! m: is an integer of intent (in).
    !    It is the row dimension of the sparse matrix
    ! n: is an integer of intent (in).
    !    It is the column dimension of the sparse matrix
    integer (kind = myint), intent (in) :: m,n

    ! MATRIX: is of the derived type {\tt ZD11\_type}
    !     with {\tt INTENT (INOUT)},
    type (zd11_type), intent (inout) :: matrix

    ! ptr: is an integer array of intent (in) of size M+1.
    !     It is the row pointer of the sparse matrix.
    ! col: is an integer array of intent (in) of size ptr(m+1)-1.
    !     It is the array of column indices.
    integer (kind = myint), INTENT (IN)  :: ptr(m+1), col(ptr(m+1)-1)


    ! info: integer scaler of INTENT (OUT).
    !       = 0 if successful
    !       = MC65_ERR_MEMORY_ALLOC if memory allocation fails.
    !       = MC65_ERR_MEMORY_deALLOC if memory deallocation fails.
    !       = MC65_ERR_RANGE_PTR if PTR(1:M+1) is not
    !         monotonically increasing
    !       = MC65_WARN_RANGE_COL: COL is not within the
    !         range of [1,n], such entries are removed
    !       = MC65_WARN_DUP_ENTRY: duplicate entries were
    !         found when cleaning the matrix and were summed
    integer (kind = myint), intent (out) :: info

    ! val: is an optional REAL array of intent (in) of
    !    size ptr(m+1)-1.
    !    It holds the entry values of the sparse matrix.
    !    If this argument is not present, the matrix is
    !    taken as pattern only.
    real (kind = myreal),  INTENT (IN), optional :: val(ptr(m+1)-1)


    ! type: an optional character array of intent (IN).
    !       If both type and A is present,
    !       the type of the MATRIX will be set as TYPE.
    !
    character (len = *), INTENT (IN), OPTIONAL :: type


    ! checking: is an optional integer scaler of INTENT (IN).
    !           If present and is
    !           1: the input data will be checked for consistency,
    !              the matrix will be cleaned by merging duplicate
    !              entries. If a column index is out of range,
    !              the entry will be excluded in the construction
    !              of the matrix.
    !              A warning will be recorded.
    !           < 0: no checking or cleaning is done
    !           Other values: current not used and has the
    !           same effect 1
    integer (kind = myint), parameter :: CHECK_RMV = 1
    integer (kind = myint), optional, INTENT (IN) :: checking

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat


    integer (kind = myint) :: nz,nzz,i,j,icol,ierr


    ! whether the matrix is pattern only
    logical :: pattern_only

    ! whether to check input data
    integer (kind = myint) :: to_check

    nzz = ptr(m+1)-1

    info = 0
    if (present(stat)) stat = 0
    to_check = 0
    if (present(checking)) then
       if (checking == CHECK_RMV) then
          to_check = CHECK_RMV
       else if (checking >= 0) then
          to_check = CHECK_RMV
       end if
    end if


    pattern_only = .true.
    if (present(val)) then
       pattern_only = .false.
       if (present(type)) then
          if (len_trim(type)>=7) then
             if (type(1:7) == "pattern".or.type(1:7) == "PATTERN") &
                  pattern_only = .true.
          end if
       end if
    end if

    if (to_check == CHECK_RMV) then
       do i = 1,m
          if (ptr(i) > ptr(i+1)) then
             info = MC65_ERR_RANGE_PTR
             return
          end if
       end do
       nz = 0
       do i = 1, nzz
          if (col(i) <= n.and.col(i) >= 1) then
             nz = nz + 1
          end if
       end do
    else
       nz = nzz
    end if


    if (pattern_only) then
       call csr_matrix_construct(matrix,m,nz,info,n = n,&
             type = "pattern",stat=ierr)
    else
       if (present(type)) then
          call csr_matrix_construct(matrix,m,nz,info,n = n,&
               type = type,stat=ierr)
       else
          call csr_matrix_construct(matrix,m,nz,info,n = n,&
               type = "general",stat=ierr)
       end if
    end if
    if (present(stat)) stat = ierr
    if (info < 0) return

    ! column index outside the range
    matrix%ptr(1) = 1
    if (nz /= nzz.and.to_check == CHECK_RMV) then
       nz = 1
       if (.not.pattern_only) then
          do i = 1, m
             do j = ptr(i), ptr(i+1)-1
                icol = col(j)
                if (icol >= 1.and.icol <= n) then
                   matrix%col(nz) = icol
                   matrix%val(nz) = val(j)
                   nz = nz + 1
                end if
             end do
             matrix%ptr(i+1) = nz
          end do
       else
          do i = 1, m
             do j = ptr(i), ptr(i+1)-1
                icol = col(j)
                if (icol >= 1.and.icol <= n) then
                   matrix%col(nz) = icol
                   nz = nz + 1
                end if
             end do
             matrix%ptr(i+1) = nz
          end do
       end if
       call csr_matrix_clean(matrix,info,stat=ierr)
       if (present(stat)) stat = ierr
       if (info < 0) return
       ! only set info for COL out of range, since info is reset when
       ! calling csr_matrix_construct or csr_matrix_clean
       if (info == MC65_WARN_DUP_ENTRY) then
          info = MC65_WARN_ENTRIES
       else
          info =  MC65_WARN_RANGE_COL
       end if
    else if (to_check == CHECK_RMV) then
       ! column index not outside the range but cleaning is needed
       matrix%ptr(1:m+1) = ptr(1:m+1)
       matrix%col(1:nz) = col(1:nz)
       if (.not.pattern_only) matrix%val(1:nz) = val(1:nz)
       call csr_matrix_clean(matrix,info,stat=ierr)
       if (present(stat)) stat = ierr
    else
       matrix%ptr(1:M+1) = ptr(1:m+1)
       matrix%col(1:nz) = col(1:nz)
       if (.not.pattern_only) matrix%val(1:nz) = val(1:nz)
    end if


  end subroutine csr_to_csr_matrix


  subroutine coo_to_csr_format(matrix,m,n,nz,irn,jcn,&
       info,val,type,checking,stat)
    !
    ! convert coordinate format (irn,jcn,val) to CSR matrix format

    ! matrix:  is of the derived type {\tt ZD11\_type}
    !          with {\tt INTENT (INOUT)}. It is
    !          the sparse matrix to be generated from the
    !          coordinated formatted
    !          input data.
    ! m: is an INTEGER scaler of INTENT (IN). It holds the row size
    ! n: is an INTEGER scaler of INTENT (IN). It holds the
    !    column size
    ! nz:  is an INTEGER scaler. It holds the number of
    !      nonzeros in the sparse matrix
    ! irn: is an INTEGER ARRAY with INTENT (IN), of size = NZ.
    !      It holds the row indices of the matrix
    ! jcn:  is an INTEGER ARRAY with INTENT (IN), of size = NZ.
    !       It holds column indices of the matrix
    ! info:  is an {\tt INTEGER} scaler of {\tt INTENT (OUT)}.
    !       INFO = 0 if the subroutine completes successful
    !       INFO = MC65_ERR_MEMORY_ALLOC if memory allocation failed.
    !       INFO = MC65_ERR_MEMORY_DEALLOC if memory
    !              deallocation failed.
    !       {\tt INFO = MC65\_WARN\_RANGE\_IRN} if {\tt IRN}
    !           is not within the
    !           range of {\tt [1,M]}. Related entries excluded
    !       {\tt INFO = MC65\_WARN\_RANGE\_JCN} if {\tt JCN}
    !           is not within the
    !           range of {\tt [1,N]}. Related entries excluded
    !       INDO = MC65_WARN_DUP_ENTRY: duplicate entries were
    !              found when cleaning the matrix
    !              and were summed
    ! val:  is an optional REAL (double precision REAL for
    !       HSL_MC65_double)
    !       ARRAY  with INTENT (IN), of size = NZ.
    !       It holds entry values of the matrix
    integer (kind=myint), intent (in) :: m, n, nz
    type (zd11_type), intent (inout) :: matrix
    integer (kind=myint), intent (in) :: irn(nz)
    integer (kind=myint), intent (in) :: jcn(nz)
    integer (kind=myint), intent (out) :: info
    real (kind = myreal), intent (in), optional :: val(nz)

    ! type: is an OPTIONAL CHARACTER array of unspecified
    !       size with INTENT (IN).
    !       If both val and type are present,
    !       the type of the MATRIX will be set as type.
    character (len=*), intent (in), optional :: type

    ! checking: is an optional integer scaler of INTENT (IN).
    !           If present and is
    !           1: the input data will be checked for consistency,
    !              the matrix will be cleaned by merging duplicate
    !              entries.
    !              If a column index is out of range, the entry will
    !              be excluded in the construction of the matrix.
    !              A warning will be recorded.
    !           < 0: no checking or cleaning is done
    !              Other values: current not used and has the same
    !              effect as 1
    integer (kind = myint), parameter :: CHECK_RMV = 1
    integer (kind = myint), optional, INTENT (IN) :: checking

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! nzz: actual number of nonzero

    integer (kind=myint) :: ierr,nzz,i
    integer (kind = myint), dimension (:), allocatable :: nrow,iia

    ! whether the matrix is pattern only
    logical :: pattern_only

    ! whether to check input data. 0 not, /= 0 yes
    integer (kind = myint) :: to_check

    logical irn_out_range,jcn_out_range

    irn_out_range = .false.
    jcn_out_range = .false.

    info = 0
    if (present(stat)) stat = 0
    to_check = 0
    if (present(checking)) then
       if (checking == CHECK_RMV) then
          to_check = CHECK_RMV
       else if (checking >= 0) then
          to_check = CHECK_RMV
       end if
    end if

    pattern_only = .true.
    if (present(val)) then
       pattern_only = .false.
       if (present(type)) then
          if (len_trim(type)>=7) then
             if (type(1:7) == "pattern".or.type(1:7) == "PATTERN") &
                  pattern_only = .true.
          end if
       end if
    end if

    if (to_check == CHECK_RMV) then
       nzz = 0
       do i = 1, nz
          if (jcn(i) <= n.and.jcn(i) >= 1) then
             if (irn(i) <= m.and.irn(i) >= 1) then
                nzz = nzz + 1
             else
                irn_out_range = .true.
             end if
          else
             jcn_out_range = .true.
          end if
       end do
    else
       nzz = nz
    end if

    if (pattern_only) then
       call csr_matrix_construct(matrix,m,nzz,info,n = n,&
            type = "pattern",stat=ierr)
    else
       if (present(type)) then
          call csr_matrix_construct(matrix,m,nzz,info,n = n,&
               type = type,stat=ierr)
       else
          call csr_matrix_construct(matrix,m,nzz,info,n = n,&
               type = "general",stat=ierr)
       end if
    end if
    if (present(stat)) stat = ierr
    if (info < 0) return

    allocate(nrow(m),iia(m+1), stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_ALLOC
       return
    end if

    if (pattern_only) then
       call coo_to_csr_private(to_check,m,n,nz,irn,jcn,matrix%ptr,matrix%col,&
            nrow,iia,pattern_only)
    else
       call coo_to_csr_private(to_check,m,n,nz,irn,jcn,matrix%ptr,matrix%col,&
            nrow,iia,pattern_only,val = val,a = matrix%val)
    end if

    deallocate(nrow, iia, stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       info = MC65_ERR_MEMORY_DEALLOC
       return
    end if

    if (to_check == CHECK_RMV) then
       call csr_matrix_clean(matrix,info,stat=ierr)
       if (present(stat)) stat = ierr
       if (info < 0) return
       if (irn_out_range .and. info==MC65_WARN_DUP_ENTRY) then
          info = MC65_WARN_ENTRIES
          return
       end if
       if (jcn_out_range .and. info==MC65_WARN_DUP_ENTRY) then
          info = MC65_WARN_ENTRIES
          return
       end if
       if (irn_out_range) info = MC65_WARN_RANGE_IRN
       if (jcn_out_range) info = MC65_WARN_RANGE_JCN
       if (irn_out_range .and. jcn_out_range) info = MC65_WARN_RANGE_BOTH
    end if

  end subroutine coo_to_csr_format


  subroutine coo_to_csr_private(to_check,m,n,nz,irn,jcn,ia,ja,&
       nrow,iia,pattern_only,val,a)
    ! wrapper routine for coo_to_csr_format for performance
    ! to_check: whether to check input data. 0 not, /= 0 yes
    ! m,n: row/column dimension
    ! ia: row pointer
    ! ja: column pointer
    ! nrow: number of entries in a row
    ! iia: row pointer
    ! a: entry values
    integer (kind = myint) :: to_check
    integer (kind = myint), intent (in) :: n,m,nz,irn(*),jcn(*)
    integer (kind = myint), intent (out) :: ia(*),ja(*),nrow(*),iia(*)
    real (kind = myreal), intent (OUT), optional :: a(*)
    real (kind = myreal), intent (IN), optional :: val(*)
    logical :: pattern_only

    integer (kind=myint) :: i
    integer (kind = myint) :: ii,jj

    if (to_check == 0) then
       nrow(1:m) = 0
       do i = 1, nz
          nrow(irn(i)) = nrow(irn(i)) + 1
       end do
       ia (1) = 1
       do i = 1, m
          ia(i+1) = ia(i) + nrow(i)
       end do

       iia(1:m+1) = ia(1:m+1)
       if (pattern_only) then
          do i = 1, nz
             ii = irn(i)
             jj = jcn (i)
             ja(iia(ii)) = jj
             iia(ii) = iia(ii) + 1
          end do
       else
          do i = 1, nz
             ii = irn(i)
             jj = jcn (i)
             ja(iia(ii)) = jj
             a(iia(ii)) = val(i)
             iia(ii) = iia(ii) + 1
          end do
       end if
    else
       ! checking range, exclude out of range entries
       nrow(1:m) = 0
       do i = 1, nz
          ii = irn(i)
          jj = jcn(i)
          if (ii >= 1.and.ii <= m.and.jj >= 1.and.jj <= n) &
               nrow(ii) = nrow(ii) + 1
       end do
       ia (1) = 1
       do i = 1, m
          ia(i+1) = ia(i) + nrow(i)
       end do

       iia(1:m+1) = ia(1:m+1)
       if (pattern_only) then
          do i = 1, nz
             ii = irn(i)
             jj = jcn (i)
             if (ii >= 1.and.ii <= m.and.jj >= 1.and.jj <= n) then
                ja(iia(ii)) = jj
                iia(ii) = iia(ii) + 1
             end if
          end do
       else
          do i = 1, nz
             ii = irn(i)
             jj = jcn (i)
             if (ii >= 1.and.ii <= m.and.jj >= 1.and.jj <= n) then
                ja(iia(ii)) = jj
                a(iia(ii)) = val(i)
                iia(ii) = iia(ii) + 1
             end if
          end do
       end if
    end if
  end subroutine coo_to_csr_private


  subroutine csr_matrix_to_coo(matrix,m,n,nz,irn,jcn,info,val,stat)
    ! convert the ZD11_type matrix format to coordinate form.
    ! Storage space for the array pointers  {\tt IRN},
    ! {\tt JCN} and {\tt VAL}
    ! are allocated inside the subroutine, and should be deallocated
    ! by the user once these array pointers are no longer used.
    ! being used. When the matrix is of type "pattern"
    ! and VAL is present, VAL will be allocated a size of 0.

    ! matrix:  is of the derived type {\tt ZD11\_TYPE} with
    !   {\tt INTENT (IN)}.
    type (zd11_type), intent (in) :: MATRIX

    ! m: is an INTEGER scaler of INTENT (OUT). It is the row
    !    dimension
    ! n: is an INTEGER scaler of INTENT (OUT). It is the column
    !    dimension
    ! nz: is an INTEGER scaler of INTENT (OUT). It is the number
    !    of nonzeros
    integer (kind = myint), intent (out) :: m,n,nz

    ! irn: is an INTEGER allocatable array. It holds the row indices
    ! jcn:  is an INTEGER allocatable array. It holds the column indices
    ! val:  is an OPTIONAL INTEGER allocatable array. When present,
    !       and the matrix is not of type "pattern", VAL holds
    !       the entry values
    !       when present and the matrix is of type "pattern"
    !       VAL will be allocated a size of 0.
    integer (kind = myint), allocatable :: irn(:)
    integer (kind = myint), allocatable :: jcn(:)
    real (kind = myreal), optional, allocatable :: val(:)

    ! info:  is an {\tt INTEGER} scaler of {\tt INTENT (OUT)}.
    !       INFO = 0 if the subroutine completes successful
    !       INFO = MC65_ERR_MEMORY_ALLOC if memory allocation failed.
    integer (kind = myint), intent (out) :: info

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    ! working integers
    integer (kind = myint) :: ierr

    ! pattern_only: whether the matrix is pattern only
    ! entry_values: whether val array is to be allocated and filled
    logical :: pattern_only,entry_values

    info = 0
    if (present(stat)) stat = 0
    pattern_only = .false.
    entry_values = .false.
    if (ZD11_get(matrix%type) == "pattern") pattern_only = .true.
    if ((.not.pattern_only).and.present(val)) entry_values = .true.

    m = matrix%m; n = matrix%n

    nz = matrix%ptr(m+1)-1

    if (nz < 0) return
    if (entry_values) then
       allocate (irn(1:nz),jcn(1:nz),val(1:nz),stat = ierr)
    else
       allocate (irn(1:nz),jcn(1:nz),stat = ierr)
       IF (PRESENT(VAL)) allocate(val(0),stat = ierr)
    end if
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       INFO = MC65_ERR_MEMORY_ALLOC
       return
    end if

    if (entry_values) then
       call csr_matrix_to_coo_private(m,nz,matrix%ptr,matrix%col,irn,jcn, &
                                      matrix%val,val)
    else
       call csr_matrix_to_coo_private(m,nz,matrix%ptr,matrix%col,irn,jcn)
    end if

  end subroutine csr_matrix_to_coo

  subroutine csr_matrix_to_coo_private(m,nz,ia,ja,irn,jcn,a,val)
    ! private wrapper routine called by csr_matrix_to_coo
    ! for efficient
    !
    ! m: row dimension
    ! nz: number of nonzeros
    ! ia: row pointer
    ! ja: column indices
    ! a: entry values, optional
    ! irn: row indices
    ! jcn: column indices
    ! val: entry values, optional
    integer (kind = myint), intent (in) :: m
    integer (kind = myint), intent (out) :: nz
    integer (kind = myint), dimension (*), intent (in) :: ia,ja
    real (kind = myreal), dimension (*), intent (in), optional :: a
    integer (kind = myint), dimension (*), intent (out) :: irn,jcn
    real (kind = myreal), dimension (*), intent (out), optional :: val

    integer (kind = myint) :: i,n1,i1,i2


    nz = 1
    if (present(a)) then
       do i=1,m
          i1 = ia(i); i2 = ia(i+1)-1
          n1 = i2-i1+1
          irn(nz:nz+n1-1) = i
          jcn(nz:nz+n1-1) = ja(i1:i2)
          val(nz:nz+n1-1) = a(i1:i2)
          nz = nz + n1
       end do
    else
       do i=1,m
          i1 = ia(i); i2 = ia(i+1)-1
          n1 = i2-i1+1
          irn(nz:nz+n1-1) = i
          jcn(nz:nz+n1-1) = ja(i1:i2)
          nz = nz + n1
       end do
    end if
    nz = nz -1
  end subroutine csr_matrix_to_coo_private

  subroutine csr_matrix_remove_diagonal(matrix)
    ! remove the diagonal entry of a matrix. This is used when
    ! the matrix represents a graph and there should not be
    ! diagonals in such a matrix

    ! matrix: is of the derived type {\tt ZD11\_TYPE} with
    !         {\tt INTENT (INOUT)}.
    !         It is the matrix whose diagonal is to be removed
    type (zd11_type), intent (inout) :: matrix

    ! m: row dimension
    integer (kind = myint) :: m

    m = matrix%m
    if (ZD11_get(matrix%type) /= "pattern" ) then
       call csr_matrix_remove_diag_private(m,matrix%ptr,matrix%col,matrix%val)
    else
       call csr_matrix_remove_diag_private(m,matrix%ptr,matrix%col)
    end if
  end subroutine csr_matrix_remove_diagonal

  subroutine csr_matrix_remove_diag_private(m,ia,ja,a)
    ! private  subroutine called by csr_matrix_remove_diag
    ! for performance reason
    !
    ! m: row dimension
    ! ia: row pointer
    ! ja: column indices
    ! a: entry values
    integer (kind = myint), intent (in) :: m
    integer (kind = myint), intent (inout) :: ia(*),ja(*)
    real (kind = myreal), optional, intent (inout) :: a(*)

    integer (kind = myint) :: i,j,l1,nz


    nz = 1
    l1 = 1
    if (present(a)) then
       do i = 1, m
          do j = l1, ia(i+1) - 1
             if (ja(j) /= i) then
                ja(nz) = ja(j)
                a(nz) = a(j)
                nz = nz + 1
             end if
          end do
          l1 = ia(i+1)
          ia(i+1) = nz
       end do
    else
       do i = 1, m
          do j = l1, ia(i+1) - 1
             if (ja(j) /= i) then
                ja(nz) = ja(j)
                nz = nz + 1
             end if
          end do
          l1 = ia(i+1)
          ia(i+1) = nz
       end do
    end if
  end subroutine csr_matrix_remove_diag_private


  subroutine csr_matrix_diagonal_first(matrix,info)
    ! subroutine csr_matrix_diagonal_first(matrix,info)
    !
    ! changes the storage arrange so that the diagonal of
    ! row i of the matrix is always
    ! at position matrix%ptr(i) of the matrix%col and matrix%val
    !    arrays

    ! matrix:  is of the derived type {\tt ZD11\_TYPE} with
    !    {\tt INTENT (INOUT)}.
    type (zd11_type), intent (inout) :: matrix

    ! INFO: an INTEGER scaler of INTENT (OUT).
    !       INFO = 0 if the subroutine completes successfully; and
    !       INFO = MC65_WARN_MV_DIAG_MISSING if some diagonal
    !       elements are missing.
    integer (kind = myint), intent (OUT) :: info

    integer (kind = myint) :: m

    info = 0

    m = matrix%m
    if (zd11_get(matrix%type) == "pattern") then
       call csr_matrix_diagonal_first_priv(m,matrix%ptr,matrix%col,info)
    else
       call csr_matrix_diagonal_first_priv(m,matrix%ptr,matrix%col,info, &
                                           matrix%val)
    end if
  end subroutine csr_matrix_diagonal_first

  subroutine csr_matrix_diagonal_first_priv(m,ia,ja,info,a)
    ! private routine for csr_matrix_diagonal_first for efficiency

    ! ia: row pointer to matrix
    ! ja: column indices
    integer (kind = myint), dimension (*), intent (inout) :: ia,ja
    ! a: entry values
    real (kind = myreal), dimension (*), intent (inout), optional :: a

    integer (kind = myint) :: m,i,j,jdiag,jf,info
    real (kind = myreal) :: tmp

    if (present(a)) then
       do i=1,m
          jdiag = -1
          jf = ia(i)
          do j = jf,ia(i+1)-1
             if (ja(j) == i) then
                jdiag = j
                exit
             end if
          end do
          if (jdiag < 0) then
             INFO = MC65_WARN_MV_DIAG_MISSING
             cycle
          end if
          ! swap the j-th and the jdiag-th entry in the row
          tmp = a(jdiag)
          a(jdiag) = a(jf)
          a(jf) = tmp
          ja(jdiag) = ja(jf)
          ja(jf) = i
       end do
    else
       do i=1,m
          jdiag = -1
          jf = ia(i)
          do j = jf,ia(i+1)-1
             if (ja(j) == i) then
                jdiag = j
                exit
             end if
          end do
          if (jdiag < 0) then
             INFO = MC65_WARN_MV_DIAG_MISSING
             cycle
          end if
          ja(jdiag) = ja(jf)
          ja(jf) = i
       end do
    end if
  end subroutine csr_matrix_diagonal_first_priv


  subroutine csr_matrix_write(matrix,file_name,info,form,row_ord,&
       col_ord,stat)
    ! write the matrix in the file <file_name> with
    ! various <form>. If row_ord or/and col_ord is present,
    ! the reordered matrix will be written
    !  (unless form = "unformatted"
    ! or form = "hypergraph", in which case no reordering will
    !  take place).

    ! matrix: is of the derived type {\tt ZD11\_TYPE} with
    !         {\tt INTENT (IN)}.
    !         It is the matrix to be written.
    type (zd11_type), intent (in) :: matrix

    ! file_name: is a CHARACTER array of assumed size with
    !   {\tt INTENT (IN)}.
    ! It is the file to write the matrix in.
    character (len = *), intent (in) :: file_name

    ! form: an OPTIONAL CHARACTER array of assumed size with
    !         {\tt INTENT (IN)}. It
    !         is the format in which the matrix is to be written.
    !         Not supplying <form>
    !         or a <form> that is not one of the following strings
    !         has the effect of
    !         form = "ija".
    !         In the follow let M, N and NZ be the row, column
    !         indices and the number of entries in the matrix.
    !         TYPE be the type of the matrix as a character array
    !         with trailing blanks removed.
    !         LEN_TYPE be the string length of TYPE.
    ! form =
    !       "gnuplot": the matrix is written in <file_name> and
    !                the gnuplot command in <file_name.com>,
    !                then in gnuplot, you can type
    !                "load 'file_name'" to view the matrix
    !       "hypergraph": the matrix is written in <file_name>
    !                 in hypergraph format.
    !                       N,M
    !                       row_indices_of_column_1
    !                       row_indices_of_column_2
    !                       row_indices_of_column_3
    !                       ...
    !
    !
    !       "unformatted": the matrix is written in <file_name>
    !            as unformatted data:
    !                         LEN_TYPE
    !                         TYPE
    !                         M N NZ
    !                         row_pointers
    !                         column_indices
    !                         entry_values (Only if the matrix is
    !                         not of type "pattern")
    !
    !       "ij": the matrix is written in <file_name> in
    !             coordinate format,
    !             with only row and column indices. That is
    !                         M N NZ
    !                         ....
    !                         a_row_index a_column_index
    !                         another_row_index another_column_index
    !                         ...
    !
    !       "ija" or "coordinate": the matrix is written in
    !              <file_name> in coordinate format.
    !              Note that entry values will be written only
    !              if the matrix is
    !              not of type "pattern"
    !               M N NZ
    !               ....
    !               row_index column_index a(row_index,column_index)
    !               ....
    !
    character (len = *), intent (in), optional :: form

    ! INFO: an INTEGER of INTENT (OUT).
    !       INFO = 0 if the subroutine completes successfully
    !       INFO = MC65_ERR_NO_VACANT_UNIT if no vacant unit has
    !              been found
    !              in MC65_MATRIX_WRITE
    !       INFO = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       INFO = MC65_ERR_MEMORY_DEALLOC if memory
    !              deallocation failed
    INTEGER (KIND = MYINT), INTENT (OUT) :: INFO

    ! row_ord: is an OPTIONAL INTEGER array of size M with
    !          INTENT (IN). When present, a row index
    !          I will become ROW_ORD(I) in the reordered
    !          matrix.
    ! col_ord:  is an OPTIONAL INTEGER array of size N with
    !          INTENT (IN). When present, a column index
    !          I will become COL_ORD(I) in the reordered
    !          matrix.
    integer (kind = myint), dimension (matrix%m), intent (in), &
         optional :: row_ord
    integer (kind = myint), dimension (matrix%n), intent (in), &
         optional :: col_ord

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat


    ! =============== local variables ==============

    ! pattern_only: whether the matrix is of type "pattern"
    logical :: pattern_only

    integer (kind = myint) :: m,n,nz,ierr

    if (zd11_get(matrix%type) == "pattern") then
       pattern_only = .true.
    else
       pattern_only = .false.
    end if

    info = 0;
    if (present(stat)) stat = 0

    m = matrix%m
    n = matrix%n
    nz = matrix%ptr(m+1)-1

    if (.not.present(form)) then
       goto 10
    else if (form == "gnuplot".or.form == "GNUPLOT") then
       if (present(row_ord).and.present(col_ord)) then
          call csr_matrix_write_gnuplot(file_name,m,n,matrix%ptr,matrix%col,&
               info,row_ord,col_ord)
       else if (present(row_ord)) then
          call csr_matrix_write_gnuplot(file_name,m,n,matrix%ptr,matrix%col,&
               info,row_ord = row_ord)
       else if (present(col_ord)) then
          call csr_matrix_write_gnuplot(file_name,m,n,matrix%ptr,matrix%col,&
               info,col_ord = col_ord)
       else
          call csr_matrix_write_gnuplot(file_name,m,n,matrix%ptr,matrix%col,&
               info)
       end if
       return
    else if (form == "hypergraph" .or. form == "HYPERGRAPH") then
       call csr_matrix_write_hypergraph(matrix,file_name,info,stat=ierr)
       if (present(stat)) stat = ierr
       return
    else if (form == "unformatted".or.form == "UNFORMATTED") then
       if (pattern_only) then
          call csr_matrix_write_unformatted(file_name,m,n,nz,&
               matrix%type,matrix%ptr,matrix%col,info)
       else
          call csr_matrix_write_unformatted(file_name,m,n,nz,&
               matrix%type,matrix%ptr,matrix%col,info,matrix%val)
       end if
       return
    else if (form == "ij".or.form == "IJ") then
       pattern_only = .true.
       goto 10
    else if (form == "ija".or.form == "coordinate".or.&
         form == "IJA".or.form == "COORDINATE") then
    end if

10  if (pattern_only) then
       if (present(row_ord).and.present(col_ord)) then
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,row_ord = row_ord,col_ord = col_ord)
       else if (present(row_ord)) then
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,row_ord = row_ord)
       else if (present(col_ord)) then
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,col_ord = col_ord)
       else
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info)
       end if
    else
       if (present(row_ord).and.present(col_ord)) then
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,matrix%val,row_ord = row_ord,col_ord = col_ord)
       else if (present(row_ord)) then
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,matrix%val,row_ord = row_ord)
       else if (present(col_ord)) then
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,matrix%val,col_ord = col_ord)
       else
          CALL CSR_MATRIX_WRITE_ija(file_name,m,n,nz,matrix%ptr,matrix%col, &
                   info,matrix%val)
       end if
    end if

  end subroutine csr_matrix_write


  subroutine CSR_MATRIX_WRITE_ija(file_name,m,n,nz,ia,ja,info,&
       a,row_ord,col_ord)
    ! file_name: is a CHARACTER array of assumed size with
    !            {\tt INTENT (IN)}.
    ! It is the file to write the matrix in.
    ! m: row dimension
    ! n: column dimension
    ! nz: number of nonzeros
    ! ia: row pointer
    ! ja: column indices
    ! info: flag to show whether successful or not
    ! a: entry values
    ! row_ord: row ordering. NewIndex = row_ord(OldIndex)
    ! col_ord: col ordering. NewIndex = col_ord(OldIndex)
    character (len = *), intent (in) :: file_name
    integer (kind = myint), intent (in) :: m,n,nz,ia(*),ja(*)
    integer (kind = myint), intent (inout) :: info
    real (kind = myreal), intent (in), optional :: a(*)
    integer (kind = myint), intent (in), optional :: row_ord(*),&
         col_ord(*)
    integer (kind = myint) :: file_unit
    integer (kind = myint) :: i,k

    file_unit = vacant_unit()
    if (file_unit < 0) then
       info = MC65_ERR_NO_VACANT_UNIT
       return
    end if
    open(unit = file_unit, file = file_name, form = "formatted")
    write(file_unit,"(i10,i10,i10)") m, n, nz

    if (present(row_ord).and.present(col_ord)) then
       if (present(a)) then
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10,g12.4)") &
                  (row_ord(i),col_ord(ja(k)),a(k),k = ia(i),ia(i+1)-1)
          end do
       else
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10)") &
                  (row_ord(i),col_ord(ja(k)),k = ia(i),ia(i+1)-1)
          end do
       end if
    else if (present(row_ord)) then
       if (present(a)) then
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10,g12.4)") &
                  (row_ord(i),ja(k),a(k),k = ia(i),ia(i+1)-1)
          end do
       else
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10)") &
                  (row_ord(i),ja(k),k = ia(i),ia(i+1)-1)
          end do
       end if
    else if (present(col_ord)) then
       if (present(a)) then
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10,g12.4)") &
                  (i,col_ord(ja(k)),a(k),k = ia(i),ia(i+1)-1)
          end do
       else
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10)") &
                  (i,col_ord(ja(k)),k = ia(i),ia(i+1)-1)
          end do
       end if
    else
       if (present(a)) then
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10,g12.4)") &
                  (i,ja(k),a(k),k = ia(i),ia(i+1)-1)
          end do
       else
          do i = 1, m
             if (ia(i)>=ia(i+1)) cycle
             write(file_unit,"(2i10)") (i,ja(k),k = ia(i),ia(i+1)-1)
          end do
       end if
    end if
    close(file_unit)

  end subroutine CSR_MATRIX_WRITE_ija

  subroutine csr_matrix_write_unformatted(file_name,m,n,nz,&
       type,ia,ja,info,a)
    ! file_name: is a CHARACTER array of assumed size with
    !  {\tt INTENT (IN)}.
    ! It is the file to write the matrix in.
    ! m: row dimension
    ! n: column dimension
    ! nz: number of nonzeros
    ! type: character pointer of matrix type
    ! ia: row pointer
    ! ja: column indices
    ! info: flag to show whether successful or not
    ! a: entry values
    character (len = *), intent (in) :: file_name
    integer (kind = myint), intent (in) :: m,n,nz,ia(*),ja(*)
    character, allocatable, dimension (:) :: type
    integer (kind = myint), intent (inout) :: info
    real (kind = myreal), intent (in), optional :: a(*)
    integer (kind = myint) :: file_unit

    file_unit = vacant_unit()
    if (file_unit < 0) then
       info = MC65_ERR_NO_VACANT_UNIT
       return
    end if
    open(unit = file_unit, file = file_name, form = "unformatted")
    write(file_unit) len_trim(zd11_get(type))
    write(file_unit) zd11_get(type)
    write(file_unit) m,n,nz
    write(file_unit) ia(1:m+1)
    if (nz >= 0) then
       write(file_unit) ja(1:nz)
       if (present(a)) write(file_unit) a(1:nz)
    end if
    close(file_unit)
  end subroutine csr_matrix_write_unformatted

  subroutine csr_matrix_write_gnuplot(file_name,m,n,ia,ja,&
       info,row_ord,col_ord)
    ! write a matrix in a gnuplot data file
    ! "file_name" and the corresponding gnuplot command
    ! to view the matrix in "file_name.com"
    ! The matrix is reordered using col_ord and row_ord
    ! before writing




    ! file_name: is a CHARACTER array of assumed size with
    !   {\tt INTENT (IN)}.
    ! It is the file to write the matrix in.
    ! m: row dimension
    ! n: column dimension
    ! ia: row pointer
    ! ja: column indices
    ! info: info tag.
    ! row_ord: row ordering. NewIndex = row_ord(OldIndex)
    ! col_ord: col ordering. NewIndex = col_ord(OldIndex)
    character (len = *), intent (in) :: file_name
    integer (kind = myint), intent (in) :: m,n,ia(*),ja(*)

    ! col_ord, row_ord: column and row ordering. if a
    !  column col_ord(i) = 1,
    ! column i will be the first column in the ordered matrix
    integer (kind = myint), dimension (*), intent (in), optional :: &
         col_ord, row_ord
    integer (kind = myint), intent (OUT) :: info

    integer (kind = myint) :: i,j
    integer (kind = myint) :: unit1,unit2
    character (len = 120):: fmt_col

    unit1 = vacant_unit()
    if (unit1 < 0) then
       info = MC65_ERR_NO_VACANT_UNIT
       return
    end if
    open(unit = unit1,file = file_name)

    unit2 = vacant_unit()
    if (unit2 < 0) then
       info = MC65_ERR_NO_VACANT_UNIT
       return
    end if
    open(unit = unit2,file = file_name/&
         &/".com")

    write(unit2,"(a,i12,a)") "set xrange [0:",n+1,"]"
    write(unit2,"(a,i12,a)") "set yrange [0:",m+1,"]"
    write(unit2,"(a,i12,a,i12,a)") 'set ytics ("1" ',&
         m,', "',m,'"  1 )'
    write(unit2,"(a)") 'set noxtics'

    ! use in formating for printing number of columns
    fmt_col = "(a,a,i"/&
         &/int2str(digit(n),digit(digit(n)))/&
         &/",a,i12,a)"
    write(unit2,fmt_col) 'set x2tics ("1" 1', &
         ', "', &
         n, &
         '" ', &
         n, &
         ' )'
    write(unit2,"(a)") 'set nolabel '
    !    write(unit2,"(a)") 'set term post portrait '
    !    write(unit2,"(a,a)") 'set out "',file_name/&
    !         &/'.ps"'
    ! to set x and y equal scale
    write(unit2,"(a)") 'set size  square'

    write(unit2,"(a,a,a)") 'plot "',file_name,&
         '" using 1:2 notitle w p ps 0.3'

    if (present(row_ord).and.present(col_ord)) then
       write(unit1,'(2i10)') ((col_ord(ja(j)), m+1-row_ord(i), &
            j=ia(i), ia(i+1)-1), i = 1, m)
    else if (present(row_ord)) then
       write(unit1,'(2i10)') ((ja(j), m+1-row_ord(i), &
            j=ia(i), ia(i+1)-1), i = 1, m)
    else if (present(col_ord)) then
       write(unit1,'(2i10)') ((col_ord(ja(j)), m+1-i, &
            j=ia(i), ia(i+1)-1), i = 1, m)
    else
       write(unit1,'(2i10)') ((ja(j), m+1-i, j=ia(i), ia(i+1)-1), &
            i = 1, m)
    end if

    close(unit1)
    close(unit2)
  end subroutine csr_matrix_write_gnuplot



  subroutine csr_matrix_write_hypergraph(matrix,file_name,info,stat)
    ! write the matrix in hypergraph format.
    ! that is,
    ! row_dimension_of_A^T column_dimension_of_A^T
    ! column_indices_of_row1_of_A^T
    ! column_indices_of_row2_of_A^T
    ! column_indices_of_row3_of_A^T
    ! .....


    ! matrix: the matrix to be written
    type (zd11_type), intent (in) :: matrix

    ! file_name: the file name to write the data into
    character (len = *), intent (in) :: file_name

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat

    integer (kind = myint) :: i,n,m,l1,l2,ierr
    integer (kind = myint) :: unit1
    type (zd11_type) :: matrix2
    character*80 fmt
    integer (kind = myint), intent (inout) :: info

    if (present(stat)) stat = 0

    call csr_matrix_transpose(matrix,matrix2,info,stat=ierr)
    if (present(stat)) stat = ierr
    if (info < 0) return

    unit1 = vacant_unit()
    if (unit1 < 0) then
       info = MC65_ERR_NO_VACANT_UNIT
       return
    end if
    open(unit = unit1,file = file_name)

    m = matrix%m
    n = matrix%n
    write(unit1,"(2i12)") n,m
    do i = 1, n
       l1 = matrix2%ptr(i); l2 = matrix2%ptr(i+1)-1
       write(fmt,"('(',i12,'i8)')") max(1,l2-l1+1)
       write(unit1,fmt) matrix2%col(l1:l2)
    end do

    call csr_matrix_destruct(matrix2,info,stat=ierr)
    if (present(stat)) stat = ierr
    if (info < 0) return
    close(unit1)
  end subroutine csr_matrix_write_hypergraph








  subroutine csr_matrix_read(matrix,file_name,info,form,stat)
    ! read the matrix as dumped by csr_matrix_write,
    !  unformatted only

    ! matrix:  is of the derived type {\tt ZD11\_TYPE}
    !  with {\tt INTENT (INOUT)}.
    type (zd11_type), intent (inout) :: matrix

    ! file_name:  is a {\tt CHARACTER} array of assumed size
    !              with {\tt INTENT (IN)}.
    !              It is the file to write the matrix in.
    character (len = *), intent (in) :: file_name

    ! INFO: is a INTEGER of INTENT (OUT).
    !       INFO = 0 if the subroutine returns successfully
    ! info = MC65_mem_alloc if memory allocation failed
    ! info = MC65_mem_dealloc if memory deallocation failed
    !       INFO = MC65_ERR_READ_WRONGFORM if in
    !              MC65_MATRIX_READ,form is present
    !              but is of an unsupported form.
    !       INFO = MC65_ERR_READ_FILE_MISS if in MC65_MATRIX_READ,
    !              the file <file_name> does not exist.
    !       INFO = MC65_ERR_READ_OPEN if in MC65_MATRIX_READ,
    !              opening of the file <file_name> returns iostat /= 0
    !       INFO = MC65_ERR_NO_VACANT_UNIT if no vacant unit
    !              has been found in MC65_MATRIX_WRITE
    integer (kind = myint), intent (out) :: INFO


    ! form:  is an {\tt OPTIONAL CHARACTER} array of
    !        assumed size with
    !       {\tt INTENT (IN)}. It is the format of the input file
    !        Only form = "unformatted" is supported. If
    !        form is not present, by default
    !        the subroutine assume that the input file is unformatted.
    !        If form is present but form != "unformatted",
    !        the subroutine will return with info =
    !        MC65_ERR_READ_WRONGFORM
    character (len = *), intent (in), optional :: form

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat


    ! file_unit: open the file in this unit
    integer (kind = myint) :: file_unit
    ! m,n: row and column dimension
    ! nz: number of nonzeros
    integer (kind = myint) :: m,n,nz,ierr

    integer (kind = myint) :: ios
    logical :: existed, opened
    integer (kind = myint), parameter :: maxlen = 2000
    character (len=maxlen) :: type
    integer (kind = myint) :: len_type
    integer (kind = myint) :: iform, UNFORMAT = 1, UNSUPPORTED = 2

    info = 0
    if (present(stat)) stat = 0

    iform = UNFORMAT
    if (present(form)) then
       if (form == "unformatted") then
          iform = UNFORMAT
          !  else if ...
          ! .... other forms comes here
       else
          iform = UNSUPPORTED
          info = MC65_ERR_READ_WRONGFORM
          return
       end if
    end if


    file_unit = vacant_unit()
    if (file_unit < 0) then
       info = MC65_ERR_NO_VACANT_UNIT
       return
    end if

    inquire(file = file_name, exist = existed, opened = opened, &
         iostat = ios)
    if (.not.existed) then
       info = MC65_ERR_READ_FILE_MISS
       return
    end if
    if (opened) then
       rewind(file_unit)
    end if
    if (ios /= 0) then
       info = MC65_ERR_READ_OPEN
    end if

    if (iform == UNFORMAT) then
       open(unit = file_unit, file = file_name, form = "unformatted")
    else
       open(unit = file_unit, file = file_name, form = "formatted")
    end if


    if (iform == UNFORMAT) then
       type = ""
       read(file_unit) len_type
       if (len_type > maxlen) then
          INFO = MC65_ERR_READ_MAXLEN
          return
       end if
       read(file_unit) type(1:len_type)
       read(file_unit) m,n,nz
       call csr_matrix_construct(matrix,m,nz,info,n=n,type=type,stat=ierr)
       if (present(stat)) stat = ierr
       if (info < 0) then
          return
       end if
       read(file_unit) matrix%ptr(1:matrix%m+1)
       if (nz >= 1) then
          read(file_unit) matrix%col(1:nz)
          if (ZD11_get(matrix%type) /= "pattern" ) &
               read(file_unit) matrix%val(1:nz)
       end if

       ! else if (...)
       ! other format will be here ...
    end if

    close(file_unit)

  end subroutine csr_matrix_read



  subroutine csr_matrix_condense(matrix,info,col_wgt_real,&
       col_wgt_int,realloc,iremove,stat)
    !   subroutine csr_matrix_condense(matrix,info[,&
    !        col_wgt_real,col_wgt_int,realloc])
    ! given a matrix, finds all columns that has
    ! <= iremove (by default 0) element and delete them, also
    ! merge columns with exactly the same pattern.
    ! it only works on matrix with pattern only.
    ! The matrix should also be cleaned using
    ! matrix_clean to agglomerate repeated
    ! entries!

    ! MATRIX is of the derived type {\tt ZD11\_TYPE} with
    !        {\tt INTENT (INOUT)}. on exit all columns that
    !        has no more than one entry are deleted;
    !        columns that have the same pattern are merged into one.
    !        The column dimension of the matrix is changed
    !        accordingly.
    type (zd11_type), intent (inout) :: matrix


    ! INFO  is an {\tt INTEGER} of {\tt INTENT (OUT)}.
    !       INFO = 0 if the subroutine completes successfully
    !       INFO = MC65_ERR_CONDENSE_NOPAT if the matrix
    !              is not of type {\tt "pattern"}
    !       INFO = MC65_ERR_MEMORY_ALLOC if memory allocation failed
    !       INFO = MC65_ERR_MEMORY_DEALLOC if memory
    !              deallocation failed.
    integer (kind = myint), intent (out) :: info

    ! col_wgt_real: is an OPTIONAL REAL array with INTENT (INOUT) of
    !          size N (the column dimension). When present,
    !          On entry, the first N elements contain the column
    !          weights of the matrix,
    !          on exit, the first NSUP elements of col_wgt
    !          holds the super-column weight of the condensed matrix.
    !          Here nsup is the new column dimension of the
    !          condensed matrix. Note that the two optional arguments
    !          col_wgt_real and col_wgt_int are supplied so that
    !          the user is free to treat column weight either as
    !          real arrays, or as integer arrays.
    !          The correct usage is therefore
    !          to supply, if appropriate, either col_wgt_real or
    !          col_wgt_int as a keyword argument
    !          ({\tt col_wgt_real = colwgt} or
    !          (\tt col_wgt_int = colwgt},
    !          depending on whether colwgt is real or integer)
    ! col_wgt_int: is an OPTIONAL INTEGER array with INTENT (INOUT)
    !          of size N  (the column dimension). When present,
    !          On entry, the first N elements contain the column
    !           weights of the matrix,
    !          on exit, the first NSUP elements of col_wgt
    !          holds the super-column weight of the condensed matrix.
    !          Here nsup is the new column dimension of the
    !          condensed matrix.
    !          Note that the two optional arguments
    !          col_wgt_real and col_wgt_int are supplied so that
    !          the user is free to treat column weight either as
    !          real arrays, or as integer arrays.
    !          The correct usage is therefore
    !          to supply, if appropriate, either col_wgt_real or
    !          col_wgt_int as a keyword argument
    !          ({\tt col_wgt_real = colwgt} or
    !          (\tt col_wgt_int = colwgt},
    !          depending on whether colwgt is real or integer).
    integer (kind = myint), dimension (matrix%n), optional, &
         intent (inout) :: col_wgt_int
    real (kind = myreal), dimension (matrix%n), optional, &
         intent (inout) :: col_wgt_real


    ! realloc:  is an optional real scaler of INTENT (IN).
    !   It is used to control the reallocation of storage.
    ! \begin{itemize}
    !  \item{} if {\tt REALLOC < 0}, no reallocation.
    !  \item{} if {\tt REALLOC == 0}, reallocation is carried out
    !  \item{} if {\tt REALLOC > 0}, reallocation iff memory
    !       saving is greater than
    !      {\tt REALLOC*100}\%. For example,
    !      if {\tt REALLOC = 0.5}, reallocation will be carried out
    !      if saving in storage is greater than 50\%
    !  \end{itemize}
    ! If {\tt REALLOC} is not present, no reallocation is carried out.
    real (kind = myreal), INTENT(IN), optional :: realloc

    ! iremove: an OPTIONAL INTEGER scalar. When present,
    !    any columns that have <= iremove
    !    entries are removed. If not present, only columns
    !    of zero entries are removed.
    integer (kind = myint), intent (in), optional :: iremove

    ! stat: is an optional integer saler of INTENT (OUT). If supplied,
    ! on exit it holds the error tag for memory allocation
    integer (kind = myint), optional, intent (out) :: stat
    ! ===================== local variables =======================

    ! col_wgt_new: working array. the new super-column weight.
    real (kind = myreal), dimension (:), allocatable :: col_wgt_new_r
    integer (kind = myint), dimension (:), allocatable :: &
         col_wgt_new_i

    ! f2c: the array which says which original column is
    !      merged into which super-column
    ! numsup: how many original columns make up the super-columns
    ! new: the new identity of a supercolumn
    ! flag: in which row the first appearance of a super-column?
    integer (kind = myint), dimension (:), allocatable :: &
         f2c,numsup,new,flag

    ! nsup: number of super-columns so far
    ! m,n: row and column size
    ! newid: new group of an old super-column
    ! i,j,jj: loop index
    ! scol: super-column id ofan column entry in the current row
    ! nz: number of nonzeros in the condensed matrix
    ! jjs: starting index of original row in the matrix
    integer (kind = myint) :: nsup,n,m,newid,i,j,jj,scol,nz,jjs,&
         nz_old,ierr

    ! number of super-columns after getting rid of columns
    ! with less than one entry
    integer (kind = myint) :: nsup_new

    ! col_count: a counter for (the number of entries - 1) in
    !            each column
    integer (kind = myint), dimension (:), allocatable :: col_count
    ! col_count_new: a counter for (the number of entries - 1)
    !                in each column in the condensed matrix
    !                (of super-columns)
    ! sup_to_sup: super-column index after condensing and that after
    ! getting rid of columns with only <= 1 entries
    integer (kind = myint), dimension (:), allocatable :: &
         col_count_new, sup_to_sup

    ! num_remove: any column of <= num_remove entries are
    !             removed. by default num_remove = 0
    integer (kind = myint) :: num_remove

    num_remove = 0
    if (present(iremove)) num_remove = iremove

    info = 0 ! by default the subroutine is successful.
    if (present(stat)) stat = 0

    n = matrix%n; m = matrix%m

    ! check that the matrix is pattern only!
    if (ZD11_get(matrix%type) /= "pattern" ) then
       info = MC65_ERR_CONDENSE_NOPAT
       return
    end if

    ! if no column no need to work further.
    if (n <= 0) return

    ! all columns assigned as super-column 1 first
    allocate(f2c(n),numsup(n),new(n),flag(n),col_count(n),stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       INFO = MC65_ERR_MEMORY_ALLOC
       return
    end if

    f2c = 1
    numsup(1) = n
    flag(1) = 0
    nsup = 1
    col_count = 0

    ! loop over each row, assign each column entry as a new
    !    super-column
    do i = 1, m
       ! column entries in this row will change their super-column
       ! number, so reduce the number of super-columns of their
       ! old group
       do jj = matrix%ptr(i),matrix%ptr(i+1)-1
          j = matrix%col(jj)
          col_count(j) = col_count(j) + 1
          scol = f2c(j)
          numsup(scol) = numsup(scol) - 1
       end do
       ! assign new groups
       do jj = matrix%ptr(i),matrix%ptr(i+1)-1
          j = matrix%col(jj)
          scol = f2c(j)
         ! has this super-column been seen before in row i?
          if (flag(scol) < i) then
             ! the first appearance of scol in row i
             flag(scol) = i
             if (numsup(scol) > 0) then
                nsup = nsup + 1
                flag(nsup) = i
                f2c(j) = nsup
                numsup(nsup) = 1
                new(scol) = nsup
             else
                ! if zero, the old group disappear anyway, there is
                ! no need to add a new group name, just use
                ! the old one
                numsup(scol) = 1
                new(scol) = scol
             end if
          else
             newid = new(scol)
             f2c(j) = newid
             numsup(newid) = numsup(newid) + 1
          end if
       end do
    end do

    ! now get rid of columns of <= num_remove entries by modifying f2c
    allocate(col_count_new(nsup),sup_to_sup(nsup),stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       INFO = MC65_ERR_MEMORY_ALLOC
       return
    end if

    ! The following array statement can't be used, as multiple columns
    ! (with the _same_ non-zero pattern) can be assigned to each supercolumn
    ! so we need to use a loop instead
    !col_count_new(f2c(1:n)) = col_count(1:n)
    do i = 1, n
       col_count_new(f2c(i)) = col_count(i)
    end do
    nsup_new = 0
    do i = 1, nsup
       if (col_count_new(i) > num_remove) then
          nsup_new = nsup_new + 1
          sup_to_sup(i) = nsup_new
       else
          sup_to_sup(i) = -1
       end if
    end do
    f2c(1:n) = sup_to_sup(f2c(1:n))
    nsup = nsup_new

    deallocate(col_count,col_count_new,sup_to_sup,stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       INFO = MC65_ERR_MEMORY_DEALLOC
       return
    end if

    ! condense the matrix
    flag = 0
    nz = 1
    do i = 1, m
       ! loop over each row
       jjs = matrix%ptr(i)
       matrix%ptr(i) = nz
       do jj = jjs,matrix%ptr(i+1)-1
          j = matrix%col(jj)
          scol = f2c(j)
          ! forget about columns with less than num_remove entries
          if (scol < 0) cycle
          ! if this the first time this group appear?
          ! if so this column will remain
          if (flag(scol) < i) then
             flag(scol) = i
             matrix%col(nz) = scol
             nz = nz + 1
          end if
       end do
    end do
    nz_old = size(matrix%col)
    matrix%ptr(m+1) = nz
    matrix%n = nsup

    nz = nz -1
    if (present(realloc)) then
       if (realloc == 0.0_myreal.or.(realloc > 0.and.&
            nz_old > (1.0_myreal + realloc)*nz)) then
          call csr_matrix_reallocate(matrix,nz,info)
          if (info < 0) return
       end if
    end if


    ! condense the column weight
    deallocate(numsup,new,flag,stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       INFO = MC65_ERR_MEMORY_DEALLOC
       return
    end if

    if (present(col_wgt_real)) then
       allocate(col_wgt_new_r(nsup),stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          INFO = MC65_ERR_MEMORY_ALLOC
          return
       end if

       col_wgt_new_r = 0
       do i = 1, n
          scol = f2c(i)
          ! forget about columns with less than num_remove entries
          if (scol < 0) cycle
          col_wgt_new_r(scol) = col_wgt_new_r(scol) + col_wgt_real(i)
       end do
       col_wgt_real(1:nsup) = col_wgt_new_r
       deallocate(col_wgt_new_r,stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          INFO = MC65_ERR_MEMORY_DEALLOC
       end if
    end if


    if (present(col_wgt_int)) then
       allocate(col_wgt_new_i(nsup),stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          INFO = MC65_ERR_MEMORY_ALLOC
          return
       end if

       col_wgt_new_i = 0
       do i = 1, n
          scol = f2c(i)
          ! forget about columns with less than num_remove entries
          if (scol < 0) cycle
          col_wgt_new_i(scol) = col_wgt_new_i(scol) + col_wgt_int(i)
       end do
       col_wgt_int(1:nsup) = col_wgt_new_i
       deallocate(col_wgt_new_i,stat = ierr)
       if (present(stat)) stat = ierr
       if (ierr /= 0) then
          INFO = MC65_ERR_MEMORY_DEALLOC
       end if
    end if


    deallocate(f2c,stat = ierr)
    if (present(stat)) stat = ierr
    if (ierr /= 0) then
       INFO = MC65_ERR_MEMORY_DEALLOC
    end if


  end subroutine csr_matrix_condense

  function csr_matrix_is_pattern(matrix) result (pattern)
    ! check is a matrix is of type "pattern"
    ! MATRIX is of the derived type {\tt ZD11\_TYPE}
    !        with {\tt INTENT (IN)}.
    !        On entry it is the matrix whose type is to
    !        be determined
    type (zd11_type), intent (in) :: matrix
    ! pattern: logical scalar. On exit
    ! is set to true if the matrix is of type "pattern"
    !  and false if not
    logical :: pattern
    if (ZD11_get(matrix%type) == "pattern")  then
       pattern = .true.
    else
       pattern = .false.
    end if
  end function csr_matrix_is_pattern

! ======================== private routines from here =============


  function vacant_unit()  result (file_unit)
    ! utility routine to find a vacant unit and
    ! returns a unit number of a unit that exists and is
    ! not connected

    ! max_unit: maximum number of units to check
    integer (kind = myint) :: max_unit
    parameter (max_unit = 500)
    ! file_unit: the vacant unit
    integer (kind = myint) :: file_unit

    logical :: existed, opened
    integer (kind = myint) :: ios

    do file_unit = 10, max_unit
       inquire (unit = file_unit, exist = existed, &
            opened = opened, iostat = ios)
       if (existed .and. .not. opened .and. ios == 0) return
    end do

    file_unit = -1

  end function vacant_unit




  subroutine expand1(p1,upbound_new1,info,ierr)
    ! reallocate memory for an integer or real allocatable array
    !     of up to 4 dimension.
    ! usage: expand(array,new_bound_1, ...,new_bound_n)
    ! a space of size "new_bound_1 X ... X new_bound_n"
    ! will be allocated and the content of array
    ! will be copied to the beginning of this
    ! memory.
    real (kind = myreal), allocatable:: p2(:),p1(:)
    integer (kind = myint) upbound_new1
    integer (kind = myint) ub(1),lb(1),upb(1)
    integer (kind = myint) info,i,ierr
    ierr=0

    upb(1)=upbound_new1

    lb=lbound(p1)
    ub=ubound(p1)


    allocate(p2(lb(1):upb(1)),stat=ierr)
    if (ierr == 0) then
       p2=0
       do i=lb(1),min(ub(1),upb(1))
          p2(i)=p1(i)
       end do
       deallocate(p1,stat = ierr)
       allocate(p1(lb(1):upb(1)),stat=ierr)
       if (ierr == 0) then
          do i=lb(1),min(ub(1),upb(1))
             p1(i)=p2(i)
          end do
       else
          info = MC65_ERR_MEMORY_DEALLOC
          return
       end if
    else
       info = MC65_ERR_MEMORY_ALLOC
    end if


  end subroutine expand1

  ! ============ for integer (kind = myint) arrays ================

  subroutine iexpand1(p1,upbound_new1,info,ierr)
    integer (kind = myint), allocatable:: p2(:),p1(:)
    integer (kind = myint) upbound_new1
    integer (kind = myint) ub(1),lb(1),upb(1)
    integer (kind = myint) info,i,ierr

    upb(1)=upbound_new1

    lb=lbound(p1)
    ub=ubound(p1)


    allocate(p2(lb(1):upb(1)),stat=ierr)
    if (ierr == 0) then
       p2=0
       do i=lb(1),min(ub(1),upb(1))
          p2(i)=p1(i)
       end do
       deallocate(p1,stat = ierr)
       allocate(p1(lb(1):upb(1)),stat=ierr)
       if (ierr == 0) then
          do i=lb(1),min(ub(1),upb(1))
             p1(i)=p2(i)
          end do
       else
          info = MC65_ERR_MEMORY_DEALLOC
          return
       end if
    else
       info = MC65_ERR_MEMORY_ALLOC
    end if

  end subroutine iexpand1



  function int2str(i,idigit) result (the_str)
    ! converting an integer to a string.
    integer (kind = myint):: i,idigit

    character (len = idigit) :: the_str

    character (len = 20) :: fmt_str

    fmt_str = "(I                 )"

    write(fmt_str(3:19),"(I6)") idigit

    write(the_str,fmt_str) i

  end function int2str

  function digit(i)
    ! return the number of digit of an integer
    integer (kind = myint) :: i,digit

    integer (kind = myint) :: absi


    absi = abs(i)

    digit = 1

    do while (absi >= 10)

       absi = absi / 10

       digit = digit + 1

    end do

    if (i < 0) digit = digit + 1

  end function digit





!!$  subroutine csr_matrix_reset_type(matrix,type,info)
!!$    ! subroutine csr_matrix_reset_type(matrix,type,info):
!!$    ! reset the type of the matrix
!!$
!!$    ! matrix: is of the derived type ZD11_type, INTENT (INOUT).
!!$    !         This is the matrix whose type is to be changed.
!!$    type (zd11_type), intent (inout) :: matrix
!!$
!!$    ! type: is a character array of unspecified length with
!!$    !       INTENT (IN).
!!$    !       It hold the new matrix type to be set.
!!$    character (len=*), intent (in) :: type
!!$
!!$    ! info: integer scaler of INTENT (OUT).
!!$    !       = 0 if the subroutine returned successfully
!!$    !       = MC65_ERR_MEMORY_ALLOC if memory allocation failed
!!$    !       = MC65_ERR_MEMORY_DEALLOC if memory deallocation failed
!!$    integer (kind = myint), intent (out) :: info
!!$
!!$    ! =============== local variables ===============
!!$
!!$    integer (kind =  myint) :: ierr
!!$
!!$    info = 0
!!$
!!$    if (allocated(matrix%type)) deallocate(matrix%type,stat = ierr)
!!$    if (ierr /= 0) then
!!$       info = MC65_ERR_MEMORY_DEALLOC
!!$       return
!!$    end if
!!$    call zd11_put(matrix%type,type,stat = ierr)
!!$    if (ierr /= 0) then
!!$       info = MC65_ERR_MEMORY_ALLOC
!!$    end if
!!$
!!$  end subroutine csr_matrix_reset_type
!!$
!!$
!!$  subroutine csr_matrix_fill(matrix,i,row_counter,fill_entry,&
!!$       fill_value)
!!$    ! fill an entry -- fill the current
!!$    ! entry with column index fill_entry and column value
!!$    ! fill_value, in row i of the matrix.
!!$    ! row_counter is a counter given
!!$    ! for the vacant position in the rows. This subroutine
!!$    ! is used when the matrix has been allocated storage, and also
!!$    ! the row pointer has already been established, and
!!$    ! one is at the stage of filling in the entry values
!!$    ! and column indices.
!!$
!!$    ! matrix: the matrix whose entries are to be filled.
!!$    type (zd11_type), intent (inout) :: matrix
!!$
!!$    ! i: the entry is to be filled at this i-th row.
!!$    integer (kind = myint), intent (in) :: i
!!$
!!$    ! row_counter: integer array, row_counter(i)
!!$    !              is number of entries filled in row i so far.
!!$    integer (kind = myint), dimension (*), intent (inout) :: &
!!$         row_counter
!!$
!!$    ! fill_entry: column index of the current entry
!!$    integer (kind = myint), intent (in) :: fill_entry
!!$    real (kind = myreal), intent (in), optional :: fill_value
!!$
!!$    integer (kind = myint), pointer, dimension (:) :: ja
!!$    real (kind = myreal), pointer, dimension (:) :: aa
!!$
!!$
!!$    ja => csr_matrix_getrow(matrix,i)
!!$    row_counter(i) = row_counter(i) + 1
!!$    ja(row_counter(i)) = fill_entry
!!$    if (present(fill_value)) then
!!$       aa => csr_matrix_getrowval(matrix,i)
!!$       aa(row_counter(i)) = fill_value
!!$    end if
!!$  end subroutine csr_matrix_fill
!!$
!!$
!!$  subroutine csr_matrix_component(matrix,root,mask,comp_nvtx,&
!!$       comp_nz,comp_vtx_list)
!!$    !
!!$    ! matrix_component: this subroutine finds the component
!!$    ! of the graph starting from the root. This component
!!$    ! is stored in the list of vertices "comp_vtx_list",
!!$    ! with "comp_nvtx" vertices and "comp_nz"/2 edges
!!$    ! (thus comp_nz nonzeros in the matrix describing
!!$    ! the component).
!!$
!!$    ! =================== arguments =========================
!!$
!!$    ! matrix: the graph to be inspected. This must be a symmetric
!!$    ! matrix with no diagonal
!!$    type (zd11_type), intent (in) :: matrix
!!$
!!$    ! root: the root of the current component
!!$    integer (kind = myint), intent (in) :: root
!!$
!!$    ! comp_nvtx: number of vertices in this component
!!$    ! comp_nz: number of edges*2 in the component
!!$    integer (kind = myint), intent (out) :: comp_nvtx,comp_nz
!!$
!!$    ! mask: has the dimension as the number of vertices.
!!$    ! It is zero if the vertex is not yet visited,
!!$    ! nonzero if visited. When a vertex is visited the first time,
!!$    ! its mask is set to the new index of this vertex in the component
!!$    ! comp_vtx_list: a list of the vertices in the current component.
!!$    ! mask and comp_vtx_list satisfies mask(comp_vtx(i)) = i,
!!$    ! i = 1,...,comp_nvtx
!!$    integer (kind = myint), dimension (:), intent (inout) :: mask, &
!!$         comp_vtx_list
!!$
!!$    ! ===================== local variables =====================
!!$    ! front_list: this array hold the list of vertices that
!!$    !             form the front in the breadth first search
!!$    !             starting from the root
!!$    integer (kind = myint), dimension (:), allocatable :: front_list
!!$
!!$    ! front_sta: the starting position in the front_list of the
!!$    !            current front.
!!$    ! front_sto: the ending position in the front_list of the current
!!$    ! front.
!!$    integer (kind = myint) :: front_sta, front_sto
!!$
!!$    ! n: number of vertices in the graph
!!$    ! ierr: error tag
!!$    integer (kind = myint) :: n,ierr
!!$    ! v: a vertex
!!$    ! u: neighbor of v
!!$    ! i: loop index
!!$    ! j: loop index
!!$    integer (kind = myint) :: i,j,u,v
!!$
!!$    n = matrix%m
!!$
!!$    ! the first front
!!$    allocate(front_list(n),stat = ierr)
!!$    if (ierr /= 0) stop "error allocating in matrix_component"
!!$    front_sta = 1
!!$    front_sto = 1
!!$    front_list(1) = root
!!$
!!$    comp_nvtx = 1
!!$    comp_nz = 0
!!$    comp_vtx_list(1) = root ! one vertex in the component so far
!!$    mask(root) = comp_nvtx ! mask the root
!!$
!!$    do while (front_sto-front_sta >= 0)
!!$       do i = front_sta, front_sto
!!$          v = front_list(i) ! pick a vertex from the front
!!$          ! count link to all neighbors as edges
!!$          comp_nz = comp_nz + matrix%ptr(v+1)-matrix%ptr(v)
!!$          do j = matrix%ptr(v),matrix%ptr(v+1)-1
!!$             u = matrix%col(j) ! pick its neighbor
!!$             if (mask(u) /= 0) cycle
!!$             comp_nvtx = comp_nvtx + 1 ! found a unmasked vertex
!!$             mask(u) = comp_nvtx ! mask this vertex
!!$             ! add this vertex to the component
!!$             comp_vtx_list(comp_nvtx) = u
!!$             front_list(comp_nvtx) = u ! also add it to the front
!!$          end do
!!$       end do
!!$       front_sta = front_sto + 1
!!$       front_sto = comp_nvtx
!!$    end do
!!$    deallocate(front_list,stat = ierr)
!!$    if (ierr /= 0) stop "error deallocating in matrix_component"
!!$
!!$
!!$  end subroutine csr_matrix_component
!!$
!!$  subroutine csr_matrix_crop(matrix,comp_nvtx,comp_vtx_list,&
!!$       mask,submatrix,comp_nz,info)
!!$    ! matrix_crop: generate a submatrix describing a component of
!!$    ! the graph (matrix).
!!$
!!$    ! =================== arguments =========================
!!$    ! matrix: the graph to be inspected. This must be a symmetric
!!$    ! matrix with no diagonal
!!$    type (zd11_type), intent (in) :: matrix
!!$
!!$    ! comp_nvtx: number of vertices in this component
!!$    ! comp_nz: number of edges*2 in the component (optional)
!!$    integer (kind = myint), intent (in) :: comp_nvtx
!!$    integer (kind = myint), intent (in), optional :: comp_nz
!!$    integer (kind = myint), intent (out) :: info
!!$
!!$    ! mask: has the dimension as the number of vertices.
!!$    ! The mask of a vertex is set to the new index
!!$    ! of this vertex in the component it belongs
!!$    ! comp_vtx_list: a list of the vertices in the current component
!!$    integer (kind = myint), dimension (:), intent (in) :: mask, &
!!$         comp_vtx_list
!!$
!!$    ! submatrix: the submatrix describing the component
!!$    type (zd11_type), intent (inout) :: submatrix
!!$
!!$    ! ===================== local variables =====================
!!$    ! i: loop index
!!$    ! j: loop index
!!$    ! nz: number of nonzeros in the submatrix
!!$    ! n: number of rows in the submatrix
!!$    ! m: number of columns in the submatrix
!!$    ! v: a vertex in its original index
!!$    ! l1: starting index
!!$    ! l2: stopping index
!!$    ! sl1: starting index
!!$    ! sl2: stopping index
!!$    integer (kind = myint) :: i,j,nz,n,m,v,l1,l2,sl1,sl2
!!$
!!$    ! ia: row pointer of the original graph
!!$    ! ja: column indices of the original graph
!!$    ! sia: row pointer of the subgraph
!!$    ! sja: column indices of the subgraph
!!$    integer (kind = myint), dimension (:), pointer :: ia,ja,sia,sja
!!$
!!$    ! a: entry values of the original graph
!!$    ! sa: entry values of the subgraph
!!$    real (kind = myreal), dimension (:), pointer :: a,sa
!!$
!!$    info = 0
!!$    ia => matrix%ptr
!!$    ja => matrix%col
!!$    if (.not.present(comp_nz)) then
!!$       nz = 0
!!$       do i = 1, comp_nvtx
!!$          v = comp_vtx_list(i)
!!$          nz = nz + ia(v+1) - ia(v)
!!$       end do
!!$    else
!!$       nz = comp_nz
!!$    end if
!!$
!!$    n = comp_nvtx; m = n;
!!$    call csr_matrix_construct(submatrix,n,nz,info,n = m,&
!!$         type = zd11_get(matrix%type))
!!$    sia => submatrix%ptr
!!$    sja => submatrix%col
!!$    if (ZD11_get(matrix%type) /= "pattern" ) then
!!$       sa => submatrix%val
!!$       a => matrix%val
!!$    end if
!!$    sia(1) = 1
!!$    do i = 1, comp_nvtx
!!$       v = comp_vtx_list(i)
!!$       l1 = ia(v); l2 = ia(v+1)-1
!!$       sia(i+1) = sia(i) + l2 - l1 + 1
!!$       sl1 = sia(i); sl2 = sia(i+1)-1
!!$       ! convert original index to new index in the component
!!$       sja(sl1:sl2) = mask(ja(l1:l2))
!!$       if (ZD11_get(matrix%type) /= "pattern" ) then
!!$          sa(sl1:sl2) = a(l1:l2)
!!$       end if
!!$    end do
!!$  end subroutine csr_matrix_crop
!!$
!!$
!!$
!!$
!!$  subroutine csr_matrix_crop_unsym(matrix,nrow,row_list,submatrix,info)
!!$    ! matrix_crop: generate a submatrix of the whole matrix
!!$    ! by taking all the nrow rows in the row_list.
!!$    ! column size of submatrix is assigned as the same as that
!!$    ! for matrix.
!!$
!!$    ! =================== arguments =========================
!!$    ! matrix: the matrix to be inspected.
!!$    type (zd11_type), intent (in) :: matrix
!!$
!!$    ! nrow: number of vertices in this component
!!$    integer (kind = myint), intent (in) :: nrow
!!$
!!$    ! row_list: a list of the vertices in the current component
!!$    integer (kind = myint), dimension (:), intent (in) :: row_list
!!$
!!$    ! submatrix: the submatrix describing the component
!!$    type (zd11_type), intent (inout) :: submatrix
!!$    integer (kind = myint), intent (out) :: info
!!$    ! ===================== local variables =====================
!!$    ! i: loop index
!!$    ! j: loop index
!!$    ! nz: number of nonzeros in the submatrix
!!$    ! n: number of rows in the submatrix
!!$    ! m: number of columns in the submatrix
!!$    ! v: a vertex in its original index
!!$    ! l1: starting index
!!$    ! l2: stopping index
!!$    ! sl1: starting index
!!$    ! sl2: stopping index
!!$    integer (kind = myint) :: i,j,nz,n,m,v,l1,l2,sl1,sl2
!!$
!!$    ! ia: row pointer of the original matrix
!!$    ! ja: column indices of the original matrix
!!$    ! sia: row pointer of the submatrix
!!$    ! sja: column indices of the submatrix
!!$    integer (kind = myint), dimension (:), pointer :: ia,ja,sia,sja
!!$
!!$    ! a: entry values of the original matrix
!!$    ! sa: entry values of the submatrix
!!$    real (kind = myreal), dimension (:), pointer :: a,sa
!!$
!!$    info = 0
!!$    ia => matrix%ptr
!!$    ja => matrix%col
!!$    nz = 0
!!$    do i = 1, nrow
!!$       v = row_list(i)
!!$       nz = nz + ia(v+1) - ia(v)
!!$    end do
!!$
!!$    n = nrow; m = matrix%n;
!!$    call csr_matrix_construct(submatrix,n,nz,info,n = m,&
!!$         type = zd11_get(matrix%type))
!!$
!!$    sia => submatrix%ptr
!!$    sja => submatrix%col
!!$    if (ZD11_get(matrix%type) /= "pattern" ) then
!!$       sa => submatrix%val
!!$       a => matrix%val
!!$    end if
!!$    sia(1) = 1
!!$    do i = 1, nrow
!!$       v = row_list(i)
!!$       l1 = ia(v); l2 = ia(v+1)-1
!!$       sia(i+1) = sia(i) + l2 - l1 + 1
!!$       sl1 = sia(i); sl2 = sia(i+1)-1
!!$       ! convert original index to new index in the component
!!$       sja(sl1:sl2) = ja(l1:l2)
!!$       if (ZD11_get(matrix%type) /= "pattern" ) then
!!$          sa(sl1:sl2) = a(l1:l2)
!!$       end if
!!$    end do
!!$  end subroutine csr_matrix_crop_unsym
!!$
!!$

end module HSL_MC65_double




! COPYRIGHT (c) 1995 Council for the Central Laboratory
!                    of the Research Councils
! Original date 23 March 2001
!  March 2001: threadsafe version of HSL_FA04
!  December 2003. Long lines shortened.

!-*-*-*-*-*-*-*-*  H S L _ F A 1 4   M O D U L E  *-*-*-*-*-*-*-*-*-*-*-*-*-

! 12th July 2004 Version 1.0.0. Version numbering added.

      MODULE HSL_FA14_DOUBLE

!  Portable random number generator by Linus Schrange, TOMS 5, 1979, pp 132-138
!  Fortran 95 version by Nick Gould and John Reid, RAL, September 1995

         IMPLICIT NONE

         PRIVATE
         PUBLIC :: FA14_RANDOM_REAL, FA14_RANDOM_INTEGER, FA14_GET_SEED,      &
                   FA14_SET_SEED, FA14_INITIALIZE

!  Define the working precision to be double

         INTEGER, PARAMETER :: working = KIND( 1.0D+0 )

         TYPE, PUBLIC :: FA14_seed
            PRIVATE
            INTEGER :: ix = 65535
         END TYPE
         INTEGER, PARAMETER :: a = 16807, b15 = 32768
         INTEGER, PARAMETER :: b16 = 65536, p = 2147483647
         INTEGER, PARAMETER :: b30 = 1073741824, q = 1073741823

      CONTAINS

!-*-*-*-*-*- F A 1 4 _ R A N D O M _ R E A L  S U B R O U T I N E  *-*-*-*-*-

         SUBROUTINE FA14_RANDOM_REAL ( seed, positive, random )

!  Real random number in the range [0, 1] ( if positive is .TRUE. )
!  or [-1, 1] ( if positive is .FALSE. )

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

         TYPE (FA14_seed), INTENT( INOUT ) :: seed
         LOGICAL, INTENT( IN ) :: positive
         REAL ( KIND = working ), INTENT( OUT ) :: random

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

         INTEGER :: fhi, k, leftlo, xalo, xhi
         REAL ( KIND = working ) :: x
         REAL ( KIND = working ), PARAMETER ::                                &
                one = 1.0_working, two = 2.0_working, rb16 = two ** 16,       &
                big  = one / ( two ** 31 - one ), big2 = two * big

         xhi = seed%ix / b16

! Get 16 lo bits of seed%ix and form lo product

         xalo = ( seed%ix - xhi * b16 ) * a

!  Get 15 hi order bits of lo product

         leftlo = xalo / b16

!  Form the 31 highest bits of full product

         fhi = xhi * a + leftlo

!  Get overflopast 31st bit of full product

         k = fhi / b15

!  Assemble all the parts and presubtract P. The parentheses are essential

         seed%ix = ( ( (xalo-leftlo*b16) - p ) + (fhi-k*b15)*b16 ) + k

!  Add p back in if neccessary

         IF ( seed%ix < 0 ) seed%ix = seed%ix + p

!  Multiply by 1/(2**31-1)

         xhi = seed%ix / b16
         x = FLOAT( xhi ) * rb16 + FLOAT( seed%ix - xhi * b16 )
         IF ( positive ) THEN
            random = x * big
         ELSE
            random = x * big2 - one
         END IF

         END SUBROUTINE FA14_RANDOM_REAL

!-*-*-*-*  F A 1 4 _ R A N D O M _ I N T E G E R   S U B R O U T I N E  *-*-

         SUBROUTINE FA14_RANDOM_INTEGER ( seed, n, random )

!  Integer random number in the range [1,n] if n > 1.
!  Otherwise, the value n is returned

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

         TYPE (FA14_seed), INTENT( INOUT ) :: seed
         INTEGER, INTENT( IN ) :: n
         INTEGER, INTENT( OUT ) :: random

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

         INTEGER :: be1, be2, c, d, f, fhi, g, k, leftlo
         INTEGER :: mhi, mlo, mu, nu, xalo, xhi, xlo

         IF ( n > 1 ) THEN

            xhi = seed%ix / b16

!  Get 16 lo bits of seed%ix and form lo product

            xalo = ( seed%ix - xhi * b16 ) * a

!  Get 15 hi order bits of lo product

            leftlo = xalo / b16

!  Form the 31 highest bits of full product

            fhi = xhi * a + leftlo

!  Get overflopast 31st bit of full product

            k = fhi / b15

!  Assemble all the parts and presubtract P. The parentheses are essential

            seed%ix = (((xalo-leftlo*b16) - p) + (fhi-k*b15) * b16) +  k

!  Add p back in if neccessary

            IF ( seed%ix < 0 ) seed%ix = seed%ix + p

!  Multiply by n and divide by 2**31-1 in integer arithmetic.
!  Split seed%ix and n into hi and lo parts

            xhi = seed%ix / b15 ; xlo = seed%ix - b15 * xhi
            mhi = n / b15 ; mlo = n - b15 * mhi

!  Calculate intermediate product and split into hi and lo parts.
!  Presubtract p

            f = ( xhi * mlo - p ) + xlo * mhi

!  f is > 0 if intermediate product would have overflowed

            IF ( f <= 0 ) THEN
               f = f + p ; be1 = f / b15 ; be2 = f - be1 * b15
            ELSE
               f = f - 1 ; be1 = f / b15 ; be2 = f - be1 * b15; be1 = be1 + b16
            ENDIF

!  Form product of lo parts and add in lo part of intermediate product
!  to get lo part of complete product

            g = b15 * be2 + xlo * mlo

!  Represent lo part of full product in base 2**30

            d = g / b30 ; c = xhi / 2

!  Calculate full product divided by 2**30

            f = (( 2 * ( c * mhi - q ) - 1) + mhi * ( xhi - 2 * c )) + d + be1

!  Get full product divided in base 2**31

            IF ( f <= 0 ) THEN
               f = f + p ; nu = f / 2 ; mu = f - nu * 2
            ELSE
               f = f - 1 ; nu = f / 2 ; mu = f - 2 * nu ; nu = nu + b30
            ENDIF

!  Calculate remainder of product divided by 2**31

            f = ( b30 * mu - p ) + nu + ( g - b30 * d )
            random = nu + 1

!  Add one if remainder is not < 2**31-1

            IF ( f >= 0 ) random = random + 1
         ELSE

!  If n is less than or equal to 1, set random to n.

            random = n
         END IF

         END SUBROUTINE FA14_RANDOM_INTEGER

!-*-*-*-*-*-*-  F A 1 4 _ G E T _ S E E D  S U B R O U T I N E  *-*-*-*-*-*-*-

         SUBROUTINE FA14_GET_SEED ( seed, value )

!  Determine the current word generator.

         TYPE (FA14_seed), INTENT( IN ) :: seed
         INTEGER, INTENT( OUT ) :: value

         value = seed%ix

         END SUBROUTINE FA14_GET_SEED

!-*-*-*-*-*-*-  F A 1 4 _ S E T _ S E E D   S U B R O U T I N E  *-*-*-*-*-*-

         SUBROUTINE FA14_SET_SEED ( seed, value )

!  Reset the word generator to value if value lies in the
!  interval [1, 2**31 - 1]. More generally, value is set
!  to ( value - 1 ) mod (2**31 -1) + 1

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

         TYPE (FA14_seed), INTENT( OUT ) :: seed
         INTEGER, INTENT( IN ) :: value

         seed%ix = MOD( value - 1, p ) + 1

         END SUBROUTINE FA14_SET_SEED

!-*-*-*-*-*-*-  F A 1 4 _ INITIALIZE  S U B R O U T I N E  *-*-*-*-*-*-*-

         SUBROUTINE FA14_INITIALIZE ( seed )

!  Set the word generator to its default value.

         TYPE (FA14_seed), INTENT( OUT ) :: seed

         seed%ix = 65535

         END SUBROUTINE FA14_INITIALIZE

      END MODULE HSL_FA14_DOUBLE

module hsl_maxflow
!  use hsl_zd11_double
!  use hsl_mc70_double
implicit none

type bandgraph
!    graph object of wide separator and its edges
  integer ::  nvtx  ! number of vertices
  integer ::  nedge ! number of edges
  integer, allocatable :: edges(:,:)  ! edge array
  integer, allocatable :: isAdjToSource(:)
!       isAdjToSource(u,1)  = 1 --> u is adjacent to the source
!       isAdjToSource(u,1) /= 1 --> u is NOT adjacent to the source
  integer, allocatable :: isAdjToSink(:)
!       isAdjToSink(u,1)  = 1 --> u is adjacent to the sink
!       isAdjToSink(u,1) /= 1 --> u is NOT adjacent to the sink
  integer, allocatable ::  vwts(:) ! vertex weights (optional),
!       if not present, then unit vertex weights assumed
end type bandgraph

type network
  integer ::  nnode ! number of nodes in the network
  integer ::  narc  ! number of arcs in the network
  integer ::  source! source node = 1
  integer ::  sink  ! sink node = nnode
  integer, allocatable :: inheads(:) ! for each node u
!                                      first incoming arc for u
  integer, allocatable :: outheads(:) ! for each node u
!                                       first outgoing arc for u
  integer, allocatable :: firsts(:) ! for each arc e
!                                     first node for arc e
  integer, allocatable :: seconds(:) ! for each arc e
!                                      second node for arc e
  integer, allocatable :: capacities(:) ! for each arc e
!                                         capacity of arc e
  integer, allocatable :: flows(:) ! for each arc e
!                                    flow through arc e
  integer, allocatable :: nextin(:) ! for each arc e
!                                     next incoming arc 
  integer, allocatable :: nextout(:) ! for each arc e
!                                      next outgoing arc
end type network

integer, parameter :: wp = kind(0.0d0)

contains

! ---------------------------------------------------
! mc70_maxflow
! ---------------------------------------------------
! Given a partition, get better partition using maxflow algorithm
subroutine mc70_maxflow(a_n,a_ne,a_ptr,a_row,a_weight,sumweight,a_n1,&
               a_n2,a_weight_1,a_weight_2,a_weight_sep,partition, &
               alpha,beta,msglvl,stats,cost)
implicit none
! Input matrix: a_n, a_ne, a_ptr, a_row
integer, intent(in) :: a_n ! order of matrix
integer, intent(in) :: a_ne ! number of entries in matrix (lower and
             ! upper triangle)
integer, allocatable, intent(in) :: a_ptr(:) ! On input, a_ptr(i) contains 
             ! position in a_row that entries for column i start. 
integer, allocatable, intent(in) :: a_row(:) ! On input, a_row contains row 
             ! indices of the nonzero entries. Diagonal entries have been 
             ! removed and the matrix expanded.
! At the moment weights are not used at all             
integer, intent(in) :: a_weight(a_n) ! On input, a_weight(i) contains 
             ! the weight of column i 
integer, intent(in) :: sumweight ! Sum of weights in a_weight
! Data on partition a_n1, a_n2, partition ... will be updated
integer, intent(inout) :: a_n1 ! Size of partition 1 (ie B)
integer, intent(inout) :: a_n2 ! Size of partition 2 (ie W)
integer, intent(inout) :: a_weight_1,a_weight_2,a_weight_sep ! Weighted 
             ! size of partitions and separator
integer, intent(inout) :: partition(a_n) !First a_n1 entries contain
             ! list of (local) indices in partition 1; next a_n2 entries  
             ! contain list of (local) entries in partition 2; entries in 
             ! separator are listed at the end. This is updated to the new 
             ! partition.
! Parameters alpha (for balance) and beta (for cost function)
real(wp), intent(in) :: alpha,beta
integer, intent(in) :: msglvl
! output --
!    stats[1] -- weight of vertices in S
!    stats[2] -- weight of vertices in B
!    stats[3] -- weight of vertices in W
!    stats[4] -- weight of edges in A_{S,S}
!    stats[5] -- weight of edges in A_{S,B}
!    stats[6] -- weight of edges in A_{S,W}
!    stats[7] -- weight of edges in A_{B,B}
!    stats[8] -- weight of edges in A_{B,W}
!    cost     -- cost of new partition
integer, intent(out) :: stats(8)
real(wp), intent(out) :: cost

!       TYPE(mc70_control), intent(in) :: control

type(bandgraph) :: bandg

integer, allocatable :: map(:),mapL(:),mapR(:),graph_S_S(:,:)
integer, allocatable :: dmapL(:), dmapR(:)
integer, allocatable :: sep_map(:)
integer :: a_ns, n_S, n_S_S, i, istart_s, j, j1, j2, jj, k, lp, &
                   nedge, statsR(9), statsL(9)
integer i1,i2,i0           
!       integer, parameter :: wp = kind(0.0d0)
real(wp) :: costR, costL

lp = 6

! Define bandgraph for maxflow solver

! Matlab code for first part
!%
!% form the wide separator graph for the max flow solver
!%
!   Swide = find(mapWide == 0) ; S = Swide ; n_S = length(S) ; 
!   Bwide = find(mapWide == 1) ; B = Bwide ; n_B = length(B) ; 
!   Wwide = find(mapWide == 2) ; W = Wwide ; n_W = length(W) ; 
!
!   A_S_S = A(S,S) ; A_S_B = A(S,B) ; A_S_W = A(S,W) ;
!   [rows_S_S, cols_S_S, ents_S_S] = find(A(S,S)) ;
!   n_S_S = length(rows_S_S) ;
!   counts = spones(A_S_B)*ones(n_B,1) ;
!   idx = find(counts > 0) ;
!   isAdjToSource = zeros(n_S,1) ;
!   isAdjToSource(idx) = 1 ;
!   counts = spones(A_S_W)*ones(n_W,1) ;
!   idx = find(counts > 0) ;
!   isAdjToSink = zeros(n_S,1) ;
!   isAdjToSink(idx) = 1 ;
!   graph_S_S.nvtx = n_S ;
!   graph_S_S.nedge = n_S_S ;
!   graph_S_S.edges = [rows_S_S cols_S_S] ;
!   graph_S_S.isAdjToSource = isAdjToSource ;
!   graph_S_S.isAdjToSink = isAdjToSink ;
!
! Number vertices in separator
a_ns = a_n - a_n1 - a_n2
n_s  = a_ns

! Set up map array to define in what partition each vertex lies
allocate (map(a_n))
do i = 1,a_n1
  map(partition(i)) = 1
enddo
do i = a_n1+1,a_n1+a_n2
  map(partition(i)) = 2
enddo
do i = a_n1+a_n2+1,a_n
  map(partition(i)) = 0
enddo

! Source is associated with partition B (size a_n1)
! Sink is associated with partition W   (size a_n2)

! Set number of vertices in bandgraph (separator)
bandg%nvtx = a_ns

! Allocate and determind mapping of global variables of matrix to separator set
allocate(sep_map(a_n))
do k = 1,bandg%nvtx
  i = partition(a_n1+a_n2+k)
  sep_map(i) = k
enddo

! Allocate maximum space for bandgraph edges
n_S_S = min(a_ns*a_ns,a_ne)
allocate(graph_S_S(n_S_S,2))

! Allocate arrays for bandgraph
allocate(bandg%isAdjToSource(a_ns),bandg%isAdjToSink(a_ns))
! Initialize arrays to say no vertex is in them
bandg%isAdjToSource = 0
bandg%isAdjToSink   = 0
nedge = 0

! Run through nodes in separator S and generate bandgraph
do k = 1,bandg%nvtx
  i = partition(a_n1+a_n2+k)
  j1 = a_ptr(i)
  if (i == a_n) then 
    j2 = a_ne  
  else 
    j2 = a_ptr(i+1)-1
  endif
! Run through vertices connected to vertex i seeing what partition they are in
  do jj = j1,j2
    j = a_row(jj)
! Find out in which partition node j lies using map array
    if (map(j) == 1) then
! If in partition B add vertex k to AdjToSource
       bandg%isAdjToSource(k) = 1
    endif
    if (map(j) == 2) then
! If in partition W add vertex k to AdjToSink
       bandg%isAdjToSink(k) = 1
    endif
    if (map(j) == 0) then
! If in separator add edge to bandg%edges accumulating number of edges
      nedge = nedge + 1
! To emulate matlab code      
      graph_S_S(nedge,2) = sep_map(i)
      graph_S_S(nedge,1) = sep_map(j)
    endif
  enddo
enddo
deallocate(sep_map)

bandg%nedge = nedge
allocate(bandg%edges(nedge,2))
bandg%edges(:,1) = graph_S_S(1:nedge,1)
bandg%edges(:,2) = graph_S_S(1:nedge,2)

deallocate(graph_S_S)

if (msglvl > 2) then
  write(lp,*) 'Bandgraph'
  write(lp,'(A,I4)') 'Number of vertices',bandg%nvtx
  write(lp,'(A,I4)') 'Number of edges   ',bandg%nedge
  write(lp,*) 'Edges'
  write(lp,'(10I4)') (bandg%edges(i,1),bandg%edges(i,2),i=1,bandg%nedge)
  write(lp,*) 'isAdjToSource'
  write(lp,'(10I4)') (bandg%isAdjToSource(i),i=1,bandg%nvtx)
  write(lp,*) 'isAdjToSink'
  write(lp,'(10I4)') (bandg%isAdjToSink(i),i=1,bandg%nvtx)
endif

! solve a max flow problem to find the two new maps
!
! [dmapL, dmapR] = solvemaxflow(graph_S_S, msglvl) ;
! dmapL and dmapR are allocated within solvemaxflow

call solvemaxflow(bandg,msglvl,dmapL,dmapR)

deallocate(bandg%isAdjToSource,bandg%isAdjToSink,bandg%edges)

if (msglvl > 2) then
  write(lp,*) 'dmapL ...'
  write(lp,'(10I4)') dmapL

  write(lp,*) 'dmapR ...'
  write(lp,'(10I4)') dmapR
endif

allocate (mapL(a_n),mapR(a_n))
mapL = map
mapR = map
istart_s = a_n1+a_n2
do i = 1,n_s
  mapL(partition(istart_S +i)) = dmapL(i)
  mapR(partition(istart_S +i)) = dmapR(i)
enddo
!  mapL = mapWide ; mapL(Swide) = dmapL ;
!  mapR = mapWide ; mapR(Swide) = dmapR ;
deallocate(dmapL,dmapR)
if (msglvl > 2) then
  write(lp,*) 'mapL ...'
  write(lp,'(10I4)') mapL
  write(lp,*) 'mapR ...'
  write(lp,'(10I4)') mapR
endif

! Use evaluation function to choose best partition from among these two
! Pass all of matrix A
call evalBSW(a_n,a_ptr,a_row, mapL, alpha, beta, statsL, costL)
call evalBSW(a_n,a_ptr,a_row, mapR, alpha, beta, statsR, costR)
!  statsL = evalBSW(A, mapL, params(1), params(2), msglvl) ;
!  statsR = evalBSW(A, mapR, params(1), params(2), msglvl) ;
if (msglvl > 0) then
  write(lp,'(A)') 'After evaluation of two partitions'
  write(lp,'(A)') 'left partition'
  if (statsL(9) == 1) then
    write(lp,'(A)') 'acceptable'
  else
    write(lp,'(A)') 'NOT acceptable'
  endif
  write(lp,'(A,G10.2)') 'cost', costL
  write(lp,'(A,I8,A,3I8,A)') '|S|', statsL(1), 'edges = [', &
            statsL(4), statsL(5), statsL(6),']'
  write(lp,'(A,I8,A,3I8,A)') '|B|', statsL(2), 'edges = [', &
            statsL(5), statsL(7), 0, ']'
  write(lp,'(A,I8,A,3I8,A)') '|W|', statsL(3), 'edges = [', &
            statsL(6), 0, statsL(8), ']'
   write(lp,'(A)') 'right partition'
   if (statsR(9) == 1) then
     write(lp,'(A)') 'acceptable'
   else
     write(lp,'(A)') 'NOT acceptable'
   endif
   write(lp,'(A,G10.2)') 'cost', costR
   write(lp,'(A,I8,A,3I8,A)') '|S|', statsR(1), 'edges = [', &
            statsR(4), statsR(5), statsR(6),']'
   write(lp,'(A,I8,A,3I8,A)') '|B|', statsR(2), 'edges = [', &
            statsR(5), statsR(7), 0, ']'
   write(lp,'(A,I8,A,3I8,A)') '|W|', statsR(3), 'edges = [', &
            statsR(6), 0, statsR(8), ']'
  endif

!
! Find the better of the two partitions
!

if (statsL(9) == 1 .AND. statsR(9) == 1) then
  if (msglvl > 0) write(lp,'(A)') 'both maps are acceptable'
  if (costL <= costR) then
    map = mapL
    stats = statsL(1:8)
    cost  = costL
    if (msglvl > 0) write(lp,'(A)') 'left map accepted'
  else
    map = mapR
    stats = statsR(1:8)
    cost  = costR
    if (msglvl > 0) write(lp,'(A)') 'right map accepted'
  endif
elseif (statsL(9) == 1) then
  map = mapL
  stats = statsL(1:8)
  cost  = costL
  if (msglvl > 0)  &
          write(lp,'(A)') 'right map NOT acceptable, left map accepted'
elseif (statsR(9) == 1) then
  map = mapR
  stats = statsR(1:8)
  cost  = costR
  if (msglvl > 0) &
       write(lp,'(A)') 'left map NOT acceptable, right map accepted'
else 
  if (msglvl > 0) write(lp,'(A)') 'NEITHER map acceptable'
  if (costL <= costR) then
    map = mapL
    stats = statsL(1:8)
    cost  = costL
    if (msglvl > 0) write(lp,'(A)') 'left map accepted'
  else
    map = mapR
    stats = statsR(1:8)
    cost  = costR
    if (msglvl > 0) write(lp,'(A)') 'right map accepted'
  endif
endif

deallocate (mapL,mapR)
!
! trim this best separator
!
!mapTrim, statsTrim] = trimpartition(A, [], mapWide, params, 1) ;
!call trimpartition(A, map, alpha, beta, msglvl, statsTrim, costTrim)

! Now update partition
a_n1 = 0
a_n2 = 0
a_ns = 0

do i =1,a_n
  if (map(i) == 1) then
    a_n1 = a_n1 + 1
    cycle
  endif
  if (map(i) == 2) then
    a_n2 = a_n2 + 1
    cycle
  endif
  if (map(i) == 0) then
    a_ns = a_ns + 1
    cycle
  endif
enddo

i1 = 1
i2 = a_n1 + 1
i0 = a_n1 + a_n2 + 1
do i =1,a_n
  if (map(i) == 1) then
    partition(i1) = i
    i1 = i1 + 1
    cycle
  endif
  if (map(i) == 2) then
    partition(i2) = i
    i2 = i2 + 1
    cycle
  endif
  if (map(i) == 0) then
    partition(i2) = i
    i0 = i0 + 1
    cycle
  endif
enddo

deallocate(map)

contains
include 'solvemaxflow.f90'
include 'make_network.f90'
include 'evalBSW.f90'
include 'findaugpath.f90'
include 'augmentpath.f90'
include 'findmaxflow.f90'
include 'findmincut.f90'

end subroutine mc70_maxflow

end module hsl_maxflow
