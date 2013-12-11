module hsl_maxflow
!  use hsl_zd11_double
!  use hsl_mc70_double
implicit none

! This is a version of maxflow with no bandgraph and assumption that there
! are no duplicated edges.

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
integer, intent(in) :: a_ptr(:) ! On input, a_ptr(i) contains 
             ! position in a_row that entries for column i start. 
integer, intent(in) :: a_row(:) ! On input, a_row contains row 
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

type(network) :: netw

! Work arrays
!    map,mapL,mapR of length a_n
!    dmapL,dmapR,vwts of length a_ns
!    mark1,mark2,pred,list of length number of nodes in network that
!       is bounded by 2*a_ns+2

! Eventually these arrays will be allocated higher up the chain and
! will be passed as parameters
integer, allocatable :: map(:),mapL(:),mapR(:)
integer, allocatable :: dmapL(:), dmapR(:)
integer, allocatable :: vwts(:)
integer, allocatable :: sedge(:,:)
integer, allocatable :: mark1(:),mark2(:),pred(:),list(:)
real(wp), allocatable :: imb(:)

! Local variables
integer :: a_ns, i, istart_s, j, j1, j2, k, lp, &
                   wtW, wtB, statsR(9), statsL(9)
integer i1,i2,i0,nedge
real(wp) :: costR, costL

lp = 6

!write(9,*) 'Entering maxflow'
!write(0,*) 'Entering maxflow'
allocate (map(a_n),mapL(a_n),mapR(a_n))

!
! Number vertices in separator
a_ns = a_n - a_n1 - a_n2

allocate (dmapL(2*a_ns),dmapR(a_ns),vwts(a_ns))

! Allocate network work arrays.  Length is upper bound
allocate (mark1(2*a_ns+2),mark2(2*a_ns+2),pred(2*a_ns+2),list(2*a_ns+2))
allocate (imb(2*a_ns))

! Set up map array to define in what partition each vertex lies
! At same time set weights for partition (can check with Sue's input)
wtB = 0
do i = 1,a_n1
  k = partition(i)
  map(k) = 1
  wtB = wtB + a_weight(k)
enddo
wtW = 0
do i = a_n1+1,a_n1+a_n2
  k = partition(i)
  map(k) = 2
  wtW = wtW + a_weight(k)
enddo
do i = a_n1+a_n2+1,a_n
  k = partition(i)
  map(k) = 0
  vwts(i-a_n1-a_n2) = a_weight(k)
enddo


! Count edges to get upper bound on size of sedge array
nedge=0
do k= 1,a_ns
  i = partition(a_n1+a_n2+k)
  j1 = a_ptr(i)
  if (i == a_n) then 
    nedge = nedge + a_ne +1 - a_ptr(i) 
  else 
    nedge = nedge + a_ptr(i+1)- a_ptr(i)
  endif
end do
allocate(sedge(nedge,2))

!
! Generate network graph.  The structure for our maxflow algorithm
!
! Work arrays dmapL and dmapR used to hold isadjsource and isadjsink
! Work array  mapL used to hold sep_map
! Source is associated with partition B (size a_n1)
! Sink is associated with partition W   (size a_n2)
!write(9,*) 'Calling mk_network'
!write(0,*) 'Calling mk_network'
call mk_network(a_n,a_ne,a_ptr,a_row,partition,map,a_ns,msglvl,netw,vwts, &
                wtB,wtW,sedge,mapL,dmapL,dmapR,pred,list,mark1,mark2,imb)
!write(9,*) 'Leaving mk_network'
!write(0,*) 'Leaving mk_network'

!
! solve a max flow problem to find the two new maps dmapL and dmapR
!
call solvemaxflow(netw,a_ns,msglvl,dmapL,dmapR,mark1,mark2,pred,list)


if (msglvl > 2) then
  write(lp,*) 'dmapL ...'
  write(lp,'(10I4)') dmapL
  write(lp,*) 'dmapR ...'
  write(lp,'(10I4)') dmapR
endif

mapL = map
mapR = map
istart_s = a_n1+a_n2
do i = 1,a_ns
  mapL(partition(istart_S +i)) = dmapL(i)
  mapR(partition(istart_S +i)) = dmapR(i)
enddo

if (msglvl > 2) then
  write(lp,*) 'mapL ...'
  write(lp,'(10I4)') mapL
  write(lp,*) 'mapR ...'
  write(lp,'(10I4)') mapR
endif

! Use evaluation function to choose best partition from among these two
! Use Sue's weighted code
call evalBSW(a_n,a_ne,a_ptr,a_row, a_weight, mapL, alpha, beta, statsL, costL)
call evalBSW(a_n,a_ne,a_ptr,a_row, a_weight, mapR, alpha, beta, statsR, costR)

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
a_weight_1 = stats(2)
a_weight_2 = stats(3)
a_weight_sep = stats(1)

! Now update partition
! First count number of vertices in each part
a_n1 = 0
a_n2 = 0
a_ns = 0
do i =1,a_n
  if (map(i) == 1) a_n1 = a_n1 + 1
  if (map(i) == 2) a_n2 = a_n2 + 1
  if (map(i) == 0) a_ns = a_ns + 1
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
    partition(i0) = i
    i0 = i0 + 1
    cycle
  endif
enddo

deallocate(map)
deallocate(dmapL,dmapR)
deallocate (mapL,mapR)
deallocate (mark1,mark2,pred,list,vwts)
deallocate (sedge)

contains
include 'solvemaxflow.f90'
include 'findpenalty.f90'
include 'make_network.f90'
include 'evalBSW.f90'
include 'findaugpath.f90'
include 'augmentpath.f90'
include 'findmaxflow.f90'
include 'findmincut.f90'

end subroutine mc70_maxflow

end module hsl_maxflow
