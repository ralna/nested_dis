subroutine mk_network(a_n,a_ne,a_ptr,a_row,partition,map,nvtx,msglvl,netw, &
                      vwts,wtB,wTW,sedge,sep_map,isAdjToSource,isAdjToSink,&
                      count,list,mark1,mark2,imb)
! Create and return a network structure

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

integer, intent(in) :: partition(a_n) !First a_n1 entries contain
             ! list of (local) indices in partition 1; next a_n2 entries  
             ! contain list of (local) entries in partition 2; entries in 
             ! separator are listed at the end.
integer, intent(in) :: map(a_n) !First a_n1 entries contain
integer, intent(in) :: msglvl,nvtx
integer, intent(in) :: vwts(:)
integer, intent(inout) :: wtB,wtW

! Note that we still do the allocations here.  Doing it further up the
! call tree would need to allocate much more storage.
type(network), intent(out) :: netw

integer i,iarc,ii,jj,lp,narc1,narc2,narc3,narc4,narc,nnode,nedge,u,v,wtS

! Work arrays sedge(nedge,2),sep_map(a_n),isAdjToSource(nvtx),
! isAdjToSink(nvtx)
integer :: sedge(:,:)
integer :: sep_map(:)
integer :: isAdjToSource(:),isAdjToSink(:)
real(wp) :: imb(:)
integer :: count(:),mark1(:),mark2(:),list(:)

lp = 0

if (msglvl .gt.0) then
   write(lp,'(/A)') '### inside mknetwork()'
   write(lp,'(A,I8,A,I8)') 'nvtx', nvtx
endif        


! Determine mapping of global variables of matrix to separator set
do k = 1,nvtx
  i = partition(a_n1+a_n2+k)
  sep_map(i) = k
enddo

! We could use a single array although the logic for generating it is 
!   marginally more complicated.
! For nodes in separator, set isadj as
! 1 if only connected to source
! 2 if only connected to sink
! 3 if connected to source and sink
! isadj = 0


! Run through nodes in separator S and generate edges
nedge = 0
do k = 1,nvtx
  i = partition(a_n1+a_n2+k)
  j1 = a_ptr(i)
  if (i == a_n) then 
    j2 = a_ne  
  else 
    j2 = a_ptr(i+1)-1
  endif

! Run through vertices connected to vertex i
  do jj = j1,j2
    j = a_row(jj)
! Find out in which partition node j lies using map array
    if (map(j) == 1) then
! If in partition B add vertex k to AdjToSource
       isAdjToSource(k) = 1
    endif
    if (map(j) == 2) then
! If in partition W add vertex k to AdjToSink
       isAdjToSink(k) = 1
    endif
    if (map(j) == 0) then
! If in separator add edge accumulating number of edges
      nedge = nedge + 1
! Edge has this orientation to emulate matlab code      
      sedge(nedge,2) = sep_map(i)
      sedge(nedge,1) = sep_map(j)
    endif
  enddo

enddo


! narc1 is number of vertices in separator and is equal to the number of
! added edges u- to u+
narc1 = nvtx
! narc2 is number of vertices in separator connected to source
narc2 = 0
! narc3 is number of vertices in separator connected to sink
narc3 = 0

!
! count the number of arcs
! Can't do easily in above loop because of multiple edges from source/sink
! to vertex in the separator set
!
do k = 1,nvtx
  if (isAdjToSource(k) == 1) narc2 = narc2 + 1
  if (isAdjToSink(k) == 1) narc3 = narc3 + 1
enddo

narc4 = 0
do ii = 1,nedge
  u = sedge(ii,1)  
  v = sedge(ii,2)
! --- ignore self edges ---
  if (u == v) cycle 
! --- omit edges with essential vertices ---
  if (isAdjToSink(u) == 1 .and. isAdjToSource(u) ==1) cycle
  if (isAdjToSink(v) == 1 .and. isAdjToSource(v) ==1) cycle
! --- omit pairs both adjacent to source ---
  if ((isAdjToSource(u) == 1) .and. (isAdjToSource(v) == 1)) cycle
! --- omit pairs both adjacent to sink ---
  if ((isAdjToSink(u) == 1) .and. (isAdjToSink(v) == 1)) cycle
! Qualifying arc found
  narc4 = narc4 + 1
enddo

nnode = 2*nvtx + 2
netw%nnode = nnode
narc = narc1 + narc2 + narc3 + narc4
netw%narc = narc
netw%source = 1
netw%sink = 2*nvtx + 2

if (msglvl > 0) then
  write(lp,'(I8,A)') narc1,' internal arcs'
  write(lp,'(I8,A)') narc2,' arcs from source'
  write(lp,'(I8,A)') narc3,' arcs from sink'
  write(lp,'(I8,A)') narc4,' edge arcs'
  write(lp,'(I8,A)') narc, ' total arcs'
endif

!
! create the arc arrays
!

!
! Allocations done here but could be easily moved up the path although
! values very dependent on separator size.
allocate(netw%firsts(narc),netw%seconds(narc),netw%capacities(narc))
allocate(netw%flows(narc))
allocate (netw%inheads(nnode),netw%outheads(nnode),netw%nextin(narc), &
          netw%nextout(narc))

netw%firsts = -1
netw%seconds = -1

!
! (u-,u+) arcs first
!
iarc = 0
do u = 1,nvtx
  iarc = iarc + 1
  netw%firsts(iarc) = 2*u
  netw%seconds(iarc) = 2*u + 1
! We set capacities after computing imbalance penalty
! netw%capacities(iarc) = vwts(u)
enddo

if (msglvl .gt.0) write(lp,'(A,I8)') 'after (u-,u+) arcs, iarc = ', iarc

!
! (source,u) arcs
!
do u = 1,nvtx
  if (isAdjToSource(u) == 1) then
    iarc = iarc + 1
    netw%firsts(iarc) = netw%source
    netw%seconds(iarc) = 2*u
    netw%capacities(iarc) = huge(1)/2
  endif
enddo

if (msglvl .gt.0) write(lp,'(A,I8)') 'after (source,u-) arcs, iarc = ', iarc

!
! (u,sink) arcs
!
do u = 1,nvtx
  if (msglvl .gt.5) &
    write(lp,'(A,I4,A,I8)') 'isAdjToSink(',u,')= ',isAdjToSink(u)
  if (isAdjToSink(u) == 1) then
    iarc = iarc + 1
    netw%firsts(iarc) = 2*u + 1
    netw%seconds(iarc) = netw%sink
    netw%capacities(iarc) = huge(1)/2
  endif
enddo

if (msglvl .gt.0) write(lp,'(A,I8)') 'after (u,sink) arcs, iarc = ',iarc

!
! (u+,v-) arcs
!
do ii = 1,nedge
  u = sedge(ii,1)
  v = sedge(ii,2)
  if ((u /= v) .and.  &
     (isAdjToSource(u) /= 1 .OR. isAdjToSource(v) /= 1) .and.  &
     (isAdjToSink(u)   /= 1 .OR. isAdjToSink(v)   /= 1) .and.  &
     (isAdjToSource(u) /= 1 .OR. isAdjToSink(u) /= 1) .and.  &
     (isAdjToSource(v) /= 1 .OR. isAdjToSink(v) /= 1)  ) then
    iarc = iarc + 1
    netw%firsts(iarc) = 2*u + 1
    netw%seconds(iarc) = 2*v
    netw%capacities(iarc) = huge(1)/2
  endif
enddo

if (msglvl .gt.0) write(lp,'(A,I8)') 'after (u+,v-) arcs, iarc = ',iarc

!
! Generate the head vectors for in/out edges 
! and the in/out link vectors for the arcs
!

netw%inheads = -1
netw%outheads = -1
netw%nextin = -1
netw%nextout = -1
do ii = narc,1,-1
  u = netw%firsts(ii)
  v = netw%seconds(ii)
  if (msglvl .gt.1) write(lp,'(A,I8,A,I8,A,I8,A)')  'ii',ii,'arc (',u,',',v,')'
  netw%nextin(ii) = netw%inheads(v)
  netw%inheads(v) = ii
  netw%nextout(ii) = netw%outheads(u)
  netw%outheads(u) = ii
enddo

! Generate wtS
wtS = 0
do i =1,nvtx
  wtS = wtS + vwts(i)
enddo
 if (msglvl .gt.1) write(lp,*) 'Calling findpenalty'
! compute balance penalties
!call findpenalty(netw,msglvl,nvtx,vwts,wtB,wtW,count,mark1,mark2,list,head, &
!                 imb,penalty)
call findpenalty(netw,msglvl,nvtx,vwts,wtB,wtW,wtS,count,mark1,mark2,list, &
                 isAdjToSource,&
                 imb,isAdjToSink)
 if (msglvl .gt.1) write(lp,*) 'Exiting findpenalty'

! Compute network capacities for separator arcs.
do iarc = 1,nvtx
  netw%capacities(iarc) = isAdjToSink(iarc)*vwts(iarc)
enddo


!ISD Would this not be initialized in routine that computes flows?
netw%flows = 0

if (msglvl .gt.0) write(lp,'(A/)') '### leaving mknetwork()'

end subroutine mk_network
