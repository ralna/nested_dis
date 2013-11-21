subroutine make_network(bandg,msglvl,netw,weights)
! Given a band graph, create and return a network structure


! Matlab call was:
! function network = mknetwork ( bandgraph, msglvl )

implicit none
type(bandgraph), intent(in) :: bandg
integer, intent(in) :: msglvl
type(network), intent(out) :: netw
integer, optional :: weights

integer i,iarc,ii,jj,lp,narc1,narc2,narc3,narc4,narc,nnode,nvtx,nedge,u,v
integer sumweight

integer, allocatable :: data(:,:)
integer, allocatable :: sort(:),index(:)

EXTERNAL KB07AI

lp = 6

if (msglvl .gt.0) write(lp,'(/A)') '### inside mknetwork()'

sumweight = sum(bandg%vwts)

nvtx = bandg%nvtx
nedge = bandg%nedge
!%edges = bandgraph.edges ;
!%isAdjToSource = bandgraph.isAdjToSource ;
!%isAdjToSink = bandgraph.isAdjToSink ;

!ISD You don't use local variable later so no need to set it here
!if (present (weights)) then
!   vwts = bandg%vwts ;
!else
!   vwts = 1 ;
!endif

! HST: If nedge=0, KB07 prints out messaged for this case so need to suppress
! the sort
if (nedge.ne. 0) then
  allocate (sort(nedge),index(nedge))
  do i=1,nedge
    sort(i) = bandg%edges(i,1)*nvtx + bandg%edges(i,2)
  enddo
  CALL KB07AI(sort,nedge,index)

  deallocate(sort)

  allocate (data(nedge,2))

!
! sort and compress edges
!
! Initialize data to edges
  do i = 1,nedge
    data(i,1) = bandg%edges(index(i),1)
    data(i,2) = bandg%edges(index(i),2)
  enddo

  deallocate(index)

! data = bandg%edges

!  write(lp,'(A)') 'data on bandgraph'
!  write(lp,'(A,I8)') 'nedges',nedge
!  write(lp,'(A,I8)') 'vertices',nvtx
!  write(lp,'(A)') 'edges'
! write(lp,'(10I8)') data
! Sort edges into ascending order
! data = sortrows(edges, [1,2]) ;
! CALL KB05AI(data(1,2),nedge)
! CALL KB05AI(data(1,1),nedge)

! Remove duplicates
  ii = 1
  jj = 1
  do ii = 2,nedge
    if (data(ii,1) .NE. data(ii-1,1) .OR. data(ii,2) .NE. data(ii-1,2)) then
      jj = jj + 1
      data(jj,1) = data(ii,1)
      data(jj,2) = data(ii,2)
    endif
  enddo

  nedge  = jj
else
  allocate (data(nedge,2))  
end if

if (msglvl .gt.0) write(lp,'(A,I8,A,I8)') 'nvtx', nvtx, ', nedge', nedge
!   write(lp,'(A)')    %d vertices adjacent to source', ... ...
!         length(find(isAdjToSource == 1))) ;
!   write(lp,'(A)')    %d vertices adjacent to sink', ... ...
!         length(find(isAdjToSink == 1))) ;

!
! count the number of arcs
!
!ISD Not needed
! For nodes in separator, set isadj as
! 1 if only connected to source
! 2 if only connected to sink
! 3 if connected to source and sink
! isadj = 0
! narc1 is number of vertices in separator and is equal to the number of
! added edges u- to u+
narc1 = nvtx
! narc2 is number of vertices in separator connected to source
narc2 = 0
! narc3 is number of vertices in separator connected to sink
narc3 = 0
!bndB = find(isAdjToSource) ;
!bndW = find(isAdjToSink) ;
!isAdj(bndB) = isAdj(bndB) + 1 ;
!isAdj(bndW) = isAdj(bndW) + 2 ;
!narc2 = length(find(isAdjToSource == 1)) ;
!narc3 = length(find(isAdjToSink == 1)) ;
do i = 1,nvtx
  if (bandg%isAdjToSource(i) == 1) then 
!   isadj(i) = isadj(i) + 1
    narc2 = narc2 + 1
  endif
  if (bandg%isAdjToSink(i) == 1)   then
!   isadj(i) = isadj(i) + 2
    narc3 = narc3 + 1
  endif
enddo

narc4 = 0
do ii = 1,nedge
  u = data(ii,1)  
  v = data(ii,2)
! --- ignore self edges ---
  if (u == v) cycle 
! --- omit edges with essential vertices ---
  if (bandg%isAdjToSink(u) == 1 .and. bandg%isAdjToSource(u) ==1) cycle
  if (bandg%isAdjToSink(v) == 1 .and. bandg%isAdjToSource(v) ==1) cycle
! --- omit pairs both adjacent to source ---
!  rc1 = bitand(isAdj(u), 1) ; rc2 = bitand(isAdj(v), 1) ;
!  rc3 = bitand(rc1, rc2) ;
!  if rc3 == 1
!     continue ;
!  end
  if ((bandg%isAdjToSource(u) == 1) .and. (bandg%isAdjToSource(v) == 1)) cycle
! --- omit pairs both adjacent to sink ---
!  rc1 = bitand(isAdj(u), 2) ; rc2 = bitand(isAdj(v), 2) ;
!  rc3 = bitand(rc1, rc2) ;
!  if rc3 == 2
!     continue ;
!  end
  if ((bandg%isAdjToSink(u) == 1) .and. (bandg%isAdjToSink(v) == 1)) cycle
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
!  fprintf('\n    %d internal arcs', narc1) ;
  write(lp,'(I8,A)') narc1,' internal arcs'
  write(lp,'(I8,A)') narc2,' arcs from source'
  write(lp,'(I8,A)') narc3,' arcs from sink'
  write(lp,'(I8,A)') narc4,' edge arcs'
  write(lp,'(I8,A)') narc, ' total arcs'
endif

!
! create the arc arrays
!
allocate(netw%firsts(narc),netw%seconds(narc),netw%capacities(narc))
allocate(netw%flows(narc))
netw%firsts = -1
netw%seconds = -1
!if isfield(bandgraph, 'vwts')
!ISD I just suppressed all this as it doesn't make sense now
!ISD Interesting that bang%vwts and ones(narc,1) are different lengths
!ISD Also you set all capacities in following loops so no need to
! initialize here
!ISD I confess to being a little confused about your setting of capacities
!if (present(weights)) then
!  netw%capacities = bandg%vwts
!else
!  capacities = ones(narc,1)
!  netw%capacities = 1
!endif

!
! (u-,u+) arcs first
!
iarc = 0
do u = 1,nvtx
  iarc = iarc + 1
  netw%firsts(iarc) = 2*u
  netw%seconds(iarc) = 2*u + 1
  netw%capacities(iarc) = bandg%vwts(u)
enddo

if (msglvl .gt.0) write(lp,'(A,I8)') 'after (u-,u+) arcs, iarc = ', iarc

!
! (source,u) arcs
!
do u = 1,nvtx
  if (bandg%isAdjToSource(u) == 1) then
    iarc = iarc + 1
    netw%firsts(iarc) = netw%source
    netw%seconds(iarc) = 2*u
  !  netw%capacities(iarc) = nvtx*nvtx HST: can overflow - max flow in network
  !                      can be bounded from above by sumweight
    netw%capacities(iarc) = sumweight
  endif
enddo

if (msglvl .gt.0) write(lp,'(A,I8)') 'after (source,u-) arcs, iarc = ', iarc

!
! (u,sink) arcs
!
do u = 1,nvtx
  if (msglvl .gt.5) &
    write(lp,'(A,I4,A,I8)') 'bandg%isAdjToSink(',u,')= ',bandg%isAdjToSink(u)
  if (bandg%isAdjToSink(u) == 1) then
    iarc = iarc + 1
    netw%firsts(iarc) = 2*u + 1
    netw%seconds(iarc) = netw%sink
  !  netw%capacities(iarc) = nvtx*nvtx HST: can overflow - max flow in network
  !                      can be bounded from above by sumweight
    netw%capacities(iarc) = sumweight
  endif
enddo

if (msglvl .gt.0) write(lp,'(A,I8)') 'after (u,sink) arcs, iarc = ',iarc

!
! (u+,v-) arcs
!
do ii = 1,nedge
  u = data(ii,1)
  v = data(ii,2)
  ! HST: BUG - logic incorrect - added last 2 lines of logic to correspond to 
  !  logic for computing narc4 
  if ((u /= v) .and.  &
     (bandg%isAdjToSource(u) /= 1 .OR. bandg%isAdjToSource(v) /= 1)   &
    .and. (bandg%isAdjToSink(u)   /= 1 .OR. bandg%isAdjToSink(v)   /= 1)  &
    .and.  (bandg%isAdjToSource(u) /= 1 .OR. bandg%isAdjToSink(u) /= 1)   &
    .and. (bandg%isAdjToSource(v) /= 1 .OR. bandg%isAdjToSink(v) /= 1)   &
      ) then
    iarc = iarc + 1
    netw%firsts(iarc) = 2*u + 1
    netw%seconds(iarc) = 2*v
  !  netw%capacities(iarc) = nvtx*nvtx HST: can overflow - max flow in network
  !                      can be bounded from above by sumweight
    netw%capacities(iarc) = sumweight
  endif
enddo

if (msglvl .gt.0) write(lp,'(A,I8)') 'after (u+,v-) arcs, iarc = ',iarc

!
! fill the heads vectors for in/out edges 
! and the in/out link vectors for the arcs
!

allocate (netw%inheads(nnode),netw%outheads(nnode),netw%nextin(narc), &
          netw%nextout(narc))
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

!
! set the object's fields
!
!netw%nnode = nnode ;
!netw%source = source ;
!netw%sink = sink ;
!netw%inheads = inheads ;
!netw%outheads = outheads ;
!netw%narc = narc ;
!netw%firsts = firsts ;
!netw%seconds = seconds ;
!netw%capacities = capacities ;
!netw%flows = flows ;
!ISD Would this not be initialized in routine that computes flows?
netw%flows = 0
!netw%nextin = nextin ;
!netw%nextout = nextout ;

if (msglvl .gt.0) write(lp,'(A/)') '### leaving mknetwork()'

end subroutine make_network
