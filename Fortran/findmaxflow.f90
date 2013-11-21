subroutine findmaxflow(netw,msglvl)
! Find a maximum flow through the network

! Matlab call
!function [network, predmtx, flowmtx] ...
!                               = network_findmaxflow ( network, msglvl )

implicit none
integer, intent(in) :: msglvl
type(network), intent(inout) :: netw

integer avail,iarc,lp,nnode,root,sink,source,stats(2),tag

integer, allocatable :: pred(:)

!ISD As list, dist, tags, and deltas are not used have suppressed these
!integer, allocatable :: list(:), dist(:), tags(:), deltas(:), pred(:)
!integer network.n_nodevisit, network.n_arcvisit

lp = 6

if (msglvl > 0) then
   write(lp,'(A)') ''
   write(lp,'(A)') '### inside network_findmaxflow()'
endif

nnode = netw%nnode
!ISD narc not used
!narc  = netw%narc
source = netw%source
sink = netw%sink
!firsts = network.firsts ; seconds = network.seconds ;
!capacities = network.capacities ; flows = network.flows ;
!inheads = network.inheads ; outheads = network.outheads ;
!nextin = network.nextin ; nextout = network.nextout ;
!ISD  These varibales/structure do not seem to be used so suppressed
!network.n_nodevisit = 0
!network.n_arcvisit = 0

!ISD Shouldn't this be where flows are initialized to zero?

!
! allocate and initialize workspace
!
!allocate (list(nnode),dist(nnode),tags(nnode),deltas(nnode),pred(nnode))
!ISD pred is allocated in findaugpath
!allocate (pred(nnode))
!list   = 0
!dist   = 0
!tags   = 0
!deltas = 0
!pred   = 0

if (msglvl > 0) then
    write(lp,'(I8,A,I8,A)') 0, 'nodes visited', 0, 'arcs visited'
endif

tag = 0
iarc = netw%outheads(source)
!ISD I have supppressed thes etwo variables as they are not used anywhere
!predmtx = [] ;
!flowmtx = [] ;
!while iarc ~= -1
do
  if (iarc == -1) exit

  root = netw%seconds(iarc)
  if (msglvl > 0) then
      write(lp,'(A)') ''
      write(lp,'(A,I8,A,I8,A,I4,A,I4,A,I8)')                  &
      'top of outer loop, tag',tag,', iarc =',iarc, ', (',    &
      netw%firsts(iarc),',',netw%seconds(iarc),'), root',netw%seconds(iarc)
  endif
  do
    avail = netw%capacities(iarc) - netw%flows(iarc)
    if (msglvl > 0) &
       write(lp,'(A,I8,A,I4,A,I4,A,I8,A,I8,A,I8)') 'iarc',iarc,   &
             ', (',netw%firsts(iarc),',',netw%seconds(iarc),'), cap',  &
             netw%capacities(iarc),', flow',netw%flows(iarc),', avail',avail
    if (avail > 0) then
! Set in findaugpath
!       pred(root) = iarc
      if (msglvl > 0) write(lp,'(A)') '    find augmenting path using this arc'
!
!        use BFS to find path
!        [avail, pred, stats] = network_findaugpath(network, iarc, ...
!                                                   tag, msglvl-1) ;
      call findaugpath(netw,iarc,msglvl,avail,pred,stats)
!
!        use DFS to find path
!        [avail, pred, stats] = network_findaugpath2(network, iarc, ...
!                                                    tag, msglvl-1) ;
!        predmtx = [ predmtx pred ] ;
!        flowmtx = [ flowmtx netw%flows ] ;
      if (msglvl > 0) &
           write(lp,'(I8,A,I8,A)') stats(1),'nodes visited,', stats(2),  &
               'arcs visited'
      if ((avail == 0) .OR. (pred(sink) == 0)) then
        if (msglvl > 0) write(lp,'(A)') 'no augmenting path found, exit'
        exit
      endif
      if (msglvl > 0) &
           write(lp,'(A,I8,A)') 'augmenting path found, avail =',avail, &
                ' add flow to network'

!       flows = network_augmentpath(network, avail, pred, msglvl-1) ;
      call augmentpath(netw,avail,pred,msglvl)

!ISD Seems a redundant statement or am I missing some matlab subtlty
!         network.flows = flows ;
!         flows = network.flows ;
! HST: add exit for avail==0 to stop infinite loop
    else
       exit
    endif
!     keyboard ;
  enddo
  if (msglvl > 0) &
    write(lp,'(A,I8)') 'unable to proceed further down arc', iarc
  iarc = netw%nextout(iarc)
  tag = tag + 1
  if (msglvl > 0) &
      write(lp,'(A,I8,A,I8)') 'moving to arc', iarc, ', tag', tag
enddo

! HST: check whether pred is allocated because it might not be
if (allocated(pred)) then
    deallocate(pred)
end if

if (msglvl > 0) then
   write(lp,'(A)') '### leaving network_findmaxflow()'
   write(lp,'(A)') ''
endif

end subroutine findmaxflow
