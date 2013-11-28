subroutine findmaxflow(netw,msglvl,pred,list,tags,deltas)
! Find a maximum flow through the network


implicit none
integer, intent(in) :: msglvl
type(network), intent(inout) :: netw

integer avail,iarc,lp,nnode,root,sink,source,stats(2),tag

! Work arrays ... all of length nnode

integer :: pred(:),list(:),tags(:),deltas(:)

lp = 6

nnode = netw%nnode
source = netw%source
sink = netw%sink

! Network flows initialized to zero
netw%flows = 0


tag = 0
iarc = netw%outheads(source)
do
  if (iarc == -1) exit

  root = netw%seconds(iarc)
  do
    avail = netw%capacities(iarc) - netw%flows(iarc)
   ! write(*,*) 'avail',avail, iarc, netw%narc
    if (avail > 0) then
! Set in findaugpath
!
!        use BFS to find path

      call findaugpath(netw,iarc,msglvl,avail,pred,stats, &
                       list,tags,deltas)

!
!        use DFS to find path
!        [avail, pred, stats] = network_findaugpath2(network, iarc, ...
!                                                    tag, msglvl-1) ;
!        predmtx = [ predmtx pred ] ;
!        flowmtx = [ flowmtx netw%flows ] ;
      if ((avail == 0) .OR. (pred(sink) == 0)) exit

!  Update flows
      call augmentpath(netw,avail,pred,msglvl)

    else
       exit
    endif
  enddo

  iarc = netw%nextout(iarc)
  tag = tag + 1

enddo

end subroutine findmaxflow
