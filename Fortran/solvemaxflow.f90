subroutine solvemaxflow(netw,nvtx,msglvl,maps1,maps2)
! SOLVEMAXFLOW -- find two partitions of a wide separator graph
!                 by solving a max flow problem

!Matlab call
!function [mapS1, mapS2] = solvemaxflow ( bandgraph, msglvl )
!

! input --
!
!    network graph

!
! output --
!
!    mapS1[n_S] -- first map from wide separator to {0,1,2} = {S,B,W}
!    mapS2[n_S] -- second map from wide separator to {0,1,2} = {S,B,W}
!


implicit none
type(network), intent(inout) :: netw
integer, intent(in) :: msglvl,nvtx
integer, allocatable, intent(out) :: mapS1(:),mapS2(:)

integer iarc,ii,lp,narc,nnode,u,i
integer, allocatable :: mark1(:),mark2(:)

lp = 6

!
! initialize the network from the graph
!

nnode = netw%nnode
narc = netw%narc

if (msglvl > 3) then
  write(lp,*) 'Network'
  write(lp,'(A,I4)') 'Nodes',nnode
  write(lp,'(A,I4)') 'Edges',narc
  write(lp,*) 'inheads'
  write(lp,'(10I4)') (netw%inheads(i),i=1,nnode)
endif

if (msglvl > 0) then
   write(lp,'(A,I8,A,I8,A)')  'network has', nnode, 'nodes and', narc, 'arcs'
endif

!
! Find maxflow through the network using the Ford-Fulkerson algorithm
![network, predmtx] = network_findmaxflow(network, msglvl) ;
call findmaxflow(netw,msglvl)

!firsts = network.firsts ;
!seconds = network.seconds ;
!capacities = network.capacities ;
!flows = network.flows ;

if (msglvl > 1) then
  write(lp,'(A)') '    first    second  capacity      flow'
  do iarc = 1,narc
    write(lp,'(4I10)') netw%firsts(iarc), netw%seconds(iarc),    &
                       netw%capacities(iarc), netw%flows(iarc)
    if (netw%flows(iarc) > netw%capacities(iarc) .OR. netw%flows(iarc) < 0) &
                                                                     then
      write(lp,'(A)')  'ERROR'
    elseif (netw%flows(iarc) == netw%capacities(iarc)) then
      write(lp,'(A)')  'SATURATED'
    endif
  enddo
endif

!
! Find the two mincuts
!
![network, mark1, mark2] = network_findmincut(network, msglvl) ;
call findmincut(netw,msglvl,mark1,mark2)

if (msglvl >3) then
  write(lp,*) 'mark1'
  write(lp,'(10I2)') mark1
  write(lp,*) 'mark2'
  write(lp,'(10I2)') mark2
endif

!
! load the maps
!

allocate (mapS1(nvtx),mapS2(nvtx))
mapS1 = -1
mapS2 = -1

! Accumulated but not used in code so removed
!w_B1 = 0
!w_S1 = 0
!w_W1 = 0
!w_B2 = 0
!w_S2 = 0
!w_W2 = 0

do ii = 2,nnode-1,2
  u = ii/2
  if (mark1(ii) == 1) then
    if (mark1(ii+1) == 1) then
       mapS1(u) = 1
!      w_B1 = w_B1 + 1
    else
      mapS1(u) = 0
!     w_S1 = w_S1 + 1
    endif
  endif
  if (mark1(ii) == 2) then
    if (mark1(ii+1) == 2) then
      mapS1(u) = 2
!     w_W1 = w_W1 + 1
    else
      mapS1(u) = 0
!     w_S1 = w_S1 + 1
    endif
  endif
enddo

do ii = 2,nnode-1,2
   u = ii/2
   if (mark2(ii) == 1) then
      if (mark2(ii+1) == 1) then
         mapS2(u) = 1
!        w_B2 = w_B2 + 1 ;
      else
         mapS2(u) = 0
!        w_S2 = w_S2 + 1 ;
      endif
   endif
   if (mark2(ii) == 2) then
      if (mark2(ii+1) == 2) then
         mapS2(u) = 2
!        w_W2 = w_W2 + 1 ;
      else
         mapS2(u) = 0
!        w_S2 = w_S2 + 1 ;
      endif
   endif
enddo
if (msglvl > 0) then
   write(lp,'(A)')  '### leaving solvemaxflow()'
   write(lp,'(A)')  ''
endif

!Use mark1 and mark2 to generate maps1 and maps2

end subroutine solvemaxflow

