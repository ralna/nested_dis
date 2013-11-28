subroutine findmincut(netw,msglvl,mark1,mark2,list)

!  Finds one or two mincuts, one nearest the source, one nearest the sink

! Input parameters

!   netw    -- network object
!   msglvl  -- message level

! Output parameters
!   mark1 --- to identify cut set nearest source
!   mark2 --- to identify cut set nearest sink
!

! Workspace
!   list --- to hold list of nodes being searched

implicit none
integer, intent(in) :: msglvl
integer, intent(out) :: mark1(:),mark2(:),list(:)
type(network), intent(inout) :: netw

! Local variables
integer iarc,last,lp,nnode,now,sink,source,x,z

lp = 6

nnode = netw%nnode
source = netw%source
sink = netw%sink

!
! breadth first traversal from source
!

mark1 = 2
mark1(source) = 1

list = 0
now = 1
last = 1
list(now) = source
!while now <= last
do
  if (now > last) exit
  x = list(now)
  now = now + 1
  iarc = netw%outheads(x)
! while iarc ~= -1
! Run through all arcs starting at node x and putting node at end of arc
! on list if there is spare capacity on arc
  do
  if (iarc == -1) exit
    z = netw%seconds(iarc)
    if (mark1(z) == 1) then
    else 
      if (netw%flows(iarc) < netw%capacities(iarc)) then
        last = last + 1
        list(last) = z
        mark1(z) = 1
      endif
    endif
    iarc = netw%nextout(iarc)
  enddo
  iarc = netw%inheads(x)
! while iarc ~= -1
! Run through all arcs terminating at node x and putting node at start of arc
! on list if there is spare capacity on arc
  do
  if (iarc == -1) exit
    z = netw%firsts(iarc)
    if (mark1(z) == 1) then
    else 
      if (netw%flows(iarc) > 0) then
        last = last + 1
        list(last) = z
        mark1(z) = 1
      endif
    endif
    iarc = netw%nextin(iarc)
  enddo
enddo
!
! breadth first traversal from sink
!

mark2 = 1
mark2(sink) = 2

list = 0
now = 1
last = 1
list(now) = sink
!while now <= last
do 
  if (now > last) exit
  x = list(now)
  now = now + 1
  iarc = netw%outheads(x)
! while iarc ~= -1
! Run through all arcs starting at node x and putting node at end of arc
! on list if there is spare capacity on arc
  do
    if (iarc == -1) exit
    z = netw%seconds(iarc)
    if (mark2(z) == 2) then
    else 
      if (netw%flows(iarc) > 0) then
        last = last + 1
        list(last) = z
        mark2(z) = 2
      endif
    endif
    iarc = netw%nextout(iarc)
  enddo
  iarc = netw%inheads(x)
! while iarc ~= -1
! Run through all arcs terminating at node x and putting node at start of arc
! on list if there is spare capacity on arc
  do
    if (iarc == -1) exit
    z = netw%firsts(iarc)
    if (mark2(z) == 2) then
    else 
      if (netw%flows(iarc) < netw%capacities(iarc)) then
        last = last + 1
        list(last) = z
        mark2(z) = 2
      endif
    endif
    iarc = netw%nextin(iarc)
  enddo
enddo


end subroutine findmincut
