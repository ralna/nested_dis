subroutine findmincut(netw,msglvl,mark1,mark2)

! findmincut -- find one or two mincuts,
!               one nearest the source,
!               one nearest the sink

! matlab call
!
!function [network, mark1, mark2] = network_findmincut ( network, msglvl )

!   netw    -- network object
!   msglvl  -- message level


implicit none
integer, intent(in) :: msglvl
integer, allocatable, intent(out) :: mark1(:),mark2(:)
type(network), intent(inout) :: netw

integer iarc,last,lp,nnode,now,sink,source,x,z

integer, allocatable :: list(:)

lp = 6

nnode = netw%nnode
!narc = netw.narc
source = netw%source
sink = netw%sink
!firsts = network.firsts ; seconds = network.seconds ;
!capacities = network.capacities ; flows = network.flows ;
!inheads = network.inheads ; outheads = network.outheads ;
!nextin = network.nextin ; nextout = network.nextout ;

if (msglvl > 0) then
   write(lp,'(A)')  ''
   write(lp,'(A)')  '### inside network_findmincut()'
endif

if (msglvl > 3) then
  write(lp,*) 'netw%outheads'
  write(lp,'(10I4)') netw%outheads
  write(lp,*) 'netw%seconds'
  write(lp,'(10I4)') netw%seconds
  write(lp,*) 'netw%flows'
  write(lp,'(10I4)') netw%flows
  write(lp,*) 'netw%capacities'
  write(lp,'(10I4)') netw%capacities
  write(lp,*) 'netw%firsts'
  write(lp,'(10I4)') netw%firsts
  write(lp,*) 'netw%inheads'
  write(lp,'(10I4)') netw%inheads
  write(lp,*) 'netw%nextout'
  write(lp,'(10I4)') netw%nextout
  write(lp,*) 'netw%nextin'
  write(lp,'(10I4)') netw%nextin
endif  

!
! breadth first traversal from source
!
if (msglvl > 0) write(lp,'(A)')  'breadth first traversal from the source'

allocate (mark1(nnode),mark2(nnode))
mark1 = 2
mark1(source) = 1

allocate(list(nnode))
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
  do
  if (iarc == -1) exit
    z = netw%seconds(iarc)
    if (msglvl > 1) write(lp,'(A,I4)') '(1) node', z
    if (mark1(z) == 1) then
       if (msglvl > 1) write(lp,'(A,I4)')  'already marked', z
    else 
      if (msglvl > 1) &
          write(lp,'(A,I4,A,I4)')  'not marked, flow', netw%flows(iarc),  &
                                   ', cap', netw%capacities(iarc)
      if (netw%flows(iarc) < netw%capacities(iarc)) then
        last = last + 1
        list(last) = z
        mark1(z) = 1
        if (msglvl > 1) write(lp,'(A)')  'marking and putting on list'
      endif
    endif
    iarc = netw%nextout(iarc)
  enddo
  iarc = netw%inheads(x)
! while iarc ~= -1
  do
  if (iarc == -1) exit
    z = netw%firsts(iarc)
    if (msglvl > 1) write(lp,'(A,I4)') '(2) node', z
    if (mark1(z) == 1) then
      if (msglvl > 1) write(lp,'(A,I4)')  'already marked', z
    else 
      if (msglvl > 1) write(lp,'(A,I4)') 'not marked, flow', netw%flows(iarc)
      if (netw%flows(iarc) > 0) then
        last = last + 1
        list(last) = z
        mark1(z) = 1
        if (msglvl > 1) write(lp,'(A)')  'marking and putting on list'
      endif
    endif
    iarc = netw%nextin(iarc)
  enddo
enddo
!
! breadth first traversal from sink
!
if (msglvl > 0) write(lp,'(A)')  'breadth first traversal from the sink'

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
  do
    if (iarc == -1) exit
    z = netw%seconds(iarc)
    if (msglvl > 1) write(lp,'(A,I4)') '(3) node', z
    if (mark2(z) == 2) then
      if (msglvl > 1) write(lp,'(A,I4)')  'already marked', z
    else 
      if (msglvl > 1)  &
          write(lp,'(A,I4,A,I4)')  'not marked, flow', netw%flows(iarc),  &
                                   ', cap', netw%capacities(iarc)
      if (netw%flows(iarc) > 0) then
        last = last + 1
        list(last) = z
        mark2(z) = 2
        if (msglvl > 1) write(lp,'(A)')  'marking and putting on list'
      endif
    endif
    iarc = netw%nextout(iarc)
  enddo
  iarc = netw%inheads(x)
! while iarc ~= -1
  do
    if (iarc == -1) exit
    z = netw%firsts(iarc)
    if (msglvl > 1) write(lp,'(A,I4)') '(4) node', z
    if (mark2(z) == 2) then
      if (msglvl > 1) write(lp,'(A,I4)')  'already marked', z
    else 
      if (msglvl > 1) write(lp,'(A,I4)')  'not marked, flow', netw%flows(iarc)
      if (netw%flows(iarc) < netw%capacities(iarc)) then
        last = last + 1
        list(last) = z
        mark2(z) = 2
        if (msglvl > 1) write(lp,'(A)')  'marking and putting on list'
      endif
    endif
    iarc = netw%nextin(iarc)
  enddo
enddo

if (msglvl > 0) write(lp,'(A/)')  '### leaving network_findmincut()'

end subroutine findmincut
