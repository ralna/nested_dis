subroutine findaugpath(netw,iarc_m,msglvl,avail,pred,stats,  &
                       list,tags,deltas)

! Find an augmenting path starting from arc iarc_m
!                here we use a breadth first search (BFS),
!                in findaugpath2() we use a depth first search (DFS)
!                in findaugpath3() we will use a max distance
!                from source to grow the tree to find a path

!
! input --
!      netw -- network object
!    iarc_m -- label for starting arc (u,v)
!    msglvl -- message level
!

! output --
!     avail -- if nonzero, available flow on augmenting path
!      pred -- tree predecessor vector, size nnode
!              source <-- pred^m(sink) <-- pred(pred(sink)) ...
!                                      <-- pred(sink) <-- sink 
!     stats -- statistics
!              stats(1) = # nodes visited in search
!              stats(2) = # arcs visited
!

! working --
!      list -- stack vector used for depth first search, size nnode
!      tags -- mark vector, size nnode
!    deltas -- increment flow vector, size nnode

implicit none
integer, intent(in) :: iarc_m,msglvl
integer, intent(out) :: avail,stats(2)
type(network), intent(in) :: netw
integer, intent(out) :: pred(:)

integer iarc,last,lp,nnode,now,root,sink,source,v,w

integer :: list(:), tags(:), deltas(:)

integer narc, u, n_nodevisit, n_arcvisit

lp = 6

! As input variable is intent(in) we set a local value that we will update
iarc = iarc_m

nnode = netw%nnode
narc = netw%narc
source = netw%source
sink = netw%sink

stats(1) = 0
stats(2) = 0

list   = 0
!
! Initial working storage and array pred
!
tags   = 0
deltas = 0
pred   = 0

!
! check that (source,u) is an edge
! that is, that iarc as input is an edge from the source node
! 
u = netw%seconds(iarc)
if (netw%firsts(iarc) /= source) then
   write(lp,'(A,I4,A)') 'u', u, 'is not adjacent to source'
   return
endif

!
! check for available capacity
!
avail = netw%capacities(iarc) - netw%flows(iarc)
if (avail == 0) return

!
! Find augmenting path using an alternating tree
!
root = u
now = 1
last = 1
list(1) = root
tags(root) = root
tags(source) = root
tags(sink) = -1
pred(sink) = -1
deltas(root) = avail
pred(root) = iarc
n_nodevisit = 0
n_arcvisit = 0

!while now <= last
do
  if (now > last) exit
  v = list(now)
  now = now + 1
  n_nodevisit = n_nodevisit + 1

  iarc = netw%outheads(v)
!  while iarc ~= -1
  do
    if (iarc == -1) exit
    w = netw%seconds(iarc)
    n_arcvisit = n_arcvisit + 1

    if (tags(w) /= root) then

      if (netw%capacities(iarc) > netw%flows(iarc) ) then
        avail = netw%capacities(iarc) - netw%flows(iarc)

        if (avail > deltas(v)) then
          avail = deltas(v)
        endif
        deltas(w) = avail
        pred(w) = iarc
        tags(w) = root

        if (w == sink) exit

        last = last + 1
        list(last) = w

      endif
    endif
    iarc = netw%nextout(iarc)
  enddo
  if (w == sink) exit


  iarc = netw%inheads(v)
!  while iarc ~= -1
  do
    if (iarc == -1) exit
    w = netw%firsts(iarc)
    n_arcvisit = n_arcvisit + 1
    if (tags(w) /= root) then
      if (netw%flows(iarc) > 0) then
        if (avail > netw%flows(iarc)) then
          avail = netw%flows(iarc)
        endif
        deltas(w) = avail
        pred(w) = iarc
        tags(w) = root
        last = last + 1
        list(last) = w
      endif
    endif
    iarc = netw%nextin(iarc)
  enddo
  if (w == sink) exit
enddo

if (w .ne. sink) avail = 0

stats(1) = n_nodevisit
stats(2) = n_arcvisit


end subroutine findaugpath
