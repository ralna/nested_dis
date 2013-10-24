subroutine findaugpath(netw,iarc_m,msglvl,avail,pred,stats)

!matlab call
!function [avail, pred, stats] ...
!                    = network_findaugpath ( network, iarc, tag, msglvl )

! Find an augmenting path starting from arc iarc
!                here we use a breadth first search (BFS),
!                in _findaugpath2() we use a depth first search (DFS)
!                in _findaugpath3() we will use a max distance
!                from source to grow the tree to find a path

!
! input --
!   network -- network object
!      iarc -- label for starting arc (u,v)
!    msglvl -- message level
!

! output --
!     avail -- if nonzero, available flow on augmenting path
!              source <-- pred^m(sink) <-- pred(pred(sink)) ...
!                                      <-- pred(sink) <-- sink 
!      pred -- tree predecessor vector, size nnode
!     stats -- statistics
!              stats(1) = # nodes visited in search
!              stats(2) = # arcs visited
!

! working --
!      list -- stack vector used for depth first search
!      tags -- mark vector, size nnode
!    deltas -- increment flow vector, size nnode
!      pred -- predecessor vector, size nnode

implicit none
integer, intent(in) :: iarc_m,msglvl
integer, intent(out) :: avail,stats(2)
type(network), intent(in) :: netw
integer, allocatable, intent(out) :: pred(:)

integer iarc,last,lp,nnode,now,root,sink,source,v,w
!integer network.n_nodevisit, network.n_arcvisit

!ISD As dist is not used have suppressed it
integer, allocatable :: list(:), tags(:), deltas(:)

integer narc, u, n_nodevisit, n_arcvisit

lp = 6

! As input variable is intent(in) we set a local value that we will update
iarc = iarc_m

if (msglvl > 0) then
   write(lp,'(A)') ''
   write(lp,'(A)') '### inside network_findaugpath()'
endif

nnode = netw%nnode
narc = netw%narc
source = netw%source
sink = netw%sink
!firsts = network.firsts ; seconds = network.seconds ;
!capacities = network.capacities ; flows = network.flows ;
!inheads = network.inheads ; outheads = network.outheads ;
!nextin = network.nextin ; nextout = network.nextout ;

stats(1) = 0
stats(2) = 0

!
! allocate working storage
!
allocate (list(nnode),tags(nnode),deltas(nnode))
allocate (pred(nnode))
list   = 0
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
if (avail == 0) then
   return
endif

!
! drop an alternating tree
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

if (msglvl > 0) &
   write(lp,'(A,I4,A,I4)') 'starting tree from arc', iarc, ', root', root

!while now <= last
do
  if (now > last) exit
  v = list(now)
  now = now + 1
  n_nodevisit = n_nodevisit + 1

  if (msglvl > 1) then
      write(lp,'(A,I4,A,I4,A,I4)') 'v =', v,', now =',now,', last',last
      write(lp,'(A,I4)') 'out-edges from v',v
  endif

  iarc = netw%outheads(v)
!  while iarc ~= -1
  do
    if (iarc == -1) exit
    w = netw%seconds(iarc)
    n_arcvisit = n_arcvisit + 1

    if (msglvl > 1) &
         write(lp,'(A,I4,A,I4,A,I4)') 'iarc',iarc,', w',w,', tag',tags(w)

    if (tags(w) /= root) then

      if (msglvl > 1)  &
         write(lp,'(A,I4,A,I4)') 'flow',netw%flows(iarc),'capacity', &
              netw%capacities(iarc)

      if (netw%capacities(iarc) > netw%flows(iarc) ) then
        avail = netw%capacities(iarc) - netw%flows(iarc)

        if (msglvl > 1) write(lp,'(A,I4,A,I4)') 'deltas(',v,') =',deltas(v)

        if (avail > deltas(v)) then
          avail = deltas(v)
          if (msglvl > 1) write(lp,'(A,I4)') 'avail now',avail
        endif
        deltas(w) = avail
        pred(w) = iarc
        tags(w) = root

        if (msglvl > 1) &
             write(lp,'(A,I4,A,I4,A,I4)') 'deltas(',w,') =',deltas(w),  &
                   ', pred(',w,') =',pred(w)

        if (w == sink) then
          if (msglvl > 1) write(lp,'(A)') 'w = sink, exit'
          exit
        endif
        last = last + 1
        list(last) = w

        if (msglvl > 1) &
               write(lp,'(A,I4,A,I4)') 'add',w,' to list at last',last

      endif
    endif
    iarc = netw%nextout(iarc)
  enddo
  if (w == sink) then
    if (msglvl > 1) write(lp,'(A)') 'exit'
    exit
  endif

  if (msglvl > 1) write(lp,'(A,I4)')  'in-edges from v',v
  iarc = netw%inheads(v)
!  while iarc ~= -1
  do
    if (iarc == -1) exit
    w = netw%firsts(iarc)
    n_arcvisit = n_arcvisit + 1
    if (msglvl > 1) &
       write(lp,'(A,I4,A,I4,A,I4)') 'iarc',iarc,', w',w,', tag',tags(w)
    if (tags(w) /= root) then
      if (msglvl > 1) write(lp,'(A,I4)') 'flow',netw%flows(iarc)
      if (netw%flows(iarc) > 0) then
        if (avail > netw%flows(iarc)) then
          avail = netw%flows(iarc)
          if (msglvl > 1) write(lp,'(A,I4)') 'avail now', avail
        endif
        deltas(w) = avail
        pred(w) = iarc
        tags(w) = root
        if (msglvl > 1) &
            write(lp,'(A,I4,A,I4,A,I4,A,I4)') 'deltas(',w,') =',deltas(w),  &
                 ', pred(',w,') =',pred(w)
        last = last + 1
        list(last) = w
        if (msglvl > 1) &
             write(lp,'(A,I4,A,I4)') 'add',w,'to list at last',last
      endif
    endif
    iarc = netw%nextin(iarc)
  enddo
  if (w == sink) exit
enddo

!ISD rc was not otherwise defined or used so code changed
!if (w == sink) then
!   rc = 1
!else
!   avail = 0
!endif

if (w .ne. sink) avail = 0

stats(1) = n_nodevisit
stats(2) = n_arcvisit
deallocate (list,tags,deltas)
! pred still allocated since is input to augmentpath

if (msglvl > 0) write(lp,'(A/)') '### leaving network_findaugpath()'

end subroutine findaugpath
