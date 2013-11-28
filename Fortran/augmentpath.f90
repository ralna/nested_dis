subroutine augmentpath(netw,delta,pred,msglvl)

!Matlab call
!function newflows = network_augmentpath ( network, delta, pred, msglvl )
! AUGMENTPATH -- augment a path starting from sink
!

!   network -- network object
!   delta   -- increment flow 
!   pred    -- tree predecessor vector, size nnode
!   msglvl  -- message level

implicit none
integer, intent(in) :: msglvl,delta
type(network), intent(inout) :: netw

integer iarc,lp,sink,source,v,w

integer, intent(in), allocatable :: pred(:)

!ISD As list, dist, tags, and deltas are not used have suppressed these
!integer, allocatable :: list(:), dist(:), tags(:), deltas(:), pred(:)
!integer network.n_nodevisit, network.n_arcvisit

lp = 6

!nnode = netw%nnode
!narc = netw%narc
source = netw%source
sink = netw%sink
!firsts = network.firsts ; seconds = network.seconds ;
!capacities = network.capacities ; flows = network.flows ;
!inheads = network.inheads ; outheads = network.outheads ;
!nextin = network.nextin ; nextout = network.nextout ;

if (msglvl > 0) then
  write(lp,'(A)') ''
  write(lp,'(A)') '### inside network_augmentpath()'
  write(lp,'(A,I4,A,I4)')  'delta', delta, ', pred[sink] = ', pred(sink)
endif

if (delta <= 0 .OR. pred(sink) <= 0) then
  write(lp,'(A,I24,A,I4)')  'ERROR : delta', delta, ', pred(sink) = ',pred(sink)
  return
endif

!
! work back from the sink
!
w = sink
!while w ~= source
do
  if (w == source) exit
  iarc = pred(w)
!  keyboard ;
  if (msglvl > 1) &
     write(lp,'(A,I4,A,I4,A,I4,A,I4,A)') 'w',w,', iarc',iarc, &
                ', arc = (',netw%firsts(iarc),',', netw%seconds(iarc), ')'
  if (netw%firsts(iarc) == w) then
    v = netw%seconds(iarc)
    netw%flows(iarc) = netw%flows(iarc) - delta
    if (msglvl > 1) &
        write(lp,'(A,I4,A,I4)') 'backward arc flows(',iarc,') now',  &
              netw%flows(iarc)
  elseif (netw%seconds(iarc) == w) then
    v = netw%firsts(iarc)
    netw%flows(iarc) = netw%flows(iarc) + delta
    if (msglvl > 1) &
        write(lp,'(A,I4,A,I4)') 'forward arc flows(',iarc,') now',  &
              netw%flows(iarc)
  endif
  w = v
enddo

!
! set flow values
!

if (msglvl > 0) then
  write(lp,'(A)') '### leaving network_augmentpath()'
  write(lp,'(A)') ''
endif

end subroutine augmentpath
