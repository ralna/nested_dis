subroutine augmentpath(netw,delta,pred,msglvl)

! Reset flows on augmenting path
!

! Input
!   delta   -- increment flow 
!   pred    -- tree predecessor vector, size nnode
!   msglvl  -- message level

! Input/output
!   network -- network object

implicit none
integer, intent(in) :: msglvl,delta
type(network), intent(inout) :: netw

integer iarc,lp,sink,source,v,w

integer, intent(in) :: pred(:)


lp = 6

source = netw%source
sink = netw%sink

! Should set an error flag
if (delta <= 0 .OR. pred(sink) <= 0) then
  write(lp,'(A,I4,A,I4)')  'ERROR : delta', delta, ', pred(sink) = ',pred(sink)
  return
endif

!
! work back from the sink resetting network flows
!
w = sink
!while w ~= source
do
  if (w == source) exit
  iarc = pred(w)
  if (netw%firsts(iarc) == w) then
    v = netw%seconds(iarc)
    netw%flows(iarc) = netw%flows(iarc) - delta
  elseif (netw%seconds(iarc) == w) then
    v = netw%firsts(iarc)
    netw%flows(iarc) = netw%flows(iarc) + delta
  endif
  w = v
enddo

end subroutine augmentpath
