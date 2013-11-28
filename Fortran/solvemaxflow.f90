subroutine solvemaxflow(netw,nvtx,msglvl,maps1,maps2,mark1,mark2, &
                        pred,list)
! Find two partitions of a wide separator graph by solving a max flow problem


! input --
!    nvtx -- number of vertices in separator
!  msglvl -- message level control

!
! output --
!
!    mapS1[n_S] -- first map from wide separator to {0,1,2} = {S,B,W}
!    mapS2[n_S] -- second map from wide separator to {0,1,2} = {S,B,W}
!

! input/output --
!
!    network graph

! work arrays
!
!    mark1, mark2, pred, list of length nnode
!

implicit none
type(network), intent(inout) :: netw
integer, intent(in) :: msglvl,nvtx
integer, intent(out) :: mapS1(:),mapS2(:)

integer :: mark1(:),mark2(:),pred(:),list(:)

! Local variables
integer iarc,ii,lp,narc,nnode,u,i

lp = 6


nnode = netw%nnode
narc = netw%narc

!
! Find maxflow through the network using the Ford-Fulkerson algorithm
!
call findmaxflow(netw,msglvl,pred,list,mark1,mark2)


!
! Find the two mincuts
!
call findmincut(netw,msglvl,mark1,mark2,list)

!
!Use mark1 and mark2 to generate maps1 and maps2
!

mapS1 = -1
mapS2 = -1

do ii = 2,nnode-1,2
  u = ii/2
  if (mark1(ii) == 1) then
    if (mark1(ii+1) == 1) then
       mapS1(u) = 1
    else
      mapS1(u) = 0
    endif
  endif
  if (mark1(ii) == 2) then
    if (mark1(ii+1) == 2) then
      mapS1(u) = 2
    else
      mapS1(u) = 0
    endif
  endif
enddo

do ii = 2,nnode-1,2
   u = ii/2
   if (mark2(ii) == 1) then
      if (mark2(ii+1) == 1) then
         mapS2(u) = 1
      else
         mapS2(u) = 0
      endif
   endif
   if (mark2(ii) == 2) then
      if (mark2(ii+1) == 2) then
         mapS2(u) = 2
      else
         mapS2(u) = 0
      endif
   endif
enddo

end subroutine solvemaxflow

