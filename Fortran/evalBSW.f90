subroutine evalBSW (a_n, a_ne, a_ptr, a_row, a_weight, map, alpha, beta, stats, stats10)
!Matlab call
!function stats = evalBSW ( A, map, alpha, beta, msglvl )
! 
! stats = EVALBSW ( A, map, alpha, beta )
!
! input ---
!    map[nvtx] -- map from vertices to region
!      map[u] == 0 --> u in S
!      map[u] == 1 --> u in B
!      map[u] == 2 --> u in W
!    alpha --- acceptability parameter
!    beta  --- imbalance penalty parameter
!
! output --
!    stats[1] -- weight of vertices in S
!    stats[2] -- weight of vertices in B
!    stats[3] -- weight of vertices in W
!    stats[4] -- weight of edges in A_{S,S}
!    stats[5] -- weight of edges in A_{S,B}
!    stats[6] -- weight of edges in A_{S,W}
!    stats[7] -- weight of edges in A_{B,B}
!    stats[8] -- weight of edges in A_{B,W}
!    stats[9] -- 1 if acceptable, 0 if not
!       acceptable --> alpha*min(|B|,|W|) >= max(|B|,|W|) 
!    stats[10] -- cost of partition
!       cost = |S|*(1 + (beta*| |B| - |W| |)/(|B|+|S|+|W|)) ;
!
! created -- 12jan12, cca
!
implicit none
integer, intent(in) :: a_n
integer, intent(in) :: a_ne
integer, intent(in) :: map(:),a_ptr(:),a_row(:),a_weight(:)
real(wp), intent(in) :: alpha, beta
integer, intent(out) :: stats(9)
real(wp), intent(out) :: stats10
integer minBW,maxBW,nSS,nSb,nSW,nBB,nWW,nvtx,nS,nB,nW
integer j,j1,j2,jj,u,v
real(wp) diffBW
nvtx = a_n
stats(1:9) = -1
nS = 0
nB = 0
nW = 0
do u = 1,nvtx
   if (map(u) == 0) then
      nS = nS + a_weight(u)
   elseif (map(u) == 1) then
      nB = nB + a_weight(u)
   elseif (map(u) == 2) then
      nW = nW + a_weight(u)
   endif
enddo
stats(1) = nS
stats(2) = nB
stats(3) = nW
minBW = min(nB,nW)
maxBW = max(nB,nW)
diffBW = real(abs(nB - nW))/real(nS + nB + nW)
if (.false.) then
nSS = 0
nSB = 0
nSW = 0
nBB = 0
nWW = 0
![rows, cols, ents] = find(A) ;
!nzA = length(rows) ;
do j = 1,a_n
  j1 = a_ptr(j)
  if (j.eq.a_n) then
   j2 = a_ne
  else
   j2 = a_ptr(j+1)-1
  end if
  v = j
  do jj = j1,j2
    u = a_row(jj) ; 
!   v = cols(ii) ;
    if (map(u) == 0) then
      if (map(v) == 0) then
        nSS = nSS + 1
      elseif (map(v) == 1) then
        nSB = nSB + 1
      elseif (map(v) == 2) then
        nSW = nSW + 1
      endif
    elseif (map(u) == 1) then
      if (map(v) == 1) then
        nBB = nBB + 1
      endif
    elseif (map(u) == 2) then
      if (map(v) == 2) then
        nWW = nWW + 1
      endif
    endif
  enddo
enddo
stats(4) = nSS
stats(5) = nSB
stats(6) = nSW
stats(7) = nBB
stats(8) = nWW
end if
!    stats[9] -- 1 if acceptable, 0 if not
!       acceptable --> alpha*min(|B|,|W|) >= max(|B|,|W|) 
!    stats[10] -- cost of partition
!       cost = |S|*(1 + (beta*| |B| - |W| |)/(|B|+|S|+|W|)) ;
!     or
!       cost = |S|/(|B||W|)
if (alpha*minBW >= maxBW) then
   stats(9) = 1
else
   stats(9) = 0
endif
!write(lp,'(A)') 'eval routine'
!write(lp,'(A,I4)') 'ns',nS
!write(lp,'(A,D12.4)') 'diffBW',diffBW
!write(lp,'(A,D12.4)') 'beta',beta
!stats10 = real(nS)*(1.0d0 + beta*diffBW)
stats10 = (real(nS)/real(nW))/real(nB)

end subroutine evalBSW
