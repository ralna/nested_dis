subroutine findpenalty(netw,msglvl,nvtx,vwts,wtB,wtW,wtS,&
                       count,mark,mark1,list,head, &
                       imb,penalty)

!  Finds one or two mincuts, one nearest the source, one nearest the sink

! Input parameters

!   netw    -- network object
!   msglvl  -- message level

! Output parameters
!   mark --- records level set of vertex in cutset and two-level set also
!   mark1 --- records level set of vertex in cutset strating from sink
!

implicit none
integer, intent(in) :: msglvl,nvtx
integer, intent(in) :: vwts(:)
integer, intent(inout) :: wtB,wtW
integer, intent(out) :: mark(:),mark1(:)
type(network), intent(inout) :: netw

! Workspace
!   list --- to hold list of nodes being searched
integer :: list(:),head(:),penalty(:),count(:)
real(wp) :: imb(:)

! Local variables
integer iarc,inode,kptr,last,lp,maxl,minl,nnode,now,sink,source,wtS,x,z

lp = 0

nnode = netw%nnode
source = netw%source
sink = netw%sink

if (msglvl > 0) then
   write(lp,'(A)')  ''
   write(lp,'(A)')  '### inside findpenalty()'
   write(lp,'(A,2I4)')  'Source and sink nodes are: ',source,sink
endif

!
! breadth first traversal from source
!
if (msglvl > 0) write(lp,'(A)')  'breadth first traversal from the source'

!k equals level set
k = 0
mark = 0
mark(source) = k
now = 1
last = 1
list(now) = source
kptr = last

! Each pass through this loop determines node in level set k
do k = 1,nnode
if (msglvl > 1) write(lp,'(A,I4,A,I4)') 'now = ',now,' kptr = ',kptr
! Run through all nodes in the previous level set
  do 
    x = list(now)
    if (now > kptr) exit
    if (msglvl > 1) write(lp,'(A,I4,A,I4)') 'Vertex ',x,' in list position',now
    now = now + 1
! Find all nodes connected to this node
    iarc = netw%outheads(x)
    do
      if (msglvl > 1) write(lp,'(A,I4)')  'iarc =',iarc
      if (iarc == -1) exit
      z = netw%seconds(iarc)
      if (msglvl > 1) write(lp,'(A,I4)') 'z = ',z
      if (z == sink) then
        iarc = netw%nextout(iarc)
        cycle
      endif
! Jump if already visited
      if (mark(z/2) /= 0) then
        iarc = netw%nextout(iarc)
        cycle
      endif
! Add z to set of nodes associated with this level for future searching
      last = last + 1
! Go immediately to node u+
      list(last) = z+1
! Node z is in level set k
      mark(z/2) = k
      if (msglvl > 1) write(lp,'(A,I4,A,I4)') 'Vertex ',z/2,' is in level ',k
      iarc = netw%nextout(iarc)
    enddo
  enddo
  if (last == kptr) exit
  kptr = last
enddo

!
! breadth first traversal from sink
!
if (msglvl > 0) write(lp,'(A)')  'breadth first traversal from the sink'

mark1 = 0
mark1(sink) = 1
maxl = -nnode
minl = nnode
list = 0
now = 1
last = 1
kptr = 1
list(now) = sink
do k = 1,nnode
!while now <= last
if (msglvl > 1) write(lp,'(A,I4,A,I4)') 'now = ',now,' kptr = ',kptr
  do 
    if (now > kptr) exit
    x = list(now)
    if (msglvl > 1) write(lp,'(A,I4,A,I4)') 'Vertex ',x,' in list position',now
    now = now + 1
    iarc = netw%inheads(x)
    do
      if (msglvl > 1) write(lp,'(A,I4)')  'iarc =',iarc
      if (iarc == -1) exit
      z = netw%firsts(iarc)
      if (msglvl > 1) write(lp,'(A,I4)') 'z = ',z
      if (z==source) then
        iarc = netw%nextin(iarc)
        cycle
      endif
      if (mark1((z-1)/2) /= 0) then
        iarc = netw%nextin(iarc)
        cycle
      endif
      last = last + 1
      list(last) = z-1
      mark1((z-1)/2) = k
      if (msglvl > 1) write(lp,'(A,I4,A,I4)') 'Vertex ',(z-1)/2,' is in level ',k
      mark((z-1)/2) = mark((z-1)/2) - k
      minl = min(minl,mark((z-1)/2))
      maxl = max(maxl,mark((z-1)/2))
      iarc = netw%nextin(iarc)
    enddo
  enddo
  if (last == kptr) exit
  kptr = last
enddo

! Compute half-level sets
if (msglvl > 1) write(lp,'(A,2I4)') 'minl, maxl ',minl,maxl
! We will count from zero to maxl-minl
count(1:maxl-minl+1) = 0
head(1:maxl-minl+1) = -1
do i = 1,nvtx
  mark(i) = mark(i) - minl + 1
  list(i) = head(mark(i))
  head(mark(i)) = i
  count(mark(i)) = count(mark(i)) + vwts(i)
enddo
if (msglvl > 1) then
write(lp,'(A)') 'Number of vertices in each half-level set'
do i = 1,maxl-minl+1
  write(lp,'(10(2I4))') i,count(i)
enddo
end if
! Run through half-level sets computing imbalances
! wtB is weight of B
! wtW is set to total weight
! Still need to weight vertices in S
wtW = wtW + wtS
if (msglvl > 1) write(lp,('(A,2I4)')) 'wtB,wtW ', wtB,wtW
if (maxl-minl == 0) then
  imb(1) = max(real(wtB)/real(wtW),real(wtW)/real(wtB))
else
  wtW = wtW - count(1)
  do k = 1,maxl-minl
    wtW = wtW - count(k+1)
    imb(k) = max(real(wtB)/real(wtW),real(wtW)/real(wtB))
    wtB = wtB + count(k)
  enddo
endif
if (msglvl > 1) then
write(lp,'(A)') 'Imbalances'
do i = 1,maxl-minl+1
  write(lp,'(10(2G12.2))') i,imb(i)
enddo
end if
! Run through nodes in level set assigning penalty to them
do k = 1,maxl-minl+1
  inode = head(k)
  do
    if (inode == -1) exit
!   write(lp,('(A,2I4)')) 'level and inode',k,inode
    if (k == 1) then 
      penalty(inode) = floor(100*imb(1))
    else
     if (k == maxl-minl+1) penalty(inode) = floor(100*imb(maxl-minl))
     if (k > 1 .and. k < maxl-minl+1)   &
                          penalty(inode) = floor(100*min(imb(k),imb(k-1)))
    end if
    inode = list(inode)
  enddo
enddo
if (msglvl > 1) then
write(lp,'(A)') 'Computed penalties'
do i = 1,nvtx
  write(lp,'(10(2I4))') i,penalty(i)
enddo
end if

if (msglvl > 0) write(lp,'(A/)')  '### leaving findpenalty()'

end subroutine findpenalty
