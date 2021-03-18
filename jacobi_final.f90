program Main
implicit none
integer, parameter :: n=2
real*8, parameter:: error=1.0e-15
real*8, dimension(n,n)  :: h, x, ht
real*8, dimension(n)  :: arr
real*8  w(n),A(n,n)
integer i, j
! matrix h
!  h=reshape((/1.0,1.0,0.5,9.5,1.0,4.0,5.0,3.0,0.5,5.0,11.0,8.5,9.5,3.0,8.5,8.0/),(/n,n/))
  CALL RANDOM_NUMBER(h)  ! create a random matrix h
! Now symmetrize it:
  ht=transpose(h)
  h=(h+ht)/2
 ! h = reshape((/1,-1,0,-1,3,4,0,4,1/),(/3,3/))
! print a header and the original matrix
  print *, "The Original matrix h is:"
  do i=1,n
     write (*,'(10f12.6)') (h(i,j),j=1,n)
  end do
  do i=1,n
    do j=1,n
       A(i,j)=h(i,j)
    end do
  end do
  call eigenvalues(A,n,w)
  call Jacobi(h,x,error,n)

! print solutions
  print *, "The Eigenvalues are:"
  write (*,'(10f12.6)') (h(i,i),i=1,n)
  print *,"The Eigenvectors are:"
  do i = 1,n
     write (*,'(10f12.6)')  (x(i,j),j=1,n)
  end do
  do i=1,n
     arr(i)=h(i,i)
  end do
!  print * , 'arr is' , arr
  call sort_pick(n,arr,x)
  print *,"The Eigenvalues sorted in ascending order:"
     write (*,'(10f12.6)')  arr
  print *,"The corresponding Eigenvectors sorted in ascending order:"
  do i = 1,n
     write (*,'(10f12.6)')  (x(i,j),j=1,n)
  end do
end program main
subroutine eigenvalues(A,kk,eigvalues)
! to calculate eigen values from dsyev
!calling list
integer             :: i,j
integer, intent(in)             :: kk
real*8, intent(inout) :: eigvalues(kk)
!real*8, intent(inout) :: eigvectors(kk,kk)
real*8, intent(in) :: A(kk,kk)

!local 
real*8,allocatable :: work(:)
integer                      :: lwork,info

lwork = max(1,3*kk-1)
allocate(work(lwork))


call dsyev('V','U',kk,A,kk,eigvalues,WORK,LWORK,info)
print *, 'According to dsyev the eigenvalues are:'
write (*,'(10f12.6)') eigvalues
print *, 'According to dsyev the eigenvectors are:'
do i = 1,kk
   write (*,'(10f12.6)')  (A(i,j),j=1,kk)
end do
RETURN
!call dsyev('V','L',kk,A,kk,eigvectors,WORK,LWORK,info)
deallocate(work)
!if(info .neq. 0) exit
end subroutine
Subroutine sort_pick(n, arr, xx)
!USE nrtype
implicit none
integer,  intent(in) :: n 

real*8, dimension(n)  :: b
real*8, dimension(n), intent(inout)  :: arr
real*8, dimension(n,n), intent(inout)  :: xx
!Sorts an array arr into ascending numerical order
integer :: i,j
!REAL :: n=size(arr)
real*8 :: a
do j=2,n
!Pick out each element in turn.
  a=arr(j)
  b=xx(:,j)
do i=j-1,1,-1
  if (arr(i) .le. a) exit
     arr(i+1)=arr(i)
     xx(:,i+1)=xx(:,i)
  end do
arr(i+1)=a   
xx(:,i+1)=b  
end do
End subroutine sort_pick
function sgn(x) result(y) 
! sgn function         
  real*8, intent(in)   :: x
  integer :: y
  if (x > 0) then
    y = 1
  else if(x < 0) then
    y = -1
  else  
    y = 0
  end if
end function 
subroutine Jacobi(h,x,error,n)
implicit none
integer i, j, k, n
real*8, dimension(n,n):: h,x
real*8 :: error, add, bar 
real*8 :: delta, t, c, s,  tau
integer :: sgn
!integer, intent(in) :: y
!1) initialize C=1
!n=3
x = 0.0
do i=1,n
  x(i,i) = 1.0   
end do
! find the sum of all off-diagonal elements (squared)
add = 0.0  
do i=1,n
  do j=1,n
    if (i.ne.j) add = add + h(i,j)**2
  end do
end do
!print *, "this is add", add
if (add <= error) return  ! converged

! average for off-diagonal elements /2
bar = 0.5*add/float(n*n)

! 2)
do while (add.gt.error)   ! not converged
  do i=2,n  ! i from 2 to n
    do j=1,i-1 ! j from 1 to i-1
      if (h(j,i)**2 <= bar) cycle  ! do not touch small elements, for faster calculation
      add = add - 2.0*h(j,i)**2
      bar = 0.5*add/float(n*n)
!a)   
      delta = (h(i,i)-h(j,j))   !hii-hjj
      t = (2*sgn(delta)*h(i,j))/(sqrt((4*h(i,j)**2)+(abs(delta))**2)+abs(delta))   
      c = 1/(sqrt(1+t**2))   ! c
      s = c*t  ! s
!b) 
      do k=1,n    ! looping over k
        tau = -s*x(k,i)+c*x(k,j)  
        x(k,i) =  c*x(k,i)+s*x(k,j)  
        x(k,j) = tau
      end do
!c) 
      do k=1,n
        if (k .ne.i .and. k .ne. j) then
           tau = -s*h(k,i)+c*h(k,j) 
           h(k,i) =  c*h(k,i)+s*h(k,j) 
           h(k,j) = tau
           h(i,k) = h(k,i)
           h(j,k) = h(k,j)
        end if
      end do
!d)      
        tau=(c**2)*h(i,i) + (s**2)*h(j,j) + (2*c*s*h(i,j))
        h(j,j)=(s**2)*h(i,i) + (c**2)*h(j,j) - (2*c*s*h(i,j)) 
        h(i,i) = tau
        h(i,j) = 0
        h(j,i) = 0
    end do  ! j
  end do ! i
end do ! while
return
! gfortran jacobi_final.f90 -llapack -lblas -o out.out -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
        
end subroutine 
