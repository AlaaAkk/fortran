PROGRAM main

    USE OVERLAP_MATRIX
    USE INPUT_H2
    USE INPUT_H2O

    IMPLICIT NONE
 !   integer :: N
    integer, parameter :: K_H2 = 2  ! Number of basis functions
   ! integer, dimension(K_H2,3) :: qm_num_H2            
   ! real*8, dimension(K_H2,3) :: Rn_H2             
   ! real*8, dimension(K_H2,3) :: lamda_H2, alpha_H2 
    real*8, dimension(K_H2,K_H2) :: S_H2 
    REAL*8, allocatable, dimension(:,:) :: Rn_H2      ! Basis functions' centers
    REAL*8, allocatable, dimension(:,:) :: lamda_H2     ! Basis functions' angular momenta
    REAL*8, allocatable, dimension(:,:) :: alpha_H2     ! Contraction exponential coefficients
    integer, allocatable, dimension(:,:) :: qm_num_H2     ! Conttaction linear coefficients
    ! -----------------
    integer, parameter :: K_H2O = 7 
!    integer, dimension(K_H2O,3) :: qm_num_H2O              
!    real*8, dimension(K_H2O,3) :: Rn_H2O               
!    real*8, dimension(K_H2O,3) :: lamda_H2O, alpha_H2O  
!    real*8, dimension(K_H2O,K_H2O) :: S_H2O 
    real*8, dimension(K_H2O,K_H2O) :: S_H2O 
    REAL*8, allocatable, dimension(:,:) :: Rn_H2O      ! Basis functions' centers
    REAL*8, allocatable, dimension(:,:) :: lamda_H2O     ! Basis functions' angular momenta
    REAL*8, allocatable, dimension(:,:) :: alpha_H2O     ! Contraction exponential coefficients
    integer, allocatable, dimension(:,:) :: qm_num_H2O     ! Conttaction linear coefficients

    integer :: i,j

    ! -----------
    ! MOLECULE H2
    ! -----------
 ! Call the subroutine here 
    CALL load_h2('geom_h2.xyz','sto-3g-h2.dat',K_H2,Rn_H2,qm_num_H2,alpha_H2,lamda_H2) 
  !  Rn_H2(1,1:3) = (/-0.73D0, 0.0D0, 0.0D0/) ! Position of first H atom
  !  Rn_H2(2,1:3) = (/+0.73D0, 0.0D0, 0.0D0/) ! Position of the second H atom

  !  qm_num_H2(1,1:3) = (/0, 0, 0/) ! (0,0) H1 1s
  !  qm_num_H2(2,1:3) = (/0, 0, 0/) ! (0,0) H2 1s

  !  lamda_H2(1,1:3) = (/0.15432897, 0.53532814, 0.44463454 /)
  !  lamda_H2(2,1:3) = (/0.15432897, 0.53532814, 0.44463454 /)

  !  alpha_H2(1,1:3) = (/3.42525091, 0.62391373, 0.1688554/)
  !  alpha_H2(2,1:3) = (/3.42525091, 0.62391373, 0.1688554/)

    CALL S_overlap(K_H2,lamda_H2,alpha_H2,qm_num_H2,Rn_H2,S_H2)

    WRITE(*,*) "Overlap matrix of H2:"
    do i = 1,K_H2
       write (*,'(10f12.6)')  (S_H2(i,j),j=1,K_H2)
    end do
    DEALLOCATE(Rn_H2,lamda_H2,alpha_H2,qm_num_H2)
  
    CALL load_h2o('geom_h2o.xyz','sto-3g-h2o.dat',K_H2O,Rn_H2O,qm_num_H2O,alpha_H2O,lamda_H2O) 
    ! ------------
    ! MOLECULE H2O
    ! ------------
 ! Call the subroutine here  
    !Position
!    Rn_H2O(1,1:3) = (/1.2D0, 1.43D0,  0.0D0/)     ! position H1 1s
!    Rn_H2O(2,1:3) = (/1.2D0, -1.43D0, 0.0D0/)     ! position H2 1s
!    Rn_H2O(3,1:3) = (/0.0D0, 0.0D0, 0.0D0/)       ! position O1 1s
!    Rn_H2O(4,1:3) = (/0.0D0, 0.0D0, 0.0D0/)       ! position O1 2s
!    Rn_H2O(5,1:3) = (/0.0D0, 0.0D0, 0.0D0/)       ! position O1 2px
!    Rn_H2O(6,1:3) = (/0.0D0, 0.0D0, 0.0D0/)       ! position O1 2py
!    Rn_H2O(7,1:3) = (/0.0D0, 0.0D0, 0.0D0/)       ! position O1 2pz
!    ! a,b,c :C_Exponent 
!    qm_num_H2O(1,1:3) = (/0, 0, 0/) !  H 1s (0,0)
!    qm_num_H2O(2,1:3) = (/0, 0, 0/) !  H 1s (0,0)
!    qm_num_H2O(3,1:3) = (/0, 0, 0/) !  O 1s (0,0)
!    qm_num_H2O(4,1:3) = (/0, 0, 0/) !  O 2s (0,0)
!    qm_num_H2O(5,1:3) = (/1, 0, 0/) !  O 2px (1,1)
!    qm_num_H2O(6,1:3) = (/0, 1, 0/) !  O 2py (1,-1)
!    qm_num_H2O(7,1:3) = (/0, 0, 1/) !  O 2pz (1,0)
!    !lamda: Coeff
!    lamda_H2O(1,1:3) = (/0.44463454D0, 0.53532814D0,0.15432897D0/)
!    lamda_H2O(2,1:3) = (/0.44463454D0, 0.53532814D0,0.15432897D0/)
!    lamda_H2O(3,1:3) = (/0.44463454D0, 0.53532814D0,0.15432897D0/)
!    lamda_H2O(4,1:3) = (/0.70011547D0, 0.39951283D0, -0.09996723D0/)
!    lamda_H2O(5,1:3) = (/0.15591627D0, 0.60768372D0, 0.39195739D0/)
!    lamda_H2O(6,1:3) = (/0.15591627D0, 0.60768372D0, 0.39195739D0/)
!    lamda_H2O(7,1:3) = (/0.15591627D0, 0.60768372D0, 0.39195739D0/)
!    ! Alpha
!    alpha_H2O(1,1:3) = (/0.1688554D0, 0.62391373D0,3.42525091D0/)
!    alpha_H2O(2,1:3) = (/0.1688554D0, 0.62391373D0,3.42525091D0/)
!    alpha_H2O(3,1:3) = (/6.4436083D0, 23.808861D0, 130.70932D0/)
!    alpha_H2O(4,1:3) = (/0.380389D0, 1.1695961D0, 5.0331513D0/)
!    alpha_H2O(5,1:3) = (/0.380389D0, 1.1695961D0, 5.0331513D0/)
!    alpha_H2O(6,1:3) = (/0.380389D0, 1.1695961D0, 5.0331513D0/)
!    alpha_H2O(7,1:3) = (/0.380389D0, 1.1695961D0, 5.0331513D0/)

    CALL S_overlap(K_H2O,lamda_H2O,alpha_H2O,qm_num_H2O,Rn_H2O,S_H2O)
    WRITE(*,*) "Overlap matrix of H2O:"
    do i = 1,K_H2O
       write (*,'(10f12.6)')  (S_H2O(i,j),j=1,K_H2O)
    end do

    DEALLOCATE(Rn_H2O,lamda_H2O,alpha_H2O,qm_num_H2O)

end PROGRAM main
