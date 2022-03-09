MODULE INPUT_H2

    IMPLICIT NONE

    CONTAINS

        SUBROUTINE load_h2(geometry,basis,N,Rn,qm_num,alpha,lamda)

            IMPLICIT NONE

            ! INPUT
            INTEGER :: Ne                       ! Basis set size
            INTEGER, intent(in) :: N           ! Input file name
            CHARACTER(len=*), intent(in) :: geometry           ! Input file name
            CHARACTER(len=*), intent(in) :: basis              ! Input file name

            ! INTERMEDIATE VARIABLES
            INTEGER :: i,j ! Loop index

            ! ALLOCATABLE
            REAL*8, allocatable, dimension(:,:) , intent(out) :: Rn        
            integer, allocatable, dimension(:,:) , intent(out) :: qm_num    
            REAL*8, allocatable, dimension(:,:) , intent(out) :: alpha     
            REAL*8, allocatable, dimension(:,:), intent(out) :: lamda      
            CHARACTER (len=2) :: species(N)           ! Atomic potisions
            CHARACTER (len=2) :: qm          ! Atomic potisions

            WRITE(*,*) "Opening Geometry: ", geometry
            OPEN(unit=100,file=geometry,form="formatted",status="old",action="read")

            READ(100,*) Ne  ! reading number of atoms
            ALLOCATE(Rn(N,3))
            ALLOCATE(qm_num(N,3))
            ALLOCATE(alpha(N,3))
            ALLOCATE(lamda(N,3))

            DO i = 1, N
            READ(100,*) species(i), Rn(i,:)  ! read the species and save and positions
            END DO
             
            CLOSE(unit=100) ! Close the file

!            WRITE(*,*) "species? : ", (species(i), i=1,N)
!            WRITE(*,*) "geo? : "
!            write (*,'(3f12.6)') (Rn(i,:),i=1,N)

            ! Here I should read the basis file
            WRITE(*,*) "Opening sto_3g : ", basis
            OPEN(unit=200,file=basis,form="formatted",status="old",action="read")

            DO i=1,13
            READ(200,*) 
            END DO

            READ(200,*) qm 
            IF (qm .eq. 'S') then
               DO i = 1, Ne
                 qm_num(i,1:3)=(/0,0,0/)  ! read the species and save and positions
               END DO
            END IF

            i=1   
            do while (i < Ne)
               DO j = 1, 3
                  READ(200,*)  alpha(i,j), lamda(i,j)
               END DO
               alpha(i+1,:)=alpha(i,:)
               lamda(i+1,:)=lamda(i,:)
               i=i+1
            END DO

            CLOSE(unit=200) ! Close the file

        END SUBROUTINE load_h2

END MODULE INPUT_H2
