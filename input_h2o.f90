MODULE INPUT_H2O

    IMPLICIT NONE

    CONTAINS

        SUBROUTINE load_h2o(geometry,basis,N,Rn,qm_num,alpha,lamda)

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
            CHARACTER (len=2) :: qm_1          ! Atomic potisions

            WRITE(*,*) "Opening Geometry: ", geometry
            OPEN(unit=100,file=geometry,form="formatted",status="old",action="read")

            READ(100,*) Ne  ! reading number of atoms
            ALLOCATE(Rn(N,3))
            ALLOCATE(qm_num(N,3))
            ALLOCATE(alpha(N,3))
            ALLOCATE(lamda(N,3))

            DO i = 1, Ne
            READ(100,*) species(i), Rn(i,:)  ! read the species and save and positions
            END DO
             
            CLOSE(unit=100) ! Close the file

          !  WRITE(*,*) "species? : ", (species(i), i=1,Ne)  ! how many H and how many O?
          !  WRITE(*,*) "geo? : "
          !  write (*,'(3f12.6)') (Rn(i,:),i=1,Ne)

            ! Here I should read the basis file
            WRITE(*,*) "Opening sto_3g : ", basis
            OPEN(unit=200,file=basis,form="formatted",status="old",action="read")

            DO i=1,13
            READ(200,*) 
            END DO

            READ(200,*) qm_1 ! read S 
            qm_num(1,1:3)=(/0,0,0/)  ! 
            qm_num(2,1:3)=(/0,0,0/)  ! 
            qm_num(3,1:3)=(/0,0,0/)  ! 
            qm_num(4,1:3)=(/0,0,0/)  ! 
            qm_num(5,1:3) = (/1, 0, 0/) !  O 2px (1,1)
            qm_num(6,1:3) = (/0, 1, 0/) !  O 2py (1,-1)
            qm_num(7,1:3) = (/0, 0, 1/) !  O 2pz (1,0)


            i=1   
            do while (i < 2)
               DO j = 1, 3
                  READ(200,*)  alpha(i,j), lamda(i,j)
               END DO
               alpha(i+1,:)=alpha(i,:)
               lamda(i+1,:)=lamda(i,:)
               i=i+1
            END DO
         !  lamda(3,1:3)=lamda(2,1:3)
           DO i=1,3
           READ(200,*) 
           END DO
           DO j = 1, 3
              READ(200,*)  alpha(3,j),lamda(3,j)
           END DO
           READ(200,*)  ! the Sp
            i=4   
               DO j = 1, 3
                  READ(200,*)  alpha(i,j), lamda(4,j), lamda(i+1,j)
               END DO
            do while (i < 6)
               alpha(i+1,:)=alpha(i,:)
               alpha(i+2,:)=alpha(i,:)
               lamda(i+2,:)=lamda(i+1,:)
               i=i+1
            END DO
            CLOSE(unit=200) ! Close the file

        END SUBROUTINE load_h2o

END MODULE INPUT_H2O
