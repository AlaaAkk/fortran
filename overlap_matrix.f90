module OVERLAP_MATRIX


    contains
        ! Factorial function
        function factorial(n) result(fact)
            IMPLICIT NONE
            integer, intent(in) :: n
            integer :: i
            integer :: fact

            fact = 1

            do i = 2, n
                fact = fact * i
            end do

        end function factorial
        ! Double Factorial
        function factorial2(n) result(fact2)
            IMPLICIT NONE
            integer, intent(in) :: n
            integer :: i
            integer :: fact2

            if (n == -1) then
                fact2 = 1
            else if (mod(n,2) == 0) then 
                fact2 = 1

                do i = 2, n, 2
                    fact2 = fact2 * i
                end do
            else                      ! odd
                fact2 = 1

                do i = 1, n, 2
                    fact2 = fact2 * i
                end do
            end if

        end function factorial2
        ! binomial (n // k)= n!/k!(n-k)!
        function binomial(n,k)
            IMPLICIT NONE

            integer, intent(in) :: n, k
            integer :: binomial

            binomial = factorial(n)/(factorial(k)*factorial(n-k))

        end function binomial
        ! Gaussian product
        subroutine gaussian_product(alphai,alphaj,Ra,Rb,rp,gp)
            IMPLICIT NONE

            real*8, intent(in) :: alphai, alphaj 
            real*8, dimension(3), intent(in) :: Ra, Rb 
            real*8, dimension(3) :: rab 
            real*8, intent(out) :: gp 
            real*8, dimension(3), intent(out) :: rp 

            rab = Ra - Rb
            gp =  - (alphai * alphaj * dot_product(rab,rab))/ (alphai + alphaj) 
            gp = exp(gp)
            rp = (alphai * Ra + alphaj * Rb) / (alphai + alphaj) 

        end subroutine gaussian_product

        function norm(a,b,c,alphai) result(N)
            IMPLICIT NONE

            real*8, parameter :: PI = 4.0 * ATAN(1.0) 
            integer, intent(in) :: a, b, c 
            real*8, intent(in) :: alphai 
            real*8 :: N ! Normalization factor

            N = (2 * alphai / PI)**(3.0D0 / 4.0D0) 
            N = N * (8.0D0 * alphai)**((a + b + c) / 2.0D0)
            N = N * SQRT(real(factorial(a) * factorial(b) * factorial(c)))
            N = N / SQRT(real(factorial(2 * a) * factorial(2 * b) * factorial(2 * c)))

        end function norm



        function Xij(ak,al,alphai,alphaj,xA,xB,xp)
            ! Eq(41)
            IMPLICIT NONE

            real*8, parameter :: PI = 4.0 * ATAN(1.0) 
            integer, intent(in) :: ak,al 
            real*8, intent(in) :: alphai, alphaj 
            real*8, intent(in) :: xA, xB 
            real*8, intent(in) :: xp 
            real*8 :: tmp
            integer :: n 
            integer :: m 
            real*8 :: Xij

            Xij = 0.0D0

            do m = 0, ak
                do n = 0, al
                    if (mod(m+n,2) == 0) then
                        tmp = binomial(ak,m) * binomial(al,n) * ((xp-xA)**(ak-m)) * ((xp-xB)**(al-n)) 
                        tmp = tmp * SQRT(2*PI)* factorial2(m + n - 1)
                        tmp = tmp / (2.0D0 * (alphai + alphaj))**((m + n +1) / 2.0D0)

                        Xij = Xij + tmp
                 !   else 
                 !       Xij=0    
                    end if
                end do
            end do
        end function Xij



        function S_tilda(ak,bk,ck,al,bl,cl,alphai,alphaj,Ra,Rb) result(S)
            ! Eq (39)
            IMPLICIT NONE
            real*8, parameter :: PI = 4. * ATAN(1.) 
            integer, intent(in) :: ak, bk, ck, al, bl, cl 
            real*8, intent(in) :: alphai, alphaj 
            real*8, dimension(3), intent(in) :: Ra, Rb 
            real*8, dimension(3) :: rp 
            real*8 :: gp 
            real*8 :: S

            CALL gaussian_product(alphai,alphaj,Ra,Rb,rp,gp) 

            S = 1
            S = S * Xij(ak,al,alphai,alphaj,Ra(1),Rb(1),rp(1))       ! X_tilda 
            S = S * Xij(bk,bl,alphai,alphaj,Ra(2),Rb(2),rp(2))       ! Y_tilda
            S = S * Xij(ck,cl,alphai,alphaj,Ra(3),Rb(3),rp(3))       ! Z_tilda
            S = S * norm(ak,bk,ck,alphai) * norm(al,bl,cl,alphaj)   ! Normalization factors of two gaussians
            S = S * gp                                     ! gausian product

        end function S_tilda


        subroutine S_overlap(Kf,lamda,alpha,qm_num,Rn,S)
            IMPLICIT NONE
            integer, parameter :: ci=3 ! 3 for H2 and H2O
            integer, intent(in) :: Kf ! Number of basis functions
            real*8, dimension(Kf,3), intent(in) :: Rn 
            integer, dimension(Kf,3), intent(in) :: qm_num 
            real*8, dimension(Kf,3), intent(in) :: lamda 
            real*8, dimension(Kf,3), intent(in) :: alpha 
            integer :: i,j,k,l
            real*8 :: tmp
            real*8, dimension(Kf,Kf), intent(out) :: S

            S(:,:) = 0.0D0

            do i = 1,Kf
                do j = 1,Kf
                    do k = 1,ci
                        do l = 1,ci
                            tmp = lamda(i,k) * lamda(j,l)
                            tmp = tmp * S_tilda(qm_num(i,1),qm_num(i,2),qm_num(i,3),&
                                    qm_num(j,1),qm_num(j,2),qm_num(j,3),alpha(i,k),alpha(j,l),&
                                    Rn(i,:),Rn(j,:))   

                            S(i,j) = S(i,j) + tmp
                        end do ! l
                    end do ! k
                end do ! j
            end do ! i

        end subroutine S_overlap






end module OVERLAP_MATRIX
