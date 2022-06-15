module matrix

    !==============================================================
    ! This code is part of FCclasses2
    !==============================================================
    !
    ! Despription
    ! -----------
    ! Subroutines to perform matrix manipulations:
    ! *diagonalization of symmetric matrices
    ! *Lowding orthogonalization
    ! *Determinant
    ! *Matrix printing
    !
    ! Notes
    ! LAPACK subroutines are used
    !==============================================================

    implicit none

    CONTAINS

    function identity_matrix(N) result(A)
        
        !===============================
        ! Description
        ! -----------
        ! Returns an identity matrix of 
        ! size N
        ! 
        !===============================

        integer,intent(in) :: N
        real(8),dimension(N,N) :: A
        
        ! Local
        integer :: i

        A(1:N,1:N) = 0.d0
        do i=1,N
            A(i,i) = 1.d0
        enddo

        return

    end function identity_matrix

    function determinant_realsym(N,A) result(det)
    
        !=============================================================
        ! Description
        ! -------------
        ! Computes the determinant by first applying a PLU-like 
        ! decomposition. Actually, the LAPACK symmetric routine.
        ! dsytrf is used performing Bunch-Kaufman diagonal pivoting 
        ! method. It performs the transformation of input A as (*):
        !  A = L * D * L^t
        ! where L is a unit lower triangular matrix multiplied by a
        ! permutation and D is a diagonal matrix with only 1x1 or 
        ! 2x2 non-zero diagonal blocks. Since L is a unit triangular
        ! matrix multiplied by a permutation:
        ! det(L)*det(L^t) = 1
        ! So det(A) = det(D)
        ! The only complications may arise with 2x2 blocks, which are 
        ! marked by negative values of ipiv(i) = ipiv(i+1). Since the 
        ! block is symmetric, only one of the non-diagonal elements is
        ! printed: the D(k+1,k).
        ! Note that it might arrise that ipiv(i+2)=ipiv(i+1), but it
        ! conforms another block. This will confuse an if statement based
        ! only on ipiv(i) and ipiv(i+1), so we explicitely indicate a 
        ! block with a logical variable
        !
        ! (*) this is valid if 'L' is used. For 'U' the actual details change
        !==============================================================

        real(8) :: det 
        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: A
        ! Local
        real(8),dimension(N,N) :: AuxArray
        real(8) :: AuxScalar
        ! LAPACK things
        double complex,dimension(:),allocatable :: work
        integer,dimension(1:N) :: ipiv
        integer :: info, lwork
        ! Other
        integer :: alloc_status
        logical :: cycle_next=.false.
        ! Counters
        integer :: i

        ! Copy input matrix
        AuxArray(1:N,1:N) = A(1:N,1:N)

        ! Initialize LAPACK (calculating optimal size first)
        call dsytrf('L', N, AuxArray, N, ipiv, AuxScalar, -1, info)
        lwork = int(AuxScalar)
        allocate(work(1:lwork),stat=alloc_status) 
        if (alloc_status /= 0) then
            write(0,*) "ERROR: Memory cannot be allocated (DETERMINANT ROUTINE)"
            stop
        endif

        ! Perform factorization
        call dsytrf('L', N, AuxArray, N, ipiv, work, lwork, info)
        if (info < 0) then
            write(0,*) "ERROR IN dsytrf (DETERMINANT ROUTINE)"
            stop
        endif

        ! And compute the determinant
        det = 1.d0
        do i=1,N
            if (cycle_next) then
                cycle_next=.false.
                cycle
            endif
            if (ipiv(i) >= 0) then
                det = det * AuxArray(i,i)
            else if (ipiv(i) < 0 .and. ipiv(i) == ipiv(i+1)) then
                ! Compute the 2x2 block
                cycle_next=.true.
                det = det * &
                      ( AuxArray(i,i)*AuxArray(i+1,i+1) - AuxArray(i+1,i)**2 )
            ! The if block is reordered to avoid evaluating ipiv(i-1) when i=1
            ! CHECK THE RESULTS ARE STILL OK!
            ! In any case, this else should never be reached using the cycle_next approach
            else !if (ipiv(i) < 0 .and. ipiv(i) == ipiv(i-1)) then
                ! Second part of a 2x2 block
                cycle
            endif
        enddo
        
        return

    end function determinant_realsym

    function determinant_realgen(N,A) result(det)
    
        !=============================================================
        ! Description
        ! -------------
        ! Computes the determinant by first applying a PLU-like 
        ! decomposition, using the LAPACK for generic matrix
        !==============================================================

        real(8) :: det 
        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: A
        ! Local
        real(8),dimension(N,N) :: AuxArray
        real(8) :: AuxScalar
        ! LAPACK things
        double complex,dimension(:),allocatable :: work
        integer,dimension(1:N) :: ipiv
        integer :: info, lwork
        ! Counters
        integer :: i

        ! Copy input matrix
        AuxArray(1:N,1:N) = A(1:N,1:N)

        ! Perform factorization
        call dgetrf(N, N, AuxArray, N, ipiv, info)
        ! if info>0, it worked, but the determinant will be zero
        if (info < 0) then
            write(0,*) "ERROR IN dgetrf (GEN DETERMINANT ROUTINE)"
            stop
        endif

        ! And compute the determinant
        det = 1.d0
        do i=1,N
            det = det*AuxArray(i,i)
            if (ipiv(i) /= i) det=-det
        enddo
        
        return

    end function determinant_realgen

    function inverse_realsym(N,A) result(Ainv)
    
        !========================================
        ! Description
        ! Computes the determinant by first 
        ! applying a PLU decomposition.
        ! Using LAPACK
        !=======================================

        real(8),dimension(N,N) :: Ainv 
        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: A
        ! Local
        real(8) :: AuxScalar
        ! LAPACK things
        double complex,dimension(:),allocatable :: work
        integer,dimension(1:N) :: ipiv
        integer :: info, lwork
        ! Other
        integer :: alloc_status
        ! Counters
        integer :: i,j

        ! Copy input matrix
        Ainv(1:N,1:N) = A(1:N,1:N)

        ! Initialize LAPACK (calculating optimal size first)
        call dsytrf('L', N, Ainv, N, ipiv, AuxScalar, -1, info)
        lwork = int(AuxScalar)
        allocate(work(1:lwork),stat=alloc_status) 
        if (alloc_status /= 0) then
            print*, "ERROR: Memory cannot be allocated (LAPACK)"
            stop
        endif

        ! Inverse
        call dsytrf('L', N, Ainv, N, ipiv, work, lwork, info)
        if (info < 0) then
            write(0,*) "ERROR IN dsytrf (INVERSION ROUTINE): illegal value"
            stop
        elseif (info > 0) then
            write(0,*) "ERROR IN dsytrf (INVERSION ROUTINE): singularity"
            stop
        endif
        call dsytri('L', N, Ainv, N, ipiv, work, info)
        if (info /= 0) then
            write(0,*) "ERROR IN dsytri (INVERSION ROUTINE)"
            stop
        endif
        do i=1,N
        do j=i+1,N
            Ainv(i,j) = Ainv(j,i)
        enddo
        enddo

        return

    end function inverse_realsym

    function inverse_realgen(N,A) result(Ainv)
    
        !========================================
        ! Description
        ! Computes the determinant by first 
        ! applying a PLU decomposition.
        ! Using LAPACK
        !=======================================

        real(8),dimension(N,N) :: Ainv 
        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: A
        ! Local
        real(8) :: AuxScalar
        ! LAPACK things
        double complex,dimension(:),allocatable :: work
        integer,dimension(1:N) :: ipiv
        integer :: info, lwork
        ! Other
        integer :: alloc_status

        ! Copy input matrix
        Ainv(1:N,1:N) = A(1:N,1:N)

        ! Perform factorization
        call dgetrf(N, N, Ainv, N, ipiv, info)
        if (info < 0) then
            write(0,*) "ERROR IN dgetrf (INVERSION ROUTINE): illegal value"
            stop
        elseif (info > 0) then
            write(0,*) "ERROR IN dgetrf (INVERSION ROUTINE): singularity"
            stop
        endif
        ! Initialize LAPACK (calculating optimal size first)
        call dgetri(N, Ainv, N, ipiv, AuxScalar, -1, info)
        lwork = int(AuxScalar)
        allocate(work(1:lwork),stat=alloc_status) 
        if (alloc_status /= 0) then
            print*, "ERROR: Memory cannot be allocated (LAPACK)"
            stop
        endif
        call dgetri(N, Ainv, N, ipiv, work, lwork, info)
        if (info /= 0) then
            write(0,*) "ERROR IN dgetri (INVERSION ROUTINE)"
            stop
        endif

        return

    end function inverse_realgen

    subroutine diagonalize_full(A,N,U,d,alg)
        !
        ! =========================================================
        ! Wrapper code: it calls a given diagonalization SR,
        !               given a lower triangular matrix as input
        !               For symmetric matices only
        ! =========================================================
        !  Input:
        !    A(:,:)  - full matrix to diagonalize. It is not destroyed
        !            2D array dimension: (N,N)
        !    N     - matrix dimension
        !    alg   - string to select the algorith
        !  Output:
        !    U(:,:)- orthogonal matrix with eigenvectors as columns
        !    d(:)  - eigenvalues
        !

        implicit none

        !Input
        character(len=*),intent(in) :: alg
        integer, intent(in) :: N
        real(8),dimension(:,:), intent(in) :: A
        !Output
        real(8),dimension(:,:),intent(out) :: U
        real(8),dimension(:),intent(out) :: d 

        !Local
        real(8), parameter :: PREC=1.d-11
        integer :: IER
        real(8),dimension(N,N) :: Aux

        !Auxiliar
        character(len=50) :: dummy_char
        integer :: i,j,k


        select case (trim(adjustl(alg)))
            !The one from the master (M1 year). This is left for testing ONLY
            case("emtccm")
              !We need to feed the whole matrix, as we have
              call JACOBI(N,PREC,A,U,d)

! #ifdef USE_LAPACK
            !Only available if compiled in duble precision and if lapack libs are required
            case("lapack")
              !We need to feed the whole matrix, as we have
              !We don't destroy the original
              Aux(1:N,1:N)=A(1:N,1:N)
              call diasym(Aux(1:N,1:N),d(1:N),N)
              U(1:N,1:N)=Aux(1:N,1:N)
! #endif

            case default
              print*, "Unsupported diagonalization algorith:"//alg

        end select   

        return      

    end subroutine diagonalize_full

! #ifdef USE_LAPACK
    subroutine diasym(a,eig,n)
    !Taken from http://physics.bu.edu/py502/lectures4/examples/diatest.f90
    ! needs lapack!
    !---------------------------------------------------------!
    !Calls the LAPACK diagonalization subroutine DSYEV        !
    !input:  a(n,n) = real symmetric matrix to be diagonalized!
    !            n  = size of a                               !
    !output: a(n,n) = orthonormal eigenvectors of a           !
    !        eig(n) = eigenvalues of a in ascending order     !
    !---------------------------------------------------------!
        implicit none 

        integer,intent(in) :: n
        real(8),dimension(N,N),intent(inout) :: a
        real(8),dimension(N),intent(out) :: eig
        !Local
        integer :: l,inf
        real(8),dimension(1:1) :: work_estimate
        real(8),dimension(:),allocatable :: work 
        integer :: i

        !Needs LAPACK
        external dsyev

!          print'(/,X,A,/)', "Entering lapack diagonalization subroutine for symmetric matrices..."

        ! First estimate of the WORK size
        call dsyev('V','U',n,a,n,eig,work_estimate,-1,inf)
        ! Now we use the value estimated by LAPACK
        ! which is stored in work_estimate(1)
        l=work_estimate(1)
        allocate(work(1:l))
        call dsyev('V','U',n,a,n,eig,work,l,inf)
        if (inf /= 0) then
            print'(2X,A,I0,A,/)', "ERROR in diagonalization. (LAPACK error code: ", inf, ")"
            stop
!          else  
!              print'(2X,A,/)', "Successfull diagonalization"
        endif

    end subroutine diasym
! #endif

    subroutine orth_Lowdin(Nr,Nc,T)

        !======================================
        ! Description
        ! ------------
        ! Performs the orthogonalization of the
        ! matrix T(NrxNc). The orthogonal eigenvectors
        ! correspond to the columns: T(:,Nc).
        !
        ! Uses some Subroutines from this module
        !=======================================

        integer,intent(in) :: Nr, Nc
        real(8),dimension(:,:),intent(inout) :: T
        !Local
        real(8),dimension(1:Nr,1:Nr) :: AA
        real(8),dimension(1:Nr,1:Nr) :: BB
        real(8),dimension(1:Nr)      :: V
        !
        integer :: i, j, k, ii, jj, kk
        real(8),dimension(1:Nr,1:Nr) :: p


        ! Compute the metric matrix of T
        do i=1,Nc 
        do j=1,Nc 
            AA(i,j) = 0.d0
            do k=1,Nr
                AA(i,j)=AA(i,j)+T(k,i)*T(k,j)
            enddo
        enddo
        enddo
        ! Diagonalize metric matrix
        call diagonalize_full(AA,Nc,BB,V,"lapack")
! c     compute s^-0.5
        do ii=1,Nc
        do jj=1,Nc
        p(ii,jj)=0.d0
        do kk=1,Nc
        p(ii,jj)=p(ii,jj)+BB(jj,kk)/dsqrt(V(kk))*BB(ii,kk)
        enddo
        enddo
        enddo
! c     compute t1*s^(-1/2)
        do ii=1,Nr
        do jj=1,Nc
        BB(ii,jj)=0.d0
        do kk=1,Nc
        BB(ii,jj)=BB(ii,jj)+T(ii,kk)*p(kk,jj)
        enddo
        enddo
        enddo 
        do ii=1,Nr
        do jj=1,Nc
        T(ii,jj)=BB(ii,jj)
        enddo
        enddo

        return

    end subroutine orth_Lowdin

    subroutine JACOBI(N,RHO,FF,V,d)
    
    !=====================================================================
    ! Description
    ! Program to diagonalize a real symmetric matrix
    ! Source: EMTCCM - M1. Pais Vasco, 2010
    !
    ! Arguaments
    ! FF  is the matrix to be diagonalized (it is not modified)
    ! V   is the matrix of eigenvectors (by column)
    ! d   is the vector of eigenvalues
    ! RHO is the convergence criterium. The process ends successfully
    !     if TE<RHO, where TE is the norm of non-diagonal elements
    !
    ! Notes
    ! There is not maximum number of iterations (infinite)
    !=====================================================================
      
        implicit none
        
        integer,intent(in)::N
        real(8),dimension(:,:),intent(in)::FF
        real(8),dimension(:,:),intent(out)::V
        real(8),dimension(:),intent(out)::d
        real(8),intent(in)::RHO
        !Local
        real(8),dimension(1:N,1:N)::F !this is an auxiliar now
        real(8) :: A, COST, TE, TEN, omega, sint, U, V1, V2, V3, Z, TEM
        integer :: iter
        integer :: i, ii, ij, j, jj
        
        logical::reduced,converged
    
        print'(/,X,A)', "Entering a not so efficient diagonalization subroutine for symmetric matrices"
        print'(2X,A,E13.6,/)', "Threshol in norm:", RHO
        print'(2X,A)', "ITERATIONS..."
    
        
        ! Copy matrix to local auxiliar
        F(1:N,1:N)=FF(1:N,1:N)
        
        !C...INICIALIZACION DE LA MATRIZ DE VECTORES PROPIOS
        do I=1,N
          do J=1,N
            V(I,J)=0.0D0
          enddo
          V(I,I)=1.0D0
        enddo
        
        !C...INICIALIZACION DE LAS VARIABLES DE LA ITERACION
        A=dfloat(N)
        ITER=0
        call RMS(TE,F,N)
        TEN=TE/A
        
            
        !C...PROCESO ITERATIVO
        !Cada ciclo de iteración reduce la norma de los elementos no diagonales de forma "equilibrada"
        write(6,'(3X,A,I0)')     ' ITERATION: ',ITER
        converged=.false. 
        do while (.not.converged)
        
            ITER=ITER+1
            write(6,'(4X,A,E13.6,/)')' ERROR:     ',TE
            write(6,'(3X,A,I0)')     ' ITERATION: ',ITER
        
            reduced=.false.
            do while (.not.reduced)
            !Cada ciclo reduce el valor de los elementos que se pasen del valor medio de la norma de los elementos fuera de la diagonal
            !el ciclo de iteración se acaba cuando hayamos reducido convenientemente el valor de todos los elementos diagonales, tendremos que hacer varias pasadas
            !la variable MA controlaba si se ha producido-->cambiada por la variable lógica reduced.
        
            reduced=.true.
        
            !Buscamos en todos los elementos II,JJ no diagonales (solo la mitad, ya que es simétrica) para reducirlos 
            !hasta que todos son menores que la media de la norma (criterio adoptado)
            do II=2,N !14
                IJ=II-1
                do JJ=1,IJ !14-2
        
                !Comprobamos si el elemento no diagonal (II,JJ) cumple el criterio adoptado: si lo cumple no se aplica el algoritmo
                !si alguno no lo cumple, se repite el ciclo hasta que todos lo cumplan.
                if (DABS(F(II,JJ)).gt.TEN) then
    
                  reduced=.false.  !obliga a repetir el ciclo de reducción
            
                  !Algoritmo de Jacobi (reducción de los elementos no diagonales):
                  V1=F(JJ,JJ)
                  V2=F(II,JJ)
                  V3=F(II,II)
                  U=.5D0*(F(JJ,JJ)-F(II,II))
      
                  if (DABS(U).lt.1.d-10) then
                        OMEGA=-1.0D0
                  else
    
                    OMEGA=-F(II,JJ)/DSQRT(F(II,JJ)*F(II,JJ)+U*U)
                    Z=1.D0
                    if(U.LT.0.D0) Z=-Z
                      OMEGA=OMEGA*Z
                    endif
        
                  SINT=OMEGA/DSQRT(2.D0*(1.D0+DSQRT(1.D0-OMEGA*OMEGA)))
                  COST=DSQRT(1.D0-SINT*SINT)
        
                  !Alteramos todos los elementos de la matriz en las filas y columnas II y JJ
                  do I=1,N !13
                      !Operamos sobre la matriz F...
                      if(I.ge.II) then
                        TEM=F(I,JJ)*COST-F(I,II)*SINT
                        F(I,II)=F(I,JJ)*SINT+F(I,II)*COST
                        F(I,JJ)=TEM
                      else
                          if(I.lt.JJ) then
                              TEM=F(JJ,I)*COST-F(II,I)*SINT
                              F(II,I)=F(JJ,I)*SINT+F(II,I)*COST
                              F(JJ,I)=TEM
                          else
                              TEM=F(I,JJ)*COST-F(II,I)*SINT
                              F(II,I)=F(I,JJ)*SINT+F(II,I)*COST
                              F(I,JJ)=TEM
                          endif
                      endif
                      !Y actualizamos los vecotores propios V...
                      TEM=V(I,JJ)*COST-V(I,II)*SINT
                      V(I,II)=V(I,JJ)*SINT+V(I,II)*COST
                      V(I,JJ)=TEM
                  enddo !13
        
                  !Valor de los elementos invloucrados en el paso (II&JJ)
                  F(JJ,JJ)=V1*COST*COST+V3*SINT*SINT-2.D0*V2*SINT*COST
                  F(II,II)=V1*SINT*SINT+V3*COST*COST+2.D0*V2*SINT*COST
                  F(II,JJ)=2.D0*U*SINT*COST+V2*(COST*COST-SINT*SINT) 
    
                endif
            
              enddo !14-2
            enddo !14
        
          enddo !do while interno (ciclos de reducción)
        
        !Actualiza los valores de la norma y la norma media
        call RMS(TE,F,N)
        TEN=TE/A
        !Comprueba si ha convergido
        IF(TE.LT.RHO) then
            write(6,'(4X,A,E13.6,/)')' ERROR:     ',TE
            print'(2X,A,/)', "Successfull diagonalization"
            converged=.true.
        endif
        
        enddo !do while externo (convergencia global)
    
        d(1:N) = (/ (F(i,i), i=1,N ) /)
        
        return
    
    end subroutine JACOBI
    
    subroutine RMS(TE,F,N)

        ! This is used by JACOBI 
        
        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: F
        real(8),intent(out) :: TE
        !Local
        integer :: i, j, k
        
        TE=0.0D0
        do I=2,N
          K=I-1
          do J=1,K
            TE=TE+2.D0*F(I,J)*F(I,J)
          enddo
        enddo
        
        TE=DSQRT(TE)
        
        return
    
    end subroutine RMS

    subroutine Log_determinant_realsym(N,A,Log_det,sign_det)
    
        !=============================================================
        ! Description
        ! -------------
        ! Logaritmic version of determinant_realsym()
        ! Returns the log of the absulute value and the sign
        !==============================================================

        real(8),intent(out) :: Log_det, sign_det
        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: A
        ! Local
        real(8),dimension(N,N) :: AuxArray
        real(8) :: AuxScalar
        ! LAPACK things
        double complex,dimension(:),allocatable :: work
        integer,dimension(1:N) :: ipiv
        integer :: info, lwork
        ! Other
        integer :: alloc_status
        logical :: cycle_next=.false.
        ! Counters
        integer :: i

        ! Copy input matrix
        AuxArray(1:N,1:N) = A(1:N,1:N)

        ! Initialize LAPACK (calculating optimal size first)
        call dsytrf('L', N, AuxArray, N, ipiv, AuxScalar, -1, info)
        lwork = int(AuxScalar)
        allocate(work(1:lwork),stat=alloc_status) 
        if (alloc_status /= 0) then
            write(0,*) "ERROR: Memory cannot be allocated (DETERMINANT ROUTINE)"
            stop
        endif

        ! Perform factorization
        call dsytrf('L', N, AuxArray, N, ipiv, work, lwork, info)
        if (info < 0) then
            write(0,*) "ERROR IN dsytrf (DETERMINANT ROUTINE)"
            stop
        endif

        ! And compute the determinant
        Log_det  = 0.d0
        sign_det = 1.d0
        do i=1,N
            if (cycle_next) then
                cycle_next=.false.
                cycle
            endif
            if (ipiv(i) < 0 .and. ipiv(i) == ipiv(i+1)) then
                ! Compute the 2x2 block
                cycle_next=.true.
                AuxScalar = AuxArray(i,i)*AuxArray(i+1,i+1) - AuxArray(i+1,i)**2
                Log_det  = Log_det + Log(abs((AuxScalar)))
                sign_det = sign_det*sign(1.d0,AuxScalar)
             else if (ipiv(i) < 0 .and. ipiv(i) == ipiv(i-1)) then
                ! Second part of a 2x2 block
                cycle
             else
                AuxScalar = AuxArray(i,i)
                Log_det  = Log_det + Log(abs((AuxScalar)))
                sign_det = sign_det*sign(1.d0,AuxScalar)
            endif
        enddo
        
        return

    end subroutine Log_determinant_realsym

    !-----------------------------------------
    ! WRAPPER FUNCTION TO BLAS MATRIX PRODUCTS
    !-----------------------------------------
    function vector_dot_product(N,V1,V2) result(pes)
        ! A Wrapper to ddot to perform the dot product 
        ! Pof vectors V1(N) and V2(N)

         
        ! this definition does not work
!         INTERFACE 
!            function ddot(N,Dx,INCX,DY,INCY)
!              integer :: N, INCX, INCY
!              real(8),dimension(:) :: Dx, Dy
!            END FUNCTION ddot
!         END INTERFACE
   
        integer,intent(in)                   :: N
        real(8),dimension(:),intent(in)      :: V1,V2
        real(8)                              :: pes

        !Needs BLAS !as weird as it is, we need to define it as a real
        real(8) :: ddot
        external :: ddot
    
        pes = ddot(N,V1,1,V2,1)

        return

    end function vector_dot_product
    
    
    function vector_dot_product_ne(N,V1,V2,G) result(pes)
    
        ! Generalized dot product for non-orthogonal spaces
        ! i.e., with co-variant metric tensor not equal to identity
        ! (G matrix)
        
        integer,intent(in)                :: N
        real(8),dimension(:),intent(in)   :: V1,V2
        real(8),dimension(:,:),intent(in) :: G
        !Result
        real(8) :: pes 
        
        !Local
        integer :: i,j
        
        pes=0.d0
        do i=1,N
        do j=1,N
            pes = pes + v1(i) * G(i,j) * v2(j)
        enddo
        enddo
        
        return
    
    end function vector_dot_product_ne

    function matrix_vector_product(M,N,A,v,tA) result(p)

        !-----------------------------------
        !Multiply A array by v vector
        ! A(M,N) (in any case)
        ! v(N) (if tA=.false.)
        ! v(M) (if tA=.true.)
        !-----------------------------------

        integer,intent(in)                   :: M,N
        real(8),dimension(:,:),intent(in)    :: A
        real(8),dimension(:),intent(in)      :: v
        logical,intent(in),optional          :: tA
        real(8),dimension(:),allocatable     :: p
        !Local
        real(8),dimension(:),allocatable     :: paux
        character :: TRANS
        integer :: LDA
        logical :: tA_local

        !Needs BLAS
        external dgemv
        
        tA_local = .false.
        if (present(tA)) then
            tA_local = tA
        endif

        !First dimension as specified in the calling program
        LDA = size(A,1)
!         allocate(paux(1:max(M,N)))
        
        if (tA_local) then
            TRANS="T"
            allocate(paux(1:M))
            allocate(p(1:M))
        else
            TRANS="N"
            allocate(p(1:N))
            allocate(paux(1:N))
        endif

        call dgemv (TRANS, M, N, 1.d0, A, LDA, v, 1, 0.d0, paux, 1)
        
        p=paux

        return

    end function matrix_vector_product


    function matrix_product(NA,NB,NK,A,B,tA,tB) result(P)

        ! A Wrapper to dgemm to multiply 
        ! matrices
        ! P(NA,NB) = A(NA,K) * B(K,NB)
        ! P(NA,NB) = A^t(K,NA) * B(K,NB)
        ! and modificaitions alike

        integer,intent(in)                   :: NA,NB,NK
        real(8),dimension(:,:),intent(in)    :: A,B
        logical,intent(in),optional          :: tA,tB
        real(8),dimension(NA,NB)             :: P
        !Local
        character :: opA,opB
        integer :: LDA, LDB
        real(8),dimension(:,:),allocatable   :: Aaux,Baux
        logical :: tA_local, tB_local

        !Needs BLAS
        external dgemm
        
        tA_local = .false.
        if (present(tA)) then
            tA_local = tA
        endif
        tB_local = .false.
        if (present(tB)) then
            tB_local = tB
        endif

        !First dimension as specified in the calling program
        ! equal to NA and NB if opA,opB = 'N'
        LDA = size(A,1)
        LDB = size(B,1)
        opA = 'N'
        opB = 'N'

        !Allcate auxiliars (to avoid runtime warning with gfortran:
        ! Fortran runtime warning: An array temporary was created
        allocate(Aaux(size(A,1),size(A,2)))
        allocate(Baux(size(B,1),size(B,2)))
        Aaux=A
        Baux=B
        
        if (tA_local) opA='T'
        if (tB_local) opB='T'

        call dgemm(opA,opB,NA,NB,NK,1.d0,Aaux,LDA,Baux,LDB,0.d0,P,NA)  
        return

!         if (opA == 'N'.and. opB == 'N') then    
!             do i=1,NA
!             do j=1,NB
!                 P(i,j) = 0.d0
!                 do k=1,NK
!                     P(i,j) = P(i,j) + A(i,k) * B(k,j)
!                 enddo
!             enddo
!             enddo
!         else if (opA == 'T'.and. opB == 'N') then  
!             do i=1,NA
!             do j=1,NB
!                 P(i,j) = 0.d0
!                 do k=1,NK
!                     P(i,j) = P(i,j) + A(k,i) * B(k,j)
!                 enddo
!             enddo
!             enddo
!         else if (opA == 'N'.and. opB == 'T') then  
!             do i=1,NA
!             do j=1,NB
!                 P(i,j) = 0.d0
!                 do k=1,NK
!                     P(i,j) = P(i,j) + A(i,k) * B(j,k)
!                 enddo
!             enddo
!             enddo
!         else if (opA == 'T'.and. opB == 'T') then  
!             do i=1,NA
!             do j=1,NB
!                 P(i,j) = 0.d0
!                 do k=1,NK
!                     P(i,j) = P(i,j) + A(k,i) * B(j,k)
!                 enddo
!             enddo
!             enddo
!         endif

        return

    end function matrix_product

    function matrix_basisrot(M,N,X,A,counter) result(P)

        !=============================================
        ! Rotation by X of basis set A: 
        ! NORMAL (counter=.false.)
        ! P = X A X^t, where X(M,N) and A(N,N)
        ! INVERSE (counter=.true.)
        ! P = X^t A X, where X(N,M) and A(N,N)
        !=============================================

        integer,intent(in)                   :: M,N
        real(8),dimension(:,:),intent(in)    :: X,A
        logical,intent(in),optional          :: counter
        real(8),dimension(M,M)               :: P
        !Local
        real(8),dimension(M,N)               :: Aux
        logical :: counter_local
        
        counter_local = .false.
        if (present(counter)) then
            counter_local=counter
        endif

        if (counter_local) then
            Aux = matrix_product(M,N,N,X,A,tA=.true.)
            P   = matrix_product(M,M,N,Aux,X)
        else
            Aux = matrix_product(M,N,N,X,A)
            P   = matrix_product(M,M,N,Aux,X,tB=.true.)
        endif

        return

    end function matrix_basisrot

    function diag_basisrot(M,N,X,a,counter) result(P)

        !=============================================
        ! Rotation by X of basis set A: 
        ! NORMAL (counter=.false.)
        ! P = X A X^t, where X(M,N) and A(N,N)
        ! INVERSE (counter=.true.)
        ! P = X^t A X, where X(N,M) and A(N,N)
        !
        ! A is diagonal, and the subroutine uses 
        ! the vector a=diag(A) as input
        !=============================================

        integer,intent(in)                   :: M,N
        real(8),dimension(:,:),intent(in)    :: X
        real(8),dimension(:),intent(in)      :: a
        logical,intent(in),optional          :: counter
        real(8),dimension(M,M)               :: P
        !Local
        real(8),dimension(M,N)               :: Aux
        integer :: j
        logical :: counter_local
        
        counter_local = .false.
        if (present(counter)) then
            counter_local=counter
        endif

        if (counter_local) then
        ! NOT TESTED
            do j=1,N
                Aux(:,j) = X(j,1:M)*a(j)
            enddo
            P   = matrix_product(M,M,N,Aux,X)
        else
            do j=1,N
                Aux(:,j) = X(1:M,j)*a(j)
            enddo
            P   = matrix_product(M,M,N,Aux,X,tB=.true.)
        endif


        return

    end function diag_basisrot




     ! SORTING ROUTINES

    !< Sort diagonal of a matrix
    subroutine sort_fdiag(N,D,V,IORD,descending)
    !=========================================================
    ! Sorts eigenvalues (from Lower to Higher) in diag(D) and
    ! stores the values in V and the permutations in IORD
    !=========================================================
    
        implicit none

        integer,intent(in) :: N
        real(8), DIMENSION(:,:),intent(in) :: D
        real(8), DIMENSION(:),intent(out)  :: V
        integer,dimension(:),intent(out),optional :: IORD
        logical,intent(in),optional :: descending
    
        real(8) :: vsel
        integer :: isel, iordsel
        integer :: i,j
        logical :: rev=.false.

        if (present(descending)) then
            rev=descending
        endif

        if (present(IORD)) then
            do i=1,N
                IORD(i) = i
            enddo
        endif
    
        do i=1,N-1
            vsel=D(i,i)
            isel=i
            do j=i+1,N
                ! Normal ordering (Low to high)
                if (D(j,j)<vsel .and. .not.rev) then
                    vsel=D(j,j)
                    isel=j
                ! Reverse ordering (High to low)
                elseif (D(j,j)>vsel .and. rev) then
                    vsel=D(j,j)
                    isel=j
                endif
            enddo
            V(i)    = vsel
            if (present(IORD)) then
                iordsel    = IORD(isel)
                IORD(isel) = IORD(i)
                IORD(i)    = iordsel
            endif
        enddo
    
        return

    END SUBROUTINE sort_fdiag


    !> Sort vector of real numbers
    !! (ascending or descending order)
    subroutine sort_fvec(N,V,IORD,descending)
    
        ! From min to max is the default but can 
        ! be controlled by the descending option 
        ! (default is .false.
        ! 
        ! Works with real arrays
    
        implicit none
    
        integer,intent(in) :: N
        real(8),dimension(:),intent(inout) :: V
        integer,dimension(:),intent(out),optional :: IORD
        logical,intent(in),optional :: descending
    
        real(8) :: vsel
        integer :: isel, iordsel
        integer :: i,j
        logical :: rev=.false.
    
        if (present(descending)) then
            rev=descending
        endif

        if (present(IORD)) then
            do i=1,N
                IORD(i) = i
            enddo
        endif
    
        do i=1,N-1
            vsel=V(i)
            isel=i
            do j=i+1,N
                ! Normal ordering (Low to high)
                if (V(j)<vsel .and. .not.rev) then
                    vsel=V(j)
                    isel=j
                ! Reverse ordering (High to low)
                elseif (V(j)>vsel .and. rev) then
                    vsel=V(j)
                    isel=j
                endif
            enddo
            V(isel) = V(i)
            V(i)    = vsel
            if (present(IORD)) then
                iordsel    = IORD(isel)
                IORD(isel) = IORD(i)
                IORD(i)    = iordsel
            endif
        enddo
    
        return
    
    end subroutine sort_fvec

    !> Sort vector of integer numbers
    !! (ascending or descending order)
    subroutine sort_ivec(N,V,IORD,descending)
    
        ! From min to max is the default but can 
        ! be controlled by the descending option 
        ! (default is .false.
        ! 
        ! Works with integer arrays
    
        implicit none
    
        integer,intent(in) :: N
        integer,dimension(:),intent(inout) :: V
        integer,dimension(:),intent(out),optional :: IORD
        logical,intent(in),optional :: descending
    
        integer :: vsel
        integer :: isel, iordsel
        integer :: i,j
        logical :: rev=.false.
    
        if (present(descending)) then
            rev=descending
        endif

        if (present(IORD)) then
            do i=1,N
                IORD(i) = i
            enddo
        endif
    
        do i=1,N-1
            vsel=V(i)
            isel=i
            do j=i+1,N
                ! Normal ordering (Low to high)
                if (V(j)<vsel .and. .not.rev) then
                    vsel=V(j)
                    isel=j
                ! Reverse ordering (High to low)
                elseif (V(j)>vsel .and. rev) then
                    vsel=V(j)
                    isel=j
                endif
            enddo
            V(isel) = V(i)
            V(i)    = vsel
            if (present(IORD)) then
                iordsel    = IORD(isel)
                IORD(isel) = IORD(i)
                IORD(i)    = iordsel
            endif
        enddo
    
        return
    
    end subroutine sort_ivec


    
    subroutine rotate_angle_3D(vx,vy,vz,tx,ty,tz,Theta) 
    
        !==============================================
        !Description:
        ! ------------
        ! Subroutine to rotate the vector (vx,vy,vz) around the axis
        ! defined by (tx,ty,tz) an angle Theta (rad).
        !==============================================
    
        real(8), intent(inout) :: vx,vy,vz
        real(8), intent(in)    :: tx,ty,tz, Theta
    
        !Local
        real(8),dimension(1:3,1:3) :: R
        real(8) :: vx_tmp, vy_tmp, tmod
        real(8) :: ux, uy, uz 
    
        ! Vector u must be unitary
        tmod = sqrt(tx**2 + ty**2 + tz**2)
        ux = tx/tmod
        uy = ty/tmod
        uz = tz/tmod
    
        ! Form 3D-rotation matrix (from Wikipedia)
        R(1,1) = cos(Theta) + ux**2*(1.0-cos(Theta))
        R(1,2) = ux*uy*(1.0-cos(Theta)) - uz*sin(Theta)
        R(1,3) = ux*uz*(1.0-cos(Theta)) + uy*sin(Theta)
        R(2,1) = ux*uy*(1.0-cos(Theta)) + uz*sin(Theta)
        R(2,2) = cos(Theta) + uy**2*(1.0-cos(Theta))
        R(2,3) = uy*uz*(1.0-cos(Theta)) - ux*sin(Theta)
        R(3,1) = ux*uz*(1.0-cos(Theta)) - uy*sin(Theta)
        R(3,2) = uy*uz*(1.0-cos(Theta)) + ux*sin(Theta)
        R(3,3) = cos(Theta) + uz**2*(1.0-cos(Theta))
    
        ! Apply rotaion
        vx_tmp = vx*R(1,1) + vy*R(1,2) + vz*R(1,3)
        vy_tmp = vx*R(2,1) + vy*R(2,2) + vz*R(2,3)
        vz =     vx*R(3,1) + vy*R(3,2) + vz*R(3,3)
        vx = vx_tmp 
        vy = vy_tmp 
    
       return
    
    end subroutine rotate_angle_3D


    function rotate3D_matrix(N,M,A,Rot,tA,tR) result(Arot)
    
        !==============================================
        ! Description
        ! ------------
        ! Rotate a NxM matrix according to:
        ! Arot = R * A
        !  where R is a 3x3 rotation matrix
        !  (actually 3N'x3N' formed by repeating
        !  the same 3x3 block)
        !  and A is (NxM), where N=3*N'
        ! Allows to transpose the input A and R matrices
        ! with an optional argument 
        !
        ! Notes
        ! -----
        ! This is not the same as rotation_3D, which
        ! rotates 3D vectors. This routine rotates
        ! 3N'xM matrices
        !=================================================
      
        integer,intent(in)                     :: N,M
        real(8),dimension(:,:),intent(inout)   :: A
        real(8),dimension(1:3,1:3),intent(in)  :: rot
        real(8),dimension(N,M)                 :: Arot
        logical,intent(in),optional            :: tA
        logical,intent(in),optional            :: tR
        !Local
        ! scalar
        integer :: i, j, k, ii, jj, kk
        logical :: tA_local, tR_local
        
        tA_local = .false.
        if (present(tA)) then
            tA_local = tA
        endif
        tR_local = .false.
        if (present(tR)) then
            tR_local = tR
        endif
      
        if (tA_local) then
            if (tR_local) then
                do i=1,N
                do j=1,M
                    Arot(i,j) = 0.d0
                    ! ii runs over the 3x3 rot matrix
                    ii=mod(i+2,3)+1
                    ! kk selects the block of the NxN matrix
                    kk = 3*((i-1)/3) + 1
                    do k=1,3
                        Arot(i,j) = Arot(i,j) + Rot(k,ii) * A(j,kk)
                        kk=kk+1 
                    enddo
                enddo
                enddo
            else
                do i=1,N
                do j=1,M
                    Arot(i,j) = 0.d0
                    ! ii runs over the 3x3 rot matrix
                    ii=mod(i+2,3)+1
                    ! kk selects the block of the NxN matrix
                    kk = 3*((i-1)/3) + 1
                    do k=1,3
                        Arot(i,j) = Arot(i,j) + Rot(ii,k) * A(j,kk)
                        kk=kk+1 
                    enddo
                enddo
                enddo
            endif
        else
            if (tR_local) then
                do i=1,N
                do j=1,M
                    Arot(i,j) = 0.d0
                    ! ii runs over the 3x3 rot matrix
                    ii=mod(i+2,3)+1
                    ! kk selects the block of the NxN matrix
                    kk = 3*((i-1)/3) + 1
                    do k=1,3
                        Arot(i,j) = Arot(i,j) + Rot(k,ii) * A(kk,j)
                        kk=kk+1 
                    enddo
                enddo
                enddo
            else
                do i=1,N
                do j=1,M
                    Arot(i,j) = 0.d0
                    ! ii runs over the 3x3 rot matrix
                    ii=mod(i+2,3)+1
                    ! kk selects the block of the NxN matrix
                    kk = 3*((i-1)/3) + 1
                    do k=1,3
                        Arot(i,j) = Arot(i,j) + Rot(ii,k) * A(kk,j)
                        kk=kk+1 
                    enddo
                enddo
                enddo
            endif
        endif
    
        return
    
    end function rotate3D_matrix

    function rotate3D_vector(N,V,Rot,tR) result(Vrot)
    
        !==============================================
        ! Description
        ! ------------
        ! Rotate a N vector according to:
        ! Vrot = R * V
        !  where R is a 3x3 rotation matrix
        !  (actually 3N'x3N' formed by repeating
        !  the same 3x3 block)
        !  and V is (N), where N=3*N'
        !
        ! Notes
        ! -----
        ! Derived from rotate3D_matrix
        !=================================================
      
        integer,intent(in)                     :: N
        real(8),dimension(:),intent(inout)     :: V
        real(8),dimension(1:3,1:3),intent(in)  :: rot
        real(8),dimension(N)                   :: Vrot
        logical,intent(in),optional            :: tR
        !Local
        ! scalar
        integer :: i, j, k, ii, jj, kk
        logical :: tR_local
        
        tR_local = .false.
        if (present(tR)) then
            tR_local = tR
        endif

        if (tR_local) then
            do i=1,N
                Vrot(i) = 0.d0
                ! ii runs over the 3x3 rot matrix
                ii=mod(i+2,3)+1
                ! kk selects the block of the NxN matrix
                kk = 3*((i-1)/3) + 1
                do k=1,3
                    Vrot(i) = Vrot(i) + Rot(k,ii) * V(kk)
                    kk=kk+1 
                enddo
            enddo
        else
            do i=1,N
                Vrot(i) = 0.d0
                ! ii runs over the 3x3 rot matrix
                ii=mod(i+2,3)+1
                ! kk selects the block of the NxN matrix
                kk = 3*((i-1)/3) + 1
                do k=1,3
                    Vrot(i) = Vrot(i) + Rot(ii,k) * V(kk)
                    kk=kk+1 
                enddo
            enddo
        endif

        return
    
    end function rotate3D_vector


    function angle(vec1,vec2,n)
        real(8),intent(in),dimension(:) :: vec1, vec2
        integer,intent(in) :: n
        !Local
        real(8) :: psc,vn1,vn2,angle
        integer :: ii
        
        vn1=0.d0
        vn2=0.d0
        psc=0.d0
        do ii=1,n
            psc=psc+vec1(ii)*vec2(ii)
            vn1=vn1+vec1(ii)*vec1(ii)
            vn2=vn2+vec2(ii)*vec2(ii)
        enddo

        if(vn1.lt.1.d-10.or.vn2.lt.1.d-10) then
            angle=0.d0
        else
            angle=dacos(psc/dsqrt(vn1*vn2))*180.d0/dacos(-1.d0)
        endif

        return
    end function angle

end module matrix
