module MatrixMod

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS
    !==============================================================
    !
    ! Despription
    ! -----------
    ! Subroutines to perform matrix manipulations
    ! This is a reduced version. Only double precision
    ! floats are used (do not preprocess is needed, so as
    ! to avoid issues with compiler that might not accept
    ! preprocessing)
    !==============================================================


!MODULO CON LA SUBRUTINA DE DIAGONALIZACIÓN DE JACOBI Y RELACIONADAS
!Todas las funciones y subrutinas están definidas en doble precisión
!Fuente: material del curso intensivo del Master en Química Teórica y Modelización Computacional (Semana 2). San Sebastián, Enero-Febrero 2010 
! Original en F77 --> Reescrito en F90. Se reescribe en f90 para poder compilarlo junto al programa principal sin tener problemas (gfortran). 
! Otra opción con este compilador hubiera sido compilar la subrutina por separado (crear un objeto). El objeto se enlaza (link) al compilar el 
! programa principal.

    CONTAINS

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
! #ifdef DOUBLE
        double precision,dimension(:,:), intent(in) :: A
        !Output
        double precision,dimension(:,:),intent(out) :: U
        double precision,dimension(:),intent(out) :: d
! #else
!         real,dimension(:,:), intent(in) :: A
!         !Output
!         real,dimension(:,:),intent(out) :: U
!         real,dimension(:),intent(out) :: d
! #endif   

        !Local
! #ifdef DOUBLE
        double precision, parameter :: PREC=1.d-11
! #else
!         real, parameter :: PREC=1.d-11
! #endif
        integer :: IER

        !Auxiliar
        character(len=50) :: dummy_char
        integer :: i,j,k

        

        select case (trim(adjustl(alg)))
            !The one from the master (M1 year)
            case("emtccm")
              !We need to feed the whole matrix, as we have
              call JACOBI_SIM(N,PREC,A,U,d)

! #ifdef DOUBLE
! #ifdef USE_LAPACK
            !Only available if compiled in duble precision and if lapack libs are required
            case("lapack")
              !We need to feed the whole matrix, as we have
              !We don't destroy the original
              U(1:N,1:N)=A(1:N,1:N)
              call diasym(U,d,N)
! #endif
! #endif

            case default
              print*, "Unsupported algorith:"//alg

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
         integer l,inf
         real*8  work(n*(3+n/2)) 

         print'(/,X,A,/)', "Entering lapack diagonalization subroutine for symmetric matrices"

         l=n*(3+n/2)
         call dsyev('V','U',n,a,n,eig,work,l,inf)

     end subroutine diasym
! #endif


!C-----------------------------------------------------------------------
SUBROUTINE JACOBI_SIM(N,RHO,FF,V,d)
!C-----------------------------------------------------------------------
!CPROGRAMA PARA DIGONALIZAR UNA MATRIZ REAL Y SIMETRICA DE ORDEN N
!C
!C  F ES LA MATRIZ A DIAGONALIZAR. (NO SE MODIFICA -- nueva versión)
!C  V ES LA MATRIZ DE VECTORES PROPIOS (CADA COLUMNA UN VECTOR)
!C  d  es el vector de valores propios
!C  RHO ES EL CRITERIO DE CONVERGENCIA. EL PROCESO FINALIZA SI TE<RHO
!C    DONDE TE ES LA NORMA DE LOS ELEMENTOS NO DIAGONALES.
!C 

    !use matmod !Para la escritura de las matrices con WriteMatrix
      
! #ifdef DOUBLE
    implicit double precision (A-H,O-Z)
! #else
!     implicit real (A-H,O-Z)
! #endif
    implicit integer (i-n)
    
    integer,intent(in)::N
! #ifdef DOUBLE
    double precision,dimension(:,:),intent(inout)::V
    double precision,dimension(:,:),intent(in)::FF
    double precision,dimension(:),intent(out)::d
    double precision,intent(in)::RHO
    double precision,dimension(N,N)::F !this is an auxiliar now
! #else
!     real,dimension(:,:),intent(inout)::V
!     real,dimension(:,:),intent(in)::FF
!     real,dimension(:),intent(out)::d
!     real,dimension(N,N)::F !this is an auxiliar now
!     real,intent(in)::RHO
! #endif
    
    logical::reduced,converged

    print'(/,X,A,/)', "Entering a not so efficient diagonalization subroutine for symmetric matrices"
    
    ! Copy matrix to local auxiliar
    F=FF
    
    !C...INICIALIZACION DE LA MATRIZ DE VECTORES PROPIOS
    do I=1,N
      do J=1,N
        V(I,J)=0.0D0
      enddo
      V(I,I)=1.0D0
    enddo
    
    !C...INICIALIZACION DE LAS VARIABLES DE LA ITERACION
    A=N
    ITER=0
    call RMS(TE,F,N)
    TEN=TE/A
    
        
    !C...PROCESO ITERATIVO
    !Cada ciclo de iteración reduce la norma de los elementos no diagonales de forma "equilibrada"
    converged=.false. 
    do while (.not.converged)
    
        ITER=ITER+1
!         write(6,*)' ITERACION ..',ITER
!         write(6,*)' ERROR ......',TE
    
        reduced=.false.
        do while (.not.reduced)
        !Cada ciclo reduce el valor de los elementos que se pasen del valor medio de la norma de los elementos fuera de la diagonal
        !el ciclo de iteración se acaba cuando hayamos reducido convenientemente el valor de todos los elementos diagonales, tendremos que hacer varias pasadas
        !la variable MA controlaba si se ha producido-->cambiada por la variable lógica reduced.
    
        reduced=.true.
    
        !Buscamos en todos los elementos II,JJ no diagonales (solo la mitad, ya que es simétrica) para hacer reducirlos 
        !hasta que todos son menores que la media de la norma (criterio adoptado)
        do II=2,N !14
            IJ=II-1
            do JJ=1,IJ !14-2
    
            !Comprobamos si el elemento no diagonal (II,JJ) cumple el criterio adoptado: si lo cumple no se aplica el algoritmo
            !si alguno no lo cumple, se repite el ciclo hasta que todos lo cumplan.
! #ifdef DOUBLE
            if (DABS(F(II,JJ)).gt.TEN) then
! #else
!             if (ABS(F(II,JJ)).gt.TEN) then
! #endif

              reduced=.false.  !obliga a repetir el ciclo de reducción
        
              !Algoritmo de Jacobi (reducción de los elementos no diagonales):
              V1=F(JJ,JJ)
              V2=F(II,JJ)
              V3=F(II,II)
              U=.5D0*(F(JJ,JJ)-F(II,II))
  
! #ifdef DOUBLE
              if (DABS(U).lt.1.d-10) then
! #else
!               if (ABS(U).lt.1.d-10) then
! #endif  
                    OMEGA=-1.0D0
              else
! #ifdef DOUBLE
                OMEGA=-F(II,JJ)/DSQRT(F(II,JJ)*F(II,JJ)+U*U)
! #else
!                 OMEGA=-F(II,JJ)/SQRT(F(II,JJ)*F(II,JJ)+U*U)
! #endif  
                Z=1.D0
                if(U.LT.0.D0) Z=-Z
                  OMEGA=OMEGA*Z
                endif
    
! #ifdef DOUBLE
              SINT=OMEGA/DSQRT(2.D0*(1.D0+DSQRT(1.D0-OMEGA*OMEGA)))
              COST=DSQRT(1.D0-SINT*SINT)
! #else
!               SINT=OMEGA/SQRT(2.D0*(1.D0+SQRT(1.D0-OMEGA*OMEGA)))
!               COST=SQRT(1.D0-SINT*SINT)
! #endif 
    
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
    IF(TE.LT.RHO) converged=.true.
    
    enddo !do while externo (convergencia global)
    
    
!     do I=2,N
!       II=I-1
!       do J=1,II
!         F(J,I)=F(I,J)
!       enddo
!     enddo

! print*, "Original"
! do i=1,N
!     print'(1000F8.3)', FF(i,1:N)
! enddo

! print*, "Diagonal"
! do i=1,N
!     print'(1000F8.3)', F(i,1:N)
! enddo

    d(1:N) = (/ (F(i,i), i=1,N ) /)
    
! print*, "Autovalores"
! print'(1000F8.3)', d(1:N)

    !write(6,*)' MATRIZ DIAGONAL'
      !call WriteMatrix(N,N,F,6)
        !write(6,*)' MATRIZ DE VECTORES PROPIOS (NORMALIZADOS)'
    !  call WriteMatrix(N,N,V,6)
    
    !WRITE(6,*)' *** FIN DE LA DIAGONALIZACION ***'
    
    RETURN

END SUBROUTINE JACOBI_SIM


!C-----------------------------------------------------------------------
SUBROUTINE RMS(TE,F,N)
!C-----------------------------------------------------------------------
! #ifdef DOUBLE
    implicit double precision (A-H,O-Z)
! #else
!     implicit real (A-H,O-Z)
! #endif
    implicit integer (i-n)
    
    integer,intent(in)::N
! #ifdef DOUBLE
    double precision,dimension(:,:),intent(in)::F
    double precision,intent(out)::TE
! #else
!     real,dimension(:,:),intent(in)::F
!     real,intent(out)::TE
! #endif
    
    TE=0.0D0
    do I=2,N
      K=I-1
      do J=1,K
        TE=TE+2.D0*F(I,J)*F(I,J)
      enddo
    enddo
    
! #ifdef DOUBLE
    TE=DSQRT(TE)
! #else
!     TE=SQRT(TE)
! #endif
    
    RETURN

END SUBROUTINE RMS

!C-----------------------------------------------------------------------
        SUBROUTINE ORDEN_DIAG(N,D,V)
!C-----------------------------------------------------------------------
        !Subroutina para ordenar los elemntos diagonales aplicar las mismas permutaciones a las columnas de la matriz de vecotores propios
        
                implicit none

                INTEGER, intent(in):: N
! #ifdef DOUBLE
                double precision, DIMENSION(:,:),intent(inout) :: D,V

                double precision :: Amin, aux
                double precision,DIMENSION(1:N) :: aux2
! #else
!                 real, DIMENSION(:,:),intent(inout) :: D,V
! 
!                 real :: Amin, aux
!                 real,DIMENSION(1:N) :: aux2
! #endif

                        
                integer :: i,j,imin
                        
                do  i=1,N-1
                  Amin=D(i,i)
                  imin=i
                  
                  do j=i+1,N
                    if (D(j,j).lt.Amin) then
                     Amin=D(j,j)
                     imin=j
                    endif
                  enddo

                  aux=D(i,i)
                  aux2=V(1:N,i)
                  D(i,i)=D(imin,imin)
                  V(1:N,i)=V(1:N,imin)
                  D(imin,imin)=aux
                  V(1:N,imin)=aux2
                enddo

                RETURN
        END SUBROUTINE ORDEN_DIAG

!C-----------------------------------------------------------------------
        SUBROUTINE REV_ORDEN_DIAG(N,D,V)
!C-----------------------------------------------------------------------
        !Subroutina para ordenar los elemntos diagonales aplicar las mismas permutaciones a las columnas de la matriz de vecotores propios
        !Reverse order: from Larger to Lower
        
                implicit none

                INTEGER, intent(in):: N

! #ifdef DOUBLE
                double precision, DIMENSION(:,:),intent(inout) :: D,V

                double precision :: Amax, aux
                double precision,DIMENSION(1:N) :: aux2
! #else
!                 real, DIMENSION(:,:),intent(inout) :: D,V
! 
!                 real :: Amax, aux
!                 real,DIMENSION(1:N) :: aux2
! #endif

                        
                integer :: i,j,imax
                        
                do  i=1,N-1
                  Amax=D(i,i)
                  imax=i
                  
                  do j=i+1,N
                    if (D(j,j).gt.Amax) then
                     Amax=D(j,j)
                     imax=j
                    endif
                  enddo

                  aux=D(i,i)
                  aux2=V(1:N,i)
                  D(i,i)=D(imax,imax)
                  V(1:N,i)=V(1:N,imax)
                  D(imax,imax)=aux
                  V(1:N,imax)=aux2
                enddo

                RETURN
        END SUBROUTINE REV_ORDEN_DIAG


end module MatrixMod
