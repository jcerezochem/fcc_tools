module fft 

    use constants

    contains

    subroutine ft_calc(unt,ndata,t,xi,spec,optrans)
    
        !=============================================================
        !DESCRIPTION:
        !------------
        !Subroutine to tool to perform BACKWARD (for absorption) and
        !FORWARD (for emission) DFT, using the direct application of
        !the definition:
        !
        !   g(k+1) = sum_m {x(m+1) * exp(i*2PI*m*k/N)}
        !
        !Units:
        ! -Correlation Functions: time(fs) vs. chi(atomi units)
        !  (chi contains also the terms related with dipole moment)
        ! -Spectrun: Energy(eV) vs. lineshape (atomic units)
        !
        !NOTES
        !-----
        !The DFT is turned into an approximation to the continous FT
        !when the appropriate phase factor is taken into account, as
        !described in this StackOverflow question:
        !http://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy
        !But note that the FFT must be done over N-1 points
        !=============================================================
    
        implicit none
    
        double complex,parameter :: Im=(0.d0,1.d0)
    
        !Subroutine I/O
        integer,intent(in) :: unt
        integer,intent(in) :: ndata
        real(8),dimension(:),intent(in) :: t
        double complex,dimension(:),intent(in) :: xi
        character(len=3),intent(in) :: optrans
        real(8),dimension(:),intent(out) :: spec
    
        integer :: i, k, ios
        double complex,dimension(1:ndata) :: datos
        real(8)    :: p,q, df, dt, t0, log2
        !FFT interface stuff
        integer :: FFT_FORWARD  =-1, &
                   FFT_BACKWARD = 1, &
                   FT_DIRECTION
    
        write(unt,'(/,X,A)') "Entering in dft_calc subroutine (not FFT)"
        write(unt,'(X,A,I0)') "Ndata", ndata
        write(unt,'(X,A,2(F10.3,2X),/)') "tini and tfin", t(1), t(ndata)
        
    
        !Read data (from stdin): three colum array
        datos(1:ndata) = xi(1:ndata)
    
        if (optrans == "EMI") then
            FT_DIRECTION=FFT_FORWARD
        else
            FT_DIRECTION=FFT_BACKWARD
        endif
        call dft(datos(1:ndata),FT_DIRECTION)
    
        !Relevant data in the time and freq domains
        t0 = t(1)
        dt = (t(ndata)-t(1))/float(ndata-1)
        !df is (normal) frequency, not angular frequency (w)
        df = 1.d0/(dt*float(ndata))
        write(unt,*) "df", df !, df*N/2*fstoev
        write(unt,*) "dt", dt !, df*N/2*fstoev
            
        !To get FT from DFT, we include a phase factor, so that, FT \appox Factor * DFT
        ! Factor(k) = exp(-i dw*k * t0)*DFT
        !where k=(0,1,...,ndata/2) and dw = 2pi*df (see NOTES)
        do k=1,ndata
            datos(k) = datos(k) * dt * exp(float(FT_DIRECTION) * Im * 2.d0*pi*df * t0 * real(k-1)) /autofs
        end do
    
        !Write spectrum
        do i=1,ndata
           spec(i) = real(datos(i))
        enddo
        
        return
    
        contains
     
        subroutine dft(x,isgn)
        
            !=================================================================================
            ! Direct application of DFT definition (very ineficient)
            !
            !NOTES:
            ! Returns the positive frequencies only (half the points)
            !=================================================================================
        
            implicit none
        
            double complex,dimension(:),intent(inout)  :: x
            integer,intent(in)                         :: isgn
            !Local
            integer                                    :: N 
            double complex,dimension(size(x)+1)        :: g
            integer     :: k,m
            double complex,parameter :: Im=(0.d0,1.d0)
            real(8),parameter :: pi= 4.d0*datan(1.d0)
        
            N = size(x)
        
            do k=0,N
                g(k+1) = 0.d0
                do m=0,N-1
                    g(k+1) = g(k+1) + x(m+1) * &
                             exp(Im*2.0d0*pi*dfloat(isgn*m*k)/dfloat(N))
                enddo
            enddo
        
            x(1:N) = g(1:N)
        
            return
        
        end subroutine dft
     
    end subroutine ft_calc

end module fft
