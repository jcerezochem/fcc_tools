module fft 

    use constants

    contains

    subroutine ft_calc(ndata,t,xi,spec,optrans)
    
        !=============================================================
        !DESCRIPTION:
        !------------
        !Subroutine to tool to perform BACKWARD (for absorption) and
        !FORWARD (for emission) FFT using the Cooley-Tukey algorith
        !adapted from RosetaCode. Only FFT is done, so it is compulsory
        !that the dimension of the array is a power of 2
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
    
        integer,       parameter :: dp=selected_real_kind(15,300)
    
        double complex,parameter :: Im=(0.d0,1.d0)
    
        !Subroutine I/O
        integer,intent(in) :: ndata
        real(8),dimension(:),intent(in) :: t
        double complex,dimension(:),intent(in) :: xi
        character(len=*),intent(in) :: optrans
        real(8),dimension(:),intent(out) :: spec
    
        integer :: i, k, ios
        double complex,dimension(1:ndata) :: datos
        real(8)    :: p,q, df, dt, t0, log2
        real(8)    :: autown,evtown,fstoev,pi
        !FFT interface stuff
        integer :: FFT_FORWARD  =-1, &
                   FFT_BACKWARD = 1, &
                   FFT_DIRECTION
    
        print'(/,X,A)', "Entering in fft_calc subroutine (rosettacode)"
        print*, "Num data (total)", ndata
        print*, "Num data for DFT", ndata
        log2 = log(dfloat(ndata))/log(2.d0)
        if ( 2**INT(log2) == ndata ) then
            print'(X,A,/)', "Num. data for DFT is a power of 2. Proceeding..."
        else
           print*, "Data for DFT is not a power of 2, but this routine only performs FFT. Aborting"
           stop
        endif
    
        !Read data (from stdin): three colum array
        datos(1:ndata) = xi(1:ndata)
    
        if (optrans == "emi") then
            FFT_DIRECTION=FFT_FORWARD
        else
            FFT_DIRECTION=FFT_BACKWARD
        endif
        call fft(datos(1:ndata),FFT_DIRECTION)
    
        !-------------------------------------------
        !What is computed (BACKWARD). The DFT as:
        ! spec(k) = \sum_{j=0}^{N-1} datos(j)*exp(+Im*2pi*(j-1)*k/N)
        !without any frefactor
        !
        !How the DFT is stored:
        ! spec(1:N/2)   : 0 and positive frequencies
        ! spec(N/2+1:N) ; negagive frequencies
        !-------------------------------------------
    
        !Relevant data in the time and freq domains
        t0 = t(1)
        dt = (t(ndata)-t(1))/(ndata-1)
        !df is (normal) frequency, not angular frequency (w)
        df = 1.d0/(dt*float(ndata))
        write(6,*) "df", df !, df*N/2*fstoev
        write(6,*) "dt", dt !, df*N/2*fstoev
    
        !To get FT from DFT, we include a phase factor, so that, FT \appox Factor * DFT
        ! Factor(k) = exp(-i dw*k * t0)*DFT
        !where k=(0,1,...,N/2) and dw = 2pi*df (see NOTES)
        !Note that since we are interested in the real part, the sign in front of Im does 
        !not matter (for the imag part it should be fixed with FFT_DIRECTION)
        !Since input data have time in fs instead of au, changes dt to au with /autofs
        do k=1,ndata/2
            datos(k) = datos(k) * dt * exp(-Im * 2.d0*pi*df * t0 * real(k-1)) /autofs
        end do
    
        !Write spectrum
        spec(1:ndata/2) = real(datos(1:ndata/2))
        
        return
    
        contains
     
        recursive subroutine fft(x,isgn)
    
            !=================================================================================
            !In place Cooley-Tukey FFT (adapted from http://rosettacode.org/wiki/FFT#Fortran) 
            !
            !NOTES:
            !The subroutine is not very stable. For instance, if the variable N is named Ndata 
            !instead (and in main we use N), a segfault arises. Maybe due to a typo but 
            !anyway, some caution should be taken
            !=================================================================================
    
            implicit none
      
            complex(8), dimension(:), intent(inout)  :: x
            integer,intent(in)                       :: isgn
            !Local
            real(8)                                  :: sgn
            complex(8)                               :: t
            integer                                  :: N
            integer                                  :: i
            complex(8), dimension(:), allocatable    :: even, odd
    
            N=size(x)
      
            if(n.le.1) return
      
            allocate(odd((N+1)/2))
            allocate(even(N/2))
      
            ! divide
            odd =x(1:N:2)
            even=x(2:N:2)
      
            ! conquer
            call fft(odd,isgn)
            call fft(even,isgn)
      
            ! combine
            sgn = real(isgn,8)
            do i=1,N/2
                t=exp(sgn*cmplx(0.0d0,-2.0d0*pi*real(i-1,8)/real(N,8),8))&
                  *even(i)
                x(i)     = odd(i) + t
                x(i+N/2) = odd(i) - t
            end do
      
            deallocate(odd)
            deallocate(even)
      
        end subroutine fft
     
    end subroutine ft_calc

end module fft
