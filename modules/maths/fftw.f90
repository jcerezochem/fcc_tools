module fft

    contains

    subroutine ft_calc(unt,N,t,xi,spec,optrans)

        !=============================================================
        !DESCRIPTION:
        !------------
        !Subroutine to tool to perform BACKWARD (for absorption) and
        !FORWARD (for emission) DFT using FFTW3 library, so as to 
        !provide the lineshape spectrum.
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

        use constants
    
        implicit none
    
        include "fftw3.f"
    
        double complex,parameter :: Im=(0.d0,1.d0)
    
        !Subroutine I/O
        integer,intent(in) :: unt
        integer,intent(in) :: N
        real(8),dimension(:),intent(in) :: t
        double complex,dimension(:),intent(in) :: xi
        character(len=*),intent(in) :: optrans
        real(8),dimension(:),intent(out) :: spec
    
        integer :: i, k, ios
        double complex,dimension(1:N) :: datos, spec_im
        real(8)    :: df, dt, t0
        !FFTW3 interface stuff
        integer(8) :: plan, FFTW_DIRECTION
    
        !Copy data to manipulate them localy
        datos(1:N) = xi(1:N)
        
        write(unt,'(/,X,A)') "Entering in fft_calc subroutine (FFTW)"
        write(unt,'(X,A,I0)') "Ndata: ", N
        write(unt,'(X,A,2(F10.3,2X),/)') "tini and tfin", t(1), t(N)
    
        if (optrans == "EMI") then
            FFTW_DIRECTION=FFTW_FORWARD
        else
            FFTW_DIRECTION=FFTW_BACKWARD
        endif
        !Complex data DFT (from FFTW3)
        call dfftw_plan_dft_1d(plan,N,datos,spec_im,FFTW_DIRECTION,FFTW_ESTIMATE)
        call dfftw_execute(plan, datos, spec_im)
        call dfftw_destroy_plan(plan)
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
        dt = (t(N)-t(1))/(N-1)
        !df is (normal) frequency, not angular frequency (w)
        df = 1.d0/dt/dfloat(N)
        write(unt,*) "df", df !, df*N/2*fstoev
        write(unt,*) "dt", dt !, df*N/2*fstoev
    
        !To get FT from DFT, we include a phase factor, so that, FT \appox Factor * DFT
        ! Factor(k) = exp(-i dw*k * t0)*DFT
        !where k=(0,1,...,N/2) and dw = 2pi*df (see NOTES)
        !Note that since we are interested in the real part, the sign in front of Im does 
        !not matter (for the imag part it should be fixed with FFT_DIRECTION)
        !Since input data have time in fs instead of au, changes dt to au with /autofs
        do k=1,N
            spec_im(k) = spec_im(k) * dt * exp(dfloat(FFTW_DIRECTION)*Im * 2.d0*pi*df * t0 * real(k-1)) /autofs
            ! This should give the same result:
            !spec_im(k) = spec_im(k) * dt * exp(-Im * 2.d0*pi*df * t0 * real(k-1)) /autofs
        end do
    
        !Write spectrum
        spec(1:N) = real(spec_im(1:N))
        
        return
    
    end subroutine ft_calc

end module fft