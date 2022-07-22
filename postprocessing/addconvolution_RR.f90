program addconvolution_RR

    use alerts
    use constants
    use line_preprocess
    
    implicit none
    
    real(8),dimension(:),allocatable :: w,spect
    integer :: N
    character(len=8) :: hwhm_char
    real(8) :: hwhm_eV=0.01d0
    !Counters
    integer :: i, j
    
    ! Input defaults
    character(len=3) :: broadfun='GAU'
    ! Cut convolution after ntimes hwhm
    real(8) :: ncut_gau = -1.d0, &
               ncut_lor = -1.d0
    
    ! Additional data
    real(8) :: wi, ws, rr, ws_
    
    real(8) :: dw
    
    
    !IO
    character(len=200) :: RRdatafile='(default)',&
                          stickfile,RRoutfile
    integer :: I_STCK=10, &
               O_SPC =20, &
               O_STCK=21
    integer :: ioflag
    character(len=200) :: line
    character :: cnull
    
    
    ! Read command_line_options
    call parse_input(RRdatafile,hwhm_eV,broadfun,ncut_gau,ncut_lor)
    
    ! Set options upper case
    call set_word_upper_case(broadfun)
    
    ! Set dagaults
    if (RRdatafile == '(default)') then
        RRdatafile = 'RR_Spectrum_2D.dat'
    endif
    write(6,*) 'Reading data from: '//trim(adjustl(RRdatafile))
    
    if (ncut_gau < 0) ncut_gau = 20.0
    if (ncut_lor < 0) ncut_lor = 100.0
    
    ! Write spectrum to files
    write(hwhm_char,'(F8.3)') hwhm_eV
    RRoutfile="RR_Spectrum_2D_hwhm"//trim(adjustl(hwhm_char))//".dat"
    open(O_SPC,file=RRoutfile,status="replace")
    write(6,*) 'Writing data to  : '//trim(adjustl(RRoutfile))
    
    ! Read inout bins (sticks) 
    open(I_STCK,file=RRdatafile,iostat=ioflag)
    ! First get number of wI points inthe grid
    ! Skip header line
    read(I_STCK,*,iostat=ioflag) line
    read(I_STCK,*,iostat=ioflag) line
    ws_ = -1.0
    N=1
    read(I_STCK,'(A)',iostat=ioflag) line
    read(line,*) ws_, wi, rr
    do 
        read(I_STCK,'(A)',iostat=ioflag) line
        if (ioflag/=0) exit
        read(line,*) ws, wi, rr
        if (ws /= ws_) exit
        N = N+1
    enddo
    
    ! Allocate arrays
    allocate(spect(N),w(N))
    
    ! Now, read again all file and convolute each excitation profile
    rewind(I_STCK)
    ! Skip header line
    read(I_STCK,'(A)',iostat=ioflag) line
    write(O_SPC,'(A)') trim(line)
    read(I_STCK,'(A)',iostat=ioflag) line
    write(O_SPC,'(A)') trim(line)
    ws_ = -1.d0
    do 
        read(I_STCK,'(A)',iostat=ioflag) line
        if (ioflag/=0) exit
        read(line,*) ws, wi, rr
!             print*, rr
        if (ws /= ws_) then
            ws_ = ws
            i=0
            print'(A,F10.4,A)', 'Convoluting ', ws, ' ...'
        endif
        
        i = i+1
        if (i>N) call alert_msg('fatal','Inconsistency in input file')
        
        w(i) = wI
        spect(i) = rr
        
        if (i == N) then
            do i=1,N
                spect(i) = spect(i) / (w(i)-ws)**4
            enddo
            call add_convolution(w,spect,broadfun,hwhm_eV,ncut_gau,ncut_lor)
            do i=1,N
                write(O_SPC,'(F13.4,X,F15.4,X,E19.5)') ws, w(i), spect(i) * (w(i)-ws)**4
            enddo
        endif
        

    enddo
    
    write(6,'(X,A,/)') 'Done'

    !
                 
    stop
    
    contains
    
    subroutine add_convolution(w,spect,broadfun,hwhm_eV,ncut_gau,ncut_lor)
    
        use constants
    
        ! Args
        real(8),dimension(:),intent(inout) :: w,spect
        character(len=*),intent(in) :: broadfun
        real(8),intent(in)          :: hwhm_eV
        real(8),intent(in)          :: ncut_gau,ncut_lor
        
        ! Local
        real(8),dimension(:),allocatable :: y
        real(8) :: dx, hwhm_cm1, xr, brd, C
        integer :: N
        
        N = size(w)
        allocate(y(N))
        
        hwhm_cm1 = hwhm_eV * evtown
        
        dx = w(2) - w(1)
        
        ! Make convolution over the grid
        do i=1,N
            ! Fill x-grid
            y(i) = 0.d0
            do j=1,N
                xr = w(j) - w(i)
                ! Skip points beyond 20*HWHM from the center of the gaussian
                ! this increases efficiency and avoids IEEE_UNDERFLOW_FLAG IEEE_DENORMAL flags
                if (broadfun == 'GAU') then
                    if (dabs(xr)>hwhm_cm1*ncut_gau) cycle
                    C = dlog(2.d0)/hwhm_cm1**2
                    brd = dsqrt(C/pi) * dexp(-C * xr**2)
                else if (broadfun == 'LOR') then
                    if (dabs(xr)>hwhm_cm1*ncut_lor) cycle
                    C = hwhm_cm1/pi !dsqrt(pi)
                    brd = C/(xr**2 + hwhm_cm1**2)
                endif
                y(i) = y(i) + brd * spect(j) * dx
            enddo
        enddo
        
        ! Update output values
        spect = y
        deallocate(y)
    
        return
    
    end subroutine add_convolution
    
    
    subroutine parse_input(RRdatafile,hwhm_eV,broadfun,ncut_gau,ncut_lor)
    
        character(len=*),intent(inout) :: RRdatafile, broadfun
        real(8),intent(inout)          :: hwhm_eV, ncut_gau, ncut_lor

        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg

        argument_retrieved=.false.
        do i=1,iargc()
            if (argument_retrieved) then
                argument_retrieved=.false.
                cycle
            endif
            call getarg(i, arg) 
            select case (adjustl(arg))
                case ("-f") 
                    call getarg(i+1, RRdatafile)
                    argument_retrieved=.true.
                    
                case ("-brd") 
                    call getarg(i+1, broadfun)
                    argument_retrieved=.true.
                    
                case ("-hwhm") 
                    call getarg(i+1, arg)
                    read(arg,*) hwhm_eV
                    argument_retrieved=.true.
                    
                case ("-cut-lor") 
                    call getarg(i+1, arg)
                    read(arg,*) ncut_lor
                    argument_retrieved=.true.
                    
                case ("-cut-gau") 
                    call getarg(i+1, arg)
                    read(arg,*) ncut_gau
                    argument_retrieved=.true.
        
                case ("-h")
                    need_help=.true.

                case default
                    print*, "Unkown command line argument: "//adjustl(arg)
                    need_help = .true.
            end select
        enddo 
    
!Print options (to stderr)
        if (need_help) then

        write(0,'(/,A)') ' addconvolution_RR '
        write(0,'(A)'  ) '-----------------'
        write(0,'(A)'  ) 'Add a convolution (over wI axis) to the'
        write(0,'(A)'  ) '2D RR spectra. Requires a fine grid as input'

        call print_version(0)

        write(0,'(/,A)') 'SYNOPSIS'
        write(0,'(A)'  ) 'addconvolution_RR [-f (RRdatafile) -brd (broadtype) -hwhm (hwhm) -h]'

        write(0,'(/,A)')   'OPTIONS'
        write(0,'(A)'  )   'Flag       Description                  Current Value'
        write(0,'(A)'  )   ' -f        Correlation function file    '//trim(adjustl(RRdatafile))
        write(0,'(A)'  )   ' -brd      Broad func (GAU|LOR)         '//trim(adjustl(broadfun))
        write(0,'(A,F8.3)')' -hwhm     HWHM (in eV)              '   ,hwhm_eV
        write(0,'(A,F8.3)')' -cut-lor  Multiple of hwhm(eV) when '   ,ncut_lor  
        write(0,'(A)'  )   '           to cut the Lor convolution'
        write(0,'(A,F8.3)')' -cut-gau  Multiple of hwhm(eV) when '   ,ncut_gau
        write(0,'(A)'  )   '           to cut the Gau convolution'
        write(0,'(A)'  )   ' -h        Print help  '
        write(0,*) ''

        stop    
        endif

        return
    end subroutine parse_input

end program addconvolution_RR
