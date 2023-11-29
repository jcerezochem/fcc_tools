program convolute_RR

    use alerts
    use constants
    use line_preprocess
    
    implicit none
    
    real(8),dimension(:),allocatable :: w,spect
    real(8),dimension(:),allocatable :: ystick,xstick
    real(8) :: factor, df, dw, &
               broad, wau, dgau, dgau2, &
               C, brd, x0, xf, dx, xr
    integer :: iexp, nbroad, N
    character(len=3) :: optrans
    character(len=8) :: wI_char
    !Counters
    integer :: i, j
    
    ! Input defaults
    character(len=3) :: broadfun='LOR'
    character(len=2) :: source='1D'
    integer :: npoints=1001
    real(8) :: resol=-1.
    real(8) :: hwhm_cm1=10.d0
    real(8) :: wi_req=0.d0
    
    ! Additional data
    real(8) :: wi_sel, dev_wi, &
               wi, ws, rr, ws_,&
               ws_max, ws_min
    integer :: ni
    logical :: include_Rayleigh = .false.
    
    !IO
    character(len=200) :: RRdatafile='(default)',&
                          stickfile,convfile
    integer :: I_STCK=10, &
               O_SPC =20, &
               O_STCK=21
    integer :: ioflag
    character(len=200) :: line
    character :: cnull
    
    
    ! Read command_line_options
    call parse_input(RRdatafile,npoints,resol,hwhm_cm1,broadfun,source,wi_req,include_Rayleigh)
    
    ! Set options upper case
    call set_word_upper_case(source)
    call set_word_upper_case(broadfun)
    
    if (RRdatafile == '(default)') then
        if (source == '2D') then
            RRdatafile = 'RR_Spectrum_2D.dat'
        else
            RRdatafile = 'RR_Spectrum_VertE.dat'
        endif
        write(6,*) 'Reading data from: '//trim(adjustl(RRdatafile))
    endif
    
    ! Read inout bins (sticks) 
    open(I_STCK,file=RRdatafile,iostat=ioflag)
    if (ioflag/=0) call alert_msg('fatal','Cannot open RR datafile: '//trim(adjustl(RRdatafile)))
    if (source == '2D') then
        ! First get N=ns and select the closest indicident freq
        ws = 0.d0
        dev_wi = 1.e10
        N = 1
        ni = 0
        ws_max = 0.d0
        ! Skip header line
        read(I_STCK,*,iostat=ioflag) line
        read(I_STCK,*,iostat=ioflag) line
        do
            ws_ = ws
            read(I_STCK,'(A)',iostat=ioflag) line
            if (ioflag/=0) exit
            read(line,*) ws, wi, rr
            if (ws > ws_max) then
                ws_max = ws
            endif
            if (ws_ /= ws ) then
                N = N+1
            endif
            ! This is only done for the first ws value
            if (N == 1) then
                ni = ni + 1
                if (dabs(wi-wi_req) < dev_wi) then
                    wi_sel = wi
                    dev_wi = dabs(wi-wi_req)
                endif
            elseif (N == 2) then
                ws_min = ws
            endif
        enddo
        rewind(I_STCK)
    else if (source == '1D') then
        read(I_STCK,*,iostat=ioflag) line
        read(I_STCK,'(A)',iostat=ioflag) line
        call split_line(line,'maximum',cnull,line)
        read(line,*) wi_sel, cnull
        read(I_STCK,*,iostat=ioflag) line
        ! Change from eV to cm-1
        wi_sel = wi_sel * evtown
        i = 0
        do
            read(I_STCK,'(A)',iostat=ioflag) line
            if (ioflag/=0) exit
            read(line,*) ws_max
            if (i==1) ws_min = ws_max
            i = i+1
        enddo
        N = i
        rewind(I_STCK)
    else
        call alert_msg('fatal','Unkown -type option')
    endif
    write(6,'(X,A,F10.2,A,/)') 'Incident frequency =', wi_sel, ' cm-1'
    
    ! Generate grid
    if (include_Rayleigh) then
        ws_min = 0.d0
    else
        N = N-1
    endif
    x0 = ws_min - hwhm_cm1*6.0
    xf = ws_max + hwhm_cm1*6.0
    if (resol > 0) then
        npoints = nint(( xf - x0 )/dx) + 1
        xf = x0 + resol*dfloat(npoints-1)
    endif
    dx = ( xf - x0 )/dfloat(npoints - 1)
    
    ! Allocate
    allocate(xstick(N),ystick(N),spect(npoints),w(npoints))
    
    ! Then read data
    if (source == '2D') then
        i=0
        read(I_STCK,*,iostat=ioflag) line
        read(I_STCK,*,iostat=ioflag) line
        do
            read(I_STCK,'(A)',iostat=ioflag) line
            if (ioflag/=0) exit
            read(line,*) ws, wi, rr
            if (.not.include_Rayleigh .and. ws==0.d0) cycle
            if (wi == wi_sel) then
                i = i+1
                xstick(i) = ws
                ystick(i) = rr
            endif
        enddo
        if (i /= N) call alert_msg('fatal','reading scattered frequencies from 2D file')
    else
        read(I_STCK,*,iostat=ioflag) line
        read(I_STCK,*,iostat=ioflag) line
        read(I_STCK,*,iostat=ioflag) line
        if (.not.include_Rayleigh) read(I_STCK,*,iostat=ioflag) line
        do i=1,N
            read(I_STCK,*,iostat=ioflag) xstick(i), ystick(i)
        enddo
    endif
    
    ! Make convolution over the grid
    do i=1,npoints
        ! Fill x-grid
        w(i) = x0 + (i-1)*dx
        spect(i) = 0.d0
        do j=1,N
            xr = xstick(j) - w(i)
            ! Skip points beyond 20*HWHM from the center of the gaussian
            ! this increases efficiency and avoids IEEE_UNDERFLOW_FLAG IEEE_DENORMAL flags
            if (broadfun == 'GAU' .or. broadfun == 'VOI') then
                if (dabs(xr)>hwhm_cm1*20.0) cycle
                C = dlog(2.d0)/hwhm_cm1**2
                brd = dsqrt(C/pi) * dexp(-C * xr**2)
            else if (broadfun == 'LOR') then
                if (dabs(xr)>hwhm_cm1*100.0) cycle
                C = hwhm_cm1/pi !dsqrt(pi)
                brd = C/(xr**2 + hwhm_cm1**2)
            endif
            spect(i) = spect(i) + brd * ystick(j)
        enddo
    enddo
        
    ! Write spectrum to files
    write(wI_char,'(I0)') int(wi_sel)
    convfile="RR_convoluted_wI"//trim(adjustl(wI_char))//".dat"
    open(O_SPC,file=convfile,status="replace")
    dw  = (w(2)-w(1))/autoev !in Au
    do i=1,npoints
        write(O_SPC,'(F15.5,3X,G18.8E3)') w(i), spect(i)
    enddo
    close(O_SPC)
    !
    
                 
    stop
    
    contains
    
    subroutine parse_input(RRdatafile,npoints,resol,hwhm_cm1,broadfun,source,wi_req,include_Rayleigh)
    
        character(len=*),intent(inout) :: RRdatafile,broadfun,source
        real(8),intent(inout)          :: hwhm_cm1, wi_req, resol
        integer,intent(inout)          :: npoints
        logical,intent(inout)          :: include_Rayleigh

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
                    
                case ("-npoints") 
                    call getarg(i+1, arg)
                    read(arg,*) npoints
                    argument_retrieved=.true.
                    
                case ("-resol") 
                    call getarg(i+1, arg)
                    read(arg,*) resol
                    argument_retrieved=.true.
                    
                case ("-hwhm") 
                    call getarg(i+1, arg)
                    read(arg,*) hwhm_cm1
                    argument_retrieved=.true.
                    
                case ("-brd") 
                    call getarg(i+1, broadfun)
                    argument_retrieved=.true.
                    
                case ("-type") 
                    call getarg(i+1, source)
                    argument_retrieved=.true.
                    
                case ("-wI") 
                    call getarg(i+1, arg)
                    read(arg,*) wi_req
                    ! If -wI requested, we need 2D source
                    source='2D'
                    argument_retrieved=.true.
                    
                case ("-Rayleigh") 
                    include_Rayleigh=.true.
                case ("-noRayleigh") 
                    include_Rayleigh=.false.
        
                case ("-h")
                    need_help=.true.

                case default
                    print*, "Unkown command line argument: "//adjustl(arg)
                    need_help = .true.
            end select
        enddo 
    
!Print options (to stderr)
        if (need_help) then

        write(0,'(/,A)') ' convolute_RR '
        write(0,'(A)'  ) '-----------------'
        write(0,'(A)'  ) 'Generates a convolution of a RR'
        write(0,'(A)'  ) 'spectrum, from 1D or 2D output files'

        call print_version(0)

        write(0,'(/,A)') 'SYNOPSIS'
        write(0,'(A)'  ) 'convolute_RR [-f (RRdatafile) -type (1D|2D) -brd (broadfun) -hwhm (hwhm_cm1) '//&
                         '-wI (wI) -h]'

        write(0,'(/,A)')   'OPTIONS'
        write(0,'(A)'  )   'Flag       Description                  Current Value'
        write(0,'(A)'  )   ' -f        Correlation function file    '//trim(adjustl(RRdatafile))
        write(0,'(A,I8)')  ' -npoints  Nr points for output plot'   ,npoints 
        write(0,'(A,F8.3)')' -resol    resolution (in cm-1)     '   ,resol
        write(0,'(A)'  )   '           (if <0, ignore and use npoints)'
        write(0,'(A)'  )   ' -brd      Broad func (GAU|LOR)         '//trim(adjustl(broadfun))
        write(0,'(A,F8.3)')' -hwhm     HWHM (in eV)              '   ,hwhm_cm1
        write(0,'(A)'  )   ' -type     Type of output file         '//trim(adjustl(source))
        write(0,'(A)'  )   '           (1D|2D)'
        write(0,'(A,F8.3)')' -wI       Incident frequency (2D only)'   ,wi_req   
        write(0,'(A)'  )   '           (The closest value from the file is used)'
        write(0,'(A,L1)'  )' -[no]Rayleigh Include Rayleigh band  ',include_Rayleigh
        write(0,'(A)'  )   ' -h        Print help  '
        write(0,*) ''

        stop    
        endif

        return
    end subroutine parse_input

end program convolute_RR
