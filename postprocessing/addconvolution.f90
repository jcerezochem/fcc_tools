program addconvolution

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
    real(8) :: wi, specti
    
    real(8) :: dw
    
    
    !IO
    character(len=200) :: input_spc='spec_LS.dat',&
                          stickfile,output_spc
    integer :: I_STCK=10, &
               O_SPC =20, &
               O_STCK=21
    integer :: ioflag
    character(len=200) :: line
    character :: cnull
    
    
    ! Read command_line_options
    call parse_input(input_spc,hwhm_eV,broadfun,ncut_gau,ncut_lor)
    
    ! Set options upper case
    call set_word_upper_case(broadfun)
    
    write(6,*) 'Reading data from: '//trim(adjustl(input_spc))
    
    if (ncut_gau < 0) ncut_gau = 20.0
    if (ncut_lor < 0) ncut_lor = 100.0
    
    ! Write spectrum to files
    write(hwhm_char,'(F8.3)') hwhm_eV
    output_spc="spec_LS_add-HWHM_"//trim(adjustl(hwhm_char))//".dat"
    open(O_SPC,file=output_spc,status="replace")
    write(6,*) 'Writing data to  : '//trim(adjustl(output_spc))
    
    ! Read inout bins (sticks) 
    open(I_STCK,file=input_spc,iostat=ioflag)
    ! Read file as E(eV) LS
    N=1
    do 
        read(I_STCK,'(A)',iostat=ioflag) line
        if (ioflag/=0) exit
        read(line,*) wi, specti
        N = N+1
    enddo
    
    ! Allocate arrays
    allocate(spect(N),w(N))
    
    ! Now, read again all file and convolute each excitation profile
    rewind(I_STCK)
    i=0
    do 
        read(I_STCK,'(A)',iostat=ioflag) line
        if (ioflag/=0) exit

        i = i+1
        if (i>N) call alert_msg('fatal','Inconsistency in input file')

        read(line,*) w(i), spect(i)

    enddo
    close(I_STCK)
        
    call add_convolution(w,spect,broadfun,hwhm_eV,ncut_gau,ncut_lor)
    ! Return to RR intensity
    do i=1,N
        write(O_SPC,'(F13.4,X,F15.4,X,E19.5)') w(i), spect(i)
    enddo
    
    write(6,'(X,A,/)') 'Done'

    close(O_SPC)

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
        real(8) :: dx, xr, brd, C
        integer :: N
        
        N = size(w)
        allocate(y(N))
        
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
                    if (dabs(xr)>hwhm_eV*ncut_gau) cycle
                    C = dlog(2.d0)/hwhm_eV**2
                    brd = dsqrt(C/pi) * dexp(-C * xr**2)
                else if (broadfun == 'LOR') then
                    if (dabs(xr)>hwhm_eV*ncut_lor) cycle
                    C = hwhm_eV/pi !dsqrt(pi)
                    brd = C/(xr**2 + hwhm_eV**2)
                endif
                y(i) = y(i) + brd * spect(j) * dx
            enddo
        enddo
        
        ! Update output values
        spect = y
        deallocate(y)
    
        return
    
    end subroutine add_convolution
    
    
    subroutine parse_input(input_spc,hwhm_eV,broadfun,ncut_gau,ncut_lor)
    
        character(len=*),intent(inout) :: input_spc, broadfun
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
                    call getarg(i+1, input_spc)
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

        write(0,'(/,A)') ' addconvolution '
        write(0,'(A)'  ) '-----------------'
        write(0,'(A)'  ) 'Add a convolution to a full LS spectrum'

        call print_version(0)

        write(0,'(/,A)') 'SYNOPSIS'
        write(0,'(A)'  ) 'addconvolution [-f (input_spc) -brd (broadtype) -hwhm (hwhm) -h]'

        write(0,'(/,A)')   'OPTIONS'
        write(0,'(A)'  )   'Flag       Description                  Current Value'
        write(0,'(A)'  )   ' -f        Input spectrum               '//trim(adjustl(input_spc))
        write(0,'(A)'  )   ' -brd      Broad func (GAU|LOR)         '//trim(adjustl(broadfun))
        write(0,'(A,F8.3)')' -hwhm     HWHM (in eV)              '   ,hwhm_eV
        write(0,'(A,F8.3)')' -cut-lor  Multiple of hwhm when to  '   ,ncut_lor  
        write(0,'(A)'  )   '           cut the Lor convolution'
        write(0,'(A,F8.3)')' -cut-gau  Multiple of hwhm when to  '   ,ncut_gau
        write(0,'(A)'  )   '           cut the Gau convolution'
        write(0,'(A)'  )   ' -h        Print help  '
        write(0,*) ''

        stop    
        endif

        return
    end subroutine parse_input

end program addconvolution
