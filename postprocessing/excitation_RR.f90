program excitation_RR

    use alerts
    use constants
    use line_preprocess
    
    implicit none
    
    real(8),dimension(:),allocatable :: w,spect
    integer :: N
    character(len=3) :: optrans
    character(len=8) :: wS_char
    !Counters
    integer :: i, j
    
    ! Input defaults
    real(8) :: ws_req=0.d0
    
    ! Additional data
    real(8) :: wi, ws, rr, ws_,&
               ws_max, ws_min, &
               ws_sel, dev_ws
    logical :: read_wi
    
    
    !IO
    character(len=200) :: RRdatafile='(default)',&
                          stickfile,excfile
    integer :: I_STCK=10, &
               O_SPC =20, &
               O_STCK=21
    integer :: ioflag
    character(len=200) :: line
    character :: cnull
    
    
    ! Read command_line_options
    call parse_input(RRdatafile,ws_req)
    
    ! Set dagaults
    if (RRdatafile == '(default)') then
        RRdatafile = 'RR_Spectrum_2D.dat'
    endif
    write(6,*) 'Reading data from: '//trim(adjustl(RRdatafile))
    
    ! Read inout bins (sticks) 
    open(I_STCK,file=RRdatafile,iostat=ioflag)
    ! First get number of wI points inthe grid
    ! Skip header line
    read(I_STCK,*,iostat=ioflag) line
    read(I_STCK,*,iostat=ioflag) line
    ws_ = -1.0
    dev_ws = 9999.
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
    ! Allocate and read the data
    allocate(spect(N),w(N))
    rewind(I_STCK)
    ! Skip header line
    read(I_STCK,*,iostat=ioflag) line
    read(I_STCK,*,iostat=ioflag) line
    ws_ = -1.0
    dev_ws = 9999.
    do 
        read(I_STCK,'(A)',iostat=ioflag) line
        if (ioflag/=0) exit
        read(line,*) ws, wi, rr
!             print*, rr
        if (ws /= ws_) then
            ws_ = ws
            read_wi = .false.
            i=0
            ! Check if we are the closest till now
            
            if (dabs(ws-ws_req) < dev_ws) then
                ws_sel = ws
                dev_ws = dabs(ws-ws_req)
                read_wi = .true.
            endif
        endif
        
        if (read_wi) then
            i = i+1
            w(i) = wI
            spect(i) = rr
        endif
    enddo
    write(6,'(X,A,F10.2,A,/)') 'Scattered frequency =', ws_sel, ' cm-1'
    
    ! Write spectrum to files
    write(wS_char,'(I0)') int(ws_sel)
    excfile="RR_excitation_wS"//trim(adjustl(wS_char))//".dat"
    open(O_SPC,file=excfile,status="replace")
    do i=1,N
        write(O_SPC,'(F15.5,3X,G18.8E3)') w(i), spect(i)
    enddo
    close(O_SPC)
    !
                 
    stop
    
    contains
    
    subroutine parse_input(RRdatafile,ws_req)
    
        character(len=*),intent(inout) :: RRdatafile
        real(8),intent(inout)          :: ws_req

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
                    
                case ("-wS") 
                    call getarg(i+1, arg)
                    read(arg,*) ws_req
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

        write(0,'(/,A)') ' excitation_RR '
        write(0,'(A)'  ) '-----------------'
        write(0,'(A)'  ) 'Extract excitation profile from'
        write(0,'(A)'  ) '2D output files'

        call print_version(0)

        write(0,'(/,A)') 'SYNOPSIS'
        write(0,'(A)'  ) 'excitation_RR [-f (RRdatafile) -wS (wS) -h]'

        write(0,'(/,A)')   'OPTIONS'
        write(0,'(A)'  )   'Flag       Description               Current Value'
        write(0,'(A)'  )   ' -f        2D RR spectrum file      '//trim(adjustl(RRdatafile))
        write(0,'(A,F8.3)')' -wS       Scattered frequency(cm-1)'   ,ws_req
        write(0,'(A)'  )   '           (The closest value from the file is used)'
        write(0,'(A)'  )   ' -h        Print help  '
        write(0,*) ''

        stop    
        endif

        return
    end subroutine parse_input

end program excitation_RR
