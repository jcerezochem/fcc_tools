program reconvolute_TI

    use alerts
    use constants
    use line_preprocess
    
    implicit none
    
    real(8),parameter :: ZERO    = 1.d-10
    
    real(8),dimension(:),allocatable :: w,spect
    real(8),dimension(:),allocatable :: ystick,xstick
    real(8) :: factor, df, dw, &
               broad, wau, dgau, dgau2, &
               C, brd, x0, xf, dx, dx0, xr
    integer :: iexp, nbroad, N
    character(len=3) :: optrans
    character(len=8) :: hwhm_char
    ! Cut convolution after ntimes hwhm
    real(8) :: ncut_gau = -1.d0, &
               ncut_lor = -1.d0
    !Counters
    integer :: i, j
    
    ! Input defaults
    character(len=3) :: broadfun='GAU'
    character(len=3) :: property='OPA'
    integer :: npoints=1001, &
               npoints_save
    real(8) :: hwhm_eV=0.01d0, &
               hwhm2_eV=0.01d0
    real(8) :: Eshift=0.d0
    real(8) :: de_eV = -9999.d0 ! for IC rate
    
    !IO
    character(len=200) :: corrfile='corr.dat',&
                          fccoutfile='fcc.out'
    character(len=200) :: stickfile='Bin_Spectrum.dat',&
                          spc1file,spc2file
    integer :: I_STCK=10, &
               I_LOG=11,  &
               O_SPC1=20, &
               O_SPC2=21
    integer :: ioflag
    character(len=200) :: line
    character :: cnull
    
    ! Aux arrays
    real(8),dimension(:),allocatable :: Vec
    
    ! IC rate interpolation
    real(8) :: R_IC
    character(len=100) :: msg
    
    
    ! Read command_line_options
    call parse_input(stickfile,npoints,hwhm_eV,hwhm2_eV,broadfun,property,Eshift,ncut_gau,ncut_lor,de_eV,&
                     fccoutfile)
    
    ! Set options upper case
    call set_word_upper_case(property)
    call set_word_upper_case(broadfun)
    
    ! Tune settings
    ! Property settings
    if (property=="OPA") then
        optrans = "ABS"
        iexp=1
        factor = facabs
    elseif (property=='EMI') then
        optrans = "EMI"
        iexp=3
        factor = facemi * autoev**2
    elseif (property=='TPA') then
        optrans = "ABS"
        iexp=2
        factor = factpa / 4.d0
    elseif (property=='ECD') then
        optrans = "ABS"
        iexp=1
        factor=facecd
    elseif (property=='MCD') then
        optrans = "ABS"
        iexp=1
        factor=facmcd
    elseif (property=='CPL') then
        optrans = "EMI"
        iexp=3
        factor=faccpl * autoev**2
    elseif (property=='IC') then
        optrans = "ABS"
        factor = 1.d15/autofs * autoev
    else
        call alert_msg('fatal','Property not supported: '//trim(adjustl(property)))
    endif
    
    ! Set cut values to defaults
    if (ncut_gau < 0.e0) ncut_gau = 20.0
    if (ncut_lor < 0.e0) then
        if (property == 'IC') then
            ncut_lor = 9999.
        else
            ncut_lor = 100.d0
        endif
    endif
    
    ! Read inout bins (sticks) 
    open(I_STCK,file=stickfile,iostat=ioflag)
    if (ioflag/=0) call alert_msg('fatal','Cannot open stickfile: '//trim(adjustl(stickfile)))
    ! First get N
    i=0
    do
        read(I_STCK,*,iostat=ioflag) brd
        if (ioflag/=0) exit
        i=i+1
    enddo
    N=i
    rewind(I_STCK)
    
    ! Tune npoints for VOI
    if (broadfun == 'VOI') then
        npoints_save = npoints
        npoints = N
    endif
    
    ! Allocate
    allocate(xstick(N),ystick(N),spect(npoints),w(npoints))
    
    ! Then read data
    do i=1,N
        ! Read sticks
        read(I_STCK,*,iostat=ioflag) xstick(i), ystick(i)
        wau=xstick(i)/autoev
        ! Transform to lineshape
        if (wau /= 0.d0) then
            ystick(i) = ystick(i)/wau**iexp
        else
            ystick(i) = 0.d0
        endif
    enddo
                 
    ! Generate grid
    if (xstick(1) < xstick(N)) then
        x0 = xstick(1) !- hwhm_eV*6.0
        xf = xstick(N) !+ hwhm_eV*6.0
    else
        x0 = xstick(N) !- hwhm_eV*6.0
        xf = xstick(1) !+ hwhm_eV*6.0
    endif
    dx = ( xf - x0 )/dfloat(npoints - 1)
    if (broadfun == 'VOI') then
        dx0 = dx
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
                if (dabs(xr)>hwhm_eV*ncut_gau) cycle
                C = dlog(2.d0)/hwhm_eV**2
                brd = dsqrt(C/pi) * dexp(-C * xr**2)
            else if (broadfun == 'LOR') then
                if (dabs(xr)>hwhm_eV*ncut_lor) cycle
                C = hwhm_eV/pi !dsqrt(pi)
                brd = C/(xr**2 + hwhm_eV**2)
            endif
            spect(i) = spect(i) + brd * ystick(j)
        enddo
    enddo
    if (broadfun == 'VOI') then
        ! We use the output data as the new sticks to be convoluted
        ! with the Lor component. This works well if the grid is
        ! fine enough. I.e., if Gau component broaded the spectrum
        ! so that the grid is not missing any spectra feature (peaks...)
        ! To ensure this condition, the first step has been carried out
        ! setting npoints = N
        ystick = spect
        xstick = w
        ! Now turn to get a spectrum with the requested npoints
        if (broadfun == 'VOI') then
            npoints = npoints_save
        endif
        deallocate(w,spect)
        allocate(w(npoints),spect(npoints))
        ! Generate grid
        if (xstick(1) < xstick(N)) then
            x0 = xstick(1) !- hwhm_eV*6.0
            xf = xstick(N) !+ hwhm_eV*6.0
        else
            x0 = xstick(N) !- hwhm_eV*6.0
            xf = xstick(1) !+ hwhm_eV*6.0
        endif
        dx = ( xf - x0 )/dfloat(npoints - 1)
        do i=1,npoints
            w(i) = x0 + (i-1)*dx
            spect(i) = 0.d0
            do j=1,N
                xr = xstick(j) - w(i)
                if (dabs(xr)>hwhm2_eV*ncut_lor) cycle
                C = hwhm2_eV/pi !dsqrt(pi)
                brd = C/(xr**2 + hwhm2_eV**2)
                ! convolute
                ! We need to multiply by dx, since we are convoluting with a stick that
                ! contains the intensity integrated over x,x+dx0
                spect(i) = spect(i) + brd * ystick(j)*dx0
            enddo
        enddo
    endif
    
    
    ! For IC, get Ead
    if (property == 'IC') then
    
        ! Reverse spectrum
        allocate(Vec(npoints))
        Vec = w(npoints:1:-1)
        w = Vec
        Vec = spect(npoints:1:-1)
        spect = Vec
        deallocate(Vec)
    
        ! Get Ead if needed
        if (de_eV<-9990.d0) then
            ! Try opening posible output files
            open(I_LOG,file=fccoutfile,status='old',iostat=ioflag)
            if (ioflag==0) then
                do
                    read(I_LOG,'(A)',iostat=ioflag) line
                    if (ioflag/=0) exit
                    if (index(line,'(for Ead =')/=0) exit
                enddo
                if (ioflag==0) then
                    read(line,*) cnull,cnull,cnull,de_eV,cnull
                endif
                write(0,'(/,X,A,F12.5,A/)') "Ead read from file ("//&
                                            trim(adjustl(fccoutfile))//&
                                            "): Ead=", de_eV, ' eV'
            endif
            ! If this does not work, set Eshift to zero
            if (ioflag/=0) then
                call alert_msg('warning','Ead not set and could not be read from file,'//&
                                        trim(adjustl(fccoutfile))//&
                                        'setting to zero')
                de_eV = 0.d0
            endif
            ! If anything worked, use Eshift = 0
            if (ioflag/=0) then
                de_eV=0.d0
                write(0,'(/,X,A,F12.5,A/)') "Cannot read from file. Setting Ead=", de_eV, ' eV'
            endif
        endif
    endif
        
    ! Write spectrum to files
    write(hwhm_char,'(F8.3)') hwhm_eV
    if (property=='IC') then
        spc2file="kic_vs_Ead_TI_hwhm"//trim(adjustl(hwhm_char))//".dat"
        open(O_SPC2,file=spc2file,status="replace")
        dw  = (w(2)-w(1))/autoev !in Au
        do i=1,npoints
            w(i) = -(w(i)+Eshift) + de_eV
            write(O_SPC2,'(F15.5,3X,G18.8E3)') w(i), spect(i)*factor*2.d0*pi
        enddo
        close(O_SPC2)
    else
        spc1file="spec_Int_TI_hwhm"//trim(adjustl(hwhm_char))//".dat"
        spc2file="spec_LS_TI_hwhm"//trim(adjustl(hwhm_char))//".dat"
        open(O_SPC1,file=spc1file,status="replace")
        open(O_SPC2,file=spc2file,status="replace")
        dw  = (w(2)-w(1))/autoev !in Au
        do i=1,npoints
            wau=(w(i)+Eshift)/autoev
            ! Issues with the format. If using G18.8, instead, values <1e-100 are written
            ! as 1-100 (missing E). Solved as:
            ! https://stackoverflow.com/questions/24004824/for-three-digit-exponents-fortran-drops-the-e-in-the-output
            write(O_SPC2,'(F15.5,3X,G18.8E3)') w(i)+Eshift, spect(i)
            write(O_SPC1,'(F15.5,3X,G18.8E3)') w(i)+Eshift, spect(i)*factor*wau**iexp
        enddo
        close(O_SPC1)
        close(O_SPC2)
    endif
    !
    
    ! For IC, get the rate
    if (property == 'IC') then
        R_IC = spect(npoints)*factor*2.d0*pi ! the rate at Ead is the last pont
        msg=' IC rate constant (s-1)'
        write(6,*) ""
        write(6,*) "========================================================"
        write(6,'(X,A,X,ES10.3)') trim(adjustl(msg)), R_IC
        write(6,'(X,A,X,F8.3,A)') "(for Ead = ", w(npoints), " eV)"
        write(6,*) "========================================================"
        write(6,*) ""
        write(6,*) ""
    endif
    
                 
    stop
    
    contains
    
    subroutine parse_input(stickfile,npoints,hwhm_eV,hwhm2_eV,broadfun,property,Eshift,ncut_gau,ncut_lor,de_eV,&
                           fccoutfile)
    
        character(len=*),intent(inout) :: stickfile,broadfun,property,fccoutfile
        real(8),intent(inout)          :: hwhm_eV, hwhm2_eV, Eshift, &
                                          ncut_gau,ncut_lor, de_eV
        integer,intent(inout)          :: npoints

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
                    call getarg(i+1, stickfile)
                    argument_retrieved=.true.
                    
                case ("-npoints") 
                    call getarg(i+1, arg)
                    read(arg,*) npoints
                    argument_retrieved=.true.
                    
                case ("-hwhm") 
                    call getarg(i+1, arg)
                    read(arg,*) hwhm_eV
                    argument_retrieved=.true.
                    
                case ("-hwhm2") 
                    call getarg(i+1, arg)
                    read(arg,*) hwhm2_eV
                    argument_retrieved=.true.
                    
                case ("-brd") 
                    call getarg(i+1, broadfun)
                    argument_retrieved=.true.
                    
                case ("-prop") 
                    call getarg(i+1, property)
                    argument_retrieved=.true.
                    
                case ("-Eshift") 
                    call getarg(i+1, arg)
                    read(arg,*) Eshift
                    argument_retrieved=.true.
                    
                case ("-Ead") 
                    call getarg(i+1, arg)
                    read(arg,*) de_eV
                    argument_retrieved=.true.
                    
                case ("-fccout") 
                    call getarg(i+1, fccoutfile)
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

        write(0,'(/,A)') ' reconvolute_TI '
        write(0,'(A)'  ) '-----------------'
        write(0,'(A)'  ) 'Generates a new "convolution" of a TI'
        write(0,'(A)'  ) 'spectrum, changing the hwhm'

        call print_version(0)

        write(0,'(/,A)') 'SYNOPSIS'
        write(0,'(A)'  ) 'reconvolute_TI -f (stickfile) [-brd (broadfun) -hwhm (hwhm_eV) '//&
                         '-prop (property) -Eshift (Eshift_eV) -cut-lor (n) -cut-gau (n) -h]'

        write(0,'(/,A)')   'OPTIONS'
        write(0,'(A)'  )   'Flag       Description                  Current Value'
        write(0,'(A)'  )   ' -f        Correlation function file    '//trim(adjustl(stickfile))
        write(0,'(A,I8)')  ' -npoints  Nr points for output plot'   ,npoints 
        write(0,'(A)'  )   ' -brd      Broad func (GAU|LOR|VOI)     '//trim(adjustl(broadfun))
        write(0,'(A,F8.3)')' -hwhm     HWHM (in eV)              '   ,hwhm_eV
        write(0,'(A,F8.3)')' -hwhm2    HWHM (in eV). For Lor term'   ,hwhm2_eV
        write(0,'(A)'  )   '           in a VOI profile '
        write(0,'(A)'  )   ' -prop     Property computed            '//trim(adjustl(property))
        write(0,'(A)'  )   '           (OPA|EMI|ECD|MCP|CPL)'
        write(0,'(A,F8.3)')' -Eshift   Energy shift (eV)         '   ,Eshift   
        if (de_eV<-9990.d0) then
        write(0,'(A)'  )   ' -Ead      Adiabatic energy (eV)        (default)'
        else
        write(0,'(A,F8.3)')' -Ead      Adiabatic energy (eV)     '   ,de_eV
        endif
        write(0,'(A)'  )   '           to get the IC rate'
        write(0,'(A)'  )   '           (default: read from file)'
        write(0,'(A)'  )   ' -fccout   Output of FCclasses job      '//trim(adjustl(fccoutfile))
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

end program reconvolute_TI
