program reconvolute_TD

    use fft
    use alerts
    use constants
    use line_preprocess
    
    implicit none
    
    real(8),parameter :: ZERO    = 1.d-10
    
    double complex,parameter :: Im = (0,1), &
                                Re = (1,0)

    real(8),dimension(:),allocatable        :: t,w,spect
    double complex,dimension(:),allocatable :: corr
    double complex :: broad, exp_w0t, exp_w0t0
    real(8) :: tt,tfin,t0,c,ci,b,bi,d,di, factor, aexp, df, dw, &
               damp, wau, dgau, dgau2, dgau3, w0_eV, w0
    integer :: iexp, nbroad, N, NN
    character(len=3) :: optrans
    character(len=8) :: hwhm_char
    !Counters
    integer :: i, ibrd
    
    ! Input defaults
    character(len=3) :: broadfun='GAU'
    character(len=3),dimension(:),allocatable ::  broadfun_list
    character(len=3) :: property='OPA'
    real(8) :: hwhm_eV=0.01d0
    real(8),dimension(:),allocatable :: hwhm_eV_list
    real(8) :: Eshift=-9999.d0
    real(8) :: de_eV = -9999.d0 ! for IC rate
    logical :: do_gibbs=.true.
    
    !IO
    character(len=200) :: corrfile='corr.dat',&
                          fccoutfile='fcc.out', &
                          spc1file,spc2file,fullcorrfile
    integer :: I_CORR=10, &
               I_LOG=11,  &
               O_SPC1=20, &
               O_SPC2=21, &
               O_CORR=22, &
               S_CORR=30
    integer :: ioflag
    character(len=200) :: line
    character :: cnull
    
    ! IC rate interpolation
    integer :: npoints
    real(8) :: freq1, freq2, RI1, RI2, R_IC
    character(len=100) :: msg
    
    ! Initialize
    w0_eV = -9999.d0
    
    ! Read command_line_options
    call parse_input(corrfile,hwhm_eV_list,broadfun_list,&
                     do_gibbs,property,Eshift,de_eV,fccoutfile,w0_eV)
                     
    ! Initialize lists
    if (.not.allocated(hwhm_eV_list)) then
        allocate(hwhm_eV_list(1))
        hwhm_eV_list(1) = hwhm_eV
    endif
    if (.not.allocated(broadfun_list)) then
        allocate(broadfun_list(1))
        broadfun_list(1)= broadfun
    endif
    
    ! Set options upper case
    call set_word_upper_case(property)
    do i=1,size(broadfun_list)
        call set_word_upper_case(broadfun_list(i))
    enddo
    
    ! Manage VOIgt
    do i=1,size(broadfun_list)
        broadfun = broadfun_list(i)
        if (broadfun == 'VOI' .and. size(broadfun_list)>1) then
            call alert_msg('fatal','VOIgt profiles do not admit additional broadenings')
        endif
    enddo
    
    broadfun = broadfun_list(1)
    if (broadfun == 'VOI') then
        deallocate(broadfun_list)
        allocate(broadfun_list(2))
        broadfun_list = (/'GAU','LOR'/)
    endif
    
    
    if (size(hwhm_eV_list) /= size(broadfun_list)) then
        print*, size(hwhm_eV_list), size(broadfun_list)
        call alert_msg('fatal','Inconsistent number of broadening functions and hwhm')
    endif
    
    !TODO: manage VOIgt
                     
    ! Get Eshift if needed
    if (Eshift==-9999.d0) then
        ! Try opening posible output files
        open(I_LOG,file=fccoutfile,status='old',iostat=ioflag)
        if (ioflag==0) then
            do
                read(I_LOG,'(A)',iostat=ioflag) line
                if (ioflag/=0) exit
                if (index(line,'Eshift=')/=0) exit
            enddo
            if (ioflag==0) then
                read(line,*) cnull,Eshift
            endif
            write(0,'(/,X,A,F12.5,A/)') "Read from file ("//&
                                        trim(adjustl(fccoutfile))//&
                                        "): Eshift=", Eshift, ' eV'
        endif
        ! If this does not work, set Eshift to zero
        if (ioflag/=0) then
            call alert_msg('warning','Eshift not set and could not be read from file '//&
                                     trim(adjustl(fccoutfile))//&
                                     ', setting to zero')
            Eshift = 0.d0
        endif
        ! If anything worked, use Eshift = 0
        if (ioflag/=0) then
            Eshift=0.d0
            write(0,'(/,X,A,F12.5,A/)') "Cannot read from file. Setting Eshift=", Eshift, ' eV'
        endif
    endif
    
    !Update w0
    if (w0_eV==-9999.d0) then
        w0_eV = Eshift
    endif
    w0_eV = w0_eV - Eshift
    
    
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
        factor = factpa / (2.d0)**iexp
    elseif (property=='TPC') then
        optrans = "ABS"
        iexp=-1
        factor = factpcd / (2.d0)**iexp
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
        factor = 1.d15/autofs
    else
        call alert_msg('fatal','Property not supported: '//trim(adjustl(property)))
    endif
    
    
    ! Read correlation function (without broadening and gibbs dumping)
    open(I_CORR,file=corrfile,iostat=ioflag)
    open(S_CORR,status='scratch')
    if (ioflag/=0) call alert_msg('fatal','Cannot open corrfile: '//trim(adjustl(corrfile)))
    ! First get N
    i=0
    do
        read(I_CORR,'(A)',iostat=ioflag) line
        if (ioflag/=0) exit
        if (INDEX(line,'#') == 1 ) cycle
        read(line,*) tfin
        write(S_CORR,*) line
        if (i==0) t0 = tfin

        i=i+1
    enddo
    N=i
    close(I_CORR)
    rewind(S_CORR)
    
    ! Manage series with initial t=0 or t>0
    if (t0 == 0.d0) then
        NN = 2*(N-1) + 1
    else
        NN = 2*N
    endif
    
    ! Allocate
    allocate(t(NN),corr(NN),spect(NN))
    
    ! Grid settings (starting point)
    w0    = w0_eV*evtown/autown
    
    ! Then read data and add new broadening and gibbs function
    do i=1,N
        ! Read correlation (without broadening/gibbs)
        read(S_CORR,*,iostat=ioflag) tt, c, ci
        
        broad = 1.d0
        do ibrd = 1,size(broadfun_list)
        
            ! Broadening settings
            dgau     = hwhm_eV_list(ibrd) *evtown/autown
            broadfun = broadfun_list(ibrd)
            
            ! Compute new broadening
            if (broadfun.eq.'GAU') then
                ! Fourier tranfomr:
                ! 1/(sigma*sqrt(2pi)) * e^(-w²/2sigma²) --> 1/sqrt(2pi) * e^(-t²*(sigma²/2))
                ! HWHM to sigma
                aexp = dgau/dsqrt(2.d0*dlog(2.d0))
                ! Exponent of time-domain Gaussian (-t²*aexp): aexp=sigma²/2
                aexp = aexp**2/2.d0
                ! Time exponent in the broadening function (Guassian)
                nbroad = 2
                ! COMPUTE
                broad = broad * dexp(-aexp*(tt/autofs)**nbroad)
            elseif (broadfun.eq.'LOR') then
                ! Fourier tranfomr:
                ! 1/pi * hwhm/(w²+hwhm²) --> 1/sqrt(2pi) * e^(-t*hwhm)
                ! Exponent of time-domain decaying exp (-t*aexp): aexp=hwhm
                aexp=dgau
                ! Time exponent in the broadening function (Exponential decay)
                nbroad=1
                ! COMPUTE
                broad = broad * dexp(-aexp*(tt/autofs)**nbroad)
            elseif (broadfun.eq.'REX') then
                ! Fourier tranform:
                ! a*e^(-aw), with w>0 --> 1/sqrt(2pi) * a/(a-it)
                ! Exponent of freq-domain decaying exp: aexp=ln2/fwhm
                aexp=dlog(2.d0)/2.d0/dgau
                ! COMPUTE
                broad = broad * aexp/(aexp - Im*tt/autofs)
            elseif (broadfun.eq.'LEX') then
                ! Fourier tranform:
                ! a*e^(aw), with w<0 --> 1/sqrt(2pi) * a/(a+it)
                ! Exponent of freq-domain decaying exp: aexp=ln2/fwhm
                aexp=dlog(2.d0)/2.d0/dgau
                ! COMPUTE
                broad = broad * aexp/(aexp + Im*tt/autofs)
            else
                call alert_msg('fatal','Unkown broadening (1) function type: '//trim(adjustl(broadfun)))
            endif
        
        enddo !list of broadenings
        
        ! Compute gibbs part
        if (do_gibbs) then
            damp = (cos(pi/2.d0*tt/tfin)**2)
        else
            damp = 1.d0
        endif
        !
        ! Change begining of the grid
        if (optrans == "EMI") then
            !-1
            exp_w0t = exp(-Im*tt/autofs*w0)
        else
            !+1
            exp_w0t = exp(+Im*tt/autofs*w0)
        endif
        !
        ! Final corr (from -T to T)
        b  = dreal(broad)
        bi = aimag(broad)
        d  = dreal(exp_w0t)
        di = aimag(exp_w0t)
        if (tt == 0.d0) then
            t(N) = tt
            corr(N) = (c + ci*Im) * damp 
        else if (t0 == 0.d0) then
            ! As t0>0, but substracting 1 for t>0 
            t(i+N-1)    = tt        
            corr(i+N-1) = (c + ci*Im) * (b + bi*Im) * (d + di*Im) * damp !* exp(-Im*Eshift/autoev*tt/autofs)
            t(N-i+1)    =-tt
            corr(N-i+1) = (c - ci*Im) * (b - bi*Im) * (d - di*Im) * damp !* exp(Im*Eshift/autoev*tt/autofs)
        else
            t(i+N)      = tt        
            corr(i+N)   = (c + ci*Im) * (b + bi*Im) * (d + di*Im) * damp !* exp(-Im*Eshift/autoev*tt/autofs)
            t(N-i+1)    =-tt
            corr(N-i+1) = (c - ci*Im) * (b - bi*Im) * (d - di*Im) * damp !* exp(Im*Eshift/autoev*tt/autofs)
        endif
    enddo
    close(S_CORR)
    
    ! Change begining of the grid
    ! Note: t0=t(NN)
    if (optrans == "EMI") then
        !-1 (but opposite sign)
        exp_w0t0 = exp(+Im*t(NN)/autofs*w0)
    else
        !+1 (but opposite sign)
        exp_w0t0 = exp(-Im*t(NN)/autofs*w0)
    endif
    ! Deactivate this option (to be checked)
    exp_w0t0 = 1.d0
    
    ! Call DFT
    write(6,'(X,A,I0,A,/)') "Calling DFT with ", NN, " points"
    call ft_calc(6,NN,t(1:NN),corr(1:NN)*exp_w0t0, &
                 spect(1:NN),optrans)
                 
    !
    ! Get the full correlation function
    write(hwhm_char,'(F8.3)') hwhm_eV_list(1)
    fullcorrfile="timecorr_hwhm"//trim(adjustl(hwhm_char))//".dat"
    open(O_CORR,file=fullcorrfile,status="replace")
    do i=1,NN
        write(O_CORR,*) t(i), dreal(corr(i)), aimag(corr(i))
    enddo
    close(O_CORR)
                 
    ! Write spectrum (LS in au)
    df = 1.d0/(t(NN) - t(1)) * ( dfloat(NN-1)/dfloat(NN) )
    deallocate(t)
    allocate(w(NN))
    do i=1,NN
        !Take into account the 1/2pi factor from the FT(delta)
        spect(i) = spect(i)/2./pi
        w(I)     = w0_eV + dfloat(i-1)*df*fstoev 
    enddo
        
    ! Write spectrum to files
    if (property=='IC') then
        spc2file="kic_vs_Ead_TD_hwhm"//trim(adjustl(hwhm_char))//".dat"
        open(O_SPC2,file=spc2file,status="replace")
        dw  = (w(2)-w(1))/autoev !in Au
        do i=1,NN
            w(i) = w(i)+Eshift
            write(O_SPC2,'(F12.5,3X,G16.8)') w(i), spect(i)*factor*2.d0*pi
        enddo
        close(O_SPC2)
    else
        spc1file="spec_Int_TD_hwhm"//trim(adjustl(hwhm_char))//".dat"
        spc2file="spec_LS_TD_hwhm"//trim(adjustl(hwhm_char))//".dat"
        open(O_SPC1,file=spc1file,status="replace")
        open(O_SPC2,file=spc2file,status="replace")
        dw  = (w(2)-w(1))/autoev !in Au
        do i=1,NN
            wau=(w(i)+Eshift)/autoev
            write(O_SPC2,'(F12.5,3X,G16.8)') w(i)+Eshift, spect(i)
            write(O_SPC1,'(F12.5,3X,G16.8)') w(i)+Eshift, spect(i)*factor*wau**iexp
        enddo
        close(O_SPC1)
        close(O_SPC2)
    endif
    
    
    ! For IC, get Ead
    if (property == 'IC') then
    
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
    
    
    ! For IC, get the rate
    if (property == 'IC') then
    
        npoints = 2*N
        
        ! Non-radiative decay
        freq1=-999.d0
        do i=2,npoints
            !Locate value for given de
            if (de_eV.ge.w(i-1).and.de_eV.le.w(i)) then
                RI1 = spect(i-1)
                RI2 = spect(i)
                freq1=w(i-1)
                freq2=w(i)
            endif
        enddo
        ! Check if we are close to the borders
        if (freq1.eq.-999.d0) then
            ! Start
            if (dabs(de_eV-w(1))<ZERO) then
                RI1 = spect(1)
                RI2 = spect(2)
                freq1=w(1)
                freq2=w(2)
            elseif (dabs(de_eV-w(npoints))<ZERO) then
                RI1 = spect(npoints-1)
                RI2 = spect(npoints)
                freq1=w(npoints-1)
                freq2=w(npoints)
            endif
        endif
                
        if (freq1.eq.-999.d0) then
            write(6,*) "Ead:", de_eV, "out of range: ",w(1)," -- ",w(npoints)
            call alert_msg('warning','IC rate could not be computed. Change settings')
        else
            ! interpolate value
            R_IC = ((de_eV-freq1)/(freq2-freq1)*(RI2-RI1)+RI1) * factor*2.d0*pi
            msg=' IC rate constant (s-1)'
            write(6,*) ""
            write(6,*) "========================================================"
            write(6,'(X,A,X,ES10.3)') trim(adjustl(msg)), R_IC
            write(6,'(X,A,X,F8.3,A)') "(for Ead = ", de_eV, " eV)"
            write(6,*) "========================================================"
            write(6,*) ""
            write(6,*) ""
        endif
    endif

    
                 
    stop
    
    contains
    
    subroutine parse_input(corrfile,hwhm_eV_list,broadfun_list,&
                           do_gibbs,property,Eshift,de_eV,fccoutfile,w0_eV)
    
        character(len=*),intent(inout) :: corrfile,property,fccoutfile
        character(len=*),dimension(:),allocatable,intent(inout) ::  broadfun_list
        logical,intent(inout)          :: do_gibbs
        real(8),intent(inout)          :: Eshift, de_eV, w0_eV
        real(8),dimension(:),allocatable,intent(inout) :: hwhm_eV_list

        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg, broadfun_string, hwhm_eV_string
        integer :: n_elem
        
        broadfun_string="GAU"
        hwhm_eV_string ="0.01"

        argument_retrieved=.false.
        do i=1,iargc()
            if (argument_retrieved) then
                argument_retrieved=.false.
                cycle
            endif
            call getarg(i, arg) 
            select case (adjustl(arg))
                case ("-f") 
                    call getarg(i+1, corrfile)
                    argument_retrieved=.true.
                    
                case ("-hwhm") 
                    call getarg(i+1, hwhm_eV_string)
                    argument_retrieved=.true.
                    ! Get vector
                    call string2vector_getsize(hwhm_eV_string,n_elem,',')
                    allocate(hwhm_eV_list(n_elem))
                    call string2rvector(hwhm_eV_string,hwhm_eV_list,n_elem,',')
                    
                case ("-damp") 
                    do_gibbs=.true.
                case ("-nodamp") 
                    do_gibbs=.false.
                    
                case ("-brd") 
                    call getarg(i+1, broadfun_string)
                    argument_retrieved=.true.
                    ! Get vector
                    call string2vector_getsize(broadfun_string,n_elem,',')
                    allocate(broadfun_list(n_elem))
                    call string2cvector(broadfun_string,broadfun_list,n_elem,',')
                    
                case ("-prop") 
                    call getarg(i+1, property)
                    argument_retrieved=.true.
                    
                case ("-Eshift") 
                    call getarg(i+1, arg)
                    read(arg,*) Eshift
                    argument_retrieved=.true.
                    
                case ("-E0") 
                    call getarg(i+1, arg)
                    read(arg,*) w0_eV
                    argument_retrieved=.true.
                    
                case ("-Ead") 
                    call getarg(i+1, arg)
                    read(arg,*) de_eV
                    argument_retrieved=.true.

                case ("-fccout") 
                    call getarg(i+1, fccoutfile)
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

        write(0,'(/,A)') ' reconvolute_TD '
        write(0,'(A)'  ) '-----------------'
        write(0,'(A)'  ) 'Generates a new "convolution" of a TD'
        write(0,'(A)'  ) 'spectrum, changing broadening funct and hwhm'

        call print_version(0)

        write(0,'(/,A)') 'SYNOPSIS'
        write(0,'(A)'  ) 'reconvolute_TD -f (corrfile) [-brd (broadfun) -hwhm (hwhm_eV) '//&
                         '-prop (property) -[no]damp -Eshift (Eshift_eV) -h]'

        write(0,'(/,A)')   'OPTIONS'
        write(0,'(A)'  )   'Flag       Description                  Current Value'
        write(0,'(A)'  )   ' -f        Correlation function file    '//trim(adjustl(corrfile))
        write(0,'(A)'  )   ' -brd      List of broad func           '//trim(adjustl(broadfun_string))
        write(0,'(A)'  )   '           (GAU|LOR|LEXP|REXP)'
        write(0,'(A)'  )   ' -hwhm     List of HWHM (in eV)         '//trim(adjustl(hwhm_eV_string))
        write(0,'(A,L1)')  ' -[no]damp Add or not damping function  ',do_gibbs
        write(0,'(A)'  )   ' -prop     Property computed            '//trim(adjustl(property))
        write(0,'(A)'  )   '           (OPA|EMI|ECD|MCP|CPL)'
        if (w0_eV<-9990.d0) then
        write(0,'(A,F8.3)')' -Eshift   Energy shift (eV)            (default)'
        else
        write(0,'(A,F8.3)')' -Eshift   Energy shift (eV)         '   ,Eshift  
        endif
        if (w0_eV<-9990.d0) then
        write(0,'(A,F8.3)')' -E0       Grid orig before shift(eV)   (default)'
        else
        write(0,'(A,F8.3)')' -E0       Grid orig before shift(eV)'   ,w0_eV  
        endif
        if (de_eV<-9990.d0) then
        write(0,'(A)'  )   ' -Ead      Adiabatic energy (eV)        (default)'
        else
        write(0,'(A,F8.3)')' -Ead      Adiabatic energy (eV)     '   ,de_eV
        endif
        write(0,'(A)'  )   '           (to compute kIC)'
        write(0,'(A)'  )   ' -fccout   Output of FCclasses job      '//trim(adjustl(fccoutfile))
        write(0,'(A)'  )   ' -h        Print help  '
        write(0,*) ''

        stop    
        endif

        return
    end subroutine parse_input

end program reconvolute_TD
