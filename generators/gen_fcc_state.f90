program gen_fcc_state

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    !Description
    ! Program to generate the state files for FCclasses from the
    ! information on the outputs of QM programs.
    !
    !==============================================================

    use constants
    use line_preprocess
    use fcc_basics
    use fcc_io
    use vibrational_analysis

    implicit none

    !Molecular info
    integer :: Nat, Nvib
    real(8)                            :: E
    real(8),dimension(:),allocatable   :: X,Y,Z,Mass
    integer,dimension(:),allocatable   :: AtNum
    real(8),dimension(:),allocatable   :: Hlt, Freq, Grad
    real(8),dimension(:,:),allocatable :: L
    real(8),dimension(:),allocatable   :: RedMass
    character(len=2)                   :: atname

    ! Harmonic model PES 
    character(len=3) :: model_PES='AH'
    logical :: force_real = .false.
    
    !Filter atoms
    character(len=800) :: filter="all"
    integer,dimension(:),allocatable :: ifilter
    integer :: Nfilt

    !Additional info to prepare the input
    real(8) :: DE, T
    logical :: is_hessian = .true., &
               is_gradient= .true., &
               is_energy= .false.
               
    ! Tunings of the output files
    logical :: write_fcc2  = .false.
    logical :: write_modes = .false.

    !Auxiliars
    character :: cnull
    integer :: error
    integer,dimension(:),allocatable   :: IAux
    real(8),dimension(:),allocatable   :: Vec
    real(8),dimension(:,:),allocatable :: Aux,Aux2
    real(8) :: aaa
    !Counters
    integer :: i,ii, j,jj, k, iat,jat

    !I/O
    character(len=100) :: strfile='default', &
                          hessfile='default',&
                          gradfile='default',&
                          enerfile='default',&
                          massfile='none',&
                          outfile='default',   &
                          outmass='default',   &
                          newoutfile='default', &
                          nmfile='compute'
    character(len=10)  :: fts='guess', &
                          fth='guess', &
                          fte='guess', &
                          ftg='guess', &
                          ftnm='guess'
    character(len=100) :: msg
    integer :: ios
    integer :: I_INP = 11, &
               I_HES = 12, &
               I_GRD = 14, &
               I_ENE = 15, &
               I_MAS = 13, &
               I_NM  = 14, &
               O_STA = 20, &
               O_FCI = 21, &
               O_HES = 22, &
               O_MAS = 23, &
               O_NEW = 24

    ! Read options
    call parse_input(strfile,fts,hessfile,fth,gradfile,ftg,enerfile,fte,massfile,nmfile,ftnm,outfile,newoutfile,outmass,&
                     model_pes,force_real,filter,write_fcc2,write_modes)
    call set_word_upper_case(model_pes)

    !Open input file
    open(I_INP,file=strfile,iostat=ios,action='read',status='old')
    if (ios /= 0) then
        print*, "Error opening "//trim(adjustl(strfile))
        stop
    endif

    !Guess the file type if not given
    if (adjustl(fts) == 'guess') then
        call split_line_back(strfile,'.',cnull,fts)
    endif

    !Read input data: natoms
    call generic_natoms_reader(I_INP,fts,Nat,error)
    if (error /= 0) then
        print*, "Error reading the number of atoms", error
        stop
    endif
    rewind(I_INP)
    
    ! Get filter (if any)
    if (adjustl(filter) == "all") then
        Nfilt=Nat
        allocate(ifilter(Nat))
        ifilter(1:Nfilt)= (/(i, i=1,Nat)/)
    else
        call selection2intlist_nel(filter,Nfilt)
        allocate(ifilter(Nfilt))
        call selection2intlist(filter,ifilter,Nfilt)
    endif
    
    if (nmfile /= 'compute' .and. Nfilt /= Nat) then
        call alert_msg('fatal','-filt cannot be used when normal modes are read, not computed')
    endif
    
    
    !Allocate input data
    allocate(X(1:3*Nat),Y(1:3*Nat),Z(1:3*Nat),Mass(1:Nat))
    allocate(AtNum(1:Nat))
    !Allocate output data
    allocate(Freq(1:3*Nfilt))
    allocate(L(1:3*Nfilt,1:3*Nfilt))
    allocate(RedMass(3*Nfilt))

    if (nmfile == 'compute') then
        ! Open new fcc3 file (will overwrite previous one)
        open(O_NEW,file=newoutfile) !,iostat=ios)
        ! Add a initial INFO section
        write(O_NEW,'(A)') 'INFO'
        write(O_NEW,'(A)') 'State file generated from file: '//trim(adjustl(strfile))//' (format: '//trim(adjustl(fts))//')'
        if (gradfile/=strfile) &
        write(O_NEW,'(A)') ' with gradient from: '//trim(adjustl(gradfile))//' (format: '//trim(adjustl(ftg))//')'
        if (hessfile/=strfile) &
        write(O_NEW,'(A)') ' with Hessian from:  '//trim(adjustl(hessfile))//' (format: '//trim(adjustl(fth))//')'
        if (massfile/='none') &
        write(O_NEW,'(A)') ' with masses from:      '//trim(adjustl(massfile))
        write(O_NEW,*) ""
    endif
    
    !Read structure
    print*, "Reading structure...", fts
    call generic_structure_reader(I_INP,fts,Nat,X,Y,Z,AtNum,Mass,error)
    if (error /= 0) then
        print*, "Error reading the geometry", error
        stop
    endif
    print'(X,A,/)', "OK"

    if (nmfile == 'compute') then
        ! Print to fcc3 file
        print*, "  and writting geom to fcc3 file..."
        write(O_NEW,'(A)') 'GEOM      UNITS=ANGS'
    !,geomunits
        if (Nfilt<int(1.e7)) then
            write(O_NEW,'(I6)') Nfilt
        else
            write(O_NEW,'(I0)') Nfilt
        endif
        write(O_NEW,'(A)') 'Geometry from '//trim(adjustl(strfile))//' in xyz format (with filter: '//trim(adjustl(filter))//')'
        do ii=1,Nfilt
            i=ifilter(ii)
            atname = atname_from_atnum(AtNum(i))
            write(O_NEW,'(A2,5X,3(F12.8,X))') atname, X(i), Y(i), Z(i)
        enddo
        write(O_NEW,*) ""
    endif

    if (adjustl(massfile) /= "none") then
        open(I_MAS,file=massfile)
        print*, "  but masses read from external massfile: "//trim(adjustl(massfile))
        do i=1,Nat
            read(I_MAS,*) Mass(i)
        enddo
        close(I_MAS)
    endif
    ! Filter masses
    allocate(Vec(Nfilt))
    do ii=1,Nfilt 
        i=ifilter(ii)
        Vec(ii) = Mass(i)
    enddo
    deallocate(Mass)
    allocate(Mass(Nfilt))
    Mass=Vec
    deallocate(Vec)
    if (write_fcc2) then
        print*, "  and writting masses to mass file..."
        open(O_MAS,file=outmass)
        do ii=1,Nfilt 
            write(O_MAS,*) Mass(ii)
        enddo
        write(O_MAS,*) ""
        close(O_MAS)
        print'(X,A,/)', "OK"
    endif

    close(I_INP)

    ! We now read the allow to read the Hessian from another file
    !Guess the file type if not given4
    if (adjustl(fte) == 'guess') then
        call split_line_back(enerfile,'.',cnull,fte)
    else if (adjustl(fte) == 'fts') then
        fte = fts
    endif
    if (adjustl(fth) == 'guess') then
        call split_line_back(hessfile,'.',cnull,fth)
    else if (adjustl(fth) == 'fts') then
        fth = fts
    endif
    if (adjustl(ftg) == 'guess') then
        call split_line_back(gradfile,'.',cnull,ftg)
    else if (adjustl(ftg) == 'fts') then
        ftg = fts
    endif
    
    if (nmfile == 'compute') then
        !Read Energy
        !Only for supported filetypes
        if (fte=='fchk' .or. fte=='engrad') then
            is_energy=.true.
            !Open gradient file
            open(I_ENE,file=enerfile,iostat=ios)
            if (ios /= 0) then
                if (adjustl(model_pes) == "VH") then
                    print*, "Error opening "//trim(adjustl(gradfile))
                    stop
                else
                    print'(/,X,A,/)', "NOTE: Energy file does not exist: skiping energy"
                    is_energy=.false.
                endif
            endif
        else
            is_energy=.false.
        endif
        if (is_energy) then
            print*, "Reading Energy...", fte
            call generic_energy_reader(I_ENE,fte,E,error)
            if (error /= 0) then
                print'(X,A,/)', "Error: Energy not present in the file."
                stop
            endif
            close(I_ENE)
        
            ! Print to fcc3 file
            print*, "  and writting energy to fcc3 file..."
            write(O_NEW,'(A)') 'ENER      UNITS=AU' !,enerunits
            write(O_NEW,'(ES16.8)') E
            write(O_NEW,*) ""
            print'(X,A,/)', "OK"
        endif

        !Read Gradient
        !Only for supported filetypes
        if (ftg=='log' .or. ftg=='fchk' .or. ftg=='molcas' .or. ftg=='cfour' .or. &
            ftg=='qchem' .or. ftg=='gmx' .or. ftg=='engrad') then ! .or. ftg=='cp2k') then
            is_gradient=.true.
            allocate (Grad(1:3*Nat))
            !Open gradient file
            open(I_GRD,file=gradfile,iostat=ios)
            if (ios /= 0) then
                if (adjustl(model_pes) == "VH") then
                    print*, "Error opening "//trim(adjustl(gradfile))
                    stop
                else
                    print'(/,X,A,/)', "NOTE: Gradient file does not exist: skiping gradient"
                    is_gradient=.false.
                endif
            endif
        else
            is_gradient=.false.
        endif
        if (is_gradient) then
            print*, "Reading Gradient...", ftg
            call generic_gradient_reader(I_GRD,ftg,Nat,Grad,error)
            if (error /= 0) then
                print'(X,A,/)', "NOTE: Gradient not present: skiping gradient"
                is_gradient=.false.
            endif
            close(I_GRD)
        endif
        
        if (is_gradient) then
        
            ! Print to fcc3 file
            print*, "  and writting grad to fcc3 file..."
            write(O_NEW,'(A)') 'GRAD      UNITS=AU' !,gradunits
            allocate(Vec(3*Nfilt))
            do ii=1,3*Nfilt,3
                iat = (ii-1)/3+1
                i   = ifilter(iat)
                Vec(ii)   = Grad(3*i-2)
                Vec(ii+1) = Grad(3*i-1)
                Vec(ii+2) = Grad(3*i)
            enddo
            deallocate(Grad)
            allocate(Grad(3*Nfilt))
            Grad=Vec
            deallocate(Vec)
            write(O_NEW,'(5ES16.8)') Grad(1:3*Nfilt)
            write(O_NEW,*) ""
            print'(X,A,/)', "OK"
        endif

        !Open hessian file
        open(I_HES,file=hessfile) !,iostat=ios)
        !Read Hessian
        allocate(Hlt(1:3*Nat*(3*Nat+1)/2))
        print*, "Reading Hessian...", fth
        call generic_Hessian_reader(I_HES,fth,Nat,Hlt,error)
        if (error /= 0) then
            print'(X,A,/)', "Hessian is not present in the file. Only valid for AS"
            is_hessian = .false.
        else

            ! Apply filter
            allocate(Aux(3*Nat,3*Nat))
            k=0
            do i=1,3*Nat
            do j=1,i
                k=k+1
                Aux(i,j) = Hlt(k)
            enddo
            enddo
            allocate(Aux2(3*Nfilt,3*Nfilt))
            do ii=1,3*Nfilt,3
            do jj=1,ii,3
                iat = (ii-1)/3+1
                i   = ifilter(iat)
                jat = (jj-1)/3+1
                j   = ifilter(jat)
                ! Fill elements
                Aux2(ii+0,jj+0) = Aux(3*i-2,3*j-2)
                Aux2(ii+1,jj+0) = Aux(3*i-1,3*j-2)
                Aux2(ii+2,jj+0) = Aux(3*i-0,3*j-2)
                
                Aux2(ii+0,jj+1) = Aux(3*i-2,3*j-1)
                Aux2(ii+1,jj+1) = Aux(3*i-1,3*j-1)
                Aux2(ii+2,jj+1) = Aux(3*i-0,3*j-1)
                
                Aux2(ii+0,jj+2) = Aux(3*i-2,3*j-0)
                Aux2(ii+1,jj+2) = Aux(3*i-1,3*j-0)
                Aux2(ii+2,jj+2) = Aux(3*i-0,3*j-0)
                
                ! Sym
                Aux2(jj+0,ii+0) = Aux2(ii+0,jj+0)
                Aux2(jj+1,ii+0) = Aux2(ii+1,jj+0)
                Aux2(jj+2,ii+0) = Aux2(ii+2,jj+0)
                
                Aux2(jj+0,ii+1) = Aux2(ii+0,jj+1)
                Aux2(jj+1,ii+1) = Aux2(ii+1,jj+1)
                Aux2(jj+2,ii+1) = Aux2(ii+2,jj+1)
                
                Aux2(jj+0,ii+2) = Aux2(ii+0,jj+2)
                Aux2(jj+1,ii+2) = Aux2(ii+1,jj+2)
                Aux2(jj+2,ii+2) = Aux2(ii+2,jj+2)
            enddo
            enddo
            allocate(Vec(3*Nfilt*(3*Nfilt+1)/2))
            k=0
            do i=1,3*Nfilt
            do j=1,i
                k=k+1
                Vec(k) = Aux2(i,j)
            enddo
            enddo
            deallocate(Aux,Aux2)
            deallocate(Hlt)
            allocate(Hlt(3*Nfilt*(3*Nfilt+1)/2))
            Hlt=Vec
            deallocate(Vec)
            
            
            ! Print to fcc3 file
            print*, "  and writting lower triangular hess to fcc3 file..."
            write(O_NEW,'(A)') 'HESS      UNITS=AU' !,hessunits
            write(O_NEW,'(5ES16.8)') Hlt(1:3*Nfilt*(3*Nfilt+1)/2)
            write(O_NEW,*) ""
            print'(X,A,/)', "OK"
        endif

        Nat = Nfilt

        !Close hessian file
        close(I_HES)
    else
        is_hessian = .false.
    endif

    if (is_hessian .and. (write_fcc2 .or. write_modes)) then
        !Perform vibrational analysis
        print*, "Diagonalizing Hessian..."
        call diag_int(Nat,X,Y,Z,Mass,Hlt,Nvib,L,Freq,error)
        if (error /= 0) then
            print*, "Error in the diagonalization"
            stop
        endif
        !Transform L to Cartesian/Normalized
        call Lmwc_to_Lcart(Nat,Nvib,Mass,L,L,error)
        call Lcart_to_LcartNrm(Nat,Nvib,L,L,RedMass,error)
        !Transform Force Constants to Freq
        do i=1,Nvib
            Freq(i) = dsign(dsqrt(dabs(Freq(i))*HARTtoJ/BOHRtoM**2/AMUtoKG)/2.d0/pi/cvel/1.d2,Freq(i))
            if (Freq(i)<0.d0 .and. force_real) then
                print'(3X,A,X,F12.3)', "Warning: an imaginary frequency turned real:", Freq(i)
                Freq(i)=dabs(Freq(i))
            endif
        enddo
        if (error /= 0) then
            print*, "Error in conversion to Cartesian"
            stop
        else
            print'(X,A,/)', "OK"
        endif
    else if (nmfile /= 'compute') then
        print*, "Reading normal modes...", ftnm
        write_fcc2 = .true.
        open(I_NM,file=nmfile,iostat=ios,action='read',status='old')
        if (ios /= 0) then
            print*, "Error opening "//trim(adjustl(nmfile))
            stop
        endif
        
        !Guess the file type if not given
        if (adjustl(ftnm) == 'guess') then
            call split_line_back(nmfile,'.',cnull,ftnm)
        endif
        
        ! Read nm
        call generic_nm_reader(I_NM,ftnm,Nat,Nvib,Freq,L,error)
        if (error < 0) then
            print*, "Error reading normal modes"
            stop
        else
            print'(X,A,/)', "OK"
        endif
        close(I_NM)
    endif

    if (write_fcc2) then
    !WRITE STATE FILE
    print*, "Writting state file..."
        open(O_STA,file=outfile,status="replace",iostat=ios)
        if (ios /= 0) then
            print*, "Cannot open "//trim(adjustl(outfile))//" to write"
            stop
        endif

        if (adjustl(model_pes) == "AH") then
            do j=1,Nat
                write(O_STA,'(E17.8)',iostat=ios) X(j),Y(j),Z(j)
            enddo
            if (is_hessian .or. nmfile /= 'compute') then
                do j=1,3*Nat
                do i=1,Nvib
                    write(O_STA,'(E17.8)',iostat=ios) L(j,i)
                enddo
                enddo
                do j=1,Nvib
                    write(O_STA,'(F10.4)',iostat=ios) Freq(j)
                enddo
            endif
        elseif (adjustl(model_pes) == "VH") then
            
            call write_fchk(O_STA,"Cartesian Gradient",'R',3*Nat,Grad,IAux,error)
            call write_fchk(O_STA,"Cartesian Force Constants",'R',3*Nat*(3*Nat+1)/2,Hlt,IAux,error)
        else
            print*, "Error: Unkown model PES: "//adjustl(model_pes)
            stop
        endif
        close(O_STA)
        if (ios /= 0) then
            print*, "Error writting state file"
            stop
        else
            print'(X,A,/)', "OK"
        endif
    endif
    
    if (write_modes) then
        ! Print discription of vibrations
        open(O_STA,file='description.dat')
        do i=1,Nvib
            write(O_STA,'(A5,I4,A1,F8.2,A)') 'Mode ', i, ':', Freq(i), ' cm-1'
            write(O_STA,'(     10X,F8.2,A)') RedMass(i), ' amu'
            do jj = 1,Nat
                j=3*jj-2
                ! Replaced
                ! call atominfo_from_atmass(Mass(i),k,atname)
                k = atnum_from_atmass(Mass(jj))
                if (k == 0) then
                    atname = 'X'
                    write(msg,'(A,F10.4,A)') 'Element cannot be set from mass', Mass(i), '. AtNum set to 0, name set to X.'
                    call alert_msg('note',msg)
                else
                    atname = atname_from_atnum(k)
                endif 
                ! Factor to transform into MWC (original L elements)
                aaa = dsqrt(Mass(jj))/RedMass(i)
                ! Keep Gaussian output
                aaa = 1.d0
                write(O_STA,'(A,3X,3F8.4)') atname, L(j,i)  *aaa,&
                                                    L(j+1,i)*aaa,&
                                                    L(j+2,i)*aaa
            enddo
            write(O_STA,*)  ''
        enddo
    endif

    if (write_fcc2) then
        !We profit to generate a first input template
        print*, "Writting input template: 'fcc2_template.inp'..."
        open(O_FCI,file="fcc2_template.inp")
        DE = -1.d0
        T  = -1.d0
        call prepare_fccinput(O_FCI,Nat,Nvib,Mass,DE,T,error)
        close(O_FCI)
        if (error /= 0) then
            print*, "Error writting input template"
            stop
        else
            print'(X,A,/)', "OK"
        endif
    endif


    print*, "** Successful end **"

    !Deallocate
    deallocate(X,Y,Z,Mass,Freq,L)
    if (nmfile == 'compute') then
        deallocate(Hlt)
        if (adjustl(model_pes) == "VH") deallocate(Grad)
    endif

    stop

    contains

    subroutine parse_input(strfile,fts,hessfile,fth,gradfile,ftg,enerfile,fte,massfile,nmfile,ftnm,outfile,newoutfile,outmass,&
                           model_pes,force_real,filter,write_fcc2,write_modes)

        character(len=*),intent(inout) :: strfile,fts,hessfile,fth,gradfile,ftg,enerfile,fte,massfile,&
                                          outfile,outmass,model_pes,newoutfile,filter,nmfile,ftnm
        logical,intent(inout)          :: force_real,write_fcc2,write_modes

        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        character(len=20)  :: prfx=""

        argument_retrieved=.false.
        do i=1,iargc()
            if (argument_retrieved) then
                argument_retrieved=.false.
                cycle
            endif
            call getarg(i, arg) 
            select case (adjustl(arg))
                case ("-i") 
                    call getarg(i+1, strfile)
                    argument_retrieved=.true.
                case ("-ft") !for backward compatibility
                    call getarg(i+1, fts)
                    argument_retrieved=.true.
                case ("-fts") 
                    call getarg(i+1, fts)
                    argument_retrieved=.true.

                case ("-ih") 
                    call getarg(i+1, hessfile)
                    argument_retrieved=.true.
                case ("-fth") 
                    call getarg(i+1, fth)
                    argument_retrieved=.true.

                case ("-inm") 
                    call getarg(i+1, nmfile)
                    argument_retrieved=.true.
                case ("-ftnm") 
                    call getarg(i+1, ftnm)
                    argument_retrieved=.true.
                    
                case ("-ig") 
                    call getarg(i+1, gradfile)
                    argument_retrieved=.true.
                case ("-ftg") 
                    call getarg(i+1, ftg)
                    argument_retrieved=.true.
                    
                case ("-ie") 
                    call getarg(i+1, enerfile)
                    argument_retrieved=.true.
                case ("-fte") 
                    call getarg(i+1, fte)
                    argument_retrieved=.true.

                case ("-im") 
                    call getarg(i+1, massfile)
                    argument_retrieved=.true.

                case ("-ofcc2") 
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.
                    
                case ("-o") 
                    call getarg(i+1, newoutfile)
                    argument_retrieved=.true.

                case ("-om") 
                    call getarg(i+1, outmass)
                    argument_retrieved=.true.

                case ("-model") 
                    call getarg(i+1, model_pes)
                    argument_retrieved=.true.
                    
                case ("-filt") 
                    call getarg(i+1, filter)
                    argument_retrieved=.true.

                case ("-force-real")
                    force_real=.true.
                case ("-noforce-real")
                    force_real=.false.
                    
                case ("-write-modes")
                    write_modes=.true.
                case ("-nowrite-modes")
                    write_modes=.false.
                    
                case ("-write-fcc2")
                    write_fcc2=.true.
                case ("-nowrite-fcc2")
                    write_fcc2=.false.
        
                case ("-h")
                    need_help=.true.

                case default
                    print*, "Unkown command line argument: "//adjustl(arg)
                    need_help = .true.
            end select
        enddo 

        ! Post-processing
        !----------------------------
        if (adjustl(strfile) == 'default') then
            if (nmfile == 'compute') then
                strfile = 'output.log'
            else
                strfile = nmfile
                fts = ftnm
            endif
        endif
        if (adjustl(hessfile) == 'default') then
            ! The try to read the hessian from strfile
            hessfile = strfile
            fth = "fts"
        endif
        if (adjustl(gradfile) == 'default') then
            ! The try to read the gradient from strfile
            gradfile = strfile
            ftg = "fts"
        endif
        if (adjustl(enerfile) == 'default') then
            ! The try to read the energy from strfile
            enerfile = strfile
            fte = "fts"
        endif
        if (adjustl(outfile) == 'default') then
            ! Get relative path (split_line_back returns the line in the first part is the splitter is not present)
            if (index(strfile,"./") /= 0) then
                call split_line_back(strfile,"./",prfx,outfile)
                prfx=trim(adjustl(prfx))//"./"
            else
                outfile=strfile
                prfx=""
            endif
            call split_line_back(outfile,".",outfile,arg)
            if (adjustl(fts) /= 'guess') arg=fts
!             outfile = trim(adjustl(prfx))//&
!                       "state_"//trim(adjustl(outfile))//'_'//trim(adjustl(arg))
! better not using the relative path (in prfx):
            outfile = "state_"//trim(adjustl(outfile))//'_'//trim(adjustl(arg))
        endif
        if (adjustl(newoutfile) == 'default') then
            ! Get relative path (split_line_back returns the line in the first part is the splitter is not present)
            if (index(strfile,"./") /= 0) then
                call split_line_back(strfile,"./",prfx,newoutfile)
                prfx=trim(adjustl(prfx))//"./"
            else
                newoutfile=strfile
                prfx=""
            endif
            call split_line_back(newoutfile,".",newoutfile,arg)
            if (adjustl(fts) /= 'guess') arg=fts
!             newoutfile = trim(adjustl(prfx))//&
            newoutfile = trim(adjustl(newoutfile))//'.fcc'
                         
        endif
        if (adjustl(outmass) == 'default') then
            call split_line_back(outfile,"state_",arg,outmass)
            if (adjustl(fts) /= 'guess') arg=fts
!             outmass = trim(adjustl(prfx))//&
            outmass = "mass_"//trim(adjustl(outmass))
        endif


       !Print options (to stderr)
        if (need_help) then

        write(0,'(/,A)') ' gen_fcc_state '
        write(0,'(A)'  ) '-----------------'
        write(0,'(A)'  ) 'Generates state_files for FCclasses from the output'
        write(0,'(A)'  ) 'files obtained with different QM or MM codes, reading the '
        write(0,'(A)'  ) 'coordinates and the Hessian from them.'
        write(0,'(A)'  ) 'Additionally, an input template is also generated'
        write(0,'(A)'  ) ''
        write(0,'(A)'  ) 'Alternatively, normal modes can be read from the QM output. In this'
        write(0,'(A)'  ) 'case, Hessian and Grad are ignored'
        write(0,'(A)'  ) ''
        write(0,'(A)'  ) 'The input structure (-i) is read from any of the supported QM output'
        write(0,'(A)'  ) 'files or through the g96 (GROMOS) structure file. The file type is'
        write(0,'(A)'  ) 'guessed based on the extension, but if this fails, it can be specified'
        write(0,'(A)'  ) 'by -fts flag. If the structure file also contains the hessian, it will'
        write(0,'(A)'  ) 'be read from that file by default. Otherwise, another file containing the'
        write(0,'(A)'  ) 'Hessian should be indicated (-ih) with filetype detected based on extension'
        write(0,'(A)'  ) 'of specified (-fth). By default, only FCC3 states files are generated. '
        write(0,'(A)'  ) 'Optionally, FCC2 files (including mass file) can be written. FCC2 output would'
        write(0,'(A)'  ) 'be named automatiacally as state_{input}_{ext}, or can be indicated (-ofcc2).'
        write(0,'(A)'  ) 'Masses to carry out the analsis are taken from an internal database but they'
        write(0,'(A)'  ) 'can be specified differently in a file (with Nat lines, each with the mass of'
        write(0,'(A)'  ) 'the atom in amu) through -im flag. The Hessian and masses used in the calculations'
        write(0,'(A)'  ) 'are writen to files, which can be specified on input (-oh -om).'

        call print_version(0)

        write(0,'(/,A)') 'SYNOPSIS'
        write(0,'(A)'  ) 'gen_fcc_state -i input_file [-fts filetype-str] [-ih hess_inp_file] [-fth filetype-hess] '//&
                         '[-o output_file] [-om mass_file] [-filter X-Y,Z] [-write-fcc2] [write-modes] [-h]'

        write(0,'(/,A)') 'OPTIONS'
        write(0,'(A)'  ) 'Flag    Description       Current Value'
        write(0,'(A)'  ) ' -i     structure_file    '//trim(adjustl(strfile))
        write(0,'(A)'  ) ' -fts   filetype(str)     '//trim(adjustl(fts))
        write(0,'(A)'  ) ' -ih    hess_input_file   '//trim(adjustl(hessfile))
        write(0,'(A)'  ) ' -fth   filetype(hess)    '//trim(adjustl(fth))
        write(0,'(A)'  ) ' -ig    grad_input_file   '//trim(adjustl(gradfile))
        write(0,'(A)'  ) ' -ftg   filetype(grad)    '//trim(adjustl(ftg))
        write(0,'(A)'  ) ' -ie    ener_input_file   '//trim(adjustl(enerfile))
        write(0,'(A)'  ) ' -fte   filetype(ener)    '//trim(adjustl(fte))
        write(0,'(A)'  ) ' -inm   structure_file    '//trim(adjustl(nmfile))
        write(0,'(A)'  ) ' -ftnm  filetype(nm)      '//trim(adjustl(ftnm))
        write(0,'(A)'  ) ' -im    mass_file         '//trim(adjustl(massfile))
        write(0,'(A)'  ) ' -o     output_file(fcc3) '//trim(adjustl(newoutfile))
        write(0,'(A)'  ) ' -ofcc2 output_file(fcc2) '//trim(adjustl(outfile))
        write(0,'(A)'  ) ' -om    mass_file         '//trim(adjustl(outmass))
!         write(0,'(A)'  ) ' -model model_pes[AH|VH]  '//trim(adjustl(model_pes))
        write(0,'(A)'  ) ' -filt  Filter atoms      '//trim(adjustl(filter))
        write(0,'(A,A)') ' -write-fcc2  Write fcc2 '
        write(0,'(A,A)') '              state files',  write_fcc2
        write(0,'(A,A)') ' -write-modes Write normal'
        write(0,'(A,A)') '              modes descr',  write_modes
        write(0,'(A,A)') ' -force-real  turn real  ',  force_real
        write(0,'(A)'  ) '        all imag freqs' 
        write(0,'(A)'  ) ' -h     print help  '
        call supported_filetype_list('freq')

        write(0,'(/,A)') 'EXAMPLES'
        write(0,'(A)'  ) ' gen_fcc_state -i output.fchk'
        write(0,'(A)'  ) '  In general simple instructions as above are required.'
        write(0,'(A)'  ) '  That one will generate: state_output_fchk and a template'
        write(0,'(A)'  ) '  input: fcc.inp from a valid gaussian fchk file.'
        write(0,'(A)'  ) ''
        write(0,'(A)'  ) ' gen_fcc_state -i output-psi4.out -fts psi4'
        write(0,'(A)'  ) '  Some output files may need to specify the filetype.'
        write(0,'(A)'  ) '  Note also that in the case of psi4, the print level'
        write(0,'(A)'  ) '  must be set to 3 in order to get the Hessian printed'
        write(0,'(A)'  ) '  in the output' 

        stop    
        endif

        return
    end subroutine parse_input

end program gen_fcc_state

