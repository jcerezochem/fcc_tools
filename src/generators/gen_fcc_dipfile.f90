program gen_fcc_dipfile


    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    !DESCRIPTION
    !-----------
    !Utility to generate FCclasses electric and magnetic dipole files from a 
    !Gaussian fchk file. For a frequency job, the derivatives of the transition
    !dipoles are also stored in the fchk, and can also be read. But note that
    !if symmetry was used in the frequency calculations, the derivatives are 
    !not computed properly, and should not be used.
    !
    !===========================================================================

    use constants
    use line_preprocess
    use fcc_basics
    use fcc_io
    
    implicit none

    !Relevant data
    integer :: Nat
    character(len=250) :: jobname
    character(len=20)  :: jobtype,method,basis
    real(8)                :: dx = -1.d0 
    real(8),dimension(1:3) :: Dip
    real(8),dimension(:),allocatable :: DipD
    real(8),dimension(:),allocatable :: nac
    logical :: derivatives=.true.
    logical :: do_nac=.false.
    !Default states are taken if =-1
    integer :: Si=-1, &
               Sf=-1
               
    !Filter atoms
    character(len=20) :: filter="all"
    integer,dimension(:),allocatable :: ifilter
    integer :: Nfilt

    !Auxiliars
    integer :: N
    character :: dtype, cnull
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: I
    integer :: error
    !Selection stuff
    integer :: Nsel 
    integer,dimension(1:1000) :: nm_select
    !Counters
    integer :: j,k, jj, ii, iat

    !I/O
    character(len=100) :: inpfile='input.log', &
                          out_eldip='default', &
                          out_magdip='default',&
                          out_nac='default'
    character(len=10)  :: ft='guess'
    integer :: ios
    integer :: I_INP=10,&
               O_DIP =20

    ! Read options
    call parse_input(inpfile,ft,out_eldip,out_magdip,out_nac,derivatives,do_nac,Si,Sf,filter)

    !Open input file
    open(I_INP,file=inpfile,status="old",iostat=ios)
    if (ios /= 0) then
        print*, "Error opening "//trim(adjustl(inpfile))
        stop
    endif

    !Guess the file type if not given
    if (adjustl(ft) == 'guess') then
        call split_line_back(inpfile,'.',cnull,ft)
    endif

    !Read input data: natoms
    call generic_natoms_reader(I_INP,ft,Nat,error)
    
    ! Get filter (if any)
    if (adjustl(filter) == "all") then
        Nfilt=Nat
        allocate(ifilter(Nat))
        ifilter(1:Nfilt)= (/(j, j=1,Nat)/)
    else
        call selection2intlist_nel(filter,Nfilt)
        allocate(ifilter(Nfilt))
        call selection2intlist(filter,ifilter,Nfilt)
    endif

    !Get eldip
    print*, "Reading transition electric dipole moment..."
    if (derivatives) then
        !Allocate output array
        allocate(DipD(1:3*3*Nat))
    endif
    call generic_dip_reader(I_INP,ft,Si,Sf,derivatives,"eldip",dx,Dip,DipD,error)
    if (error /= 0) then
        print*, "Error getting eldip. Error code:", error
        stop
    else
        print'(X,A,/)', "OK"
    endif

    !WRITE ELDIP FILE
    print*, "Writting transition electric dipole file..."
    open(O_DIP,file=out_eldip,status="replace")

    !One should replace one of these by that at the other state geom
    write(O_DIP,'(3(X,E18.9))') Dip(1:3)
    write(O_DIP,'(3(X,E18.9))') Dip(1:3)

    if (derivatives) then
        do ii=1,Nfilt
            j   = (ifilter(ii)-1)*3+1
            jj=j*3-2
            write(O_DIP,'(3(X,E18.9))',iostat=ios) DipD(jj:jj+2)
            j=j+1
            jj=j*3-2
            write(O_DIP,'(3(X,E18.9))',iostat=ios) DipD(jj:jj+2)
            j=j+1
            jj=j*3-2
            write(O_DIP,'(3(X,E18.9))',iostat=ios) DipD(jj:jj+2)
        enddo
    endif
    close(O_DIP)
    if (ios /= 0) then
        print*, "Error writting eldip file"
        stop
    else
        print'(X,A,/)', "OK"
    endif
    
    !============================================================
    !Rewind input file
    rewind(I_INP)
    !============================================================

    !Get magdip
    print*, "Reading transition magnetic dipole moment..."
    call generic_dip_reader(I_INP,ft,Si,Sf,derivatives,"magdip",dx,Dip,DipD,error)
    if (error /= 0) then
        print*, "Error getting magdip. Error code:", error
        stop
    else
        print'(X,A,/)', "OK"
    endif

    !WRITE MAGDIP FILE
    print*, "Writting transition magnetic dipole file..."
    open(O_DIP,file=out_magdip,status="replace")

    !One should replace one of these by that at the other state geom
    write(O_DIP,'(3(X,E18.9))') Dip(1:3)
    write(O_DIP,'(3(X,E18.9))') Dip(1:3)

    if (derivatives) then
        do ii=1,Nfilt
            j   = (ifilter(ii)-1)*3+1
            jj=j*3-2
            write(O_DIP,'(3(X,E18.9))',iostat=ios) DipD(jj:jj+2)
            j=j+1
            jj=j*3-2
            write(O_DIP,'(3(X,E18.9))',iostat=ios) DipD(jj:jj+2)
            j=j+1
            jj=j*3-2
            write(O_DIP,'(3(X,E18.9))',iostat=ios) DipD(jj:jj+2)
        enddo
    endif
    close(O_DIP)
    if (ios /= 0) then
        print*, "Error writting eldip file"
        stop
    else
        print'(X,A,/)', "OK"
    endif

    if (derivatives) deallocate(DipD)
    
    if (do_nac) then
    
        allocate(nac(1:3*Nat))
    
        !============================================================
        !Rewind input file
        rewind(I_INP)
        !============================================================

        !Get nac
        print*, "Reading non-adiabatic coupling elements..."
        call generic_nac_reader(I_INP,ft,Si,Sf,nac,error)
        if (error /= 0) then
            print*, "Error getting NAC. Error code:", error
            stop
        else
            print'(X,A,/)', "OK"
        endif

        !WRITE NAC FILE
        print*, "Writting non-adiabatic coupling file..."
        open(O_DIP,file=out_nac,status="replace")

        do ii=1,Nfilt
            j   = (ifilter(ii)-1)*3+1
            write(O_DIP,'(3(X,E18.9))',iostat=ios) nac(j:j+2)
        enddo
        close(O_DIP)
        deallocate(nac)
        if (ios /= 0) then
            print*, "Error writting nac file"
            stop
        else
            print'(X,A,/)', "OK"
        endif

    endif
        
    print*, "** Successful end **"

    stop

    contains

    subroutine parse_input(inpfile,ft,out_eldip,out_magdip,out_nac,derivatives,do_nac,Si,Sf,filter)

        character(len=*),intent(inout) :: inpfile,ft,out_eldip,out_magdip,out_nac,filter
        logical,intent(inout)          :: derivatives,do_nac
        integer,intent(inout)          :: Si, Sf

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
                    call getarg(i+1, inpfile)
                    argument_retrieved=.true.
                case ("-ft") 
                    call getarg(i+1, ft)
                    argument_retrieved=.true.

                case ("-oe") 
                    call getarg(i+1, out_eldip)
                    argument_retrieved=.true.

                case ("-om") 
                    call getarg(i+1, out_magdip)
                    argument_retrieved=.true.
                    
                case ("-on") 
                    call getarg(i+1, out_nac)
                    argument_retrieved=.true.

                case ("-der")
                    derivatives=.true.
                case ("-noder")
                    derivatives=.false.
                    
                case ("-nac")
                    do_nac=.true.
                case ("-nomac")
                    do_nac=.false.

                case ("-Si") 
                    call getarg(i+1, arg)
                    read(arg,*) Si 
                    argument_retrieved=.true.

                case ("-Sf") 
                    call getarg(i+1, arg)
                    read(arg,*) Sf
                    argument_retrieved=.true.
                    
                case ("-filt")
                    call getarg(i+1, filter)
                    argument_retrieved=.true.
        
                case ("-h")
                    need_help=.true.

                case default
                    print*, "Unkown command line argument: "//adjustl(arg)
                    need_help = .true.
            end select
        enddo 

        ! Post-processing
        !----------------------------
        if (adjustl(out_eldip) == 'default') then
            ! Get relative path
            if (index(inpfile,"./") /= 0) then
                call split_line_back(inpfile,"./",prfx,out_eldip)
                prfx=trim(adjustl(prfx))//"./"
            else
                out_eldip = inpfile
                prfx=""
            endif
            call split_line_back(out_eldip,".",out_eldip,arg)
            if (adjustl(ft) /= 'guess') arg=ft
            out_eldip = & !trim(adjustl(prfx))//&
                      "eldip_"//trim(adjustl(out_eldip))//'_'//trim(adjustl(arg))
        endif
        if (adjustl(out_magdip) == 'default') then
            ! Get relative path
            if (index(inpfile,"./") /= 0) then
                call split_line_back(inpfile,"./",prfx,out_magdip)
                prfx=trim(adjustl(prfx))//"./"
            else
                out_magdip = inpfile
                prfx=""
            endif
            call split_line_back(out_magdip,".",out_magdip,arg)
            if (adjustl(ft) /= 'guess') arg=ft
            out_magdip = & !trim(adjustl(prfx))//&
                      "magdip_"//trim(adjustl(out_magdip))//'_'//trim(adjustl(arg))
        endif
        if (adjustl(out_nac) == 'default') then
            ! Get relative path
            if (index(inpfile,"./") /= 0) then
                call split_line_back(inpfile,"./",prfx,out_nac)
                prfx=trim(adjustl(prfx))//"./"
            else
                out_nac = inpfile
                prfx=""
            endif
            call split_line_back(out_nac,".",out_nac,arg)
            if (adjustl(ft) /= 'guess') arg=ft
            out_nac = & !trim(adjustl(prfx))//&
                      "nac_"//trim(adjustl(out_nac))//'_'//trim(adjustl(arg))
        endif


       !Print options (to stderr)
        if (need_help) then

        write(0,'(/,A)') ' gen_fcc_dipfile '
        write(0,'(A)'  ) '-----------------'
        write(0,'(A)'  ) 'Generates dipfiles (eldip, magdip) for FCclasses'
        write(0,'(A)'  ) 'from the output of a QM program'

        call print_version()

        write(0,'(/,A)') 'SYNOPSIS'
        write(0,'(A)'  ) 'gen_fcc_dipfile -i input_file [-ft filetype] [-oe out_eldip]'//&
                         ' [-om out_magdip] [-Si <initial state>] [-Sf <final state>] [-noder] [-h]'

        write(0,'(/,A)') 'OPTIONS'
        write(0,'(A)'  ) 'Flag      Description         Current Value'
        write(0,'(A)'  ) ' -i       input_file          '//trim(adjustl(inpfile))
        write(0,'(A)'  ) ' -ft      filetype            '//trim(adjustl(ft))
        write(0,'(A,I0)')' -Si      initial state       ',Si
        write(0,'(A)')   '          (-1=default)        '
        write(0,'(A,I0)')' -Sf      final state         ',Sf
        write(0,'(A)')   '          (-1=default)        '
        write(0,'(A)'  ) ' -oe      out_eldip           '//trim(adjustl(out_eldip))
        write(0,'(A)'  ) ' -om      out_magdip          '//trim(adjustl(out_magdip))
        write(0,'(A)'  ) ' -on      out_nac             '//trim(adjustl(out_nac))
        write(0,'(A)'  ) ' -filt    Filter atoms (ders) '//trim(adjustl(filter))
        write(0,'(A,L1)'  ) ' -[no]der get derivatives  ',derivatives
        write(0,'(A,L1)'  ) ' -[no]nac get NAC          ',do_nac
        write(0,'(A)'  ) ' -h       print help  '
        call supported_filetype_list('trdip')

        write(0,'(/,A)') 'EXAMPLES'
        write(0,'(A)'  ) ' get_fcc_dipfile -i qmfile.fchk'
        write(0,'(A)'  ) '   A simple instruction as above is normally sufficient to get'
        write(0,'(A)'  ) '   eldip_qmfile_fchk and magdip_qmfile_fchk (in this example)'
        write(0,'(A)'  ) ''
        write(0,'(A)'  ) ' get_fcc_dipfile -i qmfile-psi4,out -ft psi4 -Si 2 -Sf 3'
        write(0,'(A)'  ) '   If the QM program computes them (e.g. Psi4 through EOM-CCSD)'
        write(0,'(A)'  ) '   both the initial and final states can be selected'
        write(0,'(A)'  ) ''

        stop    
        endif

        return
    end subroutine parse_input

end program gen_fcc_dipfile
