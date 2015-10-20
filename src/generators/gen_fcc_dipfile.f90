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

    !Relevant data
    integer :: Nat
    character(len=250) :: jobname
    character(len=20)  :: jobtype,method,basis
    real(8)                :: dx = -1.d0 
    real(8),dimension(1:3) :: Dip
    real(8),dimension(:),allocatable :: DipD
    logical :: derivatives=.true.
    !Default states are taken if =-1
    integer :: Si=-1, &
               Sf=-1

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
    integer :: j,k, jj

    !I/O
    character(len=100) :: inpfile='input.log', &
                          out_eldip='default', &
                          out_magdip='default'
    character(len=10)  :: ft='guess'
    integer :: ios
    integer :: I_INP=10,&
               O_DIP =20

    ! Read options
    call parse_input(inpfile,ft,out_eldip,out_magdip,derivatives,Si,Sf)

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
        do j=1,3*Nat
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
        do j=1,3*Nat
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

    print*, "** Successful end **"

    stop

    contains

    subroutine parse_input(inpfile,ft,out_eldip,out_magdip,derivatives,Si,Sf)

        character(len=*),intent(inout) :: inpfile,ft,out_eldip,out_magdip
        logical,intent(inout)          :: derivatives
        integer,intent(inout)          :: Si, Sf

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

                case ("-der")
                    derivatives=.true.
                case ("-noder")
                    derivatives=.false.

                case ("-Si") 
                    call getarg(i+1, arg)
                    read(arg,*) Si 
                    argument_retrieved=.true.

                case ("-Sf") 
                    call getarg(i+1, arg)
                    read(arg,*) Sf
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
            call split_line_back(inpfile,".",out_eldip,arg)
            if (adjustl(ft) /= 'guess') arg=ft
            out_eldip  = "eldip_"//trim(adjustl(out_eldip))//'_'//trim(adjustl(arg))
        endif
        if (adjustl(out_magdip) == 'default') then
            call split_line_back(inpfile,".",out_magdip,arg)
            if (adjustl(ft) /= 'guess') arg=ft
            out_magdip  = "magdip_"//trim(adjustl(out_magdip))//'_'//trim(adjustl(arg))
        endif


       !Print options (to stderr)
        if (need_help) then

        write(0,'(/,A)') ' gen_fcc_dipfile '
        write(0,'(A)'  ) '-----------------'
        write(0,'(A)'  ) 'Generates dipfiles (eldip, magdip) for FCclasses'
        write(0,'(A)'  ) 'from the output of a QM program'

        write(0,'(/,A)') 'SYNOPSIS'
        write(0,'(A)'  ) 'gen_fcc_dipfile -i input_file [-ft filetype] [-oe out_eldip]'//&
                         ' [-om out_magdip] [-Si <initial state>] [-Sf <final state>] [-noder] [-h]'

        write(0,'(/,A)') 'OPTIONS'
        write(0,'(A)'  ) 'Flag      Description      Current Value'
        write(0,'(A)'  ) ' -i       input_file       '//trim(adjustl(inpfile))
        write(0,'(A)'  ) ' -ft      filetype         '//trim(adjustl(ft))
        write(0,'(A,I0)')' -Si      initial state    ',Si
        write(0,'(A)')   '          (-1=default)     '
        write(0,'(A,I0)')' -Sf      final state      ',Sf
        write(0,'(A)')   '          (-1=default)     '
        write(0,'(A)'  ) ' -oe      out_eldip        '//trim(adjustl(out_eldip))
        write(0,'(A)'  ) ' -om      out_magdip       '//trim(adjustl(out_magdip))
        write(0,'(A,L1)'  ) ' -[no]der get derivatives  ',derivatives
        write(0,'(A)'  ) ' -h       print help  '
        call supported_filetype_list('trdip')

        stop    
        endif

        return
    end subroutine parse_input

end program gen_fcc_dipfile
