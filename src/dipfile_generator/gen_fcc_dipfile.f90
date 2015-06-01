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

    !Relevant data
    integer :: Nat
    character(len=250) :: jobname
    character(len=20)  :: jobtype,method,basis
    real(8),dimension(1:3) :: Dip
    real(8),dimension(:),allocatable :: DipD
    logical :: derivatives

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
    character(len=100) :: inpfile, out_eldip, out_magdip
    character(len=5)   :: ft
    integer :: ios
    integer :: I_INP=10,&
               O_DIP =20

    ! Read options
    call parse_input(inpfile,ft,outfile)

    !Open input file
    open(I_INP,file=inpfile,iostat=ios)
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
    call generic_eldip_reader(I_INP,ft,iGS,iES,Dip,error)
    if (derivatives) then
        !Allocate output array
        allocate(DipD(1:3*3*Nat)
        call generic_eldipDer_reader(I_INP,ft,iGS,iES,Nat,DipD,error)
    endif

    !WRITE ELDIP FILE
    open(O_DIP,file=out_eldip,status="replace")

    !One should replace one of these by that at the other state geom
    write(O_DIP,'(3(X,E18.9))') Dip(1:3)
    write(O_DIP,'(3(X,E18.9))') Dip(1:3)

    if (derivatives) then
        do j=1,3*Nat
            jj=j*3-2
            write(O_DIP,'(3(X,E18.9))') DipD(jj:jj+2)
        enddo
    endif
    close(O_DIP)

    if (.not.get_magdip) then
        stop
     endif

    !Get magdip
    call generic_magdip_reader(I_INP,ft,iGS,iES,Dip,error)
    if (derivatives) then
        call generic_magdipDer_reader(I_INP,ft,iGS,iES,Nat,DipD,error)
    endif

    !WRITE ELDIP FILE
    open(O_DIP,file=out_magdip,status="replace")

    !One should replace one of these by that at the other state geom
    write(O_DIP,'(3(X,E18.9))') Dip(1:3)
    write(O_DIP,'(3(X,E18.9))') Dip(1:3)

    if (derivatives) then
        do j=1,3*Nat
            jj=j*3-2
            write(O_DIP,'(3(X,E18.9))') DipD(jj:jj+2)
        enddo
    endif
    close(O_DIP)

    deallocate(DipD)

    stop

    contains

    subroutine parse_input(inpfile,ft,out_eldip,out_magdip)

        character(len=*),intent(inout) :: inpfile,ft,out_eldip,out_magdip

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
            call split_line(inpfile,".",out_eldip,arg)
            if (adjustl(ft) /= 'guess') arg=ft
            out_eldip  = "eldip_"//trim(adjustl(out_eldip))//'_'//trim(adjustl(arg))
        endif
        if (adjustl(out_magdip) == 'default') then
            call split_line(inpfile,".",out_magdip,arg)
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
        write(0,'(A)'  ) 'gen_fcc_dipfile -i input_file [-ft filetype] [-oe out_eldip] [-om out_magdip] [-h]'

        write(0,'(/,A)') 'OPTIONS'
        write(0,'(A)'  ) 'Flag      Description      Current Value'
        write(0,'(A)'  ) ' -i       input_file       '//trim(adjustl(inpfile))
        write(0,'(A)'  ) ' -ft      filetype         '//trim(adjustl(ft))
        write(0,'(A)'  ) ' -oe      out_eldip        '//trim(adjustl(out_eldip))
        write(0,'(A)'  ) ' -om      out_magdip       '//trim(adjustl(out_magdip))
        write(0,'(A)'  ) ' -[no]der get derivatives  '
        write(0,'(A)'  ) ' -h       print help  '
        call supported_filetype_list('trdip')

        stop    
        endif

        return
    end subroutine parse_input

end program gen_fcc_dipfile
