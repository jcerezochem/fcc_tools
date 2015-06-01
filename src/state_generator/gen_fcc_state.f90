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
    real(8),dimension(:),allocatable   :: X,Y,Z,Mass
    real(8),dimension(:),allocatable   :: Hlt, Freq
    real(8),dimension(:,:),allocatable :: L

    !Additional info to prepare the input
    real(8) :: DE, T

    !Auxiliars
    character :: cnull
    integer :: error
    !Counters
    integer :: i, j

    !I/O
    character(len=100) :: inpfile='input.log', &
                          outfile='default'
    character(len=10)  :: ft='guess'
    integer :: ios
    integer :: I_INP = 11, &
               O_STA = 20, &
               O_FCI = 21

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

    !Allocate input data
    allocate(X(1:3*Nat),Y(1:3*Nat),Z(1:3*Nat),Mass(1:3*Nat))
    allocate(Hlt(1:3*Nat*(3*Nat+1)/2))
    !Allocate output data
    allocate(Freq(1:3*Nat))
    allocate(L(1:3*Nat,1:3*Nat))

    !Read structure
    print*, "Reading structure..."
    call generic_structure_reader(I_INP,ft,Nat,X,Y,Z,Mass,error)
    if (error /= 0) then
        print*, "Error reading the geometry"
        stop
    else
        print'(X,A,/)', "OK"
    endif
! do i=1,Nat
! print*, X(i), Y(i), Z(i), Mass(i)
! enddo

    !Read Hessian
    print*, "Reading Hessian..."
    call generic_Hessian_reader(I_INP,ft,Nat,Hlt,error)
    if (error /= 0) then
        print*, "Error reading the Hessian"
        stop
    else
        print'(X,A,/)', "OK"
    endif
! print*, Hlt(1) , Hlt(3*Nat*(3*Nat+1)/2)

    !Close input file
    close(I_INP)

    !Perform vibrational analysis
    print*, "Diagonalizing Hessian..."
    call diag_int(Nat,X,Y,Z,Mass,Hlt,Nvib,L,Freq,error)
    if (error /= 0) then
        print*, "Error in the diagonalization"
        stop
    endif
    !Transform L to Cartesian/Normalized
    call Lmwc_to_Lcart(Nat,Nvib,Mass,L,L,error)
    call Lcart_to_LcartNrm(Nat,Nvib,L,L,error)
    !Transform Force Constants to Freq
    do i=1,Nvib
        Freq(i) = dsign(dsqrt(dabs(Freq(i))*HARTtoJ/BOHRtoM**2/UMAtoKG)/2.d0/pi/clight/1.d2,Freq(i))
    enddo
    if (error /= 0) then
        print*, "Error in conversion to Cartesian"
        stop
    else
        print'(X,A,/)', "OK"
    endif

    !WRITE STATE FILE
    print*, "Writting state file..."
    open(O_STA,file=outfile,status="replace",iostat=ios)
    if (ios /= 0) then
        print*, "Cannot open "//trim(adjustl(outfile))//" to write"
        stop
    endif

    do j=1,Nat
        write(O_STA,'(E17.8)',iostat=ios) X(j),Y(j),Z(j)
    enddo
    do j=1,3*Nat
    do i=1,Nvib
        write(O_STA,'(E17.8)',iostat=ios) L(j,i)
    enddo
    enddo
    do j=1,Nvib
        write(O_STA,'(F10.4)',iostat=ios) Freq(j)
    enddo
    close(O_STA)

    if (ios /= 0) then
        print*, "Error writting state file"
        stop
    else
        print'(X,A,/)', "OK"
    endif

    if (error == 0) then
        print*, "** Successful end **"
    else
        print*, "** Finish with erors **"
    endif

    !We profit to generate a first input template
    open(O_FCI,file="fcc_template.inp")
    DE = 2.d0
    T  = 1.d-3 
    call prepare_fccinput(O_FCI,Nat,Nvib,Mass,DE,T,error)
    close(O_FCI)

    stop

    contains

    subroutine parse_input(inpfile,ft,outfile)

        character(len=*),intent(inout) :: inpfile,ft,outfile

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

                case ("-o") 
                    call getarg(i+1, outfile)
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
        if (adjustl(outfile) == 'default') then
            call split_line(inpfile,".",outfile,arg)
            if (adjustl(ft) /= 'guess') arg=ft
            outfile = "state_"//trim(adjustl(outfile))//'_'//trim(adjustl(arg))
        endif


       !Print options (to stderr)
        if (need_help) then

        write(0,'(/,A)') ' gen_fcc_state '
        write(0,'(A)'  ) '-----------------'
        write(0,'(A)'  ) 'Generates state_files for FCclasses from the output'
        write(0,'(A)'  ) 'files obtained with different QM codes, reading the '
        write(0,'(A)'  ) 'coordinates and the Hessian from them.'
        write(0,'(A)'  ) 'Additionally, an input template is also generated'

        write(0,'(/,A)') 'SYNOPSIS'
        write(0,'(A)'  ) 'gen_fcc_state -i input_file [-ft filetype] [-o output_file] [-h]'

        write(0,'(/,A)') 'OPTIONS'
        write(0,'(A)'  ) 'Flag   Description      Current Value'
        write(0,'(A)'  ) ' -i    input_file       '//trim(adjustl(inpfile))
        write(0,'(A)'  ) ' -ft   filetype         '//trim(adjustl(ft))
        write(0,'(A)'  ) ' -o    output_file      '//trim(adjustl(outfile))
        write(0,'(A)'  ) ' -h    print help  '
        call supported_filetype_list('freq')

        stop    
        endif

        return
    end subroutine parse_input

end program gen_fcc_state

