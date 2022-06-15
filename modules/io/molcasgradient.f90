program molcasgradient

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    !Description
    ! Program to transform molcas gradient (including symmetry)
    ! into a more general (unsymmetrized) file
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
    real(8),dimension(:),allocatable   :: Hlt, Freq, Grad
    real(8),dimension(:,:),allocatable :: L

    ! Harmonic model PES 
    character(len=3) :: model_PES='AH'

    !Additional info to prepare the input
    real(8) :: DE, T
    logical :: is_hessian = .true.

    !Auxiliars
    character :: cnull
    integer :: error
    integer,dimension(:),allocatable :: IAux
    !Counters
    integer :: i, j

    !I/O
    character(len=100) :: strfile='input.log', &
                          hessfile='default',&
                          gradfile='default',&
                          massfile='none',&
                          outfile='default',   &
                          outhess='default',   &
                          outmass='default'
    character(len=10)  :: fts='guess', &
                          fth='guess', &
                          ftg='guess'
    integer :: ios
    integer :: I_INP = 11, &
               I_HES = 12, &
               I_GRD = 14, &
               I_MAS = 13, &
               O_STA = 20, &
               O_FCI = 21, &
               O_HES = 22, &
               O_MAS = 23

    ! Read options
    call parse_input(strfile,fts,hessfile,fth,gradfile,ftg,massfile,outfile,outhess,outmass,model_pes)
    call set_word_upper_case(model_pes)

    !Open input file
    open(I_INP,file=strfile,iostat=ios)
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
    rewind(I_INP)

    !Allocate input data
    allocate(X(1:3*Nat),Y(1:3*Nat),Z(1:3*Nat),Mass(1:3*Nat))
    allocate(Hlt(1:3*Nat*(3*Nat+1)/2))
    if (adjustl(model_pes) == "VH") allocate (Grad(1:3*Nat))
    !Allocate output data
    allocate(Freq(1:3*Nat))
    allocate(L(1:3*Nat,1:3*Nat))

    !Read structure
    print*, "Reading structure..."
    call generic_structure_reader(I_INP,fts,Nat,X,Y,Z,Mass,error)
    if (error /= 0) then
        print*, "Error reading the geometry", error
        stop
    endif

    if (adjustl(massfile) /= "none") then
        open(I_MAS,file=massfile)
        print*, "  but masses read from external massfile: "//trim(adjustl(massfile))
        do i=1,Nat
            read(I_MAS,*) Mass(i)
        enddo
        close(I_MAS)
    endif
    print*, "  and writting masses to file..."
    open(O_MAS,file=outmass)
    do i=1,Nat 
        write(O_MAS,*) Mass(i)
    enddo
    write(O_MAS,*) ""
    close(O_MAS)
    print'(X,A,/)', "OK"

    close(I_INP)

    ! We now read the allow to read the Hessian from another file
    !Guess the file type if not given
    if (adjustl(fth) == 'guess') then
        call split_line_back(hessfile,'.',cnull,fth)
    else if (adjustl(fth) == 'fts') then
        fth = fts
    endif
    if (adjustl(ftg) == 'guess') then
        call split_line_back(gradfile,'.',cnull,fth)
    else if (adjustl(ftg) == 'fth') then
        ftg = fth
    endif
    

    !Read Gradient for VH (for now assume same file as Hess)
    if (adjustl(model_pes) == "VH") then
        !Open gradient file
        open(I_GRD,file=gradfile,iostat=ios)
        if (ios /= 0) then
            print*, "Error opening "//trim(adjustl(gradfile))
            stop
        endif
        print*, "Reading Gradient..."
        call generic_gradient_reader(I_HES,fth,Nat,Grad,error)
        if (error /= 0) then
            print'(X,A,/)', "Error: Gradient not present in the file."
            stop
        else
            print'(X,A,/)', "OK"
        endif
    endif




    stop

    contains

    subroutine parse_input(inpfile,ft)

        character(len=*),intent(inout) :: inpfile,ft

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
                case ("-ft") !for backward compatibility
                    call getarg(i+1, ft)
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

        print*, "No manual (for the moment)"

        stop    
        endif

        return
    end subroutine parse_input

end program molcasgradient

