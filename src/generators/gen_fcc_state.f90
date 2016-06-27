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
    real(8),dimension(:),allocatable   :: Hlt, Freq, Grad
    real(8),dimension(:,:),allocatable :: L

    ! Harmonic model PES 
    character(len=3) :: model_PES='AH'
    logical :: force_real = .false.

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
    call parse_input(strfile,fts,hessfile,fth,gradfile,ftg,massfile,outfile,outhess,outmass,model_pes,force_real)
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
    print*, "Reading structure...", fts
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
        print*, "Reading Gradient...", ftg
        call generic_gradient_reader(I_GRD,ftg,Nat,Grad,error)
        if (error /= 0) then
            print'(X,A,/)', "Error: Gradient not present in the file."
            stop
        else
            print'(X,A,/)', "OK"
        endif
        close(I_GRD)
    endif

    !Open hessian file
    open(I_HES,file=hessfile,iostat=ios)
    !Read Hessian
    print*, "Reading Hessian...", fth
    call generic_Hessian_reader(I_HES,fth,Nat,Hlt,error)
    if (error /= 0) then
        print'(X,A,/)', "Hessian is not present in the file. Only valid for AS"
        is_hessian = .false.
    else
        print*, "  and writting hessian (lower triangular elements) to file..."
        open(O_HES,file=outhess)
        do i=1,3*Nat*(3*Nat+1)/2
            write(O_HES,'(G16.8)') Hlt(i)
        enddo
        write(O_HES,*) ""
        close(O_HES)
        print'(X,A,/)', "OK"
    endif

    !Close hessian file
    close(I_HES)

    if (is_hessian .and. adjustl(model_pes) /= "VH") then
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
    endif

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
        if (is_hessian) then
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

    !We profit to generate a first input template
    print*, "Writting input template: 'fcc_template.inp'..."
    open(O_FCI,file="fcc_template.inp")
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

    print*, "** Successful end **"

    !Deallocate
    deallocate(X,Y,Z,Mass,Hlt,Freq,L)
    if (adjustl(model_pes) == "VH") deallocate(Grad)

    stop

    contains

    subroutine parse_input(strfile,fts,hessfile,fth,gradfile,ftg,massfile,outfile,outhess,outmass,model_pes,force_real)

        character(len=*),intent(inout) :: strfile,fts,hessfile,fth,gradfile,ftg,massfile,&
                                          outfile,outhess,outmass,model_pes
        logical,intent(inout)          :: force_real

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

                case ("-ig") 
                    call getarg(i+1, gradfile)
                    argument_retrieved=.true.
                case ("-ftg") 
                    call getarg(i+1, ftg)
                    argument_retrieved=.true.

                case ("-im") 
                    call getarg(i+1, massfile)
                    argument_retrieved=.true.

                case ("-o") 
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.

                case ("-oh") 
                    call getarg(i+1, outhess)
                    argument_retrieved=.true.

                case ("-om") 
                    call getarg(i+1, outmass)
                    argument_retrieved=.true.

                case ("-model") 
                    call getarg(i+1, model_pes)
                    argument_retrieved=.true.

                case ("-force-real")
                    force_real=.true.
                case ("-noforce-real")
                    force_real=.false.
        
                case ("-h")
                    need_help=.true.

                case default
                    print*, "Unkown command line argument: "//adjustl(arg)
                    need_help = .true.
            end select
        enddo 

        ! Post-processing
        !----------------------------
        if (adjustl(hessfile) == 'default') then
            ! The try to read the hessian from strfile
            hessfile = strfile
            fth = "fts"
        endif
        if (adjustl(gradfile) == 'default') then
            ! The try to read the hessian from strfile
            gradfile = hessfile
            ftg = "fth"
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
            outfile = trim(adjustl(prfx))//&
                      "state_"//trim(adjustl(outfile))//'_'//trim(adjustl(arg))
        endif
        if (adjustl(outhess) == 'default') then
            call split_line_back(outfile,"state_",arg,outhess)
            if (adjustl(fts) /= 'guess') arg=fts
            outhess = trim(adjustl(prfx))//&
                      "hessian_"//trim(adjustl(outhess))
        endif
        if (adjustl(outmass) == 'default') then
            call split_line_back(outfile,"state_",arg,outmass)
            if (adjustl(fts) /= 'guess') arg=fts
            outmass = trim(adjustl(prfx))//&
                      "mass_"//trim(adjustl(outmass))
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
        write(0,'(A)'  ) 'The input structure (-i) is read from any of the supported QM output'
        write(0,'(A)'  ) 'files or through the g96 (GROMOS) structure file. The file type is'
        write(0,'(A)'  ) 'guessed based on the extension, but if this fails, it can be specified'
        write(0,'(A)'  ) 'by -fts flag. If the structure file also contains the hessian, it will'
        write(0,'(A)'  ) 'be read from that file by default. Otherwise, another file containing the'
        write(0,'(A)'  ) 'Hessian should be indicated (-ih) with filetype detected based on extension'
        write(0,'(A)'  ) 'of specified (-fth). The output is named automatiacally as state_{input}_{ext}'
        write(0,'(A)'  ) 'or can be indicated (-o). Masses are taken from an internal database but you'
        write(0,'(A)'  ) 'can specify different values in a file (with Nat lines, each with the mass of'
        write(0,'(A)'  ) 'the atom in amu) through -im flag. The Hessian and masses used in the calculations'
        write(0,'(A)'  ) 'are writen to files, which can be specified on input (-oh -om).'

        call print_version()

        write(0,'(/,A)') 'SYNOPSIS'
        write(0,'(A)'  ) 'gen_fcc_state -i input_file [-fts filetype-str] [-ih hess_inp_file] [-fth filetype-hess] '//&
                         '[-o output_file] [-oh hessian_file] [-om mass_file] [-model model_PES] [-h]'

        write(0,'(/,A)') 'OPTIONS'
        write(0,'(A)'  ) 'Flag    Description      Current Value'
        write(0,'(A)'  ) ' -i     structure_file   '//trim(adjustl(strfile))
        write(0,'(A)'  ) ' -fts   filetype(str)    '//trim(adjustl(fts))
        write(0,'(A)'  ) ' -ih    hess_input_file  '//trim(adjustl(hessfile))
        write(0,'(A)'  ) ' -fth   filetype(hess)   '//trim(adjustl(fth))
        write(0,'(A)'  ) ' -ig    grad_input_file  '//trim(adjustl(gradfile))
        write(0,'(A)'  ) ' -ftg   filetype(grad)   '//trim(adjustl(ftg))
        write(0,'(A)'  ) ' -im    mass_file        '//trim(adjustl(massfile))
        write(0,'(A)'  ) ' -o     output_file      '//trim(adjustl(outfile))
        write(0,'(A)'  ) ' -oh    hess_out_file    '//trim(adjustl(outhess))
        write(0,'(A)'  ) ' -om    mass_file        '//trim(adjustl(outmass))
        write(0,'(A)'  ) ' -model model_pes[AH|VH] '//trim(adjustl(model_pes))
        write(0,'(A)'  ) ' -force-real  turn real ',  force_real
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

