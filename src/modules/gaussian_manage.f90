module gaussian_manage
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012!

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get information from 
    !  Gaussian log files. This is a med-low level module that
    !  only requires low-level subroutines (i.e. do not use 
    !  structure_types module) 
    !
    ! Notes
    !  Only low level modules are required: line_preprocess and alerts.
    ! 
    !==============================================================

    !Common declarations:
    !===================
    use line_preprocess
    use alerts
    implicit none

    contains

    subroutine summary_parser(unt,isect,section,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! New summary parser, that simply reads in the summary section 
        ! and passes a piece of it back to the calling program as a string.
        ! So, there is no need of any special types in the subroutine.
        ! 
        ! The max length of each summary line is hardcoded to 240
        !
        !Arguments
        ! unt     (inp) int /scalar   unit for the file 
        ! isect   (inp) int /scalar   number of the section
        !                             1: General info (host,user,date and 
        !                                job info: job type, method, basis set
        !                                molecular formula)
        !                             2: Route section
        !                             3: Job title
        !                             4: Molecular geometry (and atnames)
        !                             5: Electronic properties
        !                             6: Hessian (lower triangular form)
        !                             7: Gradient
        ! section (out) char/scalar   string with the selected section
        ! error_flag (out) flag  0 : success
        !                        1 : Requested secction not present in the summary
        !                            Returns the last section read
        !                     10+i : read error in line i
        !                       -i : the section size is i characters larger then
        !                            length of the character. Additional i character 
        !                            are required. The output includes the last
        !                            part of the section that fits the length
        !
        !==============================================================


        integer,intent(in)           :: unt
        integer,intent(in)           :: isect
        character(len=*),intent(out) :: section
        integer,intent(out),optional :: error_flag

        !Local
        ! Counters
        integer :: i, isection
        ! I/O
        integer :: IOstatus
        character(len=240) :: line

        !Locate summary
        error_flag = 0
        i=0
        do
            i=i+1
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) then
                error_flag = i+10
                rewind(unt)
                return
            endif
            if ( INDEX(line,"GINC") /= 0 ) exit
        enddo

        !Read sections
        isection = 1
        do isection=1,7
            error_flag = 0
            section    = ""
            do
                i=i+1
                !Check if there is space for a line. Otherwise, make 
                ! space and continue
                if (len(section)-len_trim(section) <= len_trim(line)) then
                    error_flag = error_flag - len_trim(line)
                    section=section(len_trim(line)+1:len(section))
                endif
                !Append last line to section
                section=adjustl(trim(section))//line
                !If we reach a section end "\\", split and finish section
                if ( INDEX(section,'\\') /= 0 ) then
                    call split_line(section,'\\',section,line)
                    exit
                endif
                !Read in new line
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                if ( IOstatus < 0 ) then
                    error_flag = i+10
                    rewind(unt)
                    return
                endif
            enddo
            if (isection == isect) then
                rewind(unt)
                return
            endif
            if (line == "@") then
                error_flag = 1
                rewind(unt)
                return
            endif

        enddo

        rewind(unt)
        return

    end subroutine summary_parser


    subroutine estimate_section_length(unt,isect,length,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Estimate the length of the section. This ca be used to allocate
        ! the character to hold the section content
        !
        !Arguments
        ! unt     (inp) int /scalar   unit for the file 
        ! isect   (inp) int /scalar   number of the section
        !                             1: General info (host,user,date and 
        !                                job info: job type, method, basis set
        !                                molecular formula)
        !                             2: Route section
        !                             3: Job title
        !                             4: Molecule specifications (atname,geom)
        !                             5: Electronic properties
        !                             6: Hessian (lower triangular form)
        !                             7: Gradient
        ! section (out) int /scalar   Length of the string to store the section isect
        ! error_flag (out) flag  0 : success
        !                        1 : Requested secction not present in the summary
        !                            Returns the last section read
        !                     10+i : read error in line i
        !                       -i : the section size is i characters larger then
        !                            length of the character. Additional i character 
        !                            are required. The output includes the last
        !                            part of the section that fits the length
        !
        !==============================================================


        integer,intent(in)           :: unt
        integer,intent(in)           :: isect
        integer,intent(out)          :: length
        integer,intent(out),optional :: error_flag

        !Local
        ! Counters
        integer :: i, isection
        ! I/O
        integer :: IOstatus
        character          :: cnull
        character(len=240) :: line

        !Locate summary
        error_flag = 0
        i=0
        do
            i=i+1
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) then
                error_flag = i+10
                rewind(unt)
                return
            endif
            if ( INDEX(line,"GINC") /= 0 ) exit
        enddo

        !Read sections
        isection = 1
        do isection=1,7
            length     = 0
            error_flag = 0
            do
                !Update length
                length=length+len_trim(line)

                !If we reach a section end "\\", pass to the next section
                if ( INDEX(line,'\\') /= 0 ) then
                    call split_line(line,'\\',cnull,line)
                    exit
                endif

                !Read in new line
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                if ( IOstatus < 0 ) then
                    error_flag = i+10
                    rewind(unt)
                    return
                endif
            enddo
            if (isection == isect) then
                rewind(unt)
                return
            endif
            if (line == "@") then
                error_flag = 1
                rewind(unt)
                return
            endif

        enddo

        rewind(unt)
        return

    end subroutine estimate_section_length

    subroutine read_gausslog_property(unt,property,value,io_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Get electronic property from the summary (section=5). Returs the
        ! the value in a string (so as to avoid specifying the type)
        ! USES ROUTINES IN THIS MODULE
        !
        !Arguments
        ! unt      (inp) int /scalar   unit for the file 
        ! property (inp) char/scalar   Name of the property to get
        ! value    (out) char/scalar   Value for the property
        ! io_flag  (io ) flag          On input : length of properties section to allocate
        !                                   an auxiliar character (prop_section)
        !                              On output: error_flag
        !                                   1 : Secction 5 not present in the summary
        !                                   2 : Requested property not found
        !                                   3 : Requested property not found, and 
        !                                       'prop_section' was not allocatated properly
        !                                   4 : Warning: 'prop_section' was not allocated 
        !                                       properly and the value might be incorrect
        !                                10+i : read error in line i (in summary_parser)
        !                                  -i : 'value' length is not enough to store the value
        !                                        It should be at least len=i
        !
        !Notes:
        ! -About allocatable characters:
        !  Fortran2003 supports allocatable strings (implemented in gfortran >=4.8, ifort v12 and maybe lower)
        !  This would overcome the need of explicetily allocating the
        !  auxiliar string with io_flag (but would impose a requirement to the compiler version)
        !
        !==============================================================

        integer,intent(in)           :: unt
        character(len=*),intent(in)  :: property
        character(len=*),intent(out) :: value
        integer,intent(inout)        :: io_flag
        
        !local
        character(len=io_flag)       :: prop_section
        character(len=len(property)+len(value))    :: current_property
        character(len=5)             :: auxchar
        integer :: i,j


        call summary_parser(unt,5,prop_section,io_flag)
        if (io_flag > 0) return
        !If the error_flag<0, prop_section was not allocated properly
        !(not large enough). We try to read the property anyway
        !since it might be present (if not, an error is risen)

        !Run over all properties
        do while ( len_trim(prop_section) /= 0 )
            i=len_trim(prop_section)
            call split_line(prop_section,'\',current_property,prop_section)
            !Get the value
            call split_line(current_property,'=',current_property,value)

            !Check if it is the requested property          
            if ( adjustl(current_property) == adjustl(property) ) then
                !Check if 'value' was large enough
                j = i-len_trim(prop_section)-len_trim(property)-1
                if ( j > len(value) ) then
                    io_flag=-j
                    call alert_msg("warning","Property "//trim(adjustl(property))// &
                                   " cannot be allocated. Provide len of size "//auxchar)
                    rewind(unt)
                    return
                endif
                !Warning about summary_parser errors
                if ( (io_flag < 0) .and. len_trim(prop_section) == 0) io_flag=4
                rewind(unt)
                return
            endif
        enddo

        if (io_flag < 0) then
            io_flag = 3
        else
            io_flag = 2
        endif

        rewind(unt)
        return
        
    end subroutine read_gausslog_property


    subroutine read_gausslog_natoms(unt,natoms,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Get natoms from molecular formula. It is taken from section 1, 
        ! which contains
        !  1: flag? =1?
        !  2: flag? =1?
        !  3: GINC-$HOST
        !  4: Job Type (Opt, Freq...)
        !  5: Method (RB3LYP...)
        !  6: Basis set
        !  7: Molecular formula   <===
        !  8: $USER
        !  9: $date
        ! 10: flag? =0?
        ! USES ROUTINES IN THIS MODULE
        !
        !Arguments
        ! unt      (inp) int /scalar   unit for the file 
        ! natoms   (out) char/scalar   Number of atoms
        ! error_flag (io ) flag        Error flag
        !                                   1 : Secction 5 not present in the summary
        !                                   2 : Requested property not found
        !                                   3 : Requested property not found, and 
        !                                       'prop_section' was not allocatated properly
        !                                   4 : Warning: 'prop_section' was not allocatated 
        !                                       properly and the value might be incorrect
        !                                10+i : read error in line i
        !                                  -i : 'value' length is not enough to store the value
        !                                        It should be at least len=i
        !
        !Notes:
        ! -About allocatable characters:
        !  Fortran2003 supports allocatable strings (implemented in gfortran >=4.8, ifort v12 and maybe lower)
        !  This would overcome the need of explicetily allocating the auxiliar strings
        !  (but would impose a requirement to the compiler version)
        !
        !==============================================================
        
        integer,intent(in)           :: unt
        integer,intent(out)          :: natoms
        integer,intent(out),optional :: error_flag
        
        !local
        character(len=500)           :: info_section
        character(len=100)           :: formula
        integer                      :: i,j,n
        ! array with number of atoms per element in the formula
        !  len set the maximum number (len=3 => 999)
        !  dimension set the maximum different elements (50 seems enough)
        character(len=3),dimension(50) :: nel_char
        character :: cnull

        call summary_parser(unt,1,info_section,error_flag)

        !Get formula (7th entry)
        do i=1,7
            call split_line(info_section,'\',formula,info_section)
        enddo

        !Remove (charge,spin) tag
        call split_line(formula,"(",formula,cnull)

        do i=1,len_trim(formula)
            !Convert to blank if it is not a number (48==0 to 57==9)
            if (ichar(formula(i:i)) < 48 .or. &
                ichar(formula(i:i)) > 57) formula(i:i)=" "
        enddo

        !Use line parser to get the numbers reparatelly
        call parse_line(formula,n,nel_char)

        !Now read each one and sum them
        natoms = 0
        do i=1,n
            read(nel_char(i),*) j
            natoms = natoms + j
        enddo

        rewind(unt)
        return

    end subroutine read_gausslog_natoms


    subroutine read_gauslog_geom(unt,Nat,AtName,X,Y,Z,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Get geometry and atom names from section 4. The number of atoms
        ! is also taken (since they are needed), using read_gausslog_natoms
        ! USES ROUTINES IN THIS MODULE
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        ! AtName  (out) char/vertor    Atom names
        ! X,Y,Z   (out) real/vectors   Coordinate vectors (ANGSTRONG)
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                   1 : required Z-mat, but not implemented
        !                                   2 : Unkonwn format
        !                       +/-(100000+i) : Error in read_gausslog_natoms call: +/-i
        !                       +/-(200000+i) : Error in summary_parser call: +/-i
        !
        !Notes:
        ! TODO: Add Z-matrix support (currently rises and error)
        !
        !==============================================================

        integer,intent(in)                        :: unt
        integer,intent(out)                       :: Nat
        real(kind=8),dimension(:),intent(out)     :: X,Y,Z
        character(len=*),dimension(:),intent(out) :: AtName
        integer,intent(out),optional              :: error_flag

        !Local
        integer :: i, j, nitems
        character(len=100000) :: string 
        character(len=200)    :: geom_char 
        character,dimension(4) :: dummy_char

        ! Natoms
        call read_gausslog_natoms(unt,Nat,error_flag)
        if (error_flag < 0) then
            error_flag = error_flag - 100000
            rewind(unt)
            return
        elseif (error_flag > 0) then
            error_flag = error_flag + 100000
            rewind(unt)
            return
        endif

        !Geom
        call summary_parser(unt,4,string,error_flag)
        if (error_flag < 0) then
            error_flag = error_flag - 200000
            rewind(unt)
            return
        elseif (error_flag > 0) then
            error_flag = error_flag + 200000
            rewind(unt)
            return
        endif
        !Throw "charge mult" away
        call split_line(string,'\',geom_char,string)

        !Read the first geometry entry
        call split_line(string,'\',geom_char,string)
        !Get the number of elements
        call string2cvector(geom_char,dummy_char,nitems,',')
        !The decide the format
        if (nitems == 4) then
            !Read first entry
            read(geom_char,*) AtName(1), &
                              X(1),      &
                              Y(1),      &
                              Z(1)
            !Continue with the rest of the entries
            do i=2,Nat
                call split_line(string,'\',geom_char,string)
                read(geom_char,*) AtName(i), &
                                  X(i),      &
                                  Y(i),      &
                                  Z(i)
            enddo
        elseif (nitems == 5) then
            !Read first atom
            read(geom_char,*) AtName(1), &
                              j,         &
                              X(1),      &
                              Y(1),      &
                              Z(1)
            do i=2,Nat
                call split_line(string,'\',geom_char,string)
                read(geom_char,*) AtName(i), &
                                  j,         &
                                  X(i),      &
                                  Y(i),      &
                                  Z(i)
            enddo
        elseif (nitems == 1) then
            !Read Z-mat...
            error_flag=1
            call alert_msg("warning","Z-mat reading from G09 summary not yet implemented")
            rewind(unt)
            return
        else 
            error_flag=2
            call alert_msg("warning","Unexpected structure format in G09 summary section")
            rewind(unt)
            return
        endif

        rewind(unt)
        return

    end subroutine read_gauslog_geom


    subroutine read_fchk(unt,section,data_type,N,A,I,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        !Description
        ! Generic SR to read any section of the checkpoint
        ! Enter deallocated arrays
        !Arguments
        ! unt (int;in): unit number of the log file
        ! section(char,in): name of the section to be read
        ! data_type(char,out): Integer (I) or Real (R) data read
        ! N(int,in): Number of elements to be read
        ! A(real,dimension(:)): Real array to store real data
        ! I(integer,dimension(:)): Int array to store int data
        ! error_flag(integer,out): 0: success
        !                          1: section not found
        !==============================================================

        integer,intent(in)                                       :: unt
        character(len=*),intent(in)                              :: section
        character(len=1),intent(out)                             :: data_type
        integer,intent(out)                                      :: N
        double precision, dimension(:), allocatable, intent(out) :: A
        integer,dimension(:), allocatable, intent(out)           :: I
        integer,intent(out),optional                             :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        character(len=42)  :: section_full
        character(len=1)   :: is_array
        character(len=40)  :: cdata
        !I/O
        integer :: IOstatus
        
        
        ! Search section
        if (present(error_flag)) error_flag = 0
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    call alert_msg("warning","Section '"//trim(adjustl(section))//"' not present in the FCHK file.")
                    if (present(error_flag)) error_flag=1
                    rewind(unt)
                    return
                endif
                ! 2) Found what looked for!      
                if ( INDEX(line,trim(adjustl(section))) /= 0 ) then
                    read(line,'(A42)') section_full
                    if (adjustl(section_full) == adjustl(section)) exit
                endif
        enddo

        !Get info about section from line just read
        read(line,'(A42,X,A1,3X,A1,X,A)') section_full, data_type, is_array, cdata
        if (is_array /= "N") then
            !Is not an array
            N=1
            if ( data_type == "R" ) then 
                allocate( A(1:1) )
                read(cdata,*) A(1)
            elseif ( data_type == "I" ) then
                allocate( I(1:1) )
                read(cdata,*) I(1) 
            endif
        else
            read(cdata,*) N
            if ( data_type == "R" ) then
                allocate( A(1:N) )
                read(unt,*) A(1:N)
            elseif ( data_type == "I" ) then
                allocate( I(1:N) )
                read(unt,*) I(1:N)
            endif
        endif 

        rewind(unt)
        return

    end subroutine read_fchk

    subroutine write_fchk(unt,section,data_type,N,A,I,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Generic SR to read any section of the checkpoint
        ! Enter allocated arrays
        !Arguments
        ! unt (int;in): unit number of the log file
        ! section(char,in): name of the section to be written
        ! data_type(char,in): Integer (I) or Real (R) data write
        ! N(int,in): Number of elements to be written
        ! A(real,dimension(:)): Real array to store real data
        ! I(integer,dimension(:)): Int array to store int data
        ! error_flag(integer,out): 0: success
        !                          1: write failure
        !
        ! NOTES
        ! Should is_array be an input?
        ! Not really, we can specify that is a scalar, e.g. by setting 
        ! N=0 (so we let N=1 fro vectors with length equal 1, if this 
        ! is ever the case)
        !==============================================================

        integer,intent(in) :: unt
        character(len=*),intent(in) :: section
        character(len=1),intent(in) :: data_type
        integer,intent(in) :: N
        double precision, dimension(:), intent(in) :: A
        integer,dimension(:), allocatable, intent(in) :: I
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character(len=43) :: section_full
        !I/O
        integer :: IOstatus
        
        error_flag = 0
        section_full = adjustl(section)

        !If N=0, it is an scalar
        if (N == 0) then
            if ( data_type == "I" ) then 
                write(unt,'(A43,A,I17)') section_full, data_type, I(1)
            elseif ( data_type == "R" ) then
                write(unt,'(A43,A,ES27.15)') section_full, data_type, A(1)
            endif
        else
            write(unt,'(A43,A,A5,I12)') section_full, data_type,"   N=",N
            if ( data_type == "I" ) then 
                write(unt,'(6I12)') I(1:N)
            elseif ( data_type == "R" ) then
                write(unt,'(5ES16.8)') A(1:N)
            endif
        endif 

        return

    end subroutine write_fchk

    subroutine read_gaussfchk_dip(unt,Si,Sf,derivatives,dip_type,Dip,DipD,error_flag)

        !=====================================================
        ! THIS CODE IS PART OF FCC_TOOLS
        !=====================================================
        ! Description
        !  Read transition electric or magnetic dipole moments
        !  from G09 fchk files (corresponding to a ES calculation)
        !  The information is in "ETran ..." sections:
        !  *"ETran scalars" contains:
        !    <number of ES> <?> <?> <?> <target state> <?>
        !     0 0 0...
        !  *"ETran state values"  contains
        !    ·First the properties of each excited state (up to Nes):
        !    1,  2   , 3  , 4  , 5     , 6     , 7     , 8    , 9    , 10   , 11 , 12 , 13 , 14 , 15 , 16 ...
        !    E, {muNx,muNy,muNz,muvelNx,muvelNy,muvelNz,mmagNx,mmagNy,mmagNz,unkX,unkY,unkZ,unkX,unkY,unkZ}_N=1,Nes
        !    ·Then, the derivates of each property with respect to Cartesian coordiates only for target state
        !     For each Cartesian coordiate, all derivatives are shown:
        !     dE/dx1 dmux/dx1 dmuy/dx1 ... unkZ/dx1
        !     dE/dy1 dmux/dy1 dmuy/dy1 ... unkZ/dy1
        !     ...
        !     dE/dzN dmux/dzN dmuy/dzN ... unkZ/dzN
        !
        ! Notes
        !  Only the length-gauge reponse properties are taken
        !========================================================


        integer,intent(in)              :: unt
        integer,intent(inout)           :: Si, Sf
        logical,intent(inout)           :: derivatives
        character(len=*),intent(in)     :: dip_type
        real(8),dimension(:),intent(out):: Dip 
        real(8),dimension(:),intent(out):: DipD
        integer,intent(out),optional    :: error_flag

        !Local
        !Variables for read_fchk
        real(8),dimension(:),allocatable :: A
        integer,dimension(:),allocatable :: IA
        character(len=1)                 :: data_type
        integer                          :: N
        !Other local
        integer                          :: i,j,k, jj
        !FCHK specific (to be move to the new sr)
        integer :: Ntarget, Nes, Nat
        ! Local error flag
        integer :: error_local
        character(len=3) :: dummy_char

        ! Number of excited states computed
        call read_fchk(unt,"ETran scalars",data_type,N,A,IA,error_local)
        if (error_local == 0) then
            Nes = IA(1)
            Ntarget = IA(5)
            deallocate(IA)
        else
            if (present(error_flag)) error_flag = error_local
            return
        endif
        
        !Manage defaults and requests
        if (Si == -1) Si = 0
        if (Sf == -1) Sf = Ntarget
        if (Si /= 0) then
            write(dummy_char,'(I0)') Si
            call alert_msg("fatal","TD-DFT calcs in G09 only provide trdip from/to GS, but requested S="//dummy_char)
            if (present(error_flag)) error_flag=-1
            return
        endif
        if (Sf /= Ntarget) then
            call alert_msg("note","Retrieving trdip from a state different from the target. Derivatives not available.")
            derivatives=.false.
        endif
        
        ! Read ETran state values
        call read_fchk(unt,"ETran state values",data_type,N,A,IA,error_local)
        if (error_local /= 0) then
            if (present(error_flag)) error_flag = error_local
            print*, "ERROR: 'ETran state values' section not found"
            stop
        endif
        
        print'(2X,A,I0,A,I0,A)', "Transition dipole: <",Si,'|m|',Sf,'>' 
        if (adjustl(dip_type) == "eldip") then
            j=(Ntarget-1)*16 + 2
            Dip(1:3) = A(j:j+2)
        else if (adjustl(dip_type) == "magdip") then
            j=(Ntarget-1)*16 + 8
            Dip(1:3) = A(j:j+2)/(-2.d0)
        endif

        if (derivatives) then
            !Take Nat from allocated size of DipD
            Nat = size(DipD)/9
            if (N /= Nes*16+48+3*Nat*16) then
                call alert_msg("warning","Dipole derivatives requested but not found in fchk")
                return
            else
                print'(2X,A)', "Getting dipole derivatives"
            endif
           
            !Note the order is not the same as needed for FCclasses
            ! FCHK: dmux/dx1, dmuy/dx1, dmuz/dx1,  dmux/dy1, dmuy/dy1, dmuz/dy1,  ...
            ! FCC : dmux/dx1, dmux/dy1, dmux/dz1,  dmuy/dx1, dmuy/dy1, dmuy/dz1,  ...
            do j=1,3*Nat
                k = 16*Nes+48 + 16*(j-1)
                jj = j*3-2
                if (adjustl(dip_type) == "eldip") then
                    DipD(jj:jj+2) = A(k+2:k+4)
                else if (adjustl(dip_type) == "magdip") then
                    !Note that MagDip should be divided by -2 (see FCclasses manual)
                    DipD(jj:jj+2) = A(k+8:k+10)/(-2.d0)
                endif
        
                if (DipD(jj  ) == 0.d0 .and.&
                    DipD(jj+1) == 0.d0 .and.&
                    DipD(jj+2) == 0.d0) then
                    !Symmetry is used
                    call alert_msg("warning","Looks like the computation was done with symmetry. "//&
                                             "If so, dipole ders. are not reliable!")
                endif
            enddo
        endif
        deallocate(A)

        return

    end subroutine read_gaussfchk_dip

    subroutine read_gausslog_targestate(unt,S,error_flag)

        !=====================================================
        ! THIS CODE IS PART OF FCC_TOOLS
        !=====================================================
        ! Description
        !  Read transition electric or magnetic dipole moments
        !  from G09 log files (corresponding to a TDDFT calculation)
        !
        ! Notes
        !  Only the length-gauge reponse properties are taken
        !========================================================

        integer,intent(in)              :: unt
        integer,intent(inout)           :: S
        integer,intent(out),optional    :: error_flag

        !Local
        !Variables for reading
        character(len=240)               :: line=""
        character(len=3)                 :: auxchar
        !I/O
        integer :: IOstatus
        !Other local
        integer                          :: i,j,k, ii, jj

        ! Number of excited states computed
        ! Search section
        ii = 0
        if (present(error_flag)) error_flag = 0
        do
                ii = ii + 1
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    if (present(error_flag)) error_flag = -ii
                    rewind(unt)
                    return
                endif
                ! 2) Found a Excited State section
                if ( INDEX(line,"Excited State ") /= 0 ) then
                    call split_line(line,":",line,auxchar)
                    read(line,*) auxchar, auxchar, S
                endif
                ! 3) That was the target state
                if (INDEX(line,"This state for optimization and/or second-order correction")/=0) then
                    exit
                endif
        enddo

        return

    end subroutine read_gausslog_targestate

    subroutine read_gausslog_dip(unt,Si,Sf,dip_type,Dip,error_flag)

        !=====================================================
        ! THIS CODE IS PART OF FCC_TOOLS
        !=====================================================
        ! Description
        !  Read transition electric or magnetic dipole moments
        !  from G09 log files (corresponding to a TDDFT calculation)
        !
        ! Notes
        !  Only the length-gauge reponse properties are taken
        !========================================================

        integer,intent(in)              :: unt
        integer,intent(inout)           :: Si, Sf
        character(len=*),intent(in)     :: dip_type
        real(8),dimension(:),intent(out):: Dip 
        integer,intent(out),optional    :: error_flag

        !Local
        !Variables for reading
        character(len=240)               :: line=""
        integer                          :: N
        integer                          :: Sup, Sdw
        character(len=41)                :: section
        character(len=3)                 :: auxchar
        !I/O
        integer :: IOstatus
        !Other local
        integer                          :: i,j,k, ii, jj
        !FCHK specific (to be move to the new sr)
        integer :: Nes, Nat

        
        !Manage defaults and requests
        if (Si == -1) Si = 0
        if (Sf == -1) Sf = 1
        !To locate the data in file, we need the states ordered by value
        Sup = max(Si,Sf)
        Sdw = min(Si,Sf)
        
        if (dip_type == "eldip") then
            section="Ground to excited state transition electric dipole moments (Au):"
        else if (dip_type == "magdip") then
            section="Ground to excited state transition magnetic dipole moments (Au):"
        endif

        do
                ii = ii + 1
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    if (present(error_flag)) error_flag = -ii
                    rewind(unt)
                    return
                endif
                ! 2) Found what looked for!      
                if ( INDEX(line,section) /= 0) then
                    exit
                endif
        enddo

        !Reported values are:
        !1 (labels:    X Y Z Dip. S. Osc.)
        !2 State1 
        !3 State2
        ! ...
        !
        !Set line number to retrieve
        k = Sup+1

        !Read trdip
        print'(2X,A,I0,A,I0,A)', "Transition dipole: <",Si,'|m|',Sf,'>' 
        do i=1,k
            read(unt,'(A)') line
        enddo
        read(line,*) i, Dip(1:3)
        if (dip_type == "magdip") Dip(1:3) = Dip(1:3)/(-2.d0)

        return

    end subroutine read_gausslog_dip

    subroutine get_d2num(unt,iat,ixyz,istep,error_flag)

        integer,intent(in)           :: unt
        integer,intent(out)          :: iat,ixyz,istep
        integer,intent(out),optional :: error_flag

        !Local
        !Variables for reading
        character(len=240)               :: line=""
        character(len=3)                 :: auxchar
        !I/O
        integer :: IOstatus
        !Other local
        integer                          :: i,j,k, ii, jj

        ! Number of excited states computed
        ! Search section
        ii = 0
        if (present(error_flag)) error_flag = 0
        do
                ii = ii + 1
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    if (present(error_flag)) error_flag = -ii
                    return
                endif
                ! 3) Found
                if (INDEX(line,"D2Numr:")/=0) then
                    exit
                endif
        enddo

        !Split line to get data
        ! First remove trailing dot
        call split_line(line,".",line,auxchar)
        ! Get IAtom
        call split_line(line,"=",auxchar,line)
        read(line,*,iostat=IOstatus) iat
        if (IOstatus /= 0) then
            if (present(error_flag)) error_flag=1
            call alert_msg("warnign","Error reading derivatives")
            return   
        endif
        ! Get IXYZ
        call split_line(line,"=",auxchar,line)
        read(line,*) ixyz
        ! Get IStep
        call split_line(line,"=",auxchar,line)  
        read(line,*) istep

        return         

    end subroutine get_d2num

    subroutine read_gausslog_dipders(unt,Si,Sf,dip_type,dx,DipD,error_flag)

        !=====================================================
        ! THIS CODE IS PART OF FCC_TOOLS
        !=====================================================
        ! Description
        !  Compute dipole derivatives from a file containing the
        !  steps for numerical differenciation (3-points fit):
        !  EQ-BWD   EQ    EQ+FWD
        !  der = [mu(EQ+FWD)-mu(EQ+BWD)]/2*dx
        !  The step for num der is taken from sr arg "dx"
        !  Data are extracted with read_gausslog_dip, so that routine
        !  must not rewind
        !
        ! Notes
        !  Requirement of the input file:
        !  The input file should contain the backward and forward steps
        !  for numerical diferenciation for all Cartesian coordinates.
        !  The reader tip should be place so that the first element to
        !  get is the BWD step for coordinate 1 (NOT the equilibrium).
        !  This is achieved if the eq trdip is read before (in case it
        !  is contained in the file) 
        !
        !========================================================

        use constants

        integer,intent(in)              :: unt
        integer,intent(inout)           :: Si, Sf
        character(len=*),intent(in)     :: dip_type
        real(8),intent(inout)           :: dx
        real(8),dimension(:),intent(out):: DipD
        integer,intent(out),optional    :: error_flag

        !Local
        !Variables for reading
        real(8),dimension(1:3)           :: Dip_bwd, Dip_fwd
        integer                          :: Nat
        integer                          :: iat,ixyz,istep
        !Other local
        integer                          :: i,j,k, ii, jj
        !Local error flag
        integer                          :: error_local

        !Use detault for G09 (in au)
        if (dx == -1.d0) dx = 1.d-3/BOHRtoANGS

        !Take Nat from DipD allocation
        Nat = size(DipD)/9

        !Loop over all Cartesian
        do i=1,Nat
        do k=1,3
            j = 3*(3*(i-1)+(k-1))
            !Loop over bwd and fwd steps
            call read_gausslog_dip(unt,Si,Sf,dip_type,Dip_fwd,error_local)
            if (error_local /= 0) then
                if (present(error_flag)) error_flag = error_local
                call alert_msg("warning","Derivatives requested, but cannot be obtained. Was symmetry off?")
                return
            endif
            !Check D2Num
            call get_d2num(unt,iat,ixyz,istep,error_local)
            if (iat /= i .or. ixyz /= k .or. istep /= 1) then
                if (present(error_flag)) error_flag=1
                call alert_msg("warning","Derivatives requested, but cannot be obtained. Was symmetry off?")
                return
            endif

            call read_gausslog_dip(unt,Si,Sf,dip_type,Dip_bwd,error_local)
            if (error_local /= 0) then
                if (present(error_flag)) error_flag =  error_local
                call alert_msg("warning","Derivatives requested, but cannot be obtained. Was symmetry off?")
                return
            endif
            !Check D2Num
            call get_d2num(unt,iat,ixyz,istep,error_local)
            if (iat /= i .or. ixyz /= k .or. istep /= 2) then
                if (present(error_flag)) error_flag=1
                call alert_msg("warning","Derivatives requested, but cannot be obtained. Was symmetry off?")
                return
            endif

            if (error_local /= 0) then
                if (present(error_flag)) error_flag = error_local
                return
            endif

            DipD(j+1) = (Dip_fwd(1)-Dip_bwd(1))/2.d0/dx
            DipD(j+2) = (Dip_fwd(2)-Dip_bwd(2))/2.d0/dx
            DipD(j+3) = (Dip_fwd(3)-Dip_bwd(3))/2.d0/dx
        enddo
        enddo

        return

    end subroutine read_gausslog_dipders


end module gaussian_manage