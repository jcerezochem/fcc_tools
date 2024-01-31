module line_preprocess

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    !
    ! Despription
    ! -----------
    ! Subroutines to parse and manage lines
    !
    ! History
    ! -------
    ! 29/11/13: split_line_back: changed to store the whole thing
    !           in line_b in case there is no match (as opposed to
    !           to normal split_line SR
    ! 13/02/14: added change case subroutines
    ! 17/02/14  added "string2vector" subroutine
    ! 27/02/14  added "string2vector_char" (there also is "string2vector_int")
    !==============================================================


    implicit none

    contains

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine split_line(line,splitter,line_a,line_b)

        !Split a line from a given marker. If it is not present, it does not
        !split the line (the whole is preserved in line_a

        character(len=*),intent(in):: line,splitter
        character(len=*),intent(out):: line_a,line_b

        !local
        integer :: i,j
        !Auxiliar helps when line(input) is also one 
        !of the outputs, line_a or line_b
        character(len=(len(line_a))) :: aux_line_a

        i=INDEX(line,splitter)
        if ( i == 0 ) then
            line_a=line
            line_b=""
            return
        endif
        j=len_trim(splitter)
        
        aux_line_a=line(1:i-1)
        line_b=line(i+j:)
        line_a=aux_line_a

        return

    end subroutine split_line

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine split_line_back(line,splitter,line_a,line_b)

        !Split a line from a given marker. If it is not present, it does not
        !split the line (the whole is preserved in line_a >> change: for back
        ! this is stored in line_b (so that if no match, the whole thing is in line_b
        ! the important thing for this SR
        ! >> change again: the important thing is always what is in line_a (if there is 
        !                  no extension, we still want to know the filename

        character(len=*),intent(in):: line,splitter
        character(len=*),intent(out):: line_a,line_b

        !local
        integer :: i,j
        !Auxiliar helps when line(input) is also one 
        !of the outputs, line_a or line_b
        character(len=(len(line_a))) :: aux_line_a

        !INDEX with BACK=.true., search match from the end of the string (useful to get file extensions)
        i=INDEX(line,splitter,.true.)
        if ( i == 0 ) then
            ! We need to set line_a before overriding line_b, in case line_b and line are the same on input
            line_a=line
            line_b=""
            return
        endif
        j=len_trim(splitter)
        
        aux_line_a=line(1:i-1)
        line_b=line(i+j:)
        line_a=aux_line_a

        return

    end subroutine split_line_back

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine parse_line(line,narg,arg)

        !======================================================  
        ! Description
        !  Divides a line in the words that are separated by blank spaces
        ! Arguments
        !  line   INPUT   String to be evaluated
        !  narg   OUTPUT  Number of words
        !  arg(:) OUTPUT  String vector with all the words
        !======================================================  
    
        character(len=*),intent(in) :: line
        integer,intent(out) :: narg
        character(len=*),dimension(:),intent(out) :: arg
        !Local
        character(len=len_trim(line)) :: line_pp, args
        integer :: nchar, i, id

        line_pp = adjustl(line)
        nchar=len(line_pp)

        id=1
        args=""
        do i=1, nchar
            if ( line(i:i) == " " .and. args /= "" ) then
                arg(id)=adjustl(args)
                args=""
                id=id+1
                cycle
            elseif ( line(i:i) == " " .and. args == "" ) then
                cycle
            endif
            
            args=trim(args)//line(i:i)

        enddo
        arg(id)=adjustl(args)
        if (len_trim(args) == 0) then
            narg=0
        else
            narg=id
        endif
            
        return

    end subroutine parse_line


    subroutine read_list_int(record,nelements,int_array)

        !Read an array of integers of unknownn size stored in from a character string

        character(len=*),intent(inout) :: record
        integer,intent(out) :: nelements
        integer,dimension(:),intent(out) :: int_array

        !local
        integer :: imax, i

        nelements=0
        record=adjustl(record)
        imax=len_trim(record)
        !Eliminate extra blanks
        do i=1,imax
            if (len_trim(record(i:i)) == 0 .and. len_trim(record(i-1:i-1)) == 0) cycle
            if (len_trim(record(i:i)) == 0) nelements = nelements+1
        enddo
        nelements = nelements+1
        read(record,*) int_array(1:nelements)

        return

    end subroutine read_list_int

    subroutine set_upper_case(letter)

        character(len=1),intent(inout) :: letter

        if ( (ichar(letter) >= 97) .and. &
             (ichar(letter) <= 122) ) then
            letter = char(ichar(letter)-32)
        endif
        
        return

    end subroutine set_upper_case

    subroutine set_lower_case(letter)

        character(len=1),intent(inout) :: letter

        if ( (ichar(letter) >= 65) .and. &
             (ichar(letter) <= 90) ) then
            letter = char(ichar(letter)+32)
        endif
        
        return

    end subroutine set_lower_case

    subroutine set_word_lower_case(word)

        character(len=*),intent(inout) :: word
        !Local
        integer :: i

        do i=1,len_trim(word)
            call set_lower_case(word(i:i))
        enddo
        
        return

    end subroutine set_word_lower_case

    subroutine set_word_upper_case(word)

        character(len=*),intent(inout) :: word
        !Local
        integer :: i

        do i=1,len_trim(word)
            call set_upper_case(word(i:i))
        enddo
        
        return

    end subroutine set_word_upper_case


    subroutine string2rvector(raw_vector,array_vector,n_elem,sep)

        !Description
        ! Tranforms a string of comma sepparated values into an
        ! array of such real vaues 

        character(len=*),intent(in) :: raw_vector
        real(8),dimension(:),intent(out) :: array_vector
        integer,intent(out) :: n_elem
        character(len=*),intent(in) :: sep !separador

        !Local
        character(len=len_trim(raw_vector)) :: raw_vector_copy
        character(len=240) :: auxchar
        integer :: i
    
        
        !Copy the original vector to avoid modifying it
        raw_vector_copy = trim(adjustl(raw_vector))

        !Read unknown length vector
        i=0
        do 
            i=i+1
            if (len_trim(raw_vector_copy) == 0) then
                i=i-1
                exit
            else if ( INDEX(raw_vector_copy,sep) /= 0 ) then
                call split_line(raw_vector_copy,sep,auxchar,raw_vector_copy)
                read(auxchar,*) array_vector(i)
                ! By adjustl-ing every time we avoid double counting blank spaces
                raw_vector_copy = adjustl(raw_vector_copy)
            else 
                read(raw_vector_copy,*) array_vector(i)
                exit
            endif
        enddo  
        n_elem=i

        return

    end subroutine string2rvector

    subroutine string2ivector(raw_vector,array_vector,n_elem,sep)

        !Description
        ! Tranforms a string of comma sepparated values into an
        ! array of such integer vaues (integer version) 

        character(len=*),intent(in) :: raw_vector
        integer,dimension(:),intent(out) :: array_vector
        integer,intent(out) :: n_elem
        character(len=*),intent(in) :: sep !separador

        !Local
        character(len=len_trim(raw_vector)) :: raw_vector_copy
        character(len=240) :: auxchar
        integer :: i
    
        
        !Copy the original vector to avoid modifying it
        raw_vector_copy = trim(adjustl(raw_vector))

        !Read unknown length vector
        i=0
        do 
            i=i+1
            if (len_trim(raw_vector_copy) == 0) then
                i=i-1
                exit
            else if ( INDEX(raw_vector_copy,sep) /= 0 ) then
                call split_line(raw_vector_copy,sep,auxchar,raw_vector_copy)
                read(auxchar,*) array_vector(i)
                ! By adjustl-ing every time we avoid double counting blank spaces
                raw_vector_copy = adjustl(raw_vector_copy)
            else 
                read(raw_vector_copy,*) array_vector(i)
                exit
            endif
        enddo  
        n_elem=i

        return

    end subroutine string2ivector


    subroutine string2cvector(raw_vector,array_vector,n_elem,sep)

        !Description
        ! Tranforms a string of <sep> sepparated values into an
        ! array of such characters 
        ! This version also allows to use blank spaces as separator

        character(len=*),intent(in) :: raw_vector
        character(len=*),dimension(:),intent(out) :: array_vector
        integer,intent(out) :: n_elem
        character(len=*),intent(in) :: sep !separador

        !Local
        character(len=len_trim(raw_vector)) :: raw_vector_copy
        character(len=240) :: auxchar
        integer :: i
    
        
        !Copy the original vector to avoid modifying it
        raw_vector_copy = trim(adjustl(raw_vector))

        !Read unknown length vector
        i=0
        do 
            i=i+1
            if (len_trim(raw_vector_copy) == 0) then
                i=i-1
                exit
            else if ( INDEX(raw_vector_copy,sep) /= 0 ) then
                call split_line(raw_vector_copy,sep,auxchar,raw_vector_copy)
                ! We need to read with format 'A', not default 
                ! to catch the ","
                read(auxchar,'(A)') array_vector(i)
                ! By adjustl-ing every time we avoid double counting blank spaces
                raw_vector_copy = adjustl(raw_vector_copy)
            else 
                read(raw_vector_copy,'(A)') array_vector(i)
                exit
            endif
        enddo  
        n_elem=i

        return

    end subroutine string2cvector
    

    subroutine string2vector_getsize(raw_vector,n_elem,sep)

        !Description
        ! Get size of output vector (valid of _int, _char)

        character(len=*),intent(in) :: raw_vector
        integer,intent(out) :: n_elem
        character(len=*),intent(in) :: sep !separador

        !Local
        character(len=len_trim(raw_vector)) :: raw_vector_copy
        character(len=240) :: auxchar
        integer :: i
    
        
        !Copy the original vector to avoid modifying it
        raw_vector_copy = trim(adjustl(raw_vector))

        !Read unknown length vector
        i=0
        do 
            i=i+1
            if (len_trim(raw_vector_copy) == 0) then
                i=i-1
                exit
            else if ( INDEX(raw_vector_copy,sep) /= 0 ) then
                call split_line(raw_vector_copy,sep,auxchar,raw_vector_copy)
                ! We need to read with format 'A', not default 
                ! to catch the ","
!                 read(auxchar,'(A)') array_vector(i)
                ! By adjustl-ing every time we avoid double counting blank spaces
                raw_vector_copy = adjustl(raw_vector_copy)
            else 
!                 read(raw_vector_copy,'(A)') array_vector(i)
                exit
            endif
        enddo  
        n_elem=i

        return

    end subroutine string2vector_getsize


    ! Functions that get character from numbers

    function int20char(i,length) result(c)

        ! Converts an integer into a char of a given length
        ! If length>digits, complete with zeroes

        integer,intent(in)    :: i
        integer,intent(in)    :: length
        character(len=length) :: c
        ! Local
        character(len=10) :: fmt
        integer           :: ilength, j
        character(len=10) :: dummy_char

        !If length<number of digits, rise an error
        if (i == 0) then
            ilength = 1
        elseif (i<0) then
            ilength = int(log10(float(-i)))+2
        else
            ilength = int(log10(float(i)))+1
        endif

        if (ilength>length) then
            write(0,*) "Error in int20char: more digits in number than character size"
            stop
        endif
        ! otherwise fill head with zeroes
        dummy_char = ""
        do j=1,length-ilength
            dummy_char = trim(dummy_char)//"0"
        enddo

        ! Write format
        write(fmt,'(a,i0,a)') '(A,I',ilength,')'
        write(c,fmt) trim(dummy_char), i

        return

    end function int20char

    function int2char(i,length) result(c)

        integer,intent(in)    :: i
        integer,intent(in)    :: length
        character(len=length) :: c
        ! Local
        character(len=10) :: fmt
        integer           :: ilength

        !If length<number of digits, rise an error
        ilength = int(log10(float(i)))+1
        if (ilength>length) then
            write(0,*) "Error in int2char: more digits in number than character size"
            stop
        endif

        ! Write format 
        write(fmt,'(a,i0,a)') '(I',length,')'
        write(c,fmt) i

        return

    end function int2char

    function real2char(r,length,decimals) result(c)

        real(8),intent(in)    :: r
        integer,intent(in)    :: length, decimals
        character(len=length) :: c
        ! Local
        character(len=10) :: fmt

        ! Write format 
        write(fmt,'(a,i0,a,i0,a)') '(F',length,'.',decimals,')'
        write(c,fmt) r

        return

    end function real2char

    subroutine selection2intlist_nel(selection,Nlist)

       ! Interpret a selection to the array of integers
       ! Syntaxis of the selection:
       !  #  1 to 3   => (1,2,3)
       !  #  1,3,5-7  => (1,3,5,6,7)

        character(len=*), intent(in) :: selection
        integer, intent(out) :: Nlist
        !local 
        integer :: list
        character(len=5),dimension(:),allocatable :: selection_split
        integer :: i, j, jj
        integer :: N, range_last, range_width
        logical :: is_range
        character(len=len(selection)+10) :: selection_local

        ! Tranform "," into space 
        selection_local = ""
        j=0
        do i=1,len_trim(selection)
            j=j+1
            if (selection(i:i) == ",") then
                selection_local(j:j) = " "
            else if (selection(i:i) == "-") then
                selection_local(j:j+3) = " to "
                j=j+3
            else
                selection_local(j:j) = selection(i:i)
            endif
        enddo

        call string2vector_getsize(selection_local,N," ")
        allocate(selection_split(N))
        call string2cvector(selection_local,selection_split,N," ")

        is_range = .false.
        j = 0
        do i=1,N
            if (selection_split(i) == "to") then
                is_range =  .true.
                cycle
            endif
            ! Read number
            if (.not.is_range) then
                j = j+1
                read(selection_split(i),*) list
            else
                read(selection_split(i),*) range_last
                range_width = range_last - list
                do jj = 1, range_width
                    j = j + 1
                    list = list + 1
                enddo
                is_range = .false.
            endif
        enddo
        Nlist = j

        return

    end subroutine selection2intlist_nel

    subroutine selection2intlist(selection,list,Nlist)

       ! Interpret a selection to the array of integers
       ! Syntaxis of the selection:
       !  #  1 to 3   => (1,2,3)
       !  #  1,3,5-7  => (1,3,5,6,7)

        character(len=*), intent(in) :: selection
        integer, intent(out) :: Nlist
        integer,dimension(:) :: list
        !local 
        character(len=5),dimension(:),allocatable :: selection_split
        integer :: i, j, jj
        integer :: N, range_last, range_width
        logical :: is_range
        character(len=len(selection)+10) :: selection_local

        ! Tranform "," into space 
        selection_local = ""
        j=0
        do i=1,len_trim(selection)
            j=j+1
            if (selection(i:i) == ",") then
                selection_local(j:j) = " "
            else if (selection(i:i) == "-") then
                selection_local(j:j+3) = " to "
                j=j+3
            else
                selection_local(j:j) = selection(i:i)
            endif
        enddo

        call string2vector_getsize(selection_local,N," ")
        allocate(selection_split(N))
        call string2cvector(selection_local,selection_split,N," ")

        is_range = .false.
        j = 0
        do i=1,N
            if (selection_split(i) == "to") then
                is_range =  .true.
                cycle
            endif
            ! Read number
            if (.not.is_range) then
                j = j+1
                read(selection_split(i),*) list(j)
            else
                read(selection_split(i),*) range_last
                range_width = range_last - list(j)
                do jj = 1, range_width
                    j = j + 1
                    list(j) = list(j-1) + 1
                enddo
                is_range = .false.
            endif
        enddo
        Nlist = j

        return

    end subroutine selection2intlist
    
    subroutine selection2floatlist(selection,Nlist,list)

       ! Interpret a selection to the array of floats
       ! Requires Nlist as input if ranges are requested:
       ! Syntaxis of the selection:
       !  #  1.1 to 3.2   => (1.1,···,3.2)

        character(len=*), intent(in) :: selection
        integer, intent(in) :: Nlist
        real(8),dimension(:) :: list
        !local 
        character(len=5),dimension(:),allocatable :: selection_split
        integer :: i, j, jj
        integer :: N
        real(8) :: range_last, range_width, ranges_length
        logical :: is_range
        real(8),dimension(100) :: ranges ! can hold up to 50 ranges
        integer :: n_ranges
        character(len=len(selection)+10) :: selection_local
        
        ! Tranform "," into space 
        selection_local = ""
        j=0
        is_range = .false.
        do i=1,len_trim(selection)
            if (is_range) then
                is_range = .false.
                cycle
            endif
            j=j+1
            if (selection(i:i) == ",") then
                selection_local(j:j) = " "
            else if (selection(i:i+1) == "to") then
                selection_local(j:j+3) = " to "
                j=j+3
                is_range = .true.
            else if (selection(i:i) == "-") then
                selection_local(j:j+3) = " to "
                j=j+3 
            else
                selection_local(j:j) = selection(i:i)
            endif
        enddo

        call string2vector_getsize(selection_local,N," ")
        allocate(selection_split(N))
        call string2cvector(selection_local,selection_split,N," ")

        is_range = .false.
        j = 0
        n_ranges = 0
        ranges_length = 0.d0
        do i=1,N
            if (selection_split(i) == "to") then
                is_range =  .true.
                n_ranges = n_ranges + 1
                cycle
            endif
            ! Read number
            if (.not.is_range) then
                j = j+1
                read(selection_split(i),*) list(j)
            else
                read(selection_split(i),*) range_last
                jj = 2*n_ranges - 1
                ranges(jj)   = list(j)
                ranges(jj+1) = range_last
                ranges_length = ranges_length + ranges(jj+1) - ranges(jj)
                is_range = .false.
                j = j-1
            endif
        enddo
        if (n_ranges /= 0) then
            range_width = ranges_length / dfloat(Nlist - j - 1)
        endif
        if (n_ranges > 1) then
            print*, "ERROR: too many ranges in selection2floatlist"
            stop
        endif
        do i=1,n_ranges,2
            j = j + 1
            list(j) = ranges(i)
            N = max(1,int((ranges(i+1)-ranges(i))/range_width))
            do jj = 1,N
                j = j+1
                list(j) = min(ranges(i+1),ranges(i)+jj*range_width) 
            enddo
        enddo
        if (Nlist /= j) then
            print*, "ERROR: malformed range in selection2floatlist"
            stop
        endif

        return

    end subroutine selection2floatlist
    
    
    subroutine str_replace(line,pattern0,pattern1)

        !Repace one char by other (only works for one char for the moment)

        character(len=*),intent(inout):: line
        character(len=*),intent(in)   :: pattern0,pattern1

        !local
        integer :: n, i,j

        n = len(line)
        do i=1,n
            if (line(i:i) == pattern0 ) line(i:i) = pattern1
        enddo

        return

    end subroutine str_replace


end module line_preprocess
