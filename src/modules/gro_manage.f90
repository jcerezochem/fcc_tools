module gro_manage

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information
    !   from GROMACS files (gro, g96, Hessian mtx-dump) 
    !    
    ! Notes  
    !  All subroutines rewind the file after using it
    !==============================================================

    !Common declarations:
    !===================
    use constants
    use line_preprocess
    implicit none

    contains

    subroutine read_gro_geom(unt,Nat,AtName,X,Y,Z,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from gro files
        !
        !Description
        ! Get geometry and atom names from grom files. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        ! AtName  (out) char/vertor    Atom names
        ! X,Y,Z   (out) real/vectors   Coordinate vectors (ANGSTRONG)
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                   1 : Error
        ! 
        !==============================================================

        integer,intent(in)                          :: unt
        integer,intent(out)                         :: Nat
        character(len=*), dimension(:), intent(out) :: AtName
        real(kind=8), dimension(:), intent(out)     :: X,Y,Z
        integer,intent(out),optional                :: error_flag

        !Local variables
        !=============
        character(len=240) :: line=""
        character :: cnull
        !I/O
        integer :: IOstatus
        !Counters
        integer :: i, ii


        read(unt,'(A)') cull !title
        read(unt,*) Nat

        error_flag=0

        do i=1,Nat 

            read(unt,100,iostat=IOstatus) &
                          ii,        & !residue seq number
                          cnull,     & !residue name
                          AtName(i), & 
                          !serial
                          X(i),      &
                          Y(i),      &
                          Z(i)
            if (IOstatus /= 0) then
                error_flag=1
                rewind(unt)
                return
            endif
        enddo
        X(1:Nat)=X(1:Nat)*10.d0
        Y(1:Nat)=Y(1:Nat)*10.d0
        Z(1:Nat)=Z(1:Nat)*10.d0

        return
    !gro file format
100 format(i5,2a5,5X,3f8.3,3f8.4) !Serial is not read

    end subroutine read_gro_geom

    subroutine read_g96_geom(unt,Nat,AtName,X,Y,Z,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from gro files
        !
        !Description
        ! Get geometry and atom names from grom files. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        ! AtName  (out) char/vertor    Atom names
        ! X,Y,Z   (out) real/vectors   Coordinate vectors (ANGSTRONG)
        ! io_flag (io ) flag           Error flag:
        !                                   0 : Success
        !                                   1 : Success, but only geoms read (no AtNames)
        !                                  -1 : Error
        ! 
        !==============================================================

        integer,intent(in)  :: unt
        integer,intent(out) :: Nat
        character(len=*), dimension(:), intent(out) :: AtName
        real(kind=8), dimension(:), intent(out) :: X,Y,Z
        integer,intent(out),optional :: error_flag

        !Local variables
        !=============
        character(len=240) :: line=""
        character :: cnull
        !I/O
        integer :: IOstatus
        !Counters
        integer :: i, ii

        !The file is organized in sections. Each begining with a given name
        error_flag=-1
        do

            read(unt,*,iostat=ios) line
            if (ios /= 0) exit

            if (adjustl(line) == "POSITION" ) then
                error_flag=0
                i=0
                do 
                     read(unt,'(A)') line
                     if (adjustl(line) == "END" ) exit
                     i=i+1
                     read(line,*) ii,                    & !resseq
                                  cnull,                 & !resname
                                  AtName,                & 
                                  ii,                    & !serial
                                  X(i),                  &
                                  Y(i),                  &
                                  Z(i)
                enddo
                Nat = i
                !Transform to AA
                X(1:Nat)=X(1:Nat)*10.d0
                X(1:Nat)=X(1:Nat)*10.d0
                X(1:Nat)=X(1:Nat)*10.d0
            elseif (adjustl(line) == "POSITIONRED" ) then
            !this section only has info about coordinates (no AtNames!)
            write(0,*) "NOTE: No atomic info in g96!"
                error_flag=1
                i=0
                do 
                     read(unt,'(A)') line
                     if (adjustl(line) == "END" ) exit
                     i=i+1
                     read(line,*) X(i),      &
                                  Y(i),      &
                                  X(i)
                enddo
                Nat = i
                !Transform to AA
                X(1:Nat)=X(1:Nat)*10.d0
                X(1:Nat)=X(1:Nat)*10.d0
                X(1:Nat)=X(1:Nat)*10.d0
            else
                cycle
            endif

        enddo

        if (error_flag==-1) then
            write(0,*) "ERROR: No structure read in g96 file"
            stop
        endif

        return

    end subroutine read_g96_geom


    subroutine read_gmx_hess(unt,N,Hess,error)

        !Description
        ! Read Hessian from ascii file generated from an .mtx file (with gmxdump)
        ! 
        ! Error codes: 
        !  1: not square matrix

        integer,intent(in)::unt
        real(8),dimension(:,:),intent(out)::Hess
        integer,intent(out)::N
        integer,intent(out),optional :: error_flag
        
        !local
        character :: cnull
        !Counter
        integer::i, j
         
        read(unt,'(A)') cnull
        read(unt,*) N, i
        if (N /= i) then
            error = 1
            return
        endif

        !Read in triangular form
        j=1
        do i=1,N
            read(unt,*) Hess(i,1:N)
        enddo

        ! UNIT CONVERSION                                        ! GROMACS       --> Atomic Units  
        Hess(1:N,1:N)=Hess(1:N,1:N)/CALtoJ/HtoKCALM*BOHRtoNM**2  ! KJ/mol * nm-2 --> Hartree * bohr-2

        return

    end subroutine read_gmx_hess

end module gro_manage
