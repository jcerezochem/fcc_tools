module psi4_manage

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information
    !   from Psi4 out files
    !    
    ! Notes  
    !  All subroutines rewind the file after using it
    !==============================================================

    !Common declarations:
    !===================
    use constants
    use line_preprocess
    use alerts
    implicit none

    contains

    subroutine read_psi4_natoms(unt,Nat,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from Psi4. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from Psi4. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                  -i : Read error on line i
        ! 
        !Note
        ! Need to understand better the Psi4 output
        !==============================================================

        integer,intent(in)  :: unt
        integer,intent(out) :: Nat
        integer,intent(out) :: error_flag

        !Local variables
        !=============
        character(len=240) :: line=""
        character :: cnull
        !I/O
        integer :: IOstatus
        !Counters
        integer :: i, ii
        !Units flag
        character(len=10) :: DIST_UNITS

        ! Search section
        error_flag = 0
        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            ! 2) Found what looked for!      
            if ( adjustl(line) == "==> Geometry <==" ) then
                exit
            endif
        enddo

        ! Overpass lines until reaching the target table
        do while (INDEX(line,"Center") == 0)
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( INDEX(line,"Geometry (in Angstrom)") /= 0 ) then
                DIST_UNITS="Angstrom"
            else if ( INDEX(line,"Geometry (in Bohr)") /= 0 ) then
                DIST_UNITS="Bohr"
            endif
        enddo
        !Read Table lines
        read(unt,'(X,A)',IOSTAT=IOstatus) line
        !Start reading geometry
        i=0
        read(unt,'(X,A)',IOSTAT=IOstatus) line
        do while (len_trim(line) /= 0)
            i=i+1
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) then
                print*, "Unexpected end of file while reading Geometry"
                error_flag = 1
                rewind(unt)
                return
            endif
        enddo
        Nat = i

        rewind(unt)
        return

    end subroutine read_psi4_natoms


    subroutine read_psi4_geom(unt,Nat,AtName,X,Y,Z,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Reads coordinates and atom names from Psi4. The coordinates
        ! are retuned as a 3Nat vector

        !Description
        ! Get geometry and atom names from Psi4. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        ! AtName  (out) char/vertor    Atom names
        ! X,Y,Z   (out) real/vectors   Coordinate vectors (ANGSTRONG)
        ! io_flag  (io ) flag          Error flag:
        !                                   0 : Success
        !                                  -i : Read error on line i
        ! 
        !Note
        ! Need to understand better the Psi4 output
        !==============================================================

        integer,intent(in)  :: unt
        integer,intent(out) :: Nat
        character(len=*), dimension(:), intent(out) :: AtName
        real(kind=8), dimension(:), intent(out) :: X,Y,Z
        integer,intent(out) :: error_flag

        !Local variables
        !=============
        character(len=240) :: line=""
        character :: cnull
        !I/O
        integer :: IOstatus
        !Counters
        integer :: i, ii
        !Units flag
        character(len=10) :: DIST_UNITS

        ! Search section
        error_flag = 0
        ii = 0
        do
            ii = ii + 1
            read(unt,'(A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) then
                error_flag = -ii
                rewind(unt)
                return
            endif
            ! 2) Found what looked for!      
            if ( adjustl(line) == "==> Geometry <==" ) then
                exit
            endif
        enddo

        ! Overpass lines until reaching the target table
        do while (INDEX(line,"Center") == 0)
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( INDEX(line,"Geometry (in Angstrom)") /= 0 ) then
                DIST_UNITS="Angstrom"
            else if ( INDEX(line,"Geometry (in Bohr)") /= 0 ) then
                DIST_UNITS="Bohr"
            endif
        enddo
        !Read Table lines
        read(unt,'(X,A)',IOSTAT=IOstatus) line
        !Start reading geometry
        i=0
        read(unt,'(X,A)',IOSTAT=IOstatus) line
        do while (len_trim(line) /= 0)
            i=i+1
            read(line,*) AtName(i), &
                         X(i),      &
                         Y(i),      &
                         Z(i)
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) then
                print*, "Unexpected end of file while reading Geometry"
                error_flag = 1
                rewind(unt)
                return
            endif
        enddo
        Nat = i
        
        !Manage lenght units
        if (adjustl(DIST_UNITS) == "Bohr") then
            X(1:Nat) = X(1:Nat) * BOHRtoANGS
            Y(1:Nat) = Y(1:Nat) * BOHRtoANGS
            Z(1:Nat) = Z(1:Nat) * BOHRtoANGS
        endif

        rewind(unt)
        return

    end subroutine read_psi4_geom


    subroutine read_psi4_hess(unt,Nat,Hlt,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Read Hessian from PSI4 output. Returs the triangular part of the
        ! Hessian matrix in AU
        ! 
        !Arguments
        ! unt   (inp) scalar   unit for the file
        ! Nat   (inp) scalar   Number of atoms
        ! Hlt   (out) vector   Lower triangular part of Hessian matrix (AU)
        ! error_flag (out) scalar  error_flag :
        !                                 0 : Success
        !                                -i : Read error on line i
        !                                 2 : Wrong number of elements for Hlt
        ! Notes
        !  Using variable format as suggested in:
        !  http://stackoverflow.com/questions/9881186/fortran-output-format-dependent-on-a-variable
        !  (Answer by eriktous)
        !==============================================================

        integer,intent(in) :: unt
        integer,intent(in) :: Nat
        real(kind=8), dimension(:), intent(out) :: Hlt
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        integer :: N
        !I/O
        integer :: IOstatus
        character(len=100) :: fmt
        !Counters
        integer :: i, j, k, ii, jini, jfin, &
                   iblock, nblocks, icols
        !Auxiliar arrays
        real(kind=8),dimension(3*Nat,3*Nat) :: Hpart
        
        
        !Use N to store 3*Nat
        N = 3*Nat

        ! Search section
        ii = 0
        error_flag = 0
        do
                ii = ii + 1
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                ! 2) Found what looked for!      
                if ( INDEX(line,"Force Constants in cartesian coordinates.") /= 0 ) then
                    exit
                endif
        enddo

        !Hessian elements arranged in blocks of 5 columns each
        nblocks = N/5
        if (N /= nblocks*5) nblocks=nblocks+1
        do iblock=1,nblocks
            !Three first lines useless
            read(unt,'(A)') line
            read(unt,'(A)') line
            read(unt,'(A)') line
            jini=(iblock-1)*5+1
            jfin=min(N,iblock*5)
            do i=1,N
                read(unt,*) k, Hpart(i,jini:jfin)
            enddo
        enddo

        !Get Hlt from the whole matrix in Hpart
        k = 0
        do i=1,N
            do j=1,i
                k = k + 1
                Hlt(k) = Hpart(j,i)
            enddo
        enddo
        if (k /= (N*(N+1)/2)) then
            error_flag = 2
        endif

        rewind(unt)
        return

    end subroutine read_psi4_hess

    subroutine read_psi4_dip(unt,Si,Sf,dip_type,Dip,error_flag)

        !=====================================================
        ! THIS CODE IS PART OF FCC_TOOLS
        !=====================================================
        ! Description
        !  Read transition electric or magnetic dipole moments
        !  from Psi4 out files (corresponding to a EOM-CC calculation)
        !  The information below "Length-Gauge Rotational Strength for" sections:
        !  which can be:
        !   * From GS to ES
        !     "Length-Gauge Rotational Strength for <ES> <symm>"
        !   * From ES1 to ES2
        !     "Length-Gauge Rotational Strength for <ES1> <symm> to <ES2> <symm>"
        !  Each section contains X Y Z components for:
        !   <p|mu_e|q>   
        !   <q|mu_m|p>   
        !   <p|mu_m|q>*  
        !   <q|mu_e|p>*  
        ! note: <> and <>* should be identical (Hermitian character),
        ! but they are not within CC truncation
        !
        ! Notes
        !  Only the length-gauge reponse properties are taken
        !  Tested fot CC code only
        !========================================================

        integer,intent(in)              :: unt
        integer,intent(inout)           :: Si, Sf
        character(len=*),intent(in)    :: dip_type
        real(8),dimension(:),intent(out):: Dip 
        integer,intent(out),optional    :: error_flag

        !Local
        !Variables for reading
        character(len=240)               :: line=""
        integer                          :: N
        integer                          :: Sup, Sdw
        character(len=41)                :: from_section, to_section
        character(len=3)                 :: auxchar
        !I/O
        integer :: IOstatus
        !Other local
        integer                          :: i,j,k, ii, jj
        !FCHK specific (to be move to the new sr)
        integer :: Ntarget, Nes, Nat

        ! Number of excited states computed
        ! Search section
        ii = 0
        error_flag = 0
        do
                ii = ii + 1
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                ! 2) Found what looked for!      
                if ( INDEX(line,"Number of States") /= 0 ) then
                    exit
                endif
        enddo
        call split_line(line,"=",auxchar,line)
        read(line,*) Nes
        !The first state should be the GS, so:
        Nes = Nes-1
        
        !Manage defaults and requests
        if (Si == -1) Si = 0
        if (Sf == -1) Sf = 1
        if (Si > Nes .or. Sf > Nes) then
            call alert_msg("warning","Requested state not computed in this Psi4 run.")
            error_flag=1
            return
        endif
        !To locate the data in file, we need the states ordered by value
        Sup = max(Si,Sf)
        Sdw = min(Si,Sf)
        
        !Get transition dipole (length gauge)
        ! Set the line to find
        if (Sdw /= 0) then
            write(auxchar,'(I3)') Sdw
        else
            write(auxchar,'(I3)') Sup
        endif
        from_section="Length-Gauge Rotational Strength for "//trim(adjustl(auxchar))
        if (Sdw /= 0) then
            write(auxchar,'(I3)') Sup
            to_section=" to "//trim(adjustl(auxchar))
        else
            to_section=""
        endif
        do
                ii = ii + 1
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    error_flag = -ii
                    rewind(unt)
                    return
                endif
                ! 2) Found what looked for!      
                if ( INDEX(line,trim(adjustl(from_section))) /= 0 .and. &
                     INDEX(line,trim(adjustl(to_section)))   /= 0) then
                    exit
                endif
        enddo

        !Reported values are:
        !1 (labels:    X Y Z). Read it now:
        !2 <p|mu_e|q>   <= eldip  for abs (note)
        !3 <q|mu_m|p>   <= magdip for emi (note)
        !4 <p|mu_m|q>*  <= magdip for abs (note)
        !5 <q|mu_e|p>*  <= eldip  for emi (note)
        ! note: <> and <>* should be identical (Hermitian character),
        ! but they are not within CC truncation
        !
        !Set line number to retrieve
        if      (Sf > Si .and. dip_type == "eldip") then
            k=2
        else if (Sf > Si .and. dip_type == "magdip") then
            k=4
        else if (Sf < Si .and. dip_type == "eldip") then
            k=5
        else if (Sf > Si .and. dip_type == "magdip") then
            k=3
        endif

        !Read trdip
        print'(2X,A,I0,A,I0,A)', "Transition dipole: <",Si,'|m|',Sf,'>' 
        do i=1,k
            read(unt,'(A)') line
        enddo
        read(line,*) auxchar, Dip(1:3)
        if (dip_type == "magdip") Dip(1:3) = Dip(1:3)/(-2.d0)

        return

    end subroutine read_psi4_dip


    subroutine read_psi4_dipders(unt,Si,Sf,dip_type,dx,DipD,error_flag)

        !=====================================================
        ! THIS CODE IS PART OF FCC_TOOLS
        !=====================================================
        ! Description
        !  Compute dipole derivatives from a file containing the
        !  steps for numerical differenciation (3-points fit):
        !  EQ-BWD   EQ    EQ+FWD
        !  der = [mu(EQ+FWD)-mu(EQ+BWD)]/2*dx
        !  The step for num der is taken from sr arg "dx"
        !  Data are extracted with read_psi4_dip, so that routine
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
        !Other local
        integer                          :: i,j,k, ii, jj

        !Use detault for Psi4 (in au)
        if (dx == -1.d0) dx = 5.d-3

        !Take Nat from DipD allocation
        Nat = size(DipD)/9

        !Loop over all Cartesian
        do i=1,3*Nat
            j = 3*(i-1)
            !Loop over bwd and fwd steps
            call read_psi4_dip(unt,Si,Sf,dip_type,Dip_bwd,error_flag)
            if (error_flag /= 0) then
                call alert_msg("warning","Derivatives requested, but cannot be computed from Psi4 output")
                return
            endif
            call read_psi4_dip(unt,Si,Sf,dip_type,Dip_fwd,error_flag)
            if (error_flag /= 0) then
                call alert_msg("warning","Derivatives requested, but cannot be computed from Psi4 output")
                return
            endif
            if (error_flag /= 0) return
            DipD(j+1) = (Dip_fwd(1)-Dip_bwd(1))/2.d0/dx
            DipD(j+2) = (Dip_fwd(2)-Dip_bwd(2))/2.d0/dx
            DipD(j+3) = (Dip_fwd(3)-Dip_bwd(3))/2.d0/dx
        enddo

        return

    end subroutine read_psi4_dipders


end module psi4_manage

