module fcc_manage
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

    subroutine read_fcc_natoms(unt,Nat,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Get geometry and atom names from xyz. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        !
        !==============================================================

        integer,intent(in)  :: unt
        integer,intent(out) :: Nat
        integer,intent(out),optional :: error_flag
        !Local
        integer :: IOstatus
        character(len=200) :: line
        character(len=50) :: section='GEOM', &
                             section_full
        
        do
            read(unt,'(A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) then
                call alert_msg("warning","Section '"//trim(adjustl(section))//"' not present in the FCC file.")
                if (present(error_flag)) error_flag=1
                rewind(unt)
                return
            endif
            ! 2) Found what looked for!      
            if ( INDEX(line,trim(adjustl(section))) /= 0 ) then
                read(line,'(A10)') section_full
                if (adjustl(section_full) == adjustl(section)) exit
            endif
        enddo

        read(unt,*) Nat

        return

     end subroutine read_fcc_natoms
     
    subroutine read_fcc_geom(unt,Nat,AtName,X,Y,Z,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Get geometry and atom names from xyz. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        !
        !==============================================================

        integer,intent(in)                        :: unt
        integer,intent(out)                       :: Nat
        real(kind=8),dimension(:),intent(out)     :: X,Y,Z
        character(len=*),dimension(:),intent(out) :: AtName
        integer,intent(out),optional              :: error_flag
        !Local
        integer :: i
        integer :: IOstatus
        character(len=200) :: line
        character(len=50) :: section='GEOM', &
                             section_full
        character :: cnull
        
        do
            read(unt,'(A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) then
                call alert_msg("warning","Section '"//trim(adjustl(section))//"' not present in the FCC file.")
                if (present(error_flag)) error_flag=1
                rewind(unt)
                return
            endif
            ! 2) Found what looked for!      
            if ( INDEX(line,trim(adjustl(section))) /= 0 ) then
                read(line,'(A10)') section_full
                if (adjustl(section_full) == adjustl(section)) exit
            endif
        enddo

        read(unt,*) Nat
        read(unt,*) cnull
        do i=1,Nat
            read(unt,*) AtName(i), X(i), Y(i), Z(i)
        enddo
            
        rewind(unt)

        return

     end subroutine read_fcc_geom
     
    subroutine read_fcc_grad(unt,Nat,Grad,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Get geometry and atom names from xyz. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        !
        !==============================================================

        integer,intent(in)                        :: unt
        integer,intent(in)                        :: Nat
        real(kind=8),dimension(:),intent(out)     :: Grad
        integer,intent(out),optional              :: error_flag
        !Local
        integer :: IOstatus
        character(len=200) :: line
        character(len=50) :: section='GRAD', &
                             section_full
        character :: cnull
        
        do
            read(unt,'(A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) then
                call alert_msg("warning","Section '"//trim(adjustl(section))//"' not present in the FCC file.")
                if (present(error_flag)) error_flag=1
                rewind(unt)
                return
            endif
            ! 2) Found what looked for!      
            if ( INDEX(line,trim(adjustl(section))) /= 0 ) then
                read(line,'(A10)') section_full
                if (adjustl(section_full) == adjustl(section)) exit
            endif
        enddo

        read(unt,*) Grad(1:3*Nat)
            
        rewind(unt)

        return

     end subroutine read_fcc_grad
     
    subroutine read_fcc_hess(unt,Nat,Hlt,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Get geometry and atom names from xyz. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        !
        !==============================================================

        integer,intent(in)                        :: unt
        integer,intent(in)                        :: Nat
        real(kind=8),dimension(:),intent(out)     :: Hlt
        integer,intent(out),optional              :: error_flag
        !Local
        integer :: IOstatus
        character(len=200) :: line
        character(len=50) :: section='HESS', &
                             section_full
        character :: cnull
        
        do
            read(unt,'(A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) then
                call alert_msg("warning","Section '"//trim(adjustl(section))//"' not present in the FCC file.")
                if (present(error_flag)) error_flag=1
                rewind(unt)
                return
            endif
            ! 2) Found what looked for!      
            if ( INDEX(line,trim(adjustl(section))) /= 0 ) then
                read(line,'(A10)') section_full
                if (adjustl(section_full) == adjustl(section)) exit
            endif
        enddo

        read(unt,*) Hlt(1:3*Nat*(3*Nat+1)/2)
            
        rewind(unt)

        return

     end subroutine read_fcc_hess
     
    subroutine read_state_geom(unt,Nat,GEO,Nvib,T,G)

        !====================================
        ! Description
        ! -----------
        ! Read geometry from the state files.
        !====================================

        integer,intent(in)                          :: unt
        integer,intent(in)                          :: Nat 
        real(8),dimension(:),intent(out)            :: GEO
        integer,intent(in),optional                 :: Nvib
        real(8),dimension(:,:),intent(out),optional :: T
        real(8),dimension(:),intent(out),optional   :: G
        !Local
        integer :: i, j, IOstatu

        read(unt,*,iostat=IOstatu) (GEO(i), i=1,3*Nat)
        if (IOstatu /= 0) then
            write(0,*) "ERROR: reading state file"
            stop
        endif

        if (present(T) .and. present(G)) then 
            read(unt,*) ((T(i,j), j=1,Nvib), i=1,3*Nat)
            read(unt,*) (G(i), i=1,Nvib)
        endif

        return
 
    end subroutine read_state_geom   

end module fcc_manage