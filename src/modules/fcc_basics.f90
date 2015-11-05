module fcc_basics

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines with general tasks
    !    
    ! Notes  
    !  All subroutines rewind the file after using it
    !==============================================================

    !Common declarations:
    !===================
    use line_preprocess
    implicit none

    contains


    subroutine assign_masses(Nat,atnames,mass,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Assign masses based on atom names (elements). Masses taken
        ! from G09 outputs
        !
        !Arguments
        ! Nat     (inp)  scalar   Number of atoms
        ! atnames (int)  vector   Element names
        ! mass    (out)  vector   Masses (AMU)
        ! error_flag (out) flag   0: Success
        !                         1: Some atoms names do not coincide
        !                            with an element. Set to zero
        !
        !Notes
        !==============================================================

        integer,intent(in)                        :: Nat
        character(len=*),dimension(:),intent(in)  :: atnames
        double precision,dimension(:),intent(out) :: mass
        integer,intent(out),optional              :: error_flag

        integer :: i
        character(len=2) :: atname

        error_flag = 0
        do i=1,Nat
            atname = adjustl(atnames(i))
            select case (atname)
               case ("H")
                mass(i)= 1.0078250
               case ("He")
                mass(i)= 4.0026033
               case ("Li")
                mass(i)= 7.0160045
               case ("Be")
                mass(i)= 9.0121825
               case ("B")
                mass(i)=11.0093053
               case ("C")
                mass(i)=12.000000
               case ("N")
                mass(i)=14.0030740
               case ("O")
                mass(i)=15.9949146
               case ("F")
                mass(i)=18.9984033
               case ("Ne")
                mass(i)=19.9924391
               case ("Na")
                mass(i)=22.9897697
               case ("Mg")
                mass(i)=23.9850450
               case ("Al")
                mass(i)=26.9815413
               case ("Si")
                mass(i)=27.9769284
               case ("P")
                mass(i)=30.9737634
               case ("S")
                mass(i)=31.9720718
               case ("Cl")
                mass(i)=34.9688527
               case ("Ar")
                mass(i)=39.9623831
               case ("K")
                mass(i)=38.9637079
               case ("Ca")
                mass(i)=39.9625907
               case ("Sc")
                mass(i)=44.9559136
               case ("Ti")
                mass(i)=47.9479467
               case ("V")
                mass(i)=50.9439625
               case ("Cr")
                mass(i)=51.9405097
               case ("Mn")
                mass(i)=54.9380463
               case ("Fe")
                mass(i)=55.9349393
               case ("Co")
                mass(i)=58.9331978
               case ("Ni")
                mass(i)=57.9353471
               case ("Cu")
                mass(i)=62.9295992
               case ("Zn")
                mass(i)=63.9291454
               case ("Ga")
                mass(i)=68.9255809
               case ("Ge")
                mass(i)=73.9211788
               case ("As")
                mass(i)=74.9215955
               case ("Se")
                mass(i)=79.9165205
               case ("Br")
                mass(i)=78.9183361
               case ("Kr")
                mass(i)=83.9115064
               !Desordenados
               case ("Pd")
                mass(i)=105.9032000
               case ("Pt")
                mass(i)=194.9648000
               case ("I")
                mass(i)=126.9004000
               !Default
               case default
                error_flag = 1
                mass(i)=0.00
            end select
        enddo

        return

    end subroutine assign_masses


    subroutine prepare_fccinput(unt,Nat,Nvib,Mass,DE,T,error_flag)

    
        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Prepare input file for a FCclasses run
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! Nat     (int)  int /scalar   Number of atoms
        ! Nvib    (int)  int /scalar   Number of vibrational DFs
        ! Mass    (out)  real/vector   Atomic masses (AMU)
        ! DE      (inp)  real/scalar   Adiabatic/Vertical Energy
        ! T       (inp)  real/scalar   Temperature (K)
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        integer,intent(in)              :: unt
        integer,intent(in)              :: Nat, Nvib
        real(8),dimension(:),intent(in) :: Mass
        real(8),intent(in)              :: DE, T
        integer,intent(out),optional    :: error_flag

        !Local
        integer                          :: i
        
        error_flag = 0

        write(unt,*) Nat
        write(unt,*) Nvib
        do i=1,Nat
            write(unt,*) Mass(i)
        enddo
        if ( DE < 0.d0 ) then
            write(unt,*) " <ENERGY> Adiabatic/Vertical Energy (eV)"
        else    
            write(unt,'(F12.6,X,A)') DE, "Adiabatic/Vertical Energy (eV)"
        endif
        write(unt,*) "'abs' 'ECDNO ' 'FC'"
        if ( T < 0.d0 ) then
            write(unt,*) "0.1d0 8.d-1 ! Temp(K) / BoltzThr (for TI)"
        else
            write(unt,'(F8.3,X,A)') T, "8.d-1 ! Temp(K) / BoltzThr (for TI)"
        endif
        write(unt,*) "'D ' 2048 1000.d0 0   Algorithm"
        write(unt,*) "'state_file_1'"
        write(unt,*) "'state_file_2' 'AH'  -1.d5  1.d5"
        write(unt,*) "'eldip_file'"
        write(unt,*) "'magdip_file'"
        write(unt,*) "1 rotation to overlap S1 and S2 (1=yes)"
        write(unt,*) "25"
        write(unt,*) "20"
        write(unt,*) "1.d8"
        write(unt,*) "'Gau' 1.50d0  4.00d0 1001 0.01d0"
        write(unt,*) "'NO' 'SELECT' 0.d0"
        write(unt,*) ""
    
        close(unt)
    
        return
    
    end subroutine prepare_fccinput


end module fcc_basics

