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
    use alerts
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

        use constants

        integer,intent(in)                        :: Nat
        character(len=*),dimension(:),intent(in)  :: atnames
        double precision,dimension(:),intent(out) :: mass
        integer,intent(out),optional              :: error_flag

        !Local
        integer :: i
        character(len=len(atnames(1))) :: element
        
        do i=1,Nat
            element = element_from_AtName(atnames(i))
            mass(i) = atmass_from_atname(element) 
        enddo

        return

    end subroutine assign_masses


    function element_from_AtName(AtomName) result(element)

        character(len=*),intent(in) :: AtomName
        character(len=2) :: element

        !local
        character(len=len(AtomName)) :: atname

        !1. Process atom name 
        atname = adjustl(AtomName)

        ! Sometimes, H atoms start with a number in PDB, GRO..
        if (atname(1:1) == "1" .or. &
            atname(1:1) == "2" .or. &
            atname(1:1) == "3" .or. &
            atname(1:1) == "4" .or. &
            atname(1:1) == "5" .or. &
            atname(1:1) == "6" .or. &
            atname(1:1) == "7" .or. &
            atname(1:1) == "8" .or. &
            atname(1:1) == "9") then
            if (atname(2:2) == "H" .or. &
                atname(2:2) == "h") then
                element="H"
                return
            else
                !Then we simply remove the number and go on
                atname(1:4) = atname(2:5)
                atname(5:5) = ""
            endif
        endif

        !Set first letter to upper case
        call set_upper_case(atname(1:1))

        !First solve conflicts with one-letter elements
        select case (atname(1:1))
           !==========
            case ("H")
           !==========
            !It can be H, He, Hf (not considered the lanthanide: Ho)
                ! We consider that:
                !  HE is hidrogen labeled as "E" (strange, though)
                !  He is helium
                !  HF is hidrogen labeled as "F" (strange, though)
                !  Hf is hafnium
                select case (atname(2:2))
                    case ("e")
                        element = "He"
                        call alert_msg("warning","He taken as helium")
                    case ("f")
                        element = "Hf"
                        call alert_msg("warning","He taken as hafnium")
                    case default
                        element = "H"
                        if (adjustl(element) /= adjustl(AtomName) ) &
                         call alert_msg("note",trim(adjustl(AtomName))//" taken as hydrogen")
                end select
                return
           !==========
            case ("B")
           !==========
            !It can be B, Be, Br, Ba
                select case (atname(2:2))
                    case ("a")
                        element = "Ba"
                    case ("A")
                        element = "Ba"
                        call alert_msg("note","BA taken as barium")
                    case ("e")
                        element = "Be"
                    case ("E")
                        element = "Be"
                        call alert_msg("note","BE taken as berium")
                    case ("r")
                        element = "Br"
                    case ("R")
                        element = "Br"
                        call alert_msg("note","BR taken as bromine")
                    case default
                        element = "B"
                        if (adjustl(element) /= adjustl(AtomName) ) &
                         call alert_msg("warning",trim(adjustl(AtomName))//" taken as borium")
                end select
                return
           !==========
            case ("C")
           !==========
                !C is a nightmare... It can be C Cl Cd Ca Cr Cs (not considered the lanthanide/actinides: Ce, Cm, Cf)
                ! We consider that:
                !  CD is carbon labeled as "D"
                !  Cd is Cadmium
                !  CL and Cl are chlorine (there is not usually an "L" label)
                !  CR and Cr are chromium (WARNING: chirality label?)
                !  CS and Cs are cesium (WARNING: chirality label?)
                !  CA is carbon labeled as "A" or or calcium: use more info later
                !  Ca is calcium
                select case (atname(2:2))
                    case ("d")
                        element = "Cd"
                        call alert_msg("warning","Cd taken as cadmium")
                    case ("r")
                        element = "Cr"
                        call alert_msg("warning","Cd taken as chromium")
                    case ("R")
                        element = "Cr"
                        call alert_msg("warning","CR taken as chromium")
                    case ("s")
                        element = "Cs"
                        call alert_msg("warning","Cs taken as cesium")
                    case ("S")
                        element = "Cs"
                        call alert_msg("warning","CS taken as cesium")
                    case ("l")
                        element = "Cl"
                        call alert_msg("warning","Cl taken as chlorine")
                    case ("L")
                        element = "Cl"
                        call alert_msg("warning","CL taken as chlorine")
                    case ("a")
                        ! If it has additional labels (e.g. Ca1), 
                        ! this is probably not Ca but Carbon
                        if (len_trim(atname) > 2) then
                            element = "C"
                            call alert_msg("warning",trim(atname)//" taken as carbon")
                        else
                            element = "Ca"
                            call alert_msg("warning",trim(atname)//" taken as calcium")
                        endif
                    case ("A")
                        !This case can be either C"A" or Ca. Mark with x to check later
                        element = "C"
                        call alert_msg("note","CA taken as carbone")
                    case default
                        element = "C"
                        if (adjustl(element) /= adjustl(AtomName) ) &
                         call alert_msg("note",trim(adjustl(AtomName))//" taken as carbone")
                end select
                return
           !==========
            case ("N")
           !==========
            !It can be N, Na, Ni, Nb (not considered the lanthanide/actinides: Nd, Np, No)
                ! We consider that:
                !  NB is carbon labeled as "B"
                !  Nb is niobium
                !  Ni and NI are nickel (there is not usually an "I" label)
                !  NA is nitrogen labeled as "A" or or sodium: use more info later
                !  Na is sodium
                select case (atname(2:2))
                    case ("b")
                        element = "Nb"
                    case ("i")
                        element = "Ni"
                    case ("I")
                        element = "Ni"
                    case ("a")
                        element = "Na"
                    case ("A")
                        ! If it has additional labels (e.g. Na1), 
                        ! this is probably not Na but sodium
                        if (len_trim(atname) > 2) then
                            element = "N"
                            call alert_msg("warning",trim(atname)//" taken as nitrogen")
                        else
                            element = "Na"
                            call alert_msg("warning",trim(atname)//" taken as sodium")
                        endif
                    case default
                        element = "N"
                end select
                return
           !==========
            case ("O")
           !==========
            !It can be O, Os
                ! We consider that:
                !  OS is carbon labeled as "S" (strange, although Os is more strange)
                !  Os is osmium
                select case (atname(2:2))
                    case ("s")
                        element = "Os"
                    case default
                        element = "O"
                end select
                return
           !==========
            case ("F")
           !==========
            !It can be F, Fe
                ! We consider that:
                !  Fe and FE are iron
                select case (atname(2:2))
                    case ("e")
                        element = "Fe"
                        call alert_msg("warning","Fe taken as iron")
                    case ("E")
                        element = "Fe"
                        call alert_msg("warning","FE taken as iron")
                    case default
                        element = "F"
                        if (adjustl(element) /= adjustl(AtomName) ) &
                         call alert_msg("note",trim(adjustl(AtomName))//" taken as fluorine")
                end select
                return
           !==========
            case ("P")
           !==========
            !It can be P, Pb, Po
                ! We consider that:
                !  Pb and PB are lead
                !  Po is polonium
                !  PO is P labeled "O"
                select case (atname(2:2))
                    case ("o")
                        element = "Po"
                    case ("O")
                        element = "Po"
                    case ("t")
                        element = "Pt"
                    case ("T")
                        element = "Pt"
                    case default
                        element = "P"
                end select
                return
           !==========
            case ("S")
           !==========
            !It can be S, Sr, Se, Sn, Si
                ! We consider that:
                !  Sb is antimonium 
                !  SB sulfur labeled as "B"
                select case (atname(2:2))
                    case ("i")
                        element = "Si"
                        call alert_msg("warning","Si taken as silicon")
                    case ("I")
                        element = "Si"
                        call alert_msg("warning","SI taken as silicon")
                    case ("r")
                        element = "Sr"
                        call alert_msg("warning","Sr taken as strontium")
                    case ("R")
                        element = "Sr"
                        call alert_msg("warning","SR taken as strontium")
                    case ("n")
                        element = "Sn"
                        call alert_msg("warning","Sn taken as tin (Sn)")
                    case ("N")
                        element = "Sn"
                        call alert_msg("warning","SN taken as tin (Sn)")
                    case ("b")
                        element = "Sb"
                        call alert_msg("warning","Sb taken as antimony")
                    case default
                        element = "S"
                        if (adjustl(element) /= adjustl(AtomName) ) &
                         call alert_msg("note",trim(adjustl(AtomName))//" taken as sulfur")
                end select
                return
           !==========
            case ("K")
           !==========
            !It can be K, Kr
                select case (atname(2:2))
                    case ("r")
                        element = "Kr"
                    case ("R")
                        element = "Kr"
                    case default
                        element = "K"
                end select
                return
           !==========
            case ("V")
           !==========
            !It can only be V
                element = "V"
                return
           !==========
            case ("W")
           !==========
            !It can only be W
                element = "W"
                return
           !==========
            case ("Y")
           !==========
            !It can only be Y
                element = "Y"
                return
           !==========
            case ("U")
           !==========
            !It can only be U
                element = "U"
                return
           !==========
            case ("I")
           !==========
            !It can be I, Ir
                ! We consider that:
                !  Ir and IR are iridium
                select case (atname(2:2))
                    case ("r")
                        element = "Ir"
                    case ("R")
                        element = "Ir"
                    case default
                        element = "I"
                end select
                return
        end select

        !Once one-letter conflicts are solved, the rest are trivial
        call set_lower_case(atname(2:2))
        element = atname(1:2)

        return

    end function element_from_AtName


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
    
    subroutine prepare_fcc3input(unt,Nat,Nvib,Mass,DE,T,error_flag)

    
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
        
        write(unt,'(A)') '$$$'
        write(unt,'(A)') '; Sample fcclasses3 input. For a more detailed input template use fcclasses3 -h'
        write(unt,'(A12,X,A1,X,A)') 'PROPERTY','=', 'OPA'
        write(unt,'(A12,X,A1,X,A)') 'DE','=', '(data required!)'
        write(unt,'(A12,X,A1,X,A)') 'TEMP','=','0.0'
        write(unt,'(A12,X,A1,X,A)') 'BROADFUN','=','GAU'
        write(unt,'(A12,X,A1,X,A)') 'HWHM','=','0.01'
        write(unt,'(A12,X,A1,X,A)') 'METHOD','=','TD'
        write(unt,'(A12,X,A1,X,A)') 'STATE1_FILE','=','(data required!)'
        write(unt,'(A12,X,A1,X,A)') 'STATE2_FILE','=','(data required!)'
        write(unt,'(A12,X,A1,X,A)') 'ELDIP_FILE','=','(data required!)'
        write(unt,'(A12,X,A1,X,A)') 'MAGDIP_FILE','=','(data required!)'
        
        close(unt)
    
        return
    
    end subroutine prepare_fcc3input


end module fcc_basics

