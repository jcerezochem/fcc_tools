module constants

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        !Description
        ! This module contains mathematical and physical constants given
        ! in SI units. Taken from http://physics.nist.gov/constants
        !==============================================================

        !CONSTANTS
!#ifdef DOUBLE
        double precision, parameter :: &
!#else
!        real, parameter :: &
!#endif
                           PI      = 4.0d0*datan(1.0d0),   &
                           clight  = 2.99792458D8,    &
                           plank   = 6.62606957D-34,  &
                           plankbar= 1.054571726D-34, &
                           boltz   = 1.3806488D-23,   &
                           NAv     = 6.02214129D23,   &
                           atmass  = 1.660538921D-27  


        !CONVERSION FACTORS
!#ifdef DOUBLE
        double precision, parameter :: &
!#else
!        real, parameter :: &
!#endif
                           BOHRtoAMS = 5.2917720859D-1, &
                           BOHRtoANGS= 5.2917720859D-1, &
                           UMAtoKG   = 1.66053873d-27,  &
                           UMAtoAU   = 1.82288839d3,    &
                           AUtoKG    = 9.10938291d-31,  &
                           BOHRtoM   = 5.291772083d-11, &
                           BOHRtoNM  = 5.291772083d-2,  &
                           AMStoM    = 1.d-10,          &
                           ANGStoM   = 1.d-10,          &
                           HARTtoJ   = 4.3597482d-18,   &
                           HtoKCALM  = 627.5095d0,      &
                           CALtoJ    = 4.184,           &
                           HtoeV     = 27.2114,         &
                           autown    = 2.1947463068d5    !From FCclasses Freq from AU to cm-1



        !AtNum to element name conversion
        character(len=5),dimension(103) :: atname_from_atnum
        !This should be elsewhere (constants_mod?)
        data atname_from_atnum(1:103) &
         /'H' ,                                                                                'He',&
          'Li','Be',                                                  'B' ,'C' ,'N' ,'O' ,'F' ,'Ne',&
          'Na','Mg',                                                  'Al','Si','P' ,'S' ,'Cl','Ar',&
          'K' ,'Ca','Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',&
          'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe',&
          'Cs','Ba','La',& !Lantanides:  
!                  ---------------------------------------------------
                    'Ce','Pr','Nd','Pm','Sm','Eu','Gd',&
                    'Tb','Dy','Ho','Er','Tm','Yb','Lu',&
!                  ---------------------------------------------------
                         'Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',&
         'Fr','Ra','Ac',& !Actinides:
!                  ---------------------------------------------------
                   'Th','Pa','U' ,'Np','Pu','Am','Cm',&
                   'Bk','Cf','Es','Fm','Md','No','Lr'&
!                  ---------------------------------------------------
         /

        !AtNum to atom mass conversion
        real(8),dimension(36),parameter :: atmass_from_atnum = &
           (/1.0078250,   & !H
             4.0026033,   & !He
             7.0160045,   & !Li
             9.0121825,   & !Be
             11.0093053,  & !B
             12.000000,   & !C
             14.0030740,  & !N
             15.9949146,  & !O
             18.9984033,  & !F
             19.9924391,  & !Ne  --10
             22.9897697,  & !Na
             23.9850450,  & !Mg
             26.9815413,  & !Al
             27.9769284,  & !Si
             30.9737634,  & !P
             31.9720718,  & !S
             34.9688527,  & !Cl
             39.9623831,  & !Ar
             38.9637079,  & !K
             39.9625907,  & !Ca  --20
             44.9559136,  & !Sc
             47.9479467,  & !Ti
             50.9439625,  & !V
             51.9405097,  & !Cr
             54.9380463,  & !Mn
             55.9349393,  & !Fe
             58.9331978,  & !Co
             57.9353471,  & !Ni
             62.9295992,  & !Cu
             63.9291454,  & !Zn  --30
             68.9255809,  & !Ga
             73.9211788,  & !Ge
             74.9215955,  & !As
             79.9165205,  & !Se
             78.9183361,  & !Br
             83.9115064   & !Kr
            /)

    contains

    !Helper functions to access mass, AtNum and names

    function atnum_from_atname(name) result(Zat)

        character(len=*),intent(in) :: name
        integer                     :: Zat 
        !local  
        integer                     :: i, n

        n = size(atname_from_atnum)
        do i=1,n
            if (adjustl(name)==atname_from_atnum(i)) then
                Zat = i
                exit
            endif
        enddo

        return

    end function atnum_from_atname


    function atmass_from_atname(name) result(mass)

        character(len=*),intent(in) :: name
        real(8)                     :: mass 
        !local  
        integer                     :: i, n

        mass=0.d0
        n = size(atmass_from_atnum)
        do i=1,n
            if (adjustl(name)==atname_from_atnum(i)) then
                mass = atmass_from_atnum(i)
                exit
            endif
        enddo
        if (mass == 0.d0) then
            print*, "ERROR: mass cannot be determined for: "//name
            stop
        endif

        return

    end function atmass_from_atname

    subroutine atominfo_from_atmass(mass,Zat,name) 
        !================================================
        ! Description
        ! From mass value, get atom%name and atom%Zat
        ! So we can still use old FCclasses input format
        ! We need atomnames to use guess_connect
        !================================================

        real(8),intent(in)           :: mass
        integer,intent(out)          :: Zat
        character(len=*),intent(out) :: name
        !Local
        integer :: i, n

        name = "XX"
        n = size(atmass_from_atnum)
        do i=1,n
            if (abs(mass-atmass_from_atnum(i)) < 1.d-2) then
                name = atname_from_atnum(i)
                Zat  = i
            endif
        enddo
        if (adjustl(name) == "XX") then
            write(0,'(A,F12.6)') "ERROR: atom info could not be set from mass: ",mass
            stop
        endif

        return

    end subroutine atominfo_from_atmass

end module constants
