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
        character(len=5),dimension(103) :: atom_names_from_atnum
        !This should be elsewhere (constants_mod?)
        data atom_names_from_atnum(1:103) &
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


end module constants
