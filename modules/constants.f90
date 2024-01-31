!> Physical constantns and conversion factors
module constants

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        !Description
        ! This module contains mathematical and physical constants given
        ! in SI units. Some of them taken from http://physics.nist.gov/constants
        !==============================================================

        !CONSTANTS
        double precision, parameter :: &
                           PI      = 4.0d0*datan(1.0d0), & !
                           PIP     = PI**2,              & ! (used in FCclasses)
                           NAv     = 6.02214129D23,      & ! Avogadro number
                           AvNum   = 6.02214d23,         & ! (alias)
                           plank   = 6.62606957D-34,     & ! Planck constant
                           PH   = 6.62606957D-34,     & ! Planck constant
                           plankbar= 1.054571726D-34,    & ! Planck constant over 2PI
                           kboltz  = 1.3806488D-23,      & ! Boltzman constant
                           atmass  = 1.660538921D-27,    & ! Atomic mass
                           SL      = 2.99792458D8,       & ! Speed of light
                           cvel    = 2.99792458d8,       & ! (alias)
                           ! Non-SI
                           BoK     = 0.6950d0,           & ! Boltzman constant in cm-1
                           cvelau  = 137.036915742612      ! Speed of light in au 
                           

        !CONVERSION FACTORS
        double precision, parameter :: &
                           BOHRtoANGS= 5.2917721092D-1, & !(NIST)
                           AMUtoKG   = 1.66053873d-27,  &
                           AMUtoAU   = 1.82288839d3,    &
                           AUtoKG    = 9.10938291d-31,  & !(NIST)
                           BOHRtoM   = 5.2917721092D-11,& !(NIST)
                           BOHRtoNM  = 5.2917721092D-2, & !(NIST)
                           AMStoM    = 1.d-10,          & ! exact by definition
                           ANGStoM   = 1.d-10,          & ! exact by definition
                           HARTtoJ   = 4.35974434d-18,  & !(NIST)
                           HtoKCALM  = 627.5095d0,      &
                           CALtoJ    = 4.184,           &
                           HtoeV     = 27.21138505,     & !(NIST)
                           !FCclasses:
                           autoev    = 27.2113961d0,    &
                           autown    = 2.1947463068d5,  & ! Freq from AU to cm-1
                           evtown    = 8065.5446811132, & ! Energy from eV to cm-1
                           autofs    = 2.4189d-2,       &
                           fstoev    = 4.135667516,     &
                           autoang   = 0.5291771d0,     &
                           pmass     = 1.007825d0,      &
                           peratio   = 1836.1515d0,     &
                           autoamu   = 1.d0/(peratio/pmass),&
                           facabs    = 703.300d0,       &
                           facecd    = 20.5288d0,       &
                           facemi    = 4.d0/(3.d0*cvelau**3*autoev**3*autofs*1.d-6),&
                           factpa    = 8.35150d-4,      &
                           factpcd   = 4.87555d-5,      &
                           faccpl    = facemi*4.d0/cvelau, &
                           ! CFAC      = 4.D-22*PIP*SL*AMUtoKG/plank
                           cfac      = 1.d0/autown,     &
                           frot      = 471.443539822063d0, & ! ?? (in setHT)
                           cfacRR    = 3.42236d-47*1.d30,  &
                           facICns   = 1.d0,            &
                           facmcd    =-5.98442d-3




        !AtNum to element name conversion
        ! When initializing a character array as parameter, each element must have the declared len
        ! to have more flexibility, do not declare as parameter and use the data constructor
        character(len=2),dimension(103),parameter :: atname_from_atnum = &
        (/'H ',                                                                                'He',&
          'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne',&
          'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar',&
          'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',&
          'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe',&
          'Cs','Ba','La',& !Lantanides:  
!                  ---------------------------------------------------
                    'Ce','Pr','Nd','Pm','Sm','Eu','Gd',&
                    'Tb','Dy','Ho','Er','Tm','Yb','Lu',&
!                  ---------------------------------------------------
                         'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',&
         'Fr','Ra','Ac',& !Actinides:
!                  ---------------------------------------------------
                   'Th','Pa','U ','Np','Pu','Am','Cm',&
                   'Bk','Cf','Es','Fm','Md','No','Lr'&
!                  ---------------------------------------------------
         /)

        !AtNum to atom mass conversion
        real(8),dimension(118),parameter :: atmass_from_atnum = &
        !    Masses from QCElemental (NIST)
           (/  1.0078250  ,  & !H   -- 1  
               4.0026033  ,  & !He  -- 2  
               7.0160034  ,  & !Li  -- 3  
               9.0121831  ,  & !Be  -- 4  
              11.0093054  ,  & !B   -- 5  
              12.0000000  ,  & !C   -- 6  
              14.0030740  ,  & !N   -- 7  
              15.9949146  ,  & !O   -- 8  
              18.9984032  ,  & !F   -- 9  
              19.9924402  ,  & !Ne  -- 10 
              22.9897693  ,  & !Na  -- 11 
              23.9850417  ,  & !Mg  -- 12 
              26.9815385  ,  & !Al  -- 13 
              27.9769265  ,  & !Si  -- 14 
              30.9737620  ,  & !P   -- 15 
              31.9720712  ,  & !S   -- 16 
              34.9688527  ,  & !Cl  -- 17 
              39.9623831  ,  & !Ar  -- 18 
              38.9637065  ,  & !K   -- 19 
              39.9625909  ,  & !Ca  -- 20 
              44.9559083  ,  & !Sc  -- 21 
              47.9479420  ,  & !Ti  -- 22 
              50.9439570  ,  & !V   -- 23 
              51.9405062  ,  & !Cr  -- 24 
              54.9380439  ,  & !Mn  -- 25 
              55.9349363  ,  & !Fe  -- 26 
              58.9331943  ,  & !Co  -- 27 
              57.9353424  ,  & !Ni  -- 28 
              62.9295977  ,  & !Cu  -- 29 
              63.9291420  ,  & !Zn  -- 30 
              68.9255735  ,  & !Ga  -- 31 
              73.9211778  ,  & !Ge  -- 32 
              74.9215946  ,  & !As  -- 33 
              79.9165218  ,  & !Se  -- 34 
              78.9183376  ,  & !Br  -- 35 
              83.9114977  ,  & !Kr  -- 36 
              84.9117897  ,  & !Rb  -- 37 
              87.9056125  ,  & !Sr  -- 38 
              88.9058403  ,  & !Y   -- 39 
              89.9046977  ,  & !Zr  -- 40 
              92.9063730  ,  & !Nb  -- 41 
              97.9054048  ,  & !Mo  -- 42 
              97.9072124  ,  & !Tc  -- 43 
             101.9043441  ,  & !Ru  -- 44 
             102.9054980  ,  & !Rh  -- 45 
             105.9034804  ,  & !Pd  -- 46 
             106.9050916  ,  & !Ag  -- 47 
             113.9033651  ,  & !Cd  -- 48 
             114.9038788  ,  & !In  -- 49 
             119.9022016  ,  & !Sn  -- 50 
             120.9038120  ,  & !Sb  -- 51 
             129.9062227  ,  & !Te  -- 52 
             126.9044719  ,  & !I   -- 53 
             131.9041551  ,  & !Xe  -- 54 
             132.9054520  ,  & !Cs  -- 55 
             137.9052470  ,  & !Ba  -- 56 
             138.9063563  ,  & !La  -- 57 
             139.9054431  ,  & !Ce  -- 58 
             140.9076576  ,  & !Pr  -- 59 
             141.9077290  ,  & !Nd  -- 60 
             144.9127559  ,  & !Pm  -- 61 
             151.9197397  ,  & !Sm  -- 62 
             152.9212380  ,  & !Eu  -- 63 
             157.9241123  ,  & !Gd  -- 64 
             158.9253547  ,  & !Tb  -- 65 
             163.9291819  ,  & !Dy  -- 66 
             164.9303288  ,  & !Ho  -- 67 
             165.9302995  ,  & !Er  -- 68 
             168.9342179  ,  & !Tm  -- 69 
             173.9388664  ,  & !Yb  -- 70 
             174.9407752  ,  & !Lu  -- 71 
             179.9465570  ,  & !Hf  -- 72 
             180.9479958  ,  & !Ta  -- 73 
             183.9509309  ,  & !W   -- 74 
             186.9557501  ,  & !Re  -- 75 
             191.9614770  ,  & !Os  -- 76 
             192.9629216  ,  & !Ir  -- 77 
             194.9647917  ,  & !Pt  -- 78 
             196.9665688  ,  & !Au  -- 79 
             201.9706434  ,  & !Hg  -- 80 
             204.9744278  ,  & !Tl  -- 81 
             207.9766525  ,  & !Pb  -- 82 
             208.9803991  ,  & !Bi  -- 83 
             208.9824308  ,  & !Po  -- 84 
             209.9871479  ,  & !At  -- 85 
             222.0175782  ,  & !Rn  -- 86 
             223.0197360  ,  & !Fr  -- 87 
             226.0254103  ,  & !Ra  -- 88 
             227.0277523  ,  & !Ac  -- 89 
             232.0380558  ,  & !Th  -- 90 
             231.0358842  ,  & !Pa  -- 91 
             238.0507884  ,  & !U   -- 92 
             237.0481736  ,  & !Np  -- 93 
             244.0642053  ,  & !Pu  -- 94 
             243.0613813  ,  & !Am  -- 95 
             247.0703541  ,  & !Cm  -- 96 
             247.0703073  ,  & !Bk  -- 97 
             251.0795886  ,  & !Cf  -- 98 
             252.0829800  ,  & !Es  -- 99 
             257.0951061  ,  & !Fm  -- 100
             258.0984315  ,  & !Md  -- 101
             259.1010300  ,  & !No  -- 102
             266.1198300  ,  & !Lr  -- 103
             267.1217900  ,  & !Rf  -- 104
             268.1256700  ,  & !Db  -- 105
             271.1339300  ,  & !Sg  -- 106
             270.1333600  ,  & !Bh  -- 107
             269.1337500  ,  & !Hs  -- 108
             278.1563100  ,  & !Mt  -- 109
             281.1645100  ,  & !Ds  -- 110
             282.1691200  ,  & !Rg  -- 111
             285.1771200  ,  & !Cn  -- 112
             286.1822100  ,  & !Nh  -- 113
             289.1904200  ,  & !Fl  -- 114
             289.1936300  ,  & !Mc  -- 115
             293.2044900  ,  & !Lv  -- 116
             294.2104600  ,  & !Ts  -- 117
             294.0           & !Og  -- 118
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

        n = size(atmass_from_atnum)
        do i=1,n
            if (adjustl(name)==atname_from_atnum(i)) then
                mass = atmass_from_atnum(i)
                exit
            endif
        enddo

        return

    end function atmass_from_atname
    
    
    function atnum_from_atmass(mass) result(Zat)
        
        real(8),intent(in)           :: mass
        integer                      :: Zat
        !Local
        integer :: i, n

        Zat = 0
        n = size(atmass_from_atnum)
        do i=1,n
            if (abs(mass-atmass_from_atnum(i)) < 1.d-2) then
                Zat  = i
                exit
            endif
        enddo
        ! There is no warning here. Need to be handle from the caller (if detects Zat=0)
        
        return

    end function atnum_from_atmass


end module constants
