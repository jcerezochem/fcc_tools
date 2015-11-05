module fcc_io

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to manage output files 
    !  from different QM codes. Relies in specific modules for
    !  each package.
    !    
    !==============================================================

    !Common declarations:
    !===================
    use constants
    use line_preprocess
    use fcc_basics
    use gaussian_manage
    use gamess_manage
    use psi4_manage
    use molcas_manage
    use molpro_manage
    use gmx_manage
    implicit none

    contains

    subroutine generic_natoms_reader(unt,filetype,Nat,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic natoms reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! Nat     (out)  int /scalar   Number of atoms
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        integer,intent(in)           :: unt
        character(len=*),intent(in)  :: filetype
        integer,intent(out)          :: Nat
        integer,intent(out),optional :: error_flag

        !Local
        !Variables for read_fchk
        real(8),dimension(:),allocatable :: A
        integer,dimension(:),allocatable :: IA
        character(len=1)                 :: data_type
        integer                          :: N

        error_flag = 0
        select case (adjustl(filetype))
            case("log")
             call read_gausslog_natoms(unt,Nat,error_flag)
            case("fchk")
             call read_fchk(unt,'Number of atoms',data_type,N,A,IA,error_flag)
             Nat = IA(1)
             deallocate(IA)
            case("gms")
             call read_gamess_natoms(unt,Nat,error_flag)
            case("psi4")
             call read_psi4_natoms(unt,Nat,error_flag)
            case("molcas")
             call read_molcasUnSym_natoms(unt,Nat,error_flag)
            case("molpro")
             call read_molpro_natoms(unt,Nat,error_flag)
            case("g96")
             call read_g96_natoms(unt,Nat)
            case default
             write(0,*) "Unsupported filetype:"//trim(adjustl(filetype))
             call supported_filetype_list('freq')
             error_flag = 99
         end select

         return

    end subroutine generic_natoms_reader


    subroutine generic_structure_reader(unt,filetype,Nat,X,Y,Z,Mass,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic geometry reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! Nat     (io )  int /scalar   Number of atoms
        ! X,Y,Z   (out)  real/vectors  Coordinates
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        integer,intent(in)              :: unt
        character(len=*),intent(in)     :: filetype
        integer,intent(inout)           :: Nat
        real(8),dimension(:),intent(out):: X,Y,Z
        real(8),dimension(:),intent(out):: Mass
        integer,intent(out),optional    :: error_flag

        !Local
        !Variables for read_fchk
        real(8),dimension(:),allocatable :: A
        integer,dimension(:),allocatable :: IA
        character(len=1)                 :: data_type
        integer                          :: N
        !Other local
        character(len=2),dimension(Nat)  :: AtName
        integer                          :: i,j

        error_flag = 0
        select case (adjustl(filetype))
            case("log")
             call read_gauslog_geom(unt,Nat,AtName,X,Y,Z,error_flag)
             call assign_masses(Nat,AtName,Mass,error_flag)
            case("fchk")
             call read_fchk(unt,'Current cartesian coordinates',data_type,N,A,IA,error_flag)
             do i=1,N,3
                 j = (i-1)/3+1
                 X(j) = A(i)  *BOHRtoANGS
                 Y(j) = A(i+1)*BOHRtoANGS
                 Z(j) = A(i+2)*BOHRtoANGS
             enddo
             deallocate(A)
             call read_fchk(unt,'Real atomic weights',data_type,N,A,IA,error_flag)
             do i=1,N
                 Mass(i) = A(i)
             enddo
             deallocate(A)
            case("gms")
             call read_gamess_geom(unt,Nat,AtName,X,Y,Z,error_flag)
             call assign_masses(Nat,AtName,Mass,error_flag)
            case("psi4")
             call read_psi4_geom(unt,Nat,AtName,X,Y,Z,error_flag)
             call assign_masses(Nat,AtName,Mass,error_flag)
            case("molcas")
             call read_molcasUnSym_geom(unt,Nat,AtName,X,Y,Z,error_flag)
             call assign_masses(Nat,AtName,Mass,error_flag)
            case("molpro")
             call read_molpro_geom(unt,Nat,AtName,X,Y,Z,error_flag)
             call assign_masses(Nat,AtName,Mass,error_flag)
            case("g96")
             call read_g96_geom(unt,Nat,AtName,X,Y,Z)
             call assign_masses(Nat,AtName,Mass,error_flag)
            case default
             write(0,*) "Unsupported filetype:"//trim(adjustl(filetype))
             call supported_filetype_list('freq')
             error_flag = 99
         end select

         return

    end subroutine generic_structure_reader


    subroutine generic_gradient_reader(unt,filetype,Nat,Grad,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic Hessian reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! Nat     (int)  int /scalar   Number of atoms
        ! Grad    (out)  real/vector   Gradient (AU)
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        integer,intent(in)              :: unt
        character(len=*),intent(in)     :: filetype
        integer,intent(in)              :: Nat
        real(8),dimension(:),intent(out):: Grad
        integer,intent(out),optional    :: error_flag

        !Local
        !Variables for read_fchk
        real(8),dimension(:),allocatable :: A
        integer,dimension(:),allocatable :: IA
        character(len=1)                 :: data_type
        integer                          :: N
        !For gausslog read
        character(12*size(Grad))         :: section
        !Other auxiliar
        integer                          :: i

        error_flag = 0
        select case (adjustl(filetype))
            case("log")
             call summary_parser(unt,7,section,error_flag)
             read(section,*) Grad(1:3*Nat)
            case("fchk")
             call read_fchk(unt,'Cartesian Gradient',data_type,N,A,IA,error_flag)
             if (error_flag /= 0) return
             do i=1,N
                 Grad(i) = A(i)
             enddo
             deallocate(A)
            case default
             write(0,*) "Unsupported filetype:"//trim(adjustl(filetype))
             call supported_filetype_list('grad')
             error_flag = 99
         end select

         return

    end subroutine generic_gradient_reader

    subroutine generic_Hessian_reader(unt,filetype,Nat,Hlt,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic Hessian reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! Nat     (int)  int /scalar   Number of atoms
        ! Hlt     (out)  real/vector   Hessian (lower triangular form) (AU)
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        integer,intent(in)              :: unt
        character(len=*),intent(in)     :: filetype
        integer,intent(in)              :: Nat
        real(8),dimension(:),intent(out):: Hlt
        integer,intent(out),optional    :: error_flag

        !Local
        !Variables for read_fchk
        real(8),dimension(:),allocatable :: A
        integer,dimension(:),allocatable :: IA
        character(len=1)                 :: data_type
        integer                          :: N
        !For gausslog read
        character(12*size(Hlt))          :: section
        !Other auxiliar
        integer                          :: i

        error_flag = 0
        select case (adjustl(filetype))
            case("log")
             call summary_parser(unt,6,section,error_flag)
             read(section,*) Hlt(1:3*Nat*(3*Nat+1)/2)
            case("fchk")
             call read_fchk(unt,'Cartesian Force Constants',data_type,N,A,IA,error_flag)
             if (error_flag /= 0) return
             do i=1,N
                 Hlt(i) = A(i)
             enddo
             deallocate(A)
            case("gms")
             call read_gamess_hess(unt,Nat,Hlt,error_flag)
            case("psi4")
             call read_psi4_hess(unt,Nat,Hlt,error_flag)
            case("molcas")
             call read_molcasUnSym_hess(unt,Nat,Hlt,error_flag)
            case("molpro")
             call read_molpro_hess(unt,Nat,Hlt,error_flag)
            case("gmx")
             call read_gmx_hess(unt,Nat,Hlt,error_flag)
            case default
             write(0,*) "Unsupported filetype:"//trim(adjustl(filetype))
             call supported_filetype_list('freq')
             error_flag = 99
         end select

         return

    end subroutine generic_Hessian_reader


    subroutine generic_dip_reader(unt,filetype,Si,Sf,derivatives,dip_type,dx,Dip,DipD,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic dipole reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! Nat     (io )  int /scalar   Number of atoms
        ! X,Y,Z   (out)  real/vectors  Coordinates
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        integer,intent(in)              :: unt
        character(len=*),intent(in)     :: filetype
        integer,intent(inout)           :: Si, Sf
        logical,intent(inout)           :: derivatives
        character(len=*),intent(in)     :: dip_type
        real(8),intent(inout)           :: dx !Only for numerical diffs
        real(8),dimension(:),intent(out):: Dip 
        real(8),dimension(:),intent(out):: DipD
        integer,intent(out),optional    :: error_flag
        !Local
        integer :: S
        integer :: error_local
        character(len=3) :: dummy_char

        error_local = 0
        select case (adjustl(filetype))
            case("log")
             !Get target state
             call read_gausslog_targestate(unt,S,error_local)
             if (present(error_flag)) error_flag=error_local
             if (Si == -1) Si = 0
             if (Sf == -1) Sf = S
             if (Si /= 0) then
                 write(dummy_char,'(I0)') Si
                 call alert_msg("fatal","TD-DFT calcs in G09 only provide trdip from/to GS, but requested S="//dummy_char)
                 error_local=-1
             else
                 !Need to rewind to read the first dip
                 rewind(unt)
                 !Now read dip
                 call read_gausslog_dip(unt,Si,Sf,dip_type,Dip,error_local)
                 if (derivatives) then
                     print*, " Computing derivatives.."
                     call read_gausslog_dipders(unt,Si,Sf,dip_type,dx,DipD,error_local)
                     if (error_local /= 0) derivatives=.false.
                     !Done this, we can safely reset error_flag to 0
                     error_local = 0
                 endif
             endif
            case("fchk")
             call read_gaussfchk_dip(unt,Si,Sf,derivatives,dip_type,Dip,DipD,error_local)
            case("gms")
             call alert_msg("fatal","Filetype not supported")
            case("psi4")
             call read_psi4_dip(unt,Si,Sf,dip_type,Dip,error_local)
             if (derivatives) then
                 print*, " Computing derivatives.."
                 call read_psi4_dipders(unt,Si,Sf,dip_type,dx,DipD,error_local)
                 if (error_local /= 0) derivatives=.false.
                 !Done this, we can safely reset error_flag to 0
                 error_local = 0
             endif
            case("molcas")
             call alert_msg("fatal","Filetype not supported")
            case("molpro")
             call alert_msg("fatal","Filetype not supported")
            case default
             write(0,*) "Unsupported filetype:"//trim(adjustl(filetype))
             call supported_filetype_list('freq')
             error_local = 99
         end select
         if (present(error_flag)) error_flag=error_local

         return

    end subroutine generic_dip_reader


    subroutine supported_filetype_list(properties)

        !List the supported filetypes for a given 
        !property or set of properties
        
        character(len=*),intent(in) :: properties

        write(0,*) ""
        write(0,'(A)') "Supported filetypes:"

        if (adjustl(properties) == 'freq') then
            write(0,'(A)') " Frequencies:"
            write(0,*)     "  log    : g09 log"
            write(0,*)     "  fchk   : g09 fchk"
            write(0,*)     "  gms    : GAMESS out"
            write(0,*)     "  psi4   : Psi4 out"
            write(0,*)     "  molcas : MOLCAS (UnSym?)"
            write(0,*)     "  molpro : MOLPRO out"
            write(0,*)     "  gmx    : gromacs (g96 and dumped mtx)"
            write(0,'(A)') " Gradients (vertical models):"
            write(0,*)     "  log    : g09 log"
            write(0,*)     "  fchk   : g09 fchk"

        else if (adjustl(properties) == 'trdip') then
            write(0,'(A)') " Transition dipoles:"
            write(0,*)     "   log    : g09 log"
            write(0,*)     "   fchk   : g09 fchk"
            write(0,*)     "   psi4   : Psi4 out"
            write(0,*)     ""  
            write(0,'(A)') " Derivatives:"
            write(0,*)     "   log    : g09 log  (freq)"
            write(0,*)     "   fchk   : g09 fchk (freq)"   
            write(0,*)     "   psi4   : Psi4 out (proper input)"   

        endif

        write(0,*)     ""

        return

    end subroutine supported_filetype_list




end module fcc_io

