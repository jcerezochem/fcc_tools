AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = gen_fcc_state gen_fcc_dipfile

gen_fcc_state_SOURCES = version.f90\
                        ../modules/verbosity.f90\
                        ../modules/constants.f90\
                        ../modules/alerts.f90\
                        ../modules/MatrixMod_red.f90\
                        ../modules/line_preprocess.f90\
                        ../modules/fcc_basics.f90\
                        ../modules/gamess_manage.f90\
                        ../modules/gaussian_manage.f90\
                        ../modules/gmx_manage.f90\
                        ../modules/psi4_manage.f90\
                        ../modules/molcas_manage.f90\
                        ../modules/molpro_manage.f90\
                        ../modules/turbomol_manage.f90\
                        ../modules/orca_manage.f90\
                        ../modules/fcc_io.f90\
                        ../modules/vibrational_analysis.f90\
                        gen_fcc_state.f90

gen_fcc_state_LDFLAGS = @LIBS@

gen_fcc_dipfile_SOURCES = version.f90\
                        ../modules/verbosity.f90\
                          ../modules/constants.f90\
                          ../modules/alerts.f90\
                          ../modules/MatrixMod_red.f90\
                          ../modules/line_preprocess.f90\
                          ../modules/fcc_basics.f90\
                          ../modules/gamess_manage.f90\
                          ../modules/gaussian_manage.f90\
                          ../modules/gmx_manage.f90\
                          ../modules/psi4_manage.f90\
                          ../modules/molcas_manage.f90\
                          ../modules/molpro_manage.f90\
                          ../modules/turbomol_manage.f90\
                          ../modules/orca_manage.f90\
                          ../modules/fcc_io.f90\
                          gen_fcc_dipfile.f90

gen_fcc_dipfile_LDFLAGS = @LIBS@

CLEANFILES=*.mod

# .PHONY: version.f90
version.f90: get_git_version
get_git_version:
	@bash ./get_git_version.sh $(FCFLAGS) $(LIBS) || echo "Reusing version file"


