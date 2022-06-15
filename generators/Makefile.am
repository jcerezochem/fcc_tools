AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = gen_fcc_state gen_fcc_dipfile

gen_fcc_state_SOURCES = version.f90\
                        ../modules/verbosity.f90\
                        ../modules/constants.f90\
                        ../modules/alerts.f90\
                        ../modules/line_preprocess.f90\
                        ../modules/maths/MatrixMod_red.f90\
                        ../modules/fcc_basics.f90\
                        ../modules/io/gamess_manage.f90\
                        ../modules/io/gaussian_manage.f90\
                        ../modules/io/cfour_manage.f90\
                        ../modules/io/gmx_manage.f90\
                        ../modules/io/psi4_manage.f90\
                        ../modules/io/molcas_manage.f90\
                        ../modules/io/molden_manage.f90\
                        ../modules/io/molpro_manage.f90\
                        ../modules/io/turbomol_manage.f90\
                        ../modules/io/orca_manage.f90\
                        ../modules/io/qchem_manage.f90\
                        ../modules/io/fcc_manage.f90\
                        ../modules/io/fcc_io.f90\
                        ../modules/vibrational_analysis.f90\
                        gen_fcc_state.f90

gen_fcc_state_LDFLAGS = @LIBS@

gen_fcc_dipfile_SOURCES = version.f90\
                        ../modules/verbosity.f90\
                          ../modules/constants.f90\
                          ../modules/alerts.f90\
                          ../modules/line_preprocess.f90\
                          ../modules/maths/MatrixMod_red.f90\
                          ../modules/fcc_basics.f90\
                          ../modules/io/gamess_manage.f90\
                          ../modules/io/gaussian_manage.f90\
                          ../modules/io/cfour_manage.f90\
                          ../modules/io/gmx_manage.f90\
                          ../modules/io/psi4_manage.f90\
                          ../modules/io/molcas_manage.f90\
                          ../modules/io/molden_manage.f90\
                          ../modules/io/molpro_manage.f90\
                          ../modules/io/turbomol_manage.f90\
                          ../modules/io/orca_manage.f90\
                          ../modules/io/qchem_manage.f90\
                          ../modules/io/fcc_manage.f90\
                          ../modules/io/fcc_io.f90\
                          gen_fcc_dipfile.f90

gen_fcc_dipfile_LDFLAGS = @LIBS@

CLEANFILES=*.mod

.PHONY: ../modules/line_preprocess.f90 ../modules/alerts.f90 ../modules/constants.f90

# .PHONY: version.f90
version.f90: get_git_version
get_git_version:
	@bash ./get_git_version.sh $(FCFLAGS) $(LIBS) || echo "Reusing version file"


