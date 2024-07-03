AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = reconvolute_TD reconvolute_TI convolute_RR excitation_RR addconvolution_RR addconvolution

reconvolute_TD_SOURCES =   \
          version.f90                                         \
          ../modules/line_preprocess.f90                      \
          ../modules/alerts.f90                               \
          ../modules/constants.f90                       
      
# Now use a conditional to use FFTW of DFT
if FFTW
reconvolute_TD_SOURCES += \
          ../modules/maths/fftw.f90                           
else
reconvolute_TD_SOURCES += \
          ../modules/maths/dft.f90 
endif

reconvolute_TD_SOURCES += \
          reconvolute_TD.f90

reconvolute_TD_LDFLAGS = @LIBS@

reconvolute_TI_SOURCES =   \
          version.f90                                         \
          ../modules/line_preprocess.f90                      \
          ../modules/alerts.f90                               \
          ../modules/constants.f90                            \
          reconvolute_TI.f90
          
reconvolute_TI_LDFLAGS = @LIBS@

convolute_RR_SOURCES =   \
          version.f90                                         \
          ../modules/line_preprocess.f90                      \
          ../modules/alerts.f90                               \
          ../modules/constants.f90                            \
          convolute_RR.f90

convolute_RR_LDFLAGS = @LIBS@

excitation_RR_SOURCES =   \
          version.f90                                         \
          ../modules/line_preprocess.f90                      \
          ../modules/alerts.f90                               \
          ../modules/constants.f90                            \
          excitation_RR.f90
          
excitation_RR_LDFLAGS = @LIBS@

addconvolution_RR_SOURCES =   \
          version.f90                                         \
          ../modules/line_preprocess.f90                      \
          ../modules/alerts.f90                               \
          ../modules/constants.f90                            \
          addconvolution_RR.f90

addconvolution_SOURCES =   \
          version.f90                                         \
          ../modules/line_preprocess.f90                      \
          ../modules/alerts.f90                               \
          ../modules/constants.f90                            \
          addconvolution.f90
          
addconvolution_RR_LDFLAGS = @LIBS@

CLEANFILES=*.mod

.PHONY: ../modules/line_preprocess.f90 ../modules/alerts.f90 ../modules/constants.f90

# Generate version file only if needed. To that end the
# dependency of version.f90 is a non-file-associated 
# instruction that only changes the file if requiered
# This avoids recompiling everything every time
version.f90: get_git_version
get_git_version:
	@bash ./get_git_version.sh $(FCFLAGS) $(LIBS) || echo "Reusing version file"

