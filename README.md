# fcc_tools
Tools to facilitate input generation and analysis of FCclasses (see http://www.pi.iccom.cnr.it/en/fcclasses)

## Includes:
* Generation of inptut data:
  - `fcc_gen_state`: generate state_files and an initial template for the input file, from output of different QM programs (check `gen_fcc_state -h`)
  - `fcc_gen_dipfile`: generates eldip and magdip files from output of QM programs (check `fcc_gen_dipfile -h`)
* Post-processing tools:
  - `reconvolute_TD` and `reconvolute_TI`: regenerate the spectrum with a new convolution scheme once the TI or TD calculation is done
  - `convolute_RR`: apply a convolution to the stick RR spectra produced by FCclasses3
* Analysis GUIs:
  - `fcc_analyzer_PyQt5.py`: a GUI to analyze stick spectra from the output of a TI calculation (check documentation). Requires `python3` and packages `numpy` and `matplotlib` and `PyQt5`. The interface looks as scketched below.
  - `fcc_RRinspector_PyQt5.py`: a GUI for the inspection of vibrational RR spectra at different incident frequencies (i.e. inspect data from `RR_Spectra_2D.dat`)

![fcc_analyzer_screeshot](https://github.com/jcerezochem/fcc_tools/blob/master/doc/figs/fcc_analyzer_screeshot.png)


## Install
Download the zip or git clone this repository. Then:

1. `./fix_time_stamps.sh`

2. `./configure --bindir INSTALLPATH`

(where INSTALLPATH is the location where the binaries will be placed)

3. `make`

If you have set an appropriate PATH for installation, also:

4.  `make install`

Otherwise, the programs are installed in src/generators

