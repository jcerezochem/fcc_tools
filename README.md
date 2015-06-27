# fcc_tools
Tools to facilitate input generation and analysis of FCclasses (see http://www.pi.iccom.cnr.it/fcclasses)

Includes:
* Fortran codes:
  - `fcc_gen_state`: generate state_files and an initial template for the input file, from output of different QM programs (check `gen_fcc_state -h`)
  - `fcc_gen_dipfile`: generates eldip and magdip files from output of QM programs (check `fcc_gen_dipfile -h`)
* Python scripts:
  - `fcc_analyzer.py`: analyzes stick spectra from the output of a TI calculation (check documentation). Requires `python-2.7` and packages `numpy` and `matplotlib`. A stand-alone version generated with pyinstaller is available from the `python-standalone` branch (in `python/dist/fcc_analyzer`). Note that the file is quite heave (93MB) and it might not correespond to the latest fcc_analyzer.py version (check version on output).
