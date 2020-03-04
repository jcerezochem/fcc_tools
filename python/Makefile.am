AUTOMAKE_OPTIONS = subdir-objects
# get *py files to the distribution with "make dist"
EXTRA_DIST = fcc_analyzer.py version_tag.py sum_spectra.py

# ADDITIONAL RECURSIVE TARGETS
# (see http://www.gnu.org/software/automake/manual/html_node/Subdirectories.html)
# This has been added by first including them in configure.ac through
# AM_EXTRA_RECURSIVE_TARGETS([target1 target2...])

# Build standalone binary with pyinstaller
python-standalone-local:
	pyinstaller fcc_analyzer_PyQt4.py --onefile; continue
	pyinstaller fcc_analyzer_PyQt5.py --onefile; continue

python-update-version-local:
	@echo "  "
	@echo "Updating version_tag.py..."
	@echo "Did you first git-commit your changes?"
	@echo "  "
	@./get_git_version.sh || echo "Reusing version file"

all-local: python-update-version-local

install-exec-local:
	cp *.py $(bindir)

# if we want an additional install file
#if INSTALL_PYTHONBIN
#install-exec-local:
#	cp *.py $(prefix)/bin
#	cp dist/fcc_analyzer $(prefix)/bin
#endif
