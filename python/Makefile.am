AUTOMAKE_OPTIONS = subdir-objects
# get *py files to the distribution with "make dist"
EXTRA_DIST = fcc_analyzer_PyQt5.py fcc_analyzer_PyQt5_version.py fcc_RRinspector_PyQt5.py sum_spectra.py

fcc_analyzer_PyQt5_version.py:python-update-version-local

python-update-version-local:
	@echo "  "
	@echo "Updating version_tag.py..."
	@echo "Did you first git-commit your changes?"
	@echo "  "
	@./get_git_version.sh || echo "Reusing version file"

all-local: python-update-version-local

install-exec-local:
	cp *.py $(bindir)

