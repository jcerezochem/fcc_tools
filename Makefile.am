AUTOMAKE_OPTIONS = subdir-objects

SUBDIRS = src/generators tests/ python/ 

EXTRA_DIST = LICENSE INSTRUCTIONS AUTHORS src/modules run_test.sh python/fort* doc tests/* 

.PHONY: test version.f90
test:
	./run_test.sh


