#!/bin/bash
# This script run ac tools, fix the configure script to handle MKL test properly and runs am with proper flags

aclocal; autoconf; sed -i "s/-lmkl_intel/-lmkl_intel_lp64 -lmkl_sequential -lmkl_core/g" configure; automake -a --foreign


