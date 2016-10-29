#!/bin/sh
# This is intended to be double-clickable from the Finder in OS X 10.5 or later
# or from a File Manager in Linux.  Double-clicking on it will start Priism.
. "/home/tjl10/CUDA_SIMrecon/priism-4.4.1"/Priism_setup.sh
Priism </dev/null >/dev/null 2>&1
