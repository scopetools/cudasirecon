#!/bin/sh

# This is intended to be double-clickable from the Finder in OS X or
# a file manager in Linux.  Double-clicking on it will run the Priism
# post-installation procedure.

cd `dirname $0`
if test "`uname`" = Darwin ; then
    ./post_install.sh
else
    xterm -T "Priism post_install" -e ./post_install.sh
fi
