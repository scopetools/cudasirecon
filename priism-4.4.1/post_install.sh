#!/bin/sh

PATH=/bin:/usr/bin
export PATH

curr_dir=`pwd`

if test \! \( -r Priism_setup -a -r Priism_setup.sh \); then
    echo Could not find a readable Priism_setup and Priism_setup.sh 1>&2
    echo in the current directory \($curr_dir\).  Run post_install.sh 1>&2
    echo from the top-level directory of the Priism distribution. 1>&2
    exit 1
fi

if test $# -gt 0 -a "X-n" = "X$1"; then
    ask_questions=n
else
    ask_questions=y
fi

systype=`uname -s -m | sed -e 's@^Linux ..86@Linux/x86@' -e 's@^Linux x86_64@Linux/x86_64@' -e 's@^Darwin .*@Darwin@' -e 'y@ @_@'`

if test $systype = Linux/x86_64; then
    if test \! -d Linux/x86_64; then
        systype=Linux/x86
    fi
fi

case $systype in
    Linux/x86|Linux/x86_64|Linux/alpha)
        ranlib=
        ranlib_args=
        ;;

    Darwin)
        ranlib=/usr/bin/ranlib
        ranlib_args=
        ranlib_msg_filter="grep -Ev"
        # Filter out some ranlib warnings that are believed to be harmless:
        # ranlib issues a warning about no symbols in some parts of libive.a
        # and on a PowerPC 10.2 system ranlib complains about duplicate symbols
        # in the Intel version.
        ranlib_msg_filter_args='has no symbols|same symbol defined in more than one member|defines symbol'
        ;;

    *)
        echo This does not appear to be a system on which Priism can be 1>&2
        echo used.  uname -s reports `uname -s`.  Linux \(only for x86,
        echo or x86_64\), and Darwin \(only for Mac OS X\) are the system 1>&2
        echo types that are recognized. 1>&2
        exit 1
        ;;
esac

echo Updating Priism_setup and Priism_setup.sh ...
mv Priism_setup Priism_setup.orig
mv Priism_setup.sh Priism_setup.sh.orig

sed -e '/^[ \t]*setenv IVE_BASE/c\
setenv IVE_BASE "'"$curr_dir"'"
' Priism_setup.orig > Priism_setup
sed -e '/^[ \t]*IVE_BASE=/ c\
IVE_BASE="'"$curr_dir"'"
' Priism_setup.sh.orig > Priism_setup.sh

echo IVE_BASE has been adjusted in Priism_setup and Priism_setup.sh.
echo The old versions of those files are in Priism_setup.orig and
echo Priism_setup.sh.orig.
echo

if test \! -d $systype/LIB; then
    echo There is no $systype/LIB directory in the current directory 1>&2
    echo \($curr_dir\).  Without it, Priism will not be functional. 1>&2
    exit 1
elif test \! X$ranlib = X; then
    if test \! -x $ranlib; then
        echo Could not find the command to update archive library symbol
        echo tables on your system.  That is ok.  If at some point you
        echo install the developer tools and wish to link programs with
        echo the archive libraries included with Priism, you will need
        echo to rerun post_install.sh or use $ranlib to update the
        echo tables.
        echo
    elif test \! -w $systype/LIB; then
        echo The $systype/LIB subdirectory is not writable so the archive
        echo library symbol tables have not been updated.  That is ok.
        echo If at some point you wish to link programs with the archive
        echo libraries included with Priism, you will need to have
        echo someone with appropriate privileges rerun post_install.sh
        echo or use $ranlib to update the tables.
        echo
    else
        echo Updating achive library tables ...
        find $systype/LIB -name '*.a' \! -type l -print | xargs -n 1 $ranlib $ranlib_args 2>&1 | $ranlib_msg_filter "$ranlib_msg_filter_args"
        echo
    fi
fi

echo Setting up a double-clickable Priism.command \(requires Linux or OS X 10.5 or later\) ...
cat >Priism.command <<EOF
#!/bin/sh
# This is intended to be double-clickable from the Finder in OS X 10.5 or later
# or from a File Manager in Linux.  Double-clicking on it will start Priism.
. "$curr_dir"/Priism_setup.sh
Priism </dev/null >/dev/null 2>&1
EOF
chmod a+x Priism.command
echo

case $systype in
    Linux/x86)
	matlab_version="6.5 or later installed"
        ;;

    Linux/x86_64)
        matlab_version=7
        ;;

    Darwin)
        if test X"`uname -m`" \!= Xi386 ; then 
            matlab_version="6.5 or later installed"
        else
            matlab_version=7
        fi
        ;;

    *)
        matlab_version=
        ;;
esac

if test $ask_questions = y -a \! "X$matlab_version" = X; then
    setup_file=$systype/BIN/matlab_setup.sh
    if test -r $setup_file; then
        matlab_dir=`sed -n 's/xxx_matlab_dir="\(.*\)"/\1/p' <$setup_file`
        matlab_version=`sed -n 's/xxx_matlab_version="\(.*\)"/\1/p' <$setup_file`
        echo The MatToImg and ImgToMat tools in Priism have been configured
        echo to expect Matlab $matlab_version in $matlab_dir.
        echo Do you want to change the configuration \(y/n\)?
        read ask_matlab
    else
        echo The MatToImg and ImgToMat tools in Priism convert between
        echo Priism image files and Matlab .mat files.  Do you have
        echo Matlab $matlab_version installed \(y/n\)?
        read ask_matlab
    fi
    if test "X$ask_matlab" = Xy; then
        echo What is the directory with the Matlab files?
        read matlab_dir
        if test $systype = Darwin; then
            if test X"`uname -m`" != Xi386; then
                echo 'Is the Matlab version'
                echo '    (a) >= 6.5 and < 7'
                echo '    or (b) >= 7'
                echo '(a/b)?'
                read ask_matlab
                if test "X$ask_matlab" = Xa; then
                    matlab_version=6.5
                else
                    matlab_version=7
                fi
            else
                echo 'Are you using Matlab 7.9 (R2009b) or later'
                echo 'on a 64-bit Mac with OS X 10.5 or later (y/n)?'
                read ask_matlab
                if test "X$ask_matlab" = Xy; then
                    matlab_version=7.9
                else
                    matlab_version=7
                fi
            fi
        elif test $systype = Linux/x86 -o \( $systype = Darwin -a X"`uname -m`" != Xi386 \); then
            echo 'Is the Matlab version'
            echo '    (a) >= 6.5 and < 7'
            echo '    or (b) >= 7 and < 7.3'
            echo '    or (c) >= 7.3'
	    echo '(a/b/c)?'
            read ask_matlab
            if test "X$ask_matlab" = Xa; then
                matlab_version=6.5
            elif test "X$ask_matlab" = Xb; then
                matlab_version=7
            else
                matlab_version=7.3
            fi
        elif test $systype = Linux/x86_64; then
            echo 'Is the Matlab version'
            echo '    (a) >= 7 and < 7.3'
            echo '    or (b) >= 7.3'
            echo '(a/b)?'
            read ask_matlab
            if test "X$ask_matlab" = Xa; then
                matlab_version=7
            else
                matlab_version=7.3
            fi
        fi

        echo Updating $setup_file ...
        echo
	if test -f "$setup_file"; then
	    mv $setup_file ${setup_file}.orig
	fi

	echo \# Source this file to include the Matlab shared libraries in the >$setup_file
	echo \# library path. >>$setup_file
        echo xxx_matlab_dir=\"$matlab_dir\" >>$setup_file
	echo xxx_matlab_version=\"$matlab_version\" >>$setup_file
	case $systype in
	    Linux/x86)
                echo if test \$xxx_matlab_version = 6.5\; then >>$setup_file
		echo "    " LD_LIBRARY_PATH=\"\$xxx_matlab_dir/extern/lib/glnx86:\$xxx_matlab_dir/sys/os/glnx86\$\{LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH\}\" >>$setup_file
                echo "    " MATLAB_VERSION_STRING=-6.5 >>$setup_file
                echo else >>$setup_file
                echo "    " LD_LIBRARY_PATH=\"\$xxx_matlab_dir/bin/glnx86:\$xxx_matlab_dir/sys/os/glnx86\$\{LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH\}\" >>$setup_file
                echo "    " MATLAB_VERSION_STRING=-'"$xxx_matlab_version"' >>$setup_file
                echo fi >>$setup_file
		echo export LD_LIBRARY_PATH MATLAB_VERSION_STRING >>$setup_file
		;;

            Linux/x86_64)
                echo LD_LIBRARY_PATH=\"\$xxx_matlab_dir/bin/glnxa64:\$xxx_matlab_dir/sys/os/glnxa64\$\{LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH\}\" >>$setup_file
                echo MATLAB_VERSION_STRING=-'"$xxx_matlab_version"' >>$setup_file
		echo export LD_LIBRARY_PATH MATLAB_VERSION_STRING >>$setup_file
                if test -d Linux/x86/BIN; then
                    alt_setup_file=Linux/x86/BIN/matlab_setup.sh
                    echo Updating $alt_setup_file ...
                    echo
                    if test -f "$alt_setup_file"; then
                        mv $alt_setup_file ${alt_setup_file}.orig
                    fi

                    echo \# Source this file to include the Matlab shared libraries in the >$alt_setup_file
                    echo \# library path. >>$alt_setup_file
                    echo xxx_matlab_dir=\"$matlab_dir\" >>$alt_setup_file
	            echo xxx_matlab_version=\"$matlab_version\" >>$alt_setup_file
                    echo LD_LIBRARY_PATH=\"\$xxx_matlab_dir/bin/glnx86:\$xxx_matlab_dir/sys/os/glnx86\$\{LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH\}\" >>$alt_setup_file
                    echo MATLAB_VERSION_STRING=-'"$xxx_matlab_version"' >>$alt_setup_file
	   	    echo export LD_LIBRARY_PATH MATLAB_VERSION_STRING >>$alt_setup_file
                fi
                ;;

	    Darwin)
                echo if test \$xxx_matlab_version = 6.5\; then >>$setup_file
		echo "    " DYLD_FALLBACK_LIBRARY_PATH=\"\$xxx_matlab_dir/extern/lib/mac:\$xxx_matlab_dir/sys/os/mac\$\{DYLD_FALLBACK_LIBRARY_PATH:+:\$DYLD_FALLBACK_LIBRARY_PATH\}\" >>$setup_file
                echo else >>$setup_file
                echo "    " if test X\"\`uname -m\`\" = Xi386\; then >>$setup_file
                echo "        " DYLD_FALLBACK_LIBRARY_PATH=\"\$xxx_matlab_dir/bin/maci64:\$xxx_matlab_dir/bin/maci\$\{DYLD_FALLBACK_LIBRARY_PATH:+:\$DYLD_FALLBACK_LIBRARY_PATH\}\" >>$setup_file
                echo "    " else >>$setup_file
		echo "        " DYLD_FALLBACK_LIBRARY_PATH=\"\$xxx_matlab_dir/bin/mac:\$xxx_matlab_dir/sys/os/mac\$\{DYLD_FALLBACK_LIBRARY_PATH:+:\$DYLD_FALLBACK_LIBRARY_PATH\}\" >>$setup_file
                echo "    " fi >>$setup_file
                echo fi >>$setup_file
                echo if test \$xxx_matlab_version != 6.5 -a \$xxx_matlab_version != 7\; then >>$setup_file
                echo "    " MATLAB_VERSION_STRING=-'"$xxx_matlab_version"' >>$setup_file
                echo else >>$setup_file
                echo "    " MATLAB_VERSION_STRING=\"\" >>$setup_file
                echo fi >>$setup_file
		echo export DYLD_FALLBACK_LIBRARY_PATH MATLAB_VERSION_STRING >>$setup_file
		if test -d Darwin64/BIN -a x7.9 = "x$matlab_version" ; then
		    cp $setup_file Darwin64/BIN/matlab_setup.sh
                fi
		;;
	esac
    fi
fi

if test -r grecsrv/bin/grecsrv ; then
    echo Updating grecsrv/bin/grecsrv
    cp grecsrv/bin/grecsrv grecsrv/bin/grecsrv.orig
    sed -e '/^[ \t]*IVE_BASE=/ c\
IVE_BASE="'"$curr_dir"'"
' grecsrv/bin/grecsrv.orig > grecsrv/bin/grecsrv
fi

echo post_install.sh complete.
