#
# A file meant to be sourced by users of sh (or compatible shells) to set
# the environment for the most recent release of Priism version 4 and
# related programs (EMCAT).  Definitions that you shall need to modify are
# IVE_BASE and IVE_ENV_SETUP.  There are heuristics to set IVE_WORKING_SET
# and IVE_SIZE.  You can either modify those heuristics or set IVE_WORKING_SET
# and IVE_SIZE after sourcing this file.  The defaults for IVE_OLD_LOOK,
# IVE_WORKING_UNIT, IVE_PGTHRESH, and IVE_SHMDIR may also need to be modified.
#
# Setting any of the following environment variables prior to sourcing
# this file affects the definitions set:
#
# IVE_NO_INTERNAL_LQT
#     If set on x86 or x86_64 Linux systems or 64-bit Intel Mac systems when
#     using the 64-bit version of Priism, the value of LIBQUICKTIME_PLUGIN_DIR
#     is not set to point to the internal version of libquicktime.  If you do
#     that, you may also want to make modifications so that the Priism
#     executables pick up the alternate version of the libquicktime shared
#     library.  On Linux, two ways to do that are to set the environment
#     variable, LD_LIBRARY_PATH, to include the path to the alternate
#     libquicktime library or to alter the links in Linux/x86/LIB and
#     Linux/x86_64/LIB to point to that alternate library.  On the Mac,
#     two ways to do that are to set DYLD_LIBRARY_PATH to include the path
#     to the alternative libquicktime library or to alter the link in
#     Darwin64/LIB to point to the alternative libquicktime library.
#
# IVE_PGTHRESH
#     Overrides the default of 1024 pixels.
#
# IVE_PREFER_32
#     If on systems where either 32-bit or 64-bit operation is possible
#     (currently only x86_64 Linux systems or Intel Macs with 64-bit
#     processors), use the 32-bit executables and libraries rather than
#     the 64-bit ones.
#
# IVE_WORKING_UNIT
#     Overrides the default value of the smaller of the value of
#     IVE_WORKING_SET or 128 megabytes.
#

if test -x /bin/sed; then
  xxx_sed=/bin/sed
elif test -x /usr/bin/sed; then
  xxx_sed=/usr/bin/sed
else
  unset xxx_sed
fi
if test -x /bin/uname; then
  xxx_uname=/bin/uname
elif test -x /usr/bin/uname; then
  xxx_uname=/usr/bin/uname
else
  unset xxx_uname
fi
if test -x /bin/expr; then
  xxx_expr=/bin/expr
elif test -x /usr/bin/expr; then
  xxx_expr=/usr/bin/expr
else
  unset xxx_expr
fi
if test -x /bin/cat; then
  xxx_cat=/bin/cat
elif test -x /usr/bin/cat; then
  xxx_cat=/usr/bin/cat
else
  unset xxx_cat
fi

# For OS dependent information, query the system type and version number.
unset xxx_os_vers
if test -n "$xxx_uname" -a -n "$xxx_sed"; then
  xxx_systype=`$xxx_uname -s -m | $xxx_sed -e 's@^Linux ..86@Linux/x86@' -e 's@^Linux x86_64@Linux/x86_64@' -e 's@^Darwin .*@Darwin@' -e 'y@ @_@'`
  xxx_os_vers=`$xxx_uname -r`
else
  xxx_systype=unknown
fi

# IVE_BASE points to the top level directory for the Priism distribution.
IVE_BASE="/home/tjl10/CUDA_SIMrecon/priism-4.4.1"
export IVE_BASE

case $- in
  *i*) if test -r "$IVE_BASE"/CONFIG/Version.def -a -n "$xxx_cat"; then
         xxx_vers_name=`$xxx_cat "$IVE_BASE"/CONFIG/Version.def`
       else
         xxx_vers_name=4.0
       fi
       ;;
  *)   ;;
esac

if test $xxx_systype = Linux/x86_64; then
    if test \! -d "$IVE_BASE"/Linux/x86_64 -o -n "${IVE_PREFER_32+X}"; then
        xxx_systype=Linux/x86
    fi
elif test $xxx_systype = Darwin -a -n "$xxx_sed" -a -n "$xxx_uname"; then
    if test -z "${IVE_PREFER_32+X}" -a \( x"`$xxx_uname -m`" = xi386 -o \
        x"`$xxx_uname -m`" = xx86_64 \) -a \
        x`sysctl -n hw.optional.x86_64 2>/dev/null` = x1 -a \
        x`sw_vers | $xxx_sed -n 's/^ProductVersion:[^0123456789]\{1,\}\([0123456789]\{1,\}\.[0123456789]\{1,\}\).\{1,\}/\1/p'` \!= x10.4; then
        xxx_systype=Darwin64
    fi
fi

if test $xxx_systype = Linux/x86 -o \
        $xxx_systype = Linux/x86_64 -o \
        $xxx_systype = Darwin -o \
        $xxx_systype = Darwin64; then

    # If the environment in which batch jobs run differs from that for
    # the interactive session, then you can set IVE_ENV_SETUP to hold a
    # sh command to set the appropriate environment for the shell command.
    # If the environments are the same, then you will have no need to set
    # IVE_ENV_SETUP.  The example checks for a single configuration file and
    # exits if it could not be found.
    IVE_ENV_SETUP="{ test -r '${IVE_BASE}/Priism_setup.sh' && . '${IVE_BASE}/Priism_setup.sh' ; } || exit 1"
    export IVE_ENV_SETUP

    # Setting this variable to yes (case insensitive) will cause IVE to
    # mimic the old color and font scheme.
    #IVE_OLD_LOOK=yes; export IVE_OLD_LOOK

    # IVE_SIZE sets the maximum size of the shared memory file in megabytes
    # (Priism begins to flail badly when this size is exceeded).  If not set,
    # the default value is 300 megabytes.  On 32-bit or 64-bit systems with a
    # Priism that has been compiled with 32-bit support, the absolute maximum
    # possible is 4095; on 64-bit systems with a Priism that has been compiled
    # with 64-bit support the absolute maximum possible is 1048575 (this limit
    # was chosen because IRIX 64 limits the user address space to 2^40 bytes).
    # On IRIX and Mac OS X, the Priism interface with 32-bit support will not
    # start for values larger than 1000.  On a 32-bit Linux, the 32-bit Priism
    # interface will not start for values larger than 2000.  Large values may
    # unduly limit how much memory Priism applications are able to allocate
    # from the heap.
    #
    # IVE_WORKING_SET is the amount of shared memory, in megabytes, to use
    # before either discarding parts that can be recalculated or reread
    # (scaled images, data available on disk) or writing out data to remove it
    # from shared memory.  IVE_WORKING_SET must be less than IVE_SIZE to leave
    # rooom for data not allocated out of the working set (bookkeeping
    # information; image headers; scratch window data).  As a starting guess
    # you could set IVE_WORKING_SET to be IVE_SIZE minus 64 megabytes though
    # the minimum difference will depend on the kinds and amount of data
    # loaded.
    IVE_WORKING_SET=96
    unset IVE_SIZE

    unset xxx_memavailable
    case $xxx_systype in
      Linux/x86|Linux/x86_64)
        if test -r /proc/meminfo -a -n "$xxx_sed"; then
          xxx_memavailable=`$xxx_sed -n 's/^MemTotal: *\([0123456789]\{1,\}\)[0123456789]\{3\} kB/\1/p' /proc/meminfo`
        fi
        ;;

      Darwin|Darwin64)
        if test -x /usr/bin/hostinfo -a -n "$xxx_sed"; then
          # To avoid any more complexity, use 1000 (instead of the correct
          # value, 1024) when converting from megabytes to gigabytes.
          xxx_memavailable=`/usr/bin/hostinfo | $xxx_sed -n -e 's/\.[0123456789]\{0,\} \{1,\}megabytes//' -e 's/\.\([0123456789]\{2,2\}\)[0123456789]* \{1,\}gigabytes/\10\./' -e 's/^Primary memory available: \{1,\}\([0123456789]\{1,\}\)\..*/\1/p'`
        fi
        ;;
    esac

    if test -n "$xxx_memavailable" -a -n "$xxx_expr"; then
      # These are just rough guesses to try to leave more memory to the system
      # on machines with more memory since they are likely to be doing
      # something else.
      if test `$xxx_expr X$xxx_memavailable : 'X[0123456789][0123456789]*$'` -gt 1; then
        if test $xxx_memavailable -le 64; then
          IVE_WORKING_SET=`$xxx_expr $xxx_memavailable - 16`
        elif test $xxx_memavailable -le 128; then
          IVE_WORKING_SET=`$xxx_expr $xxx_memavailable - 32`
        elif test $xxx_memavailable -le 320; then
          IVE_WORKING_SET=`$xxx_expr $xxx_memavailable - 64`
        else
          if test $xxx_systype = Linux/x86_64 -o $xxx_systype = Darwin64; then
            IVE_WORKING_SET=`$xxx_expr $xxx_memavailable / 2 + $xxx_memavailable / 4`
            IVE_SIZE=`$xxx_expr $IVE_WORKING_SET + $IVE_WORKING_SET / 8`
            if test $IVE_SIZE -gt 1048575; then
              IVE_SIZE=1048575
              IVE_WORKING_SET=917500
            fi
            export IVE_SIZE
          else
            IVE_WORKING_SET=256
          fi
        fi
      fi
      unset xxx_memavailable
    fi
    export IVE_WORKING_SET

    # IVE_WORKING_UNIT is the amount, in megabytes, that is allocated when
    # the shared memory file must be extended to accomodate new images.  
    # IVE_WORKING_UNIT should be at least as large as a single section of
    # data (if not, Priism can fail to load new data whose section size is
    # larger than the IVE_WORKING_UNIT when you have already filled the
    # shared memory with data whose section size is less than or equal to
    # the IVE_WORKING_UNIT).
    if test -z "$IVE_WORKING_UNIT"; then
      IVE_WORKING_UNIT=128; export IVE_WORKING_UNIT
      if test $IVE_WORKING_UNIT -gt "$IVE_WORKING_SET"; then
          IVE_WORKING_UNIT="$IVE_WORKING_SET"
      fi
    fi

    # For reads from windows that involve less than IVE_PGTHRESH pixels
    # and for which the data is not available in memory, reading back in
    # a whole section will be bypassed and only the lines containing the
    # needed data will be read.  This has been set so an application reading
    # a single pixel or whole line for typical data sets will not cause
    # a whole section to be read and so that reads on whole sections
    # (bigger than 16 x 16) will cause the data to be swapped in as needed.
    # A value of 0 for IVE_PGTHRESH disables this feature.
    if test -z "$IVE_PGTHRESH"; then
      IVE_PGTHRESH=1024; export IVE_PGTHRESH
    fi

    # This specifies the directory for the shared memory file.  Because it
    # can be large (100s of megabytes), this directory should have plenty
    # of space.  It is also beneficial if access to the directory is fast.
    IVE_SHMDIR="${IVE_SHMDIR:-${TMPDIR:-/var/tmp}}"
    export IVE_SHMDIR

    #
    # Settings below this line are independent of the system.
    #
    IVE_EXE="${IVE_BASE}/${xxx_systype}/BIN"; export IVE_EXE

    # $IVE_EXE as first in the path after taking them out if they are
    # already present.
    if test -n "$xxx_sed"; then
      PATH="${IVE_EXE}"${PATH:+`echo ":${PATH}:" | $xxx_sed -e "s%:${IVE_EXE}:%:%g" -e 's%:$%%'`}
    else
      PATH="${IVE_EXE}${PATH:+:${PATH}}"
    fi
    export PATH

    # Set up path to the plugins in the internal version of libquicktime on
    # x86 or x86_64 Linux systems unless IVE_NO_INTERNAL_LQT is set.
    if test -z "${IVE_NO_INTERNAL_LQT+X}" -a \( $xxx_systype = Linux/x86 -o $xxx_systype = Linux/x86_64 -o $xxx_systype = Darwin64 \); then
        if test $xxx_systype = Linux/x86; then
          LIBQUICKTIME_PLUGIN_DIR="${IVE_BASE}"/libquicktime/x86/lib/libquicktime
        elif test $xxx_systype = Linux/x86_64; then
          LIBQUICKTIME_PLUGIN_DIR="${IVE_BASE}"/libquicktime/x86_64/lib/libquicktime
        else
          LIBQUICKTIME_PLUGIN_DIR="${IVE_BASE}"/libquicktime/Darwin/lib/libquicktime
        fi
        export LIBQUICKTIME_PLUGIN_DIR
    fi

    # Put $IVE_BASE/$xxx_systype/LIB first in the library search path after
    # removing it if it is already there.
    case $xxx_systype in
      Darwin|Darwin64)
        # If DYLD_FALLBACK_LIBRARY_PATH is not set or is empty, set it to its
        # default value ($HOME/lib:/usr/local/lib:/lib:/usr/lib according to
        # the dyld man page).
        DYLD_FALLBACK_LIBRARY_PATH="${DYLD_FALLBACK_LIBRARY_PATH:-${HOME:+${HOME}/lib:}/usr/local/lib:/lib:/usr/lib}"
        if test -n "$xxx_sed"; then
          # If available, include the 32-bit versions of the Priism libraries
          # for backwards compatibility with programs not packaged with Priism
          # which use the Priism libraries.
          if test $xxx_systype = Darwin64 -a -d "${IVE_BASE}"/Darwin/LIB; then
            DYLD_FALLBACK_LIBRARY_PATH="${IVE_BASE}/Darwin/LIB"`echo ":${DYLD_FALLBACK_LIBRARY_PATH}:" | $xxx_sed -e "s%:${IVE_BASE}/Darwin/LIB:%:%g" -e 's%:$%%'`
          fi
          DYLD_FALLBACK_LIBRARY_PATH="${IVE_BASE}/${xxx_systype}/LIB"`echo ":${DYLD_FALLBACK_LIBRARY_PATH}:" | $xxx_sed -e "s%:${IVE_BASE}/${xxx_systype}/LIB:%:%g" -e 's%:$%%'`
        else
          if test $xxx_systype = Darwin64 -a -d "${IVE_BASE}"/Darwin/LIB; then
            DYLD_FALLBACK_LIBRARY_PATH="${IVE_BASE}/Darwin/LIB:${DYLD_FALLBACK_LIBRARY_PATH}"
          fi
          DYLD_FALLBACK_LIBRARY_PATH="${IVE_BASE}/${xxx_systype}/LIB:${DYLD_FALLBACK_LIBRARY_PATH}"
        fi
        export DYLD_FALLBACK_LIBRARY_PATH
        ;;
      *)
        xxx_misc_tiff_lib=
        xxx_misc_vp1000_lib=
        if test -n "$xxx_sed"; then
          if test -n "$xxx_misc_tiff_lib"; then
            LD_LIBRARY_PATH="${xxx_misc_tiff_lib}"${LD_LIBRARY_PATH:+`echo ":${LD_LIBRARY_PATH}:" | $xxx_sed -e "s%:${xxx_misc_tiff_lib}:%:%g" -e 's%:$%%'`}
          fi
          if test -n "$xxx_misc_vp1000_lib"; then
            LD_LIBRARY_PATH="${xxx_misc_vp1000_lib}"${LD_LIBRARY_PATH:+`echo ":${LD_LIBRARY_PATH}:" | $xxx_sed -e "s%:${xxx_misc_vp1000_lib}:%:%g" -e 's%:$%%'`}
          fi
          # If available, include the 32-bit versions of the Priism libraries
          # for backwards compatibility with programs not packaged with Priism
          # which use Priism libraries.
          if test $xxx_systype = Linux/x86_64; then
              if test -d "${IVE_BASE}"/Linux/x86/LIB; then
                  LD_LIBRARY_PATH="${IVE_BASE}/Linux/x86/LIB"${LD_LIBRARY_PATH:+`echo ":${LD_LIBRARY_PATH}:" | $xxx_sed -e "s%:${IVE_BASE}/Linux/x86/LIB:%:%g" -e 's%:$%%'`}
              fi
          fi
          LD_LIBRARY_PATH="${IVE_BASE}/${xxx_systype}/LIB"${LD_LIBRARY_PATH:+`echo ":${LD_LIBRARY_PATH}:" | $xxx_sed -e "s%:${IVE_BASE}/${xxx_systype}/LIB:%:%g" -e 's%:$%%'`}
        else
          LD_LIBRARY_PATH="${IVE_BASE}/${xxx_systype}/LIB"
          # If available, include the 32-bit versions of the Priism libraries
          # for backwards compatibility with programs not packaged with Priism
          # which use Priism libraries.
          if test $xxx_systype = Linux/x86_64; then
              if test -d "${IVE_BASE}"/Linux/x86/LIB; then
                  LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${IVE_BASE}/Linux/x86/LIB"
               fi
          fi
          LD_LIBRARY_PATH="${LD_LIBRARY_PATH}${xxx_misc_vp1000_lib:+:${xxx_misc_vp1000_lib}}${xxx_misc_tiff_lib:+:${xxx_misc_tiff_lib}}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
        fi
        unset xxx_misc_tiff_lib
        unset xxx_misc_vp1000_lib
        export LD_LIBRARY_PATH
        ;;
    esac

    if test $xxx_systype = Darwin -o $xxx_systype = Darwin64; then
	# XFree86 4.2 includes a libXt.dylib linked as a two-level namespace
	# That causes Motif-based applications to fail unless they are linked
        # with -force_flat_namespace or -bind_at_load or run with 
        # DYLD_FORCE_FLAT_NAMESPACE set.  Priism used the
        # -force_flat_namespace hack for versions prior to 4.1.4, but that
        # causes unacceptably slow startup times for image windows when used
        # with Apple's X11.  Therefore, use the workaround below for XFree86
        # 4.2.  XFree86 4.2.0.1 has libXt.dylib linked as a flat namespace.
	if test -x "${IVE_EXE}/istwolevel"; then
            if test "X`"${IVE_EXE}/istwolevel" /usr/X11R6/lib/libXt.dylib`" = X1; then
	        set DYLD_FORCE_FLAT_NAMESPACE
                export DYLD_FORCE_FLAT_NAMESPACE
            fi
        fi
    elif test $xxx_systype = Linux/x86; then
        # This helps reduce the number of 1022 rendering errors from the
        # Volume Viewer plugin which uses the Volume Pro 1000 board.
        VLIResizeImage=true
        export VLIResizeImage
    fi

    case $- in
      *i*) echo "Priism $xxx_vers_name set up; type Priism to run it";;
      *) ;;
    esac

else

    case $- in
      *i*) if test -n "$xxx_uname"; then
	     echo "Sorry, Priism $xxx_vers_name not known to run with" `$xxx_uname - sr`
           else
             echo "Sorry, unknown system type; can not run Priism $xxx_vers_name"
           fi
           ;;
      *)   ;;
    esac

fi

unset xxx_sed
unset xxx_uname
unset xxx_cat
unset xxx_expr
unset xxx_systype
unset xxx_os_vers
case $- in
  *i*)   unset xxx_vers_name;;
  *) ;;
esac
