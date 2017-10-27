# Source this file to include the Matlab shared libraries in the
# library path.
xxx_matlab_dir="/opt/matlab"
xxx_matlab_version="7.3"
LD_LIBRARY_PATH="$xxx_matlab_dir/bin/glnxa64:$xxx_matlab_dir/sys/os/glnxa64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
MATLAB_VERSION_STRING=-"$xxx_matlab_version"
export LD_LIBRARY_PATH MATLAB_VERSION_STRING
