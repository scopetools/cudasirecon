#!/bin/bash
# makeotf psf.dv otf.otf -angle -0.804300 -ls 0.2035 -na 1.4 -nimm 1.515 -fixorigin 3 20

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
APP="$parent_path/../cmake_build/cudaSirecon/cudaSireconDriver"

$APP $parent_path raw --otf-file $parent_path/psf.tif --xyres=0.08 --zres=0.125 --zresPSF=0.125 -c $parent_path/config-tiff
