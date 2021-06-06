#!/bin/bash
# makeotf psf.dv otf.otf -angle -0.804300 -ls 0.2035 -na 1.4 -nimm 1.515 -fixorigin 3 20

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

cd ../"$parent_path"

APP=./cmake_build/cudaSirecon/cudaSireconDriver

# test TIFF
$APP $parent_path/raw.tif $parent_path/proc.tif $parent_path/otf.otf --xyres=0.08 --zres=0.125 -c $parent_path/config

# test MRC
# $APP $parent_path/raw.dv $parent_path/proc.dv $parent_path/otf.otf -c $parent_path/config
