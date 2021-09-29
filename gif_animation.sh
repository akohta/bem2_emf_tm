#!/bin/bash

# settings
Delay=10
Loop=0

if [ $# -ne 1 ]; then
 echo "Usage : $0 image_directory"
 exit
fi

if [ ! -d $1 ]; then
 echo "$1 does not exist or is not a directory."
 exit
fi

convert -delay $Delay -loop $Loop $1/xy_Hz_*.png xy_Hz.gif &
convert -delay $Delay -loop $Loop $1/xy_Ex_*.png xy_Ex.gif &
convert -delay $Delay -loop $Loop $1/xy_Ey_*.png xy_Ey.gif &

wait
