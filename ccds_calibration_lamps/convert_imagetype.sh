#!/bin/sh

########################
## using “convert” from ImageMagick to do ps convert into PNG
#########################

echo PS convert to PNG, please wait the process
for INP in *.ps
do
newname=`basename $INP .ps`
convert -density 150 -geometry 100% $INP $newname.png
echo ” convert $INP to $newname.png completely”
done
echo ” process ended, please check your graphical files”
