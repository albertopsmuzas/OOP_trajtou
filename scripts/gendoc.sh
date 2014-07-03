#!/bin/bash
if [ $OOP_TRAJTOU_PATH ];then
	echo OOP_TRAJTOU_PATH is an environmental variable
else
	echo OOP_TRAJTOU_PATH isn\'t an environmental variable
	exit 1
fi
# Erase previous documentation
rm -r $OOP_TRAJTOU_PATH/doc/latex
if [ "$?" == 0 ];then echo LaTeX old documentation was erased;fi
rm -r $OOP_TRAJTOU_PATH/doc/html
if [ "$?" == 0 ];then echo HTML old documentation was erased;fi
rm -r $OOP_TRAJTOU_PATH/doc/man
if [ "$?" == 0 ];then echo Manual old documentation was erased;fi
# If this was executed in OOP_TRAJTOU_PATH, generate documentation
if [ "$PWD" != $OOP_TRAJTOU_PATH ];then
	echo This script should be executed in $OOP_TRAJTOU_PATH
	echo It was executed here: $PWD
	exit 1
fi
conf=$OOP_TRAJTOU_PATH/doc/doxygenconf/all.conf
echo Cobnfiguration file in $conf
doxygen $conf
cd $OOP_TRAJTOU_PATH/doc/latex
make pdf
mv refman.pdf ../.
echo LaTeX documentation created inside $OOP_TRAJTOU_PATH/doc/latex
echo MANual documentation created inside $OOP_TRAJTOU_PATH/doc/man
echo HTML documentation created inside $OOP_TRAJTOU_PATH/doc/html
echo PDF documentation generated in $OOP_TRAJTOU_PATH/doc/refman.pdf
echo log file generated inside $OOP_TRAJTOU_PATH/doc/logs/error.log
echo 
#
