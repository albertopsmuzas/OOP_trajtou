#!/bin/bash
#//////////////////////////////////////////////////////////////////////////
# CUSTOMIZABLE BLOCK
#//////////////////////////////////////////////////////////////////////////
# Add commads to invoke the correct compiler. It should exist
# make.inc.<compilerName> and src/make.inc.<compilerName> with the correct
# compile specifications
selections="gfortran ifort"
defaultPath=${PWD}
#//////////////////////////////////////////////////////////////////////////
# NOT SO CUSTOMIZABLE BLOCK
#//////////////////////////////////////////////////////////////////////////
if [[ -f version.txt ]];then 
	echo version.txt file exists, which means that the program may be already compiled.
	echo remove it and re-run install.sh script, but you may be breaking an already installed version of the program
	exit 1
fi
# Chose debugging option ---------
echo "Install debug version?"
select debugmode in no yes Abort;do
	case $debugmode in
		yes)
			debugmode=on
			break;;
		no)
			debugmode=off
			break;;
		Abort)
			echo Installation aborted.
			exit 2;;
		*)
			echo Incorrect option. Chose another time:;;
	esac
done
# Control -------
if [[ $debugmode == "" ]];then echo For some unexpected reason, debugmode is empty. Re-run install.sh script;exit 1;fi
# Chose compiler and create convenient make.inc and src/make.inc files -------
echo "Select a compiler:"
select compiler in $selections Abort;do
	case $compiler in
		Abort)
			echo Installation aborted
			exit 2;;
		*)
			echo Checking $compiler existence...
			which $compiler
			if [[ $? != 0 ]];then echo $compiler is not in your path;exit 1;fi
			if [[ ! -f make.inc.$compiler ]];then echo make.inc.$compiler does not exist;exit 1;fi
			if [[ ! -f src/make.inc.$compiler ]];then echo src/make.inc.$compiler does not exist;exit 1;fi
			echo Adding debugmode flag to make.inc files:
			echo "# Debug option" > src/make.inc
			echo debugmode=$debugmode >> src/make.inc
			echo Creating new make.inc file
			ln -s make.inc.$compiler make.inc
			echo Creating a new src/make.inc file
			cat src/make.inc.$compiler >> src/make.inc
			break;;
		"")
			echo Incorrect option. Chose another time:;;
		esac
done
# Chose a version name or not
echo Select a name for this version \(flavour\). This may be convenient if you have different versions of the
echo program compiled with the same options, e.g., when a new PES has been implemented, a new branch is being tested, etc.
echo If you choose custom, you will have to type a version alias. Avoid special or blanck characters.
echo Choose your flavour:
select flavour in default custom Abort;do
	case $flavour in
		default)
			echo Flavour is default.
			flavour=default
			break;;
		custom)
			echo Type custom flavour \(avoid special or blanck characters\)
			read flavour
			echo Custom flavour is: $flavour
			break;;	
		Abort)
			echo Installation aborted
			exit 2;;
		*)
	esac
done
# Chose default path for PES files
echo Select a default path for PES files:
select pesPath in default custom Abort;do
	case $pesPath in
		default)
			pesPath='~/lib/trajtouPES'
			echo Default path is $pesPath
			break;;
		custom)
			echo Type custom path. It should exists or an error will occur.:
			read pesPath
			echo Custom default PES path is: $pesPath
			break;;
		Abort)
			echo Installation aborted
			exit 2;;
		*)
	esac
done
# Store installed version ---------------
versionFlag=debug_${debugmode}_compiler_${compiler}_flavour_${flavour}
# Check environmental modules -----------
echo Checking module existence:
type module
if [[ $? == 0 ]];then 
	echo WARNING: environmental modules found. Proceeding to create a private module.
	module=yes
else
	echo WARNING: environmental modules was not found. Source file will be created instead.
	module=no
fi 
# Install private module or source file --------
case $module in
	yes)
		filename=$HOME/privatemodules/ooptrajtou/${versionFlag}
		if [[ ! -f $HOME/privatemodules ]];then
			echo $HOME/privatemodules was not found. A new folder will be created.
			echo Remember to load use.own module to access modules in this folder.
			mkdir $HOME/privatemodules
		fi
		mkdir $HOME/privatemodules/ooptrajtou
		echo A new private module will be installed at $filename
#------- creating file -----------------------
cat << EOF > $filename
#%Module 1.0
# Module for OOPtrajtou program
#
append-path    PATH             $defaultPath/bin
append-path    PATH             $defaultPath/scripts
append-path    LD_LIBRARY_PATH  $defaultPath/lib
append-path		MANPATH          $defaultPath/doc/man
setenv         OOPTRAJTOUPATH   $defaultPath
setenv			OOPTRAJTOUPES	  $pesPath
setenv			LUA_PATH			  $defaultPath/lib/lua/?.lua
module-whatis  "Sets OOPTRAJTOUPATH environmental variable and adds paths so that OOPTRAJTOU libraries and binaries can be used"
EOF
#------- end file -----------------------------
		;;
	no)
	filename=${HOME}/.OOPtrajtou_${versionFlag}_rc
	echo A new source file will be created at $filename
#------- creating file -----------------------
cat << EOF > $filename
#!/bin/bash
# Source file to use OOPtrajtou program
export OOPTRAJTOUPATH=$defaultPath
export PATH=$PATH:$OOPTRAJTOUPATH/bin:$OOPTRAJTOUPATH/scripts
export OOPTRAJTOUPES=$pesPath
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OOPTRAJTOUPATH/lib
export LUA_PATH=$LUA_PATH:$defaultPath/lib/lua/?.lua
EOF
#------- end file -----------------------------
		echo In order to use source file, type:
		echo source $filename
		echo or:
		echo . $filename
		echo To source this file in every shell, add one of the previous commands to your $HOME/.bashrc file
		;;	
esac
# Store installed version info -------------
echo Version: $versionFlag > version.txt
echo Source: $filename >> version.txt
#/////// Install complements ////////////////
# Install harald-aotus lib
cd src/utils/haraldkl-aotus-707f4da5293b
./waf configure build
cd build
mv *.mod $defaultPath/mod/.
mv libflu.a $defaultPath/lib/.
mv liblualib.a $defaultPath/lib/.
mv libaotus.a $defaultPath/lib/.
mv lua $defaultPath/bin
cd $defaultPath
# Compiling section ---------------------
echo Installing src objects:
cd src; make;cd ..
echo Link executable jobs with library:
make
echo OOPtrajtou was installed successfully!!!
echo Enjoy your time!
