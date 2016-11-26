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
			ln -s debugVersion debugVersion
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
			if [[ ! -f include/make.${compiler}.inc ]];then echo include/make.${compiler}.inc does not exist;exit 1;fi
			echo Creating new make.inc file
			ln -s include/make.${compiler}.inc make.compiler.inc
			
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
		if [[ ! -d $HOME/privatemodules ]];then
			echo $HOME/privatemodules was not found. A new folder will be created.
			echo Remember to load use.own module to access modules in this folder.
			mkdir $HOME/privatemodules
		fi
		if [[ ! -d $HOME/privatemodules/ooptrajtou ]];then
			echo $HOME/privatemodules/ooptrajtou was not found. A new folder will be created.
			mkdir $HOME/privatemodules/ooptrajtou
			echo A new private module will be installed at $filename
		fi
#------- creating file -----------------------
cat << EOF > $filename
#%Module 1.0
# Module for OOPtrajtou program
#
set                       INSTALLPATH       $defaultPath
set                       PESPATH           $pesPath
prepend-path --delim ":"  PATH              \$INSTALLPATH/bin
prepend-path --delim ":"  PATH              \$INSTALLPATH/scripts
prepend-path --delim ":"  LD_LIBRARY_PATH   \$INSTALLPATH/lib
prepend-path --delim ":"  MANPATH           \$INSTALLPATH/doc/man
prepend-path --delim ";"  LUA_PATH          \$INSTALLPATH/lib/lua/?.lua
setenv         OOPTRAJTOUPATH     \$INSTALLPATH
setenv         OOPTRAJTOUPES      \$PESPATH
setenv         OOPTRAJTOU_DEBUG   $debugmode
setenv         OOPTRAJTOUCOMPILER $compiler
module-whatis  "OOP_trajtou program installed"
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
export PATH=$PATH:$defaultPath/bin
export PATH=$PATH:$defaultPath/scripts
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$defaultPath/lib
export MANPATH=$MANPATH:$defaultPath/doc/man
export OOPTRAJTOUPATH=$defaultPath
export OOPTRAJTOUPES=$pesPath
export OOPTRAJTOU_DEBUG=$debugmode
export OOPTRAJTOUCOMPILER=$compiler
export LUA_PATH=$LUA_PATH;$defaultPath/lib/lua/?.lua
EOF
#------- end file -----------------------------
		echo In order to use source file, type:
		echo source $filename
		echo or:
		echo . $filename
		echo To source this file in every shell, add one of the previous commands to your $HOME/.bashrc file
		;;	
esac
#/////// Install complements ////////////////
make "OOPTRAJTOUPATH=$defaultPath" "OOPTRAJTOU_DEBUG=$debugmode" trajtouJobs
if [[ $? == 0 ]]; then
	echo Version: $versionFlag > version.txt
	echo Source: $filename >> version.txt
	echo OOPtrajtou was installed successfully!!!
	echo Enjoy your time!
else
	echo There may be errors in the compilation. Good Luck
fi
