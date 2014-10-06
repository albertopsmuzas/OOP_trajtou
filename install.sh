#!/bin/bash
#//////////////////////////////////////////////////////////////////////////
# CUSTOMIZABLE BLOCK
#//////////////////////////////////////////////////////////////////////////
# Add commads to invoke the correct compiler. It should exist
# make.inc.<compilerName> and src/make.inc.<compilerName> with the correct
# compile specifications
selections="ifort gfortran"
defaultPath=${PWD}
#//////////////////////////////////////////////////////////////////////////
# NOT SO CUSTOMIZABLE BLOCK
#//////////////////////////////////////////////////////////////////////////
# Chose debugging option ---------
echo "Install debug version?"
select debugmode in yes no Abort;do
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
# Check environmental modules -----------
echo Checking module existence:
which module
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
		filename=$HOME/privatemodules/ooptrajtou/debug_${debugmode}_compiler_${compiler}
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
module-whatis  "Sets OOPTRAJTOUPATH environmental variable and adds paths so that OOPTRAJTOU libraries and binaries can be used"
EOF
#------- end file -----------------------------
		;;
	no)
	filename=${HOME}/.OOPtrajtou_debug_${debugmode}_compiler_${compiler}_rc
	echo A new source file will be created at $filename
#------- creating file -----------------------
cat << EOF > $filename
#!/bin/bash
# Source file to use OOPtrajtou program
export OOPTRAJTOUPATH=$defaultPath
export PATH=$PATH:$OOPTRAJTOUPATH/bin:$OOPTRAJTOUPATH/scripts
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OOPTRAJTOUPATH/lib
EOF
#------- end file -----------------------------
		echo In order to use source file, type:
		echo source $filename
		echo or:
		echo . $filename
		echo To source this file in every shell, add one of the previous commands to your $HOME/.bashrc file
		;;	
esac
cd src; make
make
echo OOPtrajtou was installed successfully!!!
echo Enjoy your time!
