#!/bin/bash
env | grep OOPTRAJTOUPATH >/dev/null
if [[ $? != 0 ]]; then echo 'OOPTRAJTOUPATH is not an environmental variable';exit 1;fi

echo $* | grep .x >/dev/null
if [[ $? != 0 ]]; then
	echo 'Type the fortran file name, substituting the fortran extension by ".x"'
	echo 'Example:'
	echo 'To compile "test.f90" fortran file, you must type: trajtouApplyOn.sh test.x'
	exit 1
else
	make -f ${OOPTRAJTOUPATH}/MakeExternalFile $*
fi
