Compiler Support {#compiler_support}
=========================

Aotus relies on the ISO-C-Binding of the Fortran 2003 standard, and therefore
requires relatively new compiler versions.
However a wide range of different compilers are known to compile the library.

## Cray
Known to work with version 7.4.0 and newer.

## GNU gfortran
Known to work with version 4.4.5 and newer.

## IBM
Known to work with version 11.1.

## Intel
Known to work with version 12.0 and newer.

## NAG
Known to work with version 5.2.

## PGI
Does not correctly resolve the generic aot_get_val interface!
