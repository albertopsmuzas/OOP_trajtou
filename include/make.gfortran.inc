export F = gfortran
export CF = $(F) -cpp -J $(OOPTRAJTOUPATH)/mod -I $(OOPTRAJTOUPATH)/mod -c -Wall -Wtabs -O2 -ffree-line-length-none -fcheck=all -g
export LF = $(F) -cpp -I $(OOPTRAJTOUPATH)/mod -Wall -Wtabs -O2 -ffree-line-length-none -fcheck=all -g
