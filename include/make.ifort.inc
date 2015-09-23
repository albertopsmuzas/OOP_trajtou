export F = ifort
export CF = $(F) -cpp -I $(OOPTRAJTOUPATH)/mod -module $(OOPTRAJTOUPATH)/mod -gen-interfaces -c -O2 -traceback -CB -warn truncated_source -g -no-wrap-margin -DINTEL
export LF = $(F) -I $(OOPTRAJTOUPATH)/mod -O2 -g -no-wrap-margin -DINTEL
