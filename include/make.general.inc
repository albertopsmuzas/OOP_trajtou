vpath %.o $(OOPTRAJTOUPATH)/bin
vpath %.a $(OOPTRAJTOUPATH)/lib
vpath %.x $(OOPTRAJTOUPATH)/bin
ifeq ($(OOPTRAJTOU_DEBUGMODE),on)
	CF = $(CF) -DDEBUG
	LF = $(LF) -DDEBUG
endif
%.o : %.f90
	$(CF) $< -o ../bin/$@
%.o : %.f
	$(CF) $< -o ../bin/$@
