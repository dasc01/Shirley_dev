# Makefile for converters to UPF format

include ../make.sys

TLDEPS = mods libiotk libs

MODULES = \
../Modules/io_files.o \
../Modules/io_global.o

MODS = \
upf.o \
splines.o \
pointers.o \
fileio.o \
atomic_waves.o \
$(MODULES)

ADDLOCALPROJ_OBJ= \
addlocalproj_upf.o \
intradialprod.o

VALENCECOREPAIR_OBJ= \
valence_corerepair.o \
renormalize_ae_from_ps.o \
momentum_matrix.o

COREVALENCE_POSITION_OBJ= \
corevalence_position.o \
position_matrix.o

all : tldeps addlocalproj_upf.x valence_corerepair.x corevalence_position.x

addlocalproj_upf.x : $(ADDLOCALPROJ_OBJ) $(MODS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		$(ADDLOCALPROJ_OBJ) $(MODS) $(LIBOBJS)

valence_corerepair.x : $(VALENCECOREPAIR_OBJ) $(MODS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		$(VALENCECOREPAIR_OBJ) $(MODS) $(LIBOBJS)

corevalence_position.x : $(COREVALENCE_POSITION_OBJ) $(MODS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		$(COREVALENCE_POSITION_OBJ) $(MODS) $(LIBOBJS)

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f  *.x *.o *~ *.F90 *.mod *.d *.i work.pc

include make.depend
