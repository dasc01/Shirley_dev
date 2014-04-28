.SUFFIXES:
.SUFFIXES: .f .f90 .o

F90=ifort
F90FLAGS=-O3

F77=$(F90)
F77FLAGS=$(F90FLAGS)

SRC= main.f90 hybrd1.f
OBJ1=$(SRC:.f=.o)
OBJ=$(OBJ1:.f90=.o)

default: $(OBJ)
	$(F90) -o hybrd $(OBJ)

clean:
	rm *.o hybrd

.f.o:
	$(F77) $(F77FLAGS) -c $<

.f90.o:
	$(F90) $(F90FLAGS) -c $<

