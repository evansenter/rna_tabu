CFLAGS =  -O2 -g
LFLAGS = -O2 -lm
OBJ = fold_vars.o pair_mat.o read_epars.o energy_par.o params.o utils.o foldMod.o get_barrier.o
EXE=get_barrier

all: get_barrier

get_barrier:get_barrier.o fold_vars.o foldMod.o utils.o energy_par.o \
 read_epars.o params.o pair_mat.o
	gcc $(LFLAGS) -o get_barrier get_barrier.o fold_vars.o foldMod.o utils.o energy_par.o \
	read_epars.o params.o pair_mat.o

fold_vars.o: fold_vars.c fold_vars.h
	gcc $(CFLAGS) -c fold_vars.c

get_barrier.o: get_barrier.c get_barrier.h utils.h
	gcc $(CFLAGS) -c get_barrier.c

foldMod.o: foldMod.c utils.h energy_par.h fold_vars.h pair_mat.h params.h
	gcc $(CFLAGS) -c foldMod.c

utils.o: utils.c config.h
	gcc $(CFLAGS) -c utils.c

energy_par.o: energy_par.c energy_const.h intloops.h
	gcc $(CFLAGS) -c energy_par.c
clean:	
	rm -f *~ $(OBJ)
small:
	make clean
	rm $(EXE)
