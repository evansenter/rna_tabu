
CFLAGS =  -O2 -g
LFLAGS = -O2 -lm
OBJ = backtrack.o exhaustive2.o fold_vars.o  RNAlocopt.o     read_epars.o convert_Vienna.o  exhaustive.o   HP.o         pair_mat.o  RNAeval.o energy_par.o      fold.o         IL.o         params.o    utils.o foldMod.o get_barrier_tabu.o
EXE=exhaustive exhaustive2 RNAlocopt RNAeval get_barrier_tabu

all: RNAeval exhaustive exhaustive2 RNAlocopt get_barrier_tabu

get_barrier_tabu:get_barrier_tabu.o fold_vars.o foldMod.o utils.o energy_par.o \
 read_epars.o params.o pair_mat.o
	cc $(LFLAGS) -o get_barrier_tabu get_barrier_tabu.o fold_vars.o foldMod.o utils.o energy_par.o \
	read_epars.o params.o pair_mat.o


RNAeval: RNAeval.o fold_vars.o fold.o utils.o energy_par.o \
 read_epars.o params.o pair_mat.o
	cc $(LFLAGS) -o RNAeval RNAeval.o fold_vars.o fold.o utils.o \
	energy_par.o \
	read_epars.o params.o pair_mat.o

RNAlocopt: RNAlocopt.o IL.o convert_Vienna.o fold_vars.o fold.o utils.o \
	energy_par.o read_epars.o params.o pair_mat.o HP.o backtrack.o
	cc $(LFLAGS) -o RNAlocopt RNAlocopt.o backtrack.o \
	convert_Vienna.o IL.o HP.o fold_vars.o fold.o \
	utils.o energy_par.o read_epars.o params.o pair_mat.o 

RNAlocopt.o: RNAlocopt.c fold.h RNAlocopt.h
	cc $(CFLAGS) -c RNAlocopt.c 

backtrack.o: backtrack.c backtrack.h IL.h HP.h RNAlocopt.h convert_Vienna.h
	cc $(CFLAGS) -c backtrack.c

convert_Vienna.o: convert_Vienna.c convert_Vienna.h fold.h energy_const.h \
	fold_vars.h pair_mat.h params.h 
	cc $(CFLAGS) -c convert_Vienna.c

HP.o: HP.c HP.h convert_Vienna.h pair_mat.h
	cc $(CFLAGS) -c HP.c

IL.o: IL.c IL.h
	cc $(CFLAGS) -c IL.c

exhaustive2: exhaustive2.o fold_vars.o fold.o utils.o energy_par.o \
 read_epars.o params.o pair_mat.o
	cc $(LFLAGS) -o exhaustive2 exhaustive2.o fold_vars.o fold.o utils.o energy_par.o \
	read_epars.o params.o pair_mat.o

exhaustive2.o: exhaustive2.c fold.h
	cc $(CFLAGS) -c exhaustive2.c

exhaustive: exhaustive.o fold_vars.o fold.o utils.o energy_par.o \
 read_epars.o params.o pair_mat.o
	cc $(LFLAGS) -o exhaustive exhaustive.o fold_vars.o fold.o utils.o energy_par.o \
	read_epars.o params.o pair_mat.o

exhaustive.o: exhaustive.c fold.h
	cc $(CFLAGS) -c exhaustive.c

RNAeval.o: RNAeval.c fold_vars.h fold.h utils.h
	cc $(CFLAGS) -c RNAeval.c

pair_mat.o: pair_mat.c pair_mat.h fold_vars.h
	cc $(CFLAGS) -c pair_mat.c

fold_vars.o: fold_vars.c fold_vars.h
	cc $(CFLAGS) -c fold_vars.c

get_barrier_tabu.o: get_barrier_tabu.c get_barrier_tabu.h utils.h
	cc $(CFLAGS) -c get_barrier_tabu.c

foldMod.o: foldMod.c utils.h energy_par.h fold_vars.h pair_mat.h params.h
	cc $(CFLAGS) -c foldMod.c

fold.o: fold.c utils.h energy_par.h fold_vars.h pair_mat.h params.h
	cc $(CFLAGS) -c fold.c

utils.o: utils.c config.h
	cc $(CFLAGS) -c utils.c

energy_par.o: energy_par.c energy_const.h intloops.h
	cc $(CFLAGS) -c energy_par.c

read_epars.o: read_epars.c utils.h energy_const.h energy_par.h
	cc $(CFLAGS) -c read_epars.c

params.o: params.c config.h energy_par.h fold_vars.h utils.h params.h
	cc $(CFLAGS) -c params.c

clean:	
	rm -f *~ $(OBJ) a.out get_barrier_tabu
small:
	make clean
	rm $(EXE)