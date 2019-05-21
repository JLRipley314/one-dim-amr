CC = gcc -I/home/jripley/one-dim-amr

CFLAGS = -g -fmax-errors=5 -Wall -O2

LIBS = -lm -lamr -lbbhutil

LDFLAGS = -L/home/jripley/one-dim-amr

OBJECTS = main.o \
	evolution_routines.o 

DEPS_MAIN = evolution_routines.h

DEPS_EVOLUTION_ROUTINES = # no deps for now 

DEPS_CHECKS_DIAGNOSTICS = # no deps for now 

run: $(OBJECTS)
	$(CC) -o run $(OBJECTS) $(LIBS) $(LDFLAGS) $(CFLAGS)

main.o: main.c $(DEPS_MAIN)
	$(CC) -c main.c $(LIBS) $(LDFLAGS) $(CFLAGS)

evolution_routines.o: evolution_routines.c $(DEPS_EVOLUTION_ROUTINES)
	$(CC) -c evolution_routines.c $(LIBS) $(LDFLAGS) $(CFLAGS)

general_checks_diagnostics.o: general_checks_diagnostics.c $(DEPS_CHECKS_DIAGNOSTICS)
	$(CC) -c general_checks_diagnostics.c $(LIBS) $(LDFLAGS) $(CFLAGS)

clean_o:
	rm *.o
clean_output:
	rm output/*
