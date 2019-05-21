CC = gcc -I/home/jripley/one-dim-amr

CFLAGS = -g -fmax-errors=5 -Wall -O2

LIBS = -lm -lamr -lbbhutil

LDFLAGS = -L/home/jripley/one-dim-amr

OBJECTS = main.o \
	evolution_routines_GR.o \
	checks_diagnostics_general.o \
	stereographic_routines.o

DEPS_MAIN = evolution_routines_GR.h

DEPS_EVOLUTION_ROUTINES_GR = evolution_routines_GR.h \
	stereographic_routines.h

DEPS_CHECKS_DIAGNOSTICS = checks_diagnostics_general.h \
	stereographic_routines.h

DEPS_STEREOGRAPHIC = stereographic_routines.h 

run: $(OBJECTS)
	$(CC) -o run $(OBJECTS) $(LIBS) $(LDFLAGS) $(CFLAGS)

main.o: main.c $(DEPS_MAIN)
	$(CC) -c main.c $(LIBS) $(LDFLAGS) $(CFLAGS)

evolution_routines_GR.o: evolution_routines_GR.c $(DEPS_EVOLUTION_ROUTINES_GR)
	$(CC) -c evolution_routines_GR.c $(LIBS) $(LDFLAGS) $(CFLAGS)

checks_diagnostics_general.o: checks_diagnostics_general.c $(DEPS_CHECKS_DIAGNOSTICS)
	$(CC) -c checks_diagnostics_general.c $(LIBS) $(LDFLAGS) $(CFLAGS)

stereographic_routines.o: stereographic_routines.c $(DEPS_STEREOGRAPHIC)
	$(CC) -c stereographic_routines.c $(LIBS) $(LDFLAGS) $(CFLAGS)

clean_o:
	rm *.o
clean_output:
	rm output/*
