CC = gcc -I/home/jripley/one-dim-amr

CFLAGS = -g -fmax-errors=5 -Wall -O2

LIBS = -lm -lamr -lbbhutil

LDFLAGS = -L/home/jripley/one-dim-amr

OBJECTS = main.o \
	evolution_routines_GR.o \
	diagnostics_GR.o \
	diagnostics_general.o \
	file_io.o \
	stereographic_routines.o

DEPS_MAIN = evolution_routines_GR.h \
	file_io.h

DEPS_EVOLUTION_ROUTINES_GR = evolution_routines_GR.h \
	stereographic_routines.h

DEPS_DIAGNOSTICS_GENERAL = diagnostics_general.h \
	stereographic_routines.h

DEPS_DIAGNOSTICS_GR = diagnostics_GR.h \
	stereographic_routines.h

DEPS_FILE_IO = file_io.h

DEPS_STEREOGRAPHIC = stereographic_routines.h 

sim: $(OBJECTS)
	$(CC) -o sim $(OBJECTS) $(LIBS) $(LDFLAGS) $(CFLAGS)

main.o: main.c $(DEPS_MAIN)
	$(CC) -c main.c $(LIBS) $(LDFLAGS) $(CFLAGS)

evolution_routines_GR.o: evolution_routines_GR.c $(DEPS_EVOLUTION_ROUTINES_GR)
	$(CC) -c evolution_routines_GR.c $(LIBS) $(LDFLAGS) $(CFLAGS)

diagnostics_general.o: diagnostics_general.c $(DEPS_DIAGNOSTICS_GENERAL)
	$(CC) -c diagnostics_general.c $(LIBS) $(LDFLAGS) $(CFLAGS)

diagnostics_GR.o: diagnostics_GR.c $(DEPS_DIAGNOSTICS_GR)
	$(CC) -c diagnostics_GR.c $(LIBS) $(LDFLAGS) $(CFLAGS)

stereographic_routines.o: stereographic_routines.c $(DEPS_STEREOGRAPHIC)
	$(CC) -c stereographic_routines.c $(LIBS) $(LDFLAGS) $(CFLAGS)

file_io.o: file_io.c $(DEPS_FILE_IO)
	$(CC) -c file_io.c $(LIBS) $(LDFLAGS) $(CFLAGS)

clean_o:
	rm *.o
clean_output:
	rm output/*
