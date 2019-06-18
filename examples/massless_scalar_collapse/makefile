CC = gcc -I/home/jripley/one-dim-amr

CFLAGS = -g -fmax-errors=5 -Wall -O3

LIBS = -lm -lamr -lbbhutil

LDFLAGS = -L/home/jripley/one-dim-amr

OBJECTS = main.o \
	evolution_routines_GR.o \
	evolution_routines_EdGB.o \
	free_initial_data.o \
	diagnostics_GR.o \
	diagnostics_EdGB.o \
	diagnostics_general.o \
	file_io.o \
	stencils.o \
	basic_matrix_computations.o

DEPS_MAIN = evolution_routines_GR.h \
	evolution_routines_EdGB.h \
	diagnostics_GR.h \
	diagnostics_EdGB.h \
	free_initial_data.h \
	file_io.h

DEPS_EVOLUTION_ROUTINES_GR = evolution_routines_GR.h \
	stencils.h

DEPS_EVOLUTION_ROUTINES_EDGB = evolution_routines_EdGB.h \
	stencils.h

DEPS_FREE_ID = free_initial_data.h \
	stencils.h

DEPS_DIAGNOSTICS_GENERAL = diagnostics_general.h \
	stencils.h

DEPS_DIAGNOSTICS_GR = diagnostics_GR.h \
	stencils.h

DEPS_DIAGNOSTICS_EDGB = diagnostics_EdGB.h \
	stencils.h \
	basic_matrix_computations.h

DEPS_BASIC_MATRIX_COMPUTATIONS = basic_matrix_computations.h 

DEPS_FILE_IO = file_io.h

DEPS_STENCILS = stencils.h 

sim: $(OBJECTS)
	$(CC) -o sim $(OBJECTS) $(LIBS) $(LDFLAGS) $(CFLAGS)

compileTest: $(OBJECTS)
	$(CC) -o compileTest $(OBJECTS) $(LIBS) $(LDFLAGS) $(CFLAGS)

main.o: main.c $(DEPS_MAIN)
	$(CC) -c main.c $(LIBS) $(LDFLAGS) $(CFLAGS)

evolution_routines_GR.o: evolution_routines_GR.c $(DEPS_EVOLUTION_ROUTINES_GR)
	$(CC) -c evolution_routines_GR.c $(LIBS) $(LDFLAGS) $(CFLAGS)

evolution_routines_EDGB.o: evolution_routines_EdGB.c $(DEPS_EVOLUTION_ROUTINES_EDGB)
	$(CC) -c evolution_routines_EdGB.c $(LIBS) $(LDFLAGS) $(CFLAGS)

free_initial_data.o: free_initial_data.c $(DEPS_FREE_ID)
	$(CC) -c free_initial_data.c $(LIBS) $(LDFLAGS) $(CFLAGS)

diagnostics_general.o: diagnostics_general.c $(DEPS_DIAGNOSTICS_GENERAL)
	$(CC) -c diagnostics_general.c $(LIBS) $(LDFLAGS) $(CFLAGS)

diagnostics_GR.o: diagnostics_GR.c $(DEPS_DIAGNOSTICS_GR)
	$(CC) -c diagnostics_GR.c $(LIBS) $(LDFLAGS) $(CFLAGS)

diagnostics_EdGB.o: diagnostics_EdGB.c $(DEPS_DIAGNOSTICS_EDGB)
	$(CC) -c diagnostics_EdGB.c $(LIBS) $(LDFLAGS) $(CFLAGS)

basic_matrix_computations.o: basic_matrix_computations.c $(DEPS_BASIC_MATRIX_COMPUTATIONS)
	$(CC) -c basic_matrix_computations.c $(LIBS) $(LDFLAGS) $(CFLAGS)

stencils.o: stencils.c $(DEPS_STENCILS)
	$(CC) -c stencils.c $(LIBS) $(LDFLAGS) $(CFLAGS)

file_io.o: file_io.c $(DEPS_FILE_IO)
	$(CC) -c file_io.c $(LIBS) $(LDFLAGS) $(CFLAGS)

clean_o:
	rm *.o
