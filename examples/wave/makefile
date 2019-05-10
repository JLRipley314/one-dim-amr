CC = gcc

CFLAGS = -g -fmax-errors=5 -Wall -01

LDFLAGS = -lm -amr

OBJECTS = main.o \
	initial_data.o \
	PDE_solvers.o \

run: $(OBJECTS)
	$(CC) -o run $(OBJECTS) $(LDFLAGS) $(CFLAGS)

main.o: main.c $(DEPS_MAIN)
	$(CC) -o main.c $(LDFLAGS) $(CFLAGS)

initial_data.o: initial_data.c $(DEPS_ID)
	$(CC) -o initial_data.c $(LDFLAGS) $(CFLAGS)

PDE_solver.o: PDE_solvers.c $(DEPS_ID)
	$(CC) -o PDE_solvers.c $(LDFLAGS) $(CFLAGS)

clean_o:
	rm *.o
