# link together to make a static linked library

CC=gcc

CFLAGS= -Wall -fmax-errors=5 -std=gnu11 -lm

DEPS_GRID_HIER = amr_grid_hierarchy.h

DEPS_EVOLVE = amr_evolve.h amr_grid_hierarchy.h

LIB_OBJECTS = amr_evolve.o amr_grid_hierarchy.o

grid_hierarchy.o: amr_grid_hierarchy.c $(DEPS_GRID_HIER)
	$(CC) -c grid_hierarchy.c $(CFLAGS)

libamr.a: $(OBJECTS)
	ar rcs libamr.a $(LIB_OBJECTS)

clean_o:
	rm *.o 
