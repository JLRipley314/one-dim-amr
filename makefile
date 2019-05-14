# link together to make a static linked library

CC=gcc

CFLAGS= -g -Wall -fmax-errors=5 -std=gnu11 -lm

DEPS_GRID_HIER = amr_grid_hierarchy.h

DEPS_EVOLVE = amr_evolve.h \
	amr_grid_hierarchy.h

LIB_OBJECTS = amr_evolve.o \
	amr_grid_hierarchy.o

libamr.a: $(LIB_OBJECTS)
	ar rcs libamr.a $(LIB_OBJECTS)

amr_evolve.o: amr_evolve.c $(DEPS_EVOLVE)
	$(CC) -c amr_evolve.c $(CFLAGS)

grid_hierarchy.o: amr_grid_hierarchy.c $(DEPS_GRID_HIER)
	$(CC) -c amr_grid_hierarchy.c $(CFLAGS)


clean_o:
	rm *.o 
