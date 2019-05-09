# link together to make .a file

CC=gcc

CFLAGS= -Wall -fmax-errors=5 -std=gnu11 -lm

DEPS_GRID_HIER = grid_hierarchy.h

LIB_OBJECTS = grid_hierarchy.o

grid_hierarchy.o: grid_hierarchy.c $(DEPS_GRID_HIER)
	$(CC) -c grid_hierarchy.c $(CFLAGS)

libamr.a: $(OBJECTS)
	ar rcs libamr.a $(LIB_OBJECTS)

clean_o:
	rm *.o 
