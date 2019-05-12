CC = gcc -I/home/jripley/one-dim-amr

CFLAGS = -g -fmax-errors=5 -Wall -O1

LIBS = -lm -lamr

LDFLAGS = -L/home/jripley/one-dim-amr

OBJECTS = main.o \
	evolution_routines.o 

DEPS_MAIN = evolution_routines.h

DEPS_ER = 

run: $(OBJECTS)
	$(CC) -o run $(OBJECTS) $(LIBS) $(LDFLAGS) $(CFLAGS)

main.o: main.c $(DEPS_MAIN)
	$(CC) -c main.c $(LIBS) $(LDFLAGS) $(CFLAGS)

evolution_routines.o: evolution_routines.c $(DEPS_ER)
	$(CC) -c evolution_routines.c $(LIBS) $(LDFLAGS) $(CFLAGS)

clean_o:
	rm *.o
