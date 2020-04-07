#==========================================================================
TOP= ./
INC= $(TOP)include/
SRC= $(TOP)src/
BIN= $(TOP)bin/
OBJDIR= $(TOP)obj/

vpath %.c $(SRC)
vpath %.h $(INC)
vpath %.o $(OBJDIR)
#==========================================================================
CC=gcc
CFLAGS= -Wall -Wextra -fmax-errors=5 -std=gnu11 -lm -O2 -g  
SYSLIB= -lm  
#==========================================================================
OBJ= $(addprefix $(OBJDIR), \
	amr_evolve.o \
	amr_grid_hierarchy.o \
	)
DEPS= amr_evolve.h amr_grid_hierarchy.h 
#==========================================================================
LIBAMR= $(BIN)libamr.a
.PHONY: all
all: $(LIBAMR)
#==========================================================================
LIBAMR: $(OBJ)
	ar rcs $@ $^  
#==========================================================================
$(BIN)libamr.a: $(OBJ)
	ar rcs $@ $^ 
#==========================================================================
$(OBJDIR)%.o: $(SRC)%.c $(DEPS)
	$(CC) $(CFLAGS) -I$(INC) $(SYSLIB) -c -o $@ $< 
#==========================================================================
.PHONY: clean_o
clean_o:
	@rm $(OBJDIR)*.o 
