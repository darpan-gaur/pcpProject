#################################################################
#### The first block are Macros ##################
CC       = gcc         # The C compiler to use
CFLAGS   = -g          # compilation flags
RM       = rm -rf      # a command for cleaning up
INCLUDES = -I ./       # include paths
LIBS     = -lm         # libraries required at the linking stage

OBJS     = utilities.o functions.o linsolvers.o main.o # all object files
#### End the block of Macros ####################

#### The second block defines "target - prerequisities - recipe" combinations
#### The syntax is as below:
# target : prerequisities
#	recipe      (there must be a tab before the recipe)

# default target is all
default: all

all: ${OBJS} 
	@echo "Building..."
	${CC} ${CFLAGS} ${INCLUDES} ${LIBS} $^ -o execfile

# Shorthand for lines 20--25; 
%.o: %.c
	${CC} ${CFLAGS} -c $< ${INCLUDES}

#### End the block of targets ####################

#### The last block are phony targets ##################
# targets not associated with files, but only recipes
list:
	@echo $(shell ls)

clean:
	@echo "Cleaning..."
	${RM} *.o
	${RM} execfile
#### End the block of phony targets ####################
