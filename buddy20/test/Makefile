# ==============================================================
# Makefile for SubDDD test1
# ==============================================================

CFLAGS =  -W -Wtraditional -Wmissing-prototypes -Wall

LIBDIR = ../src

INCDIR = ../src


# --- full object list
OBJ =  main.o 
CCFILES =  main.cc


# --------------------------------------------------------------
# Code generation
# --------------------------------------------------------------

.SUFFIXES: .cc 

.cc.o:
	g++ -I$(INCDIR) -c -g $<



# --------------------------------------------------------------
# The primary targets.
# --------------------------------------------------------------

main:	$(OBJ)
	g++ $(CFLAGS) -o main $(OBJ) -L$(LIBDIR) -lbdd -lfl -lm
	chmod u+x main

clean:
	rm -f *.o core *~
	rm -f main




