# ==============================================================
# Makefile for BDDASTAR
# ==============================================================
CC = g++-6 -Wno-deprecated
CFLAGS =  -Wtraditional -Wmissing-prototypes

LIBDIR = ./buddy20/src

INCDIR = ./buddy20/src


# --- full object list
OBJ =     BDDASTAR.o 
         

CCFILES = BDDASTAR.cc 


# --------------------------------------------------------------
# Code generation
# --------------------------------------------------------------

.SUFFIXES: .cc .c

.cc.o:
	$(CC) -I$(INCDIR)   -g  -c  $<


# --------------------------------------------------------------
# The primary targets.
# --------------------------------------------------------------

BDDASTAR:	$(OBJ)
	$(CC)  -g   $(CFLAGS) -o BDDASTAR $(OBJ) -L$(LIBDIR) -lbdd -lm
	chmod u+x BDDASTAR

clean:

	rm -f *.o core *~
	rm -f BDDASTAR


