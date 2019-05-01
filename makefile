#makefile to build BZ-arguments library
CC = mpiicpc
CPPFLAGS =
libdir = lib/
incdir = include/
srcdir = src/


OUTFILE = $(libdir)libBZ.a
#.PHONY all
all: $(OUTFILE)

Index.o: $(srcdir)Index.cpp $(incdir)Index.h
	$(CC) $(CPPFLAGS) -c $< -o $@ -I$(incdir) 

$(OUTFILE): Index.o
	ar r $@ $^


