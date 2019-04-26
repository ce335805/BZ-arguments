#makefile to build BZ-arguments library
CC = mpiicpc
CPPFLAGS =
libdir = ./
incdir = ./
srcdir = ./


OUTFILE = $(libdir)libBZ.a
#.PHONY all
all: $(OUTFILE)

FF.o: $(srcdir)Index.cpp $(incdir)Index.h
	$(CC) $(CPPFLAGS) -c $< -o $@ -I$(incdir) 

util.o: $(srcdir)util.cpp $(incdir)util.h $(incdir)globals.h
	$(CC) $(CPPFLAGS) -c $< -o $@ -I$(incdir) 

$(OUTFILE): FF.o util.o 
	ar r $@ $^


