# Makefile
# GNU makefile for example Cahn-Hilliard benchmark code
# Questions/comments to trevor.keller@gmail.com (Trevor Keller)

# includes
incdir = ../mmsp/include

# compilers/flags
compiler = g++
flags = -O3 -Wall -I $(incdir)
pcompiler = mpic++
pflags = -O3 -Wall -I $(incdir) -include mpi.h

# dependencies
core = ../main.cpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp

# the program
cahn-hilliard: cahn-hilliard.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

parallel: cahn-hilliard.cpp $(core)
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f cahn-hilliard parallel
