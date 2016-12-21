# Makefile
# GNU makefile for Cahn-Hilliard benchmark code
# Questions/comments to trevor.keller@gmail.com (Trevor Keller)

# includes
incdir = ../mmsp/include

# compilers/flags
compiler = g++
flags = -O3 -Wall -I $(incdir)
pcompiler = /usr/bin/mpic++.openmpi
pflags = -O3 -Wall -I $(incdir) -include mpi.h

# the program
cahn-hilliard: cahn-hilliard.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: cahn-hilliard.cpp
	$(pcompiler) $(pflags) $< -o $@ -lz

clean:
	rm -f cahn-hilliard parallel
