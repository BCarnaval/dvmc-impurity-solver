#this is a very basic makefile.
#important: must use -fPIC option for so libraries.

all:
	gcc src/dvmc_gnuplot.c -o dvmc_gnuplot
	gcc -O3 -fPIC -shared -o libdvmc_speedup.so src/dvmc_speedup.c 

