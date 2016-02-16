CC1=gcc-5
CC2=gcc

all:
	$(CC2) -o openmp.out ./../tools.c openmp.c -fopenmp -lm

linux-debug:
	$(CC2) -o openmp.out ./../tools.c openmp.c -fopenmp -lm -g

mac:
	$(CC1) -o openmp.out ./../tools.c openmp.c -fopenmp -lm

mac-debug:
	$(CC1) -o openmp.out ./../tools.c openmp.c -fopenmp -lm -g

clean:
	rm *.out

