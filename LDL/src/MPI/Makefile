CC1=mpicc

all:
	$(CC1) -o mpi.out ./../tools.c mpi.c -lm
debug:
	$(CC1) -o mpi.out ./../tools.c mpi.c -lm -fsanitize=address -g

clean:
	rm *.out

