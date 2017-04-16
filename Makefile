CC=mpicc
FLAGS=-O3
EXECS=mpi_solved1 mpi_solved2 mpi_solved3 mpi_solved4 mpi_solved5 mpi_solved6 mpi_solved7 jacobi2D-mpi ssort-mpi

all: ${EXECS}

mpi_solved1: mpi_solved1.c
	${CC} ${FLAGS} $^ -o mpi_solved1

mpi_solved2: mpi_solved2.c
	${CC} ${FLAGS} $^ -o mpi_solved2

mpi_solved3: mpi_solved3.c
	${CC} ${FLAGS} $^ -o mpi_solved3

mpi_solved4: mpi_solved4.c
	${CC} ${FLAGS} $^ -o mpi_solved4

mpi_solved5: mpi_solved5.c
	${CC} ${FLAGS} $^ -o mpi_solved5

mpi_solved6: mpi_solved6.c
	${CC} ${FLAGS} $^ -o mpi_solved6

mpi_solved7: mpi_solved7.c
	${CC} ${FLAGS} $^ -o mpi_solved7

jacobi2D-mpi:jacobi2D-mpi.c
	${CC} ${FLAGS} $^ -o jacobi2D-mpi

ssort-mpi:ssort-mpi.c
	${CC} ${FLAGS} $^ -o ssort-mpi


