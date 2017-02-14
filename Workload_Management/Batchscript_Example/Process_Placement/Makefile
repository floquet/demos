
MPICC=mpicc

all: show_placement.exe

show_placement.exe: show_placement.o
	$(MPICC) -o show_placement.exe show_placement.o

show_placement.o: show_placement.c
	$(MPICC) -c show_placement.c

clean:
	/bin/rm show_placement.exe show_placement.o

