
all: prog.out a.out b.out c.out d.out

prog.out: prog.c
	icc -o prog.out prog.c
a.out:
	ln -s prog.out a.out
b.out:
	ln -s prog.out b.out
c.out:
	ln -s prog.out c.out
d.out:
	ln -s prog.out d.out

clean:
	/bin/rm prog.out a.out b.out c.out d.out
