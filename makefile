SHELL = /bin/sh

CFLAGS = -mcmodel=medium
ALL_CFLAGS = -std=c11 $(CFLAGS)

mesosim: Mesosim.o Atoms.o FileIO.o Random.o Simulation\ Aux.o Simulation.o Vector.o
	gcc $(ALL_CFLAGS) -v -o mesosim Mesosim.o Atoms.o FileIO.o Random.o Simulation\ Aux.o Simulation.o Vector.o -lm

Mesosim.o: Mesosim.c stdafx.h.gch Defs.h Geometry.h Prototypes.h Vector.h Globals.h Simulation\ Globals.h
	gcc $(ALL_CFLAGS) -c Mesosim.c

Atoms.o: Atoms.c stdafx.h.gch Defs.h Geometry.h Prototypes.h Vector.h Global\ Externs.h Simulation\ Global\ Externs.h
	gcc $(ALL_CFLAGS) -c Atoms.c

FileIO.o: FileIO.c stdafx.h.gch Defs.h Geometry.h Prototypes.h Vector.h Global\ Externs.h Simulation\ Global\ Externs.h
	gcc $(ALL_CFLAGS) -c FileIO.c

Random.o: Random.c stdafx.h.gch Defs.h Geometry.h Prototypes.h Vector.h Global\ Externs.h Simulation\ Global\ Externs.h
	gcc $(ALL_CFLAGS) -c Random.c

Simulation\ Aux.o: Simulation\ Aux.c stdafx.h.gch Defs.h Geometry.h Prototypes.h Vector.h Global\ Externs.h Simulation\ Global\ Externs.h
	gcc $(ALL_CFLAGS) -c Simulation\ Aux.c

Simulation.o: Simulation.c stdafx.h.gch Defs.h Geometry.h Prototypes.h Vector.h Global\ Externs.h Simulation\ Global\ Externs.h
	gcc $(ALL_CFLAGS) -c Simulation.c

Vector.o: Vector.c stdafx.h.gch Vector.h
	gcc $(ALL_CFLAGS) -c Vector.c

stdafx.h.gch: stdafx.h
	gcc -x c-header stdafx.h -o stdafx.h.gch

clean:
	rm *.o *.i *.s mesosim.exe
