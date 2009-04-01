all: low_level_functions.oct

%.oct: %.cc
	mkoctfile $<

clean:
	-rm *.o core octave-core *.oct *~
