OCTFILES = bspderiv.oct nrb_srf_basisfun_der__.oct\
basisfun.oct bspeval.oct basisfunder.oct nrb_srf_basisfun__.oct

all: $(OCTFILES)

low_level_functions.o: low_level_functions.cc
	mkoctfile -c $<

%.oct:  %.cc low_level_functions.o
	mkoctfile $< low_level_functions.o

clean:
	-rm -f *.o core octave-core *.oct *~
