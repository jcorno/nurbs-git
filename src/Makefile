mex=mkoctfile --mex 

all: bspeval.mex bspderiv.mex bspkntins.mex bspdegelev.mex 

bspeval.mex: bspeval.c basisfun.c findspan.c mexmat.c
	$(mex) bspeval.c basisfun.c findspan.c mexmat.c

bspderiv.mex: bspderiv.c mexmat.c
	$(mex) bspderiv.c mexmat.c

bspkntins.mex: bspkntins.c findspan.c mexmat.c
	$(mex) bspkntins.c findspan.c mexmat.c

bspdegelev.mex: bspdegelev.c mexmat.c
	$(mex) bspdegelev.c mexmat.c

clean:
	rm -f bspeval.mex bspeval.o bspderiv.mex \
	bspderiv.o bspkntins.mex bspkntins.o bspdegelev.mex bspdegelev.o

