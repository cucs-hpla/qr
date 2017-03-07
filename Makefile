CFLAGS += -Wall -std=c99 -O3 -march=native -g
LDLIBS = -llapack -lblas -lm

qr-omp : CFLAGS += -fopenmp

all : qr-omp

clean:
	rm -rf qr-omp
