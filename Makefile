CFLAGS += -Wall -std=c99 -O3
LDLIBS = -llapack -lblas -lm

qr-omp : CFLAGS += -fopenmp

all : qr-omp
