CFLAGS += -Wall -std=c99
LDLIBS = -llapack -lblas -lm

qr-omp_2 : CFLAGS += -fopenmp

all : qr-omp_2
