CFLAGS += -Wall -std=c99
LDLIBS = -llapack -lblas -lm

qr-omp : CFLAGS += -fopenmp

all : qr-omp
