CC = gcc
CFLAGS = -W -Wall -lm
LDFLAGS = -lm
OBJ = print.o matgen.o sgetrf_nopiv.o convert.o blas.o gmres.o timer.o
NVCC = nvcc
NVCFLAGS =
LDLIBS = -lcublas -lcuda

%.o: %.cu
		$(NVCC) -c $< -o $@ $(NVCFLAGS)

%.o: %.c
		$(CC) -c -o $@ $< $(CFLAGS)


driver: $(OBJ) driver.c
		$(NVCC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

all: driver

clean:
	rm *.o driver
