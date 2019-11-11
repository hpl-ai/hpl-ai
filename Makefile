CC = gcc
CFLAGS = -W -Wall -lm
LDFLAGS = -lm
OBJ = print.o matgen.o sgetrf_nopiv.o convert.o blas.o gmres.o timer.o

%.o: %.c
		$(CC) -c -o $@ $< $(CFLAGS)

driver: $(OBJ) driver.c
		$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

all: driver

clean:
	rm *.o driver
