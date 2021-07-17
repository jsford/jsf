CC=gcc
CFLAGS=-I../ -O3
DEPS = jsf_dubins.h
LIBS = -lm

all: dubins-example

dubins-example:
	$(CC) dubins-example.c $(CFLAGS) $(LIBS) -o dubins

.PHONY: clean
clean:
	rm dubins
