CC=g++
CFLAGS=-g -Wall -O2 -Wno-unused-function
HEADERS=kstring.h
OBJECTS=$(HEADERS:.h=.o)
LIBS=-lz -lm

all:ttestFilter

ttestFilter: TtestFilter.c $(HEADERS) $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $< -o $@ $(LIBS)

test: test.c
	g++ $< -o $@
