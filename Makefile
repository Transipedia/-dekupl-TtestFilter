CC=g++
CFLAGS=-g -Wall -O2 -Wno-unused-function
HEADERS=kstring.h
OBJECTS=
LIBS=-lz -lm

all:TtestFilter

TtestFilter: TtestFilter.c $(HEADERS) $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $< -o $@ $(LIBS)

test: test.c
	g++ $< -o $@
