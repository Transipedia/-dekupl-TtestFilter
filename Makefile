CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
HEADERS=
OBJECTS=$(HEADERS:.h=.o)
LIBS=-lz -lm

all:ttestFilter

ttestFilter: TtestFilter.c $(HEADERS) $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $< -o $@ $(LIBS)
