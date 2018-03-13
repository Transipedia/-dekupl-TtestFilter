CC=g++
CFLAGS=-g -Wall -O2 -Wno-unused-function
HEADERS=kstring.h
OBJECTS=
LIBS=-lz -lm

.PHONY: test

all:TtestFilter

TtestFilter: TtestFilter.c $(HEADERS) $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $< -o $@ $(LIBS)

test:
	echo "T-test"
	./TtestFilter -p 0.2 test/bigger-counts.tsv.gz test/sample_conditions_full.tsv A B
	echo "Poisson"
	./TtestFilter -p 0.2 -l test/bigger-counts.tsv.gz test/sample_conditions_full.tsv A B
