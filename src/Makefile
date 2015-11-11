# Project: GenomeTester v4.0

VERSION = 4.0
CXX  = gcc

# C Files

LISTMAKER_SOURCES = \
	glistmaker.c \
	wordtable.c \
	wordmap.c \
	fasta.c \
	buffer.c \
	sequence.c \
	queue.c queue.h \
	common.c \
	utils.c

LISTMAKER2_SOURCES = \
	glistmaker2.c \
	wordtable.c \
	wordmap.c \
	fasta.c \
	buffer.c \
	sequence.c \
	queue.c \
	common.c \
	utils.c \
	trie.c

LISTQUERY_SOURCES = \
	glistquery.c \
	wordtable.c \
	wordmap.c \
	fasta.c \
	buffer.c \
	sequence.c \
	common.c \
	utils.c

LISTCOMPARE_SOURCES = \
	glistcompare.c \
	wordtable.c \
	wordmap.c \
	fasta.c \
	buffer.c \
	sequence.c \
	common.c \
	utils.c
	
RELEASEFLAGS = -O3
DEBUGFLAGS = -O0 -g
LIBS = -lm -lpthread
INCS = -I.
BINS  = glistmaker glistquery glistcompare

#CXXFLAGS = $(INCS) $(DEBUGFLAGS) -Wall 
CXXFLAGS = $(INCS) $(RELEASEFLAGS) -Wall 

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BINS) all-after

glistmaker: $(LISTMAKER_SOURCES)
	$(CXX) $(LISTMAKER_SOURCES) -o glistmaker $(LIBS) $(CXXFLAGS)

glistmaker2: $(LISTMAKER2_SOURCES)
	$(CXX) $(LISTMAKER2_SOURCES) -o glistmaker2 $(LIBS) $(DEBUGFLAGS) -Wall

glistquery: $(LISTQUERY_SOURCES)
	$(CXX) $(LISTQUERY_SOURCES) -o glistquery $(LIBS) $(CXXFLAGS)

glistcompare: $(LISTCOMPARE_SOURCES)
	$(CXX) $(LISTCOMPARE_SOURCES) -o glistcompare $(LIBS) $(CXXFLAGS)

clean: clean-custom
	rm -f *.o $(BINS)

depend:
	$(CC) $(CFLAGS) -M *.c > .depend

include .depend
