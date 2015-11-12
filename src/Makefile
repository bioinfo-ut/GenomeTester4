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
	trie.c trie.h \
	utils.c

LISTQUERY_SOURCES = \
	glistquery.c \
	wordtable.c \
	wordmap.c \
	fasta.c \
	buffer.c \
	sequence.c \
	common.c \
	queue.c queue.h \
	utils.c

LISTCOMPARE_SOURCES = \
	glistcompare.c \
	wordtable.c \
	wordmap.c \
	fasta.c \
	buffer.c \
	sequence.c \
	common.c \
	queue.c queue.h \
	utils.c

SNPCALLER_SOURCES = \
	snpcaller.c \
	trie.c trie.h \
	wordmap.c wordmap.h \
	sequence.c sequence.h \
	wordtable.c wordtable.h \
	fasta.c fasta.h \
	thread-pool.c thread-pool.h \
	queue.c queue.h \
	utils.c utils.h

NIPT_SOURCES = \
	binomial.c binomial.h \
	simplex.c simplex.h \
	utils.c utils.h \
	nipt.c

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

snpcaller: $(SNPCALLER_SOURCES)
	$(CXX) $(SNPCALLER_SOURCES) -o snpcaller $(LIBS) $(DEBUGFLAGS) -Wall

nipt: $(NIPT_SOURCES)
	$(CXX) $(NIPT_SOURCES) -o nipt $(LIBS) $(DEBUGFLAGS) -Wall
	
clean: clean-custom
	rm -f *.o $(BINS)

depend:
	$(CC) $(CFLAGS) -M *.c > .depend

include .depend
