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

GDISTRIBUTION_SOURCES = \
	gdistribution.c \
	wordmap.c wordmap.h \
	wordtable.c wordtable.h \
	utils.c utils.h \
	queue.c queue.h \
	fasta.c fasta.h \
	sequence.c sequence.h

GMERCOUNTER_SOURCES = \
	gmercounter.c \
	trie.c trie.h \
	wordmap.c wordmap.h \
	sequence.c sequence.h \
	wordtable.c wordtable.h \
	fasta.c fasta.h \
	thread-pool.c thread-pool.h \
	queue.c queue.h \
	utils.c utils.h \
	database.c database.h

GMER_CALLER_SOURCES = \
	binomial.c binomial.h \
	genotypes.c genotypes.h \
	simplex.c simplex.h \
	thread-pool.c thread-pool.h \
	utils.c utils.h \
	gmer_caller.c

TEST_SOURCES = \
	binomial.c binomial.h \
	fasta.c fasta.h \
	queue.c queue.h \
	sequence.c sequence.h \
	simplex.c simplex.h \
	utils.c utils.h \
	wordmap.c wordmap.h \
	wordtable.c wordtable.h \
	distro.c

RELEASEFLAGS = -O3
DEBUGFLAGS = -O0 -g
LIBS = -lm -lpthread -lrt
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

gdistribution: $(GDISTRIBUTION_SOURCES)
	$(CXX) $(GDISTRIBUTION_SOURCES) -o gdistribution $(LIBS) $(DEBUGFLAGS) -Wall
	
gmer_counter: $(GMERCOUNTER_SOURCES)
	$(CXX) $(GMERCOUNTER_SOURCES) -o gmer_counter $(LIBS) $(CXXFLAGS) -Wall

gmer_caller: $(GMER_CALLER_SOURCES)
	$(CXX) $(GMER_CALLER_SOURCES) -o gmer_caller $(LIBS) $(CXXFLAGS) -Wall

distro: $(TEST_SOURCES)
	$(CXX) $(TEST_SOURCES) -o distro $(LIBS) $(DEBUGFLAGS) -Wall
	
clean: clean-custom
	rm -f *.o $(BINS)

depend:
	$(CC) $(CFLAGS) -M *.c > .depend

include .depend
