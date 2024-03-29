# Project: GenomeTester v4.0

VERSION = 4.0
CXX  = gcc

# C Files

LISTMAKER_SOURCES = \
	glistmaker.c \
	listmaker-queue.c listmaker-queue.h \
	queue.c queue.h \
	bloom.c bloom.h \
	fasta.c fasta.h \
	file-array.c file-array.h \
	buffer.c buffer.h \
	sequence.c sequence.h \
	sequence-block.c sequence-block.h \
	sequence-collection.c sequence-collection.h \
	sequence-file.c sequence-file.h \
	sequence-source.c sequence-source.h \
	sequence-stream.c sequence-stream.h \
	sequence-zstream.c sequence-zstream.h \
	set-operations.c set-operations.h \
	common.c common.h \
	utils.c utils.h \
	word-array-sorted.c word-array-sorted.h \
	word-dict.c word-dict.h \
	word-list.c word-list.h \
	word-list-sorted.c word-list-sorted.h \
	word-list-stream.c word-list-stream.h \
	word-map.c word-map.h \
	word-table.c word-table.h

LISTQUERY_SOURCES = \
	glistquery.c \
	listmaker-queue.c listmaker-queue.h \
	fasta.c fasta.h \
	bloom.c bloom.h \
	buffer.c buffer.h \
	file-array.c file-array.h \
	index-map.c index-map.h \
	sequence.c sequence.h \
	sequence-collection.c sequence-collection.h \
	sequence-block.c sequence-block.h \
	sequence-file.c sequence-file.h \
	sequence-source.c sequence-source.h \
	sequence-stream.c sequence-stream.h \
	sequence-zstream.c sequence-zstream.h \
	set-operations.c set-operations.h \
	common.c common.h \
	queue.c queue.h \
	utils.c utils.h \
	word-array-sorted.c word-array-sorted.h \
	word-dict.c word-dict.h \
	word-index.c word-index.h \
	word-list.c word-list.h \
	word-list-sorted.c word-list-sorted.h \
        word-list-stream.c word-list-stream.h \
        word-map.c word-map.h \
	word-table.c word-table.h

LISTCOMPARE_SOURCES = \
	glistcompare.c \
	listmaker-queue.c listmaker-queue.h \
	bloom.c bloom.h \
	fasta.c \
	buffer.c \
	file-array.c file-array.h \
	index-map.c index-map.h \
	sequence.c \
	sequence-collection.c sequence-collection.h \
	sequence-block.c sequence-block.h \
	sequence-file.c sequence-file.h \
	sequence-source.c sequence-source.h \
	sequence-stream.c sequence-stream.h \
	sequence-zstream.c sequence-zstream.h \
	common.c \
	queue.c queue.h \
	utils.c \
	word-array-sorted.c word-array-sorted.h \
	word-dict.c word-dict.h \
	word-index.c word-index.h \
	word-list.c word-list.h \
	word-list-sorted.c word-list-sorted.h \
	word-list-stream.c word-list-stream.h \
	word-map.c word-map.h \
	word-table.c word-table.h

GMER_COUNTER_SOURCES = \
	gmer_counter.c \
	database.c database.h \
	index.c index.h \
	listmaker-queue.c listmaker-queue.h \
	trie.c trie.h \
	sequence.c sequence.h \
	sequence-collection.c sequence-collection.h \
	sequence-block.c sequence-block.h \
	sequence-file.c sequence-file.h \
	sequence-source.c sequence-source.h \
	sequence-stream.c sequence-stream.h \
	sequence-zstream.c sequence-zstream.h \
	fasta.c fasta.h \
	thread-pool.c thread-pool.h \
	queue.c queue.h \
	utils.c utils.h \
	word-list.c word-list.h \
	word-table.c word-table.h

GMER_CALLER_SOURCES = \
	binomial.c binomial.h \
	genotypes.c genotypes.h \
	simplex.c simplex.h \
	thread-pool.c thread-pool.h \
	utils.c utils.h \
	gmer_caller.c

GINDEXER_SOURCES = \
	gindexer.c \
	fasta.c fasta.h \
	sequence.c sequence.h \
	sequence-source.c sequence-source.h \
	sequence-stream.c sequence-stream.h

GASSEMBLER_SOURCES = \
	binomial.c binomial.h \
	database.c database.h \
	fasta.c fasta.h \
	index.c index.h \
	matrix.c matrix.h \
	queue.c queue.h \
	sequence.c sequence.h \
	sequence-block.c sequence-block.h \
	sequence-collection.c sequence-collection.h \
	sequence-file.c sequence-file.h \
	sequence-source.c sequence-source.h \
	trie.c trie.h \
	utils.c utils.h \
	gassembler.c

ALEQ_SOURCES = \
	binomial.c binomial.h \
	sequence.c sequence.h \
	simplex.c simplex.h \
	utils.c utils.h \
	aleq.c

KMER_PREDICTOR_SOURCES = \
	utils.c utils.h \
	word-list-stream.c word-list-stream.h \
	word-array-sorted.c word-array-sorted.h \
	word-list.c word-list.h \
	kmer-predictor.c

ARIKKEI_SOURCES = \
	libarikkei/arikkei-strlib.c libarikkei/arikkei-strlib.h \
	libarikkei/arikkei-utils.c libarikkei/arikkei-utils.h

AZ_SOURCES = \
	az/base.h \
	az/class.c az/class.h \
	az/interface.c az/interface.h \
	az/reference.c az/reference.h \
	az/object.c az/object.h \
	az/primitives.c az/primitives.h \
	az/serialization.c az/serialization.h \
	az/types.c az/types.h

RELEASEFLAGS = -O3
DEBUGFLAGS = -O0 -g
AZ_FLAGS = -DAZ_NO_PROPERTIES -DAZ_NO_STRING -DAZ_NO_VALUE
LIBS = -lm -lpthread -lrt -lz
INCS = -I.
BINS  = glistmaker glistquery glistcompare

#CXXFLAGS = $(INCS) $(DEBUGFLAGS) -Wall 
CXXFLAGS = $(INCS) $(RELEASEFLAGS) -Wall -std=c99 -D_DEFAULT_SOURCE -D_GNU_SOURCE -D_POSIX_C_SOURCE=200809L -D_XOPEN_SOURCE

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BINS) all-after

glistmaker: $(LISTMAKER_SOURCES) $(AZ_SOURCES) $(ARIKKEI_SOURCES)
	$(CXX) $(LISTMAKER_SOURCES) $(AZ_SOURCES) $(ARIKKEI_SOURCES) -o glistmaker $(LIBS) $(CXXFLAGS) $(AZ_FLAGS)

glistquery: $(LISTQUERY_SOURCES) $(AZ_SOURCES) $(ARIKKEI_SOURCES)
	$(CXX) $(LISTQUERY_SOURCES) $(AZ_SOURCES) $(ARIKKEI_SOURCES) -o glistquery $(LIBS) $(CXXFLAGS) $(AZ_FLAGS)

glistcompare: $(LISTCOMPARE_SOURCES) $(AZ_SOURCES) $(ARIKKEI_SOURCES)
	$(CXX) $(LISTCOMPARE_SOURCES) $(AZ_SOURCES) $(ARIKKEI_SOURCES) -o glistcompare $(LIBS) $(CXXFLAGS) $(AZ_FLAGS)

gindexer: $(GINDEXER_SOURCES) $(AZ_SOURCES) $(ARIKKEI_SOURCES)
	$(CXX) $(GINDEXER_SOURCES) $(AZ_SOURCES) $(ARIKKEI_SOURCES) -o gindexer $(LIBS) $(CXXFLAGS) $(AZ_FLAGS)

gmer_counter: $(GMER_COUNTER_SOURCES) $(AZ_SOURCES) $(ARIKKEI_SOURCES)
	$(CXX) $(GMER_COUNTER_SOURCES) $(AZ_SOURCES) $(ARIKKEI_SOURCES) -o gmer_counter $(LIBS) $(CXXFLAGS) $(AZ_FLAGS)

gmer_caller: $(GMER_CALLER_SOURCES)
	$(CXX) $(GMER_CALLER_SOURCES) -o gmer_caller $(LIBS) $(CXXFLAGS) -Wall

gassembler: $(GASSEMBLER_SOURCES) $(AZ_SOURCES) $(ARIKKEI_SOURCES)
	$(CXX) $(GASSEMBLER_SOURCES) $(AZ_SOURCES) $(ARIKKEI_SOURCES) -o gassembler $(LIBS) $(CXXFLAGS) $(AZ_FLAGS)

aleq: $(ALEQ_SOURCES)
	$(CXX) $(ALEQ_SOURCES) -o aleq $(LIBS) $(CXXFLAGS) -Wall

kmer_predictor: $(KMER_PREDICTOR_SOURCES)
	$(CXX) $(KMER_PREDICTOR_SOURCES)  $(AZ_SOURCES) $(ARIKKEI_SOURCES) -o kmer_predictor $(LIBS) $(AZ_FLAGS) $(CXXFLAGS)
	
clean: clean-custom
	rm -f *.o $(BINS)

depend:
	$(CC) $(CFLAGS) -M *.c > .depend

include .depend
