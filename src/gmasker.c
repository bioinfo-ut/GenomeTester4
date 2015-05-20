/*
 * GenomeTester4
 *
 * A toolkit for creating and manipulating k-mer lists from biological sequences
 * 
 * Copyright (C) 2014 University of Tartu
 *
 * Authors: Maarja Lepamets and Lauris Kaplinski
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "wordmap.h"
#include "fasta.h"
#include "sequence.h"
#include "buffer.h"
#include "utils.h"
#include "common.h"

#define MAX_LISTS 100
#define DEFAULT_COEF 1.0

typedef struct _MaskerData {
	FILE *maskedfasta;
	wordmap *map;
	parameters *p;
	int *maskbeginning;
	unsigned int cutoff;
	unsigned int nmask5p;
	unsigned int nmask3p;
} MaskerData;

static void print_help (int exit_value);


int 
main (int argc, const char *argv[]) 
{
	int argidx, inputfirst, inputlast;
	int dosoftmasking = 0;
	unsigned int nlists = 0, listidx;
	unsigned int listpos[MAX_LISTS], listcomponents[MAX_LISTS];
	double pc = 0.95;
	char *end;
	size_t size;
	wordmap **maps;
	
	/* parsing and checking the commandline arguments */
	for (argidx = 1; argidx < argc; argidx++) {
  		if (!strcmp(argv[argidx], "-v") || !strcmp(argv[argidx], "--version")) {
			fprintf (stdout, "gmasker v%d.%d\n", VERSION_MAJOR, VERSION_MINOR);
			return 0;
		} else if (!strcmp (argv[argidx], "-h") || !strcmp (argv[argidx], "--help") || !strcmp (argv[argidx], "-?")) {
			print_help (0);

		} else if (argidx == 1) {
			/* get the locations of the input fasta files */
			if (argv[argidx][0] == '-') {
				fprintf (stderr, "Error: No FastA/FastQ file specified!\n");
				print_help (1);
			}
			inputfirst = argidx;
			while (argv[argidx + 1] && argv[argidx + 1][0] != '-') {
				argidx += 1;
			}
			inputlast = argidx;
		} else if (!strcmp (argv[argidx], "-l") || !strcmp (argv[argidx], "--list")) {
			if (nlists == MAX_LISTS) {
				fprintf (stderr, "Error: Maximum number of lists that can be specified is %u!\n", MAX_LISTS);
				return (1);
			}
			/* get the positions of list files */
			listpos[nlists] = argidx;
			while (argv[argidx + 1]) {
				if (argv[argidx + 1][0] == '-') {
					if (!argv[argidx + 1][1]) break;
					strtod (argv[argidx + 1] + 1, &end);
					if (*end != 0) break;
				}
				argidx += 1;
			}
			listcomponents[nlists] = argidx - listpos[nlists];
			nlists += 1;
		} else if (!strcmp (argv[argidx], "-p") || !strcmp (argv[argidx], "--prob_cutoff")) {
			if (!argv[argidx + 1]) {
				fprintf(stderr, "Warning: No probability cutoff specified! Using the default value: %.3f.\n", pc);
				argidx += 1;
				continue;
			}
			pc = strtod (argv[argidx + 1], &end);
			if (*end != 0 || pc < 0 || pc > 1) {
				fprintf(stderr, "Error: Invalid probability cutoff: %s! Must be a floating point number between 0 and 1.\n", argv[argidx + 1]);
				print_help (1);
			}
			argidx += 1;
		} else if (!strcmp (argv[argidx], "-s") || !strcmp (argv[argidx], "--soft_mask")) {
			dosoftmasking = 1;
		} else {
			fprintf (stderr, "Error: Unknown argument: %s!\n", argv[argidx]);
			print_help (1);
		}
	}
	
	/* parsing the list files */
	if (nlists == 0) {
		fprintf (stderr, "Error: No list files specified!\n");
		print_help (1);
	}
	
	maps = (wordmap **) malloc (nlists * sizeof (wordmap *));
	
	for (listidx = 0; listidx < nlists; listidx++) {
		unsigned int pos = listpos[listidx] + 1;
		parameters p = { 0 };
		
		maps[listidx] = wordmap_new (argv[pos]);
		if (!maps[listidx]) {
			fprintf(stderr, "Error: List file not found at position %u!\n", pos);
			print_help (1);	
		}
		p.wordlength = (maps[listidx])->header->wordlength;
		p.coef = DEFAULT_COEF;
		
		if (listcomponents[listidx] > 1) {
			p.coef = (argv[pos + 1][0] == '-') ? strtod (argv[pos + 1] + 1, &end) * (-1.0) : strtod (argv[pos + 1], &end);
			if (*end != 0) {
				fprintf (stderr, "Error: Invalid coefficient: %s! Must be a floating point number.\n", argv[pos + 1]);
				print_help (1);
			}
		}
		
		if (listcomponents[listidx] > 2) {
			p.nmm = strtol (argv[pos + 2], &end, 10);
			if (*end != 0 || p.nmm > p.wordlength) {
				fprintf (stderr, "Error: Invalid number of mismatches: %s! Must be a positive integer smaller than word length %u.\n", argv[pos + 2], p.wordlength);
				print_help (1);
			}
		}
		
		if (listcomponents[listidx] > 3) {
			p.pm3 = strtod (argv[pos + 3], &end);
			if (*end != 0 || p.pm3 > p.wordlength) {
				fprintf (stderr, "Error: Invalid number of 3' perfect matches: %s! Must be a positive integer smaller than word length %u.\n", argv[pos + 3], p.wordlength);
				print_help (1);
			}
			if (p.nmm + p.pm3 > p.wordlength) {
				fprintf (stderr, "Error: Too many mismatches or 3' perfect matches for current word length %u.\n", p.wordlength);
				return 1;
			}
		}

		
	}
	  
	/*map = wordmap_new (list);*/
	
	/* run masker with user-specified parameters */
	/*for (argidx = inputlast; argidx >= inputfirst; argidx--) {
		const char *fasta = mmap_by_filename (argv[argidx], &size);
		run_masker (fasta, size, map, 2, 1);  
	}*/
	
	
	
	
	
	return 0;
	
}

void 
flush_sequence_and_mask_n_from_end (FILE *out, SequenceBuffer *sbuffer, unsigned int n, int dosoftmasking, char cmask)
{
	char c;
	unsigned int i = sequence_buffer_get_index_gap (sbuffer);
	while (i > 0) {
		c = sequence_buffer_get_next_char (sbuffer);
		if (i > n + sbuffer->nwhitespace) fprintf (out, "%c", c);
		else if (c <= ' ') fprintf (out, "%c", c);
		else fprintf (out, "%c", cmask);
		i -= 1;
	}	
	return;
}

void
flush_sequence_and_mask_n_from_beginning (FILE *out, SequenceBuffer *sbuffer, unsigned int n, unsigned int unprinted, int dosoftmasking, char cmask)
{
	char c;
	unsigned int i = sequence_buffer_get_index_gap (sbuffer);
	while (i > unprinted) {
		c = sequence_buffer_get_next_char (sbuffer);
		if (n > 0 && c > ' ') {
			fprintf (out, "%c", cmask);
			n -= 1;
		}
		else fprintf (out, "%c", c);
		i -= 1;
	}	
	return;	
}

void 
print_masked_sequence (FILE *out, SequenceBuffer *sbuffer, unsigned int n, unsigned int unprinted, int dosoftmasking, char cmask)
{
	char c;
	unsigned int i = sequence_buffer_get_index_gap (sbuffer);
	while (i > unprinted) {
		c = sequence_buffer_get_next_char (sbuffer);
		if (n > 0 && c > ' ') {
			fprintf (out, "%c", cmask);
			n -= 1;
		}
		else fprintf (out, "%c", c);
		i -= 1;
	}	
	return;	
}

int
mask_sequence (FastaReader *reader, void *data)
{
	MaskerData *mdata = (MaskerData *) data;
	FILE *maskedfasta = mdata->maskedfasta;
	wordmap *map = mdata->map;
	parameters *p = mdata->p;
	int *maskbeginning = mdata->maskbeginning;
	unsigned int nmask5p = mdata->nmask5p;
	unsigned int nmask3p = mdata->nmask3p;
	unsigned int cutoff = mdata->cutoff;
	
	SequenceBuffer *sbuffer = reader->sbuffer;
	unsigned long long word;
	unsigned int freq;
	int mask = 0;
	int v = 0;
	
	if (reader->noword) {
		flush_sequence_and_mask_n_from_beginning (maskedfasta, reader->sbuffer, *maskbeginning ? nmask3p : 0, nmask5p + reader->sbuffer->nwhitespace, 0, 'N');
		return v;
	}
	
	word = reader->currentword;
	//printf("\notsin sõna: %s\n", word_to_string (word, reader->wordlength));
	
	freq = wordmap_search_query (map, word, p, 0, 0, 0, NULL);
	//printf("leitud sagedus: %u\n", freq);
	
	// FIXME: lineaarkombinatsioon
	if (freq >= cutoff) mask = 1;
	
	
	if (mask) {
		flush_sequence_and_mask_n_from_beginning (maskedfasta, reader->sbuffer, *maskbeginning ? nmask3p : 0, nmask5p + reader->sbuffer->nwhitespace, 0, 'N');
		flush_sequence_and_mask_n_from_end (maskedfasta, reader->sbuffer, nmask5p, 0, 'N');
		*maskbeginning = 1;
	} else if (sequence_buffer_is_full (reader->sbuffer)) {
		flush_sequence_and_mask_n_from_beginning (maskedfasta, reader->sbuffer, *maskbeginning ? nmask3p : 0, nmask5p + reader->sbuffer->nwhitespace, 0, 'N');
	}
	
	
	return v;
}


int 
print_header (FastaReader *reader, void *data)
{
	MaskerData *mdata = (MaskerData *) data;
	fprintf (mdata->maskedfasta, ">%s\n", reader->seqname);
	return 0;
}

int 
flush_sequence_end (FastaReader *reader, void *data)
{
	MaskerData *mdata = (MaskerData *) data;
	FILE *maskedfasta = mdata->maskedfasta;
	int *maskbeginning = mdata->maskbeginning;
	unsigned int nmask3p = mdata->nmask3p;
	flush_sequence_and_mask_n_from_beginning (maskedfasta, reader->sbuffer, *maskbeginning ? nmask3p : 0, 0, 0, 'N');
	*maskbeginning = 0;
	mdata->maskbeginning = maskbeginning;
	return 0;
}


int
run_masker (const char *fasta, size_t fsize, wordmap *map, unsigned int nmask5p, unsigned int nmask3p)
{
	MaskerData mdata = { 0 };
	parameters p = { 0 };
	int maskbeginning = 0;
	
	int v;
	//FILE *maskedfasta;
	unsigned int cutoff = 1;
	
	/* FIXME: hardcoded */
	//maskedfasta = fopen ("masked_seq.fasta", "w");
	mdata.maskedfasta = stdout;
	mdata.map = map;
	mdata.cutoff = cutoff;
	mdata.nmask5p = nmask5p;
	mdata.nmask3p = nmask3p;
	mdata.maskbeginning = &maskbeginning;
	
	p.wordlength = map->header->wordlength;
	mdata.p = &p;
	
	
	v = fasta_reader_read (fasta, fsize, map->header->wordlength, (void *) &mdata, 0, 1, print_header, flush_sequence_end, mask_sequence);
	
	
	//fclose (maskedfasta);
	
	if (v) return v;
	return 0;	
	
}


static void
print_help (int exit_value)
{
	fprintf (stdout, "Usage: gmasker <INPUTFILES> [OPTIONS]\n");
	fprintf (stdout, "Options:\n");
	fprintf (stdout, "    -v, --version            - print version information and exit\n");
	fprintf (stdout, "    -h, --help               - print this usage screen and exit\n");

	
	fprintf (stdout, "    -D                       - increase debug level\n");
	exit (exit_value);
}
