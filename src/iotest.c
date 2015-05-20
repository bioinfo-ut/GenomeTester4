
#include <string.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>

static double
get_time (void)
{
	struct timeval tv;
	gettimeofday (&tv, NULL);
	double dtval = (tv.tv_sec + tv.tv_usec / 1000000.0);
	return dtval;
}

int
main (int argc, const char *argv[])
{
  const char *filenames[2];
  int nfiles, i;
  int streamin = 0;
  int streamout = 0;
  int blocksize = 8;
  unsigned long long filesize = 10000000000;
  for (i = 1; i < argc; i++) {
    if (!strcmp (argv[i], "-blocksize")) {
      blocksize = atoi (argv[++i]);
    } else if (!strcmp (argv[i], "-filesize")) {
      filesize = strtoll (argv[++i], NULL, 10);
    } else if (!strcmp (argv[i], "-streamin")) {
      streamin = 1;
    } else if (!strcmp (argv[i], "-streamout")) {
      streamout = 1;
    } else {
      if (nfiles < 2) filenames[nfiles++] = argv[i];
    }
  }

  fprintf (stderr, "File size is %lld block size is %d\n", filesize, blocksize);
  
  if (streamin) {
    FILE *ifs;
    double start, end;
    char *b;
    unsigned long long size;
    fprintf (stdout, "Stream reading (block = %d)\n", blocksize);
    start = get_time ();
    b = malloc (blocksize);
    ifs = fopen (filenames[0], "r");
    size = 0;
    while (!feof (ifs)) {
      fread (b, blocksize, 1, ifs);
      size += blocksize;
    }
    fclose (ifs);
    end = get_time ();
    fprintf (stdout, "Stream reading (size = %lld, block = %d) %.2f\n", size, blocksize, end - start);
  }

  if (streamout) {
    FILE *ofs;
    double start, end;
    char *b;
    unsigned long long i;
    unsigned int j;
    fprintf (stdout, "Stream writing (block = %d)\n", blocksize);
    start = get_time ();
    b = malloc (blocksize);
    for (i = 0; i < blocksize; i++) b[i] = i % 256;
    ofs = fopen (filenames[0], "w");
    for (i = 0; i < filesize; i += blocksize) {
      for (j = 0; j < blocksize; j++) b[j] = j % 256;
      fwrite (b, blocksize, 1, ofs);
    }
    fclose (ofs);
    end = get_time ();
    fprintf (stdout, "Stream writing (size = %lld, block = %d) %.2f\n", filesize, blocksize, end - start);
  }
}
