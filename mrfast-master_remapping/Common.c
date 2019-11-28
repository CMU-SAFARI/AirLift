/*
 * Copyright (c) <2008 - 2020>, University of Washington, Simon Fraser University, 
 * Bilkent University and Carnegie Mellon University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the names of the University of Washington, Simon Fraser University, 
 *   Bilkent University, Carnegie Mellon University,
 *   nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  Authors: 
  Farhad Hormozdiari
	  farhadh AT uw DOT edu
  Faraz Hach
	  fhach AT cs DOT sfu DOT ca
  Can Alkan
	  calkan AT gmail DOT com
  Hongyi Xin
	  gohongyi AT gmail DOT com
  Donghyuk Lee
	  bleups AT gmail DOT com
*/


#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <zlib.h>
#include <string.h>
#include "Common.h"


unsigned short 			SEQ_LENGTH = 0;
long long		       	memUsage = 0;
/**********************************************/
FILE *fileOpen(char *fileName, char *mode)
{
  FILE *fp;
  fp = fopen (fileName, mode);
  if (fp == NULL)
    {
      fprintf(stderr, "Error: Cannot Open the file %s\n", fileName);
      fflush(stderr);
      exit(0);
    }
  return fp;
}
/**********************************************/
gzFile fileOpenGZ(char *fileName, char *mode)
{
  gzFile gzfp;
  gzfp = gzopen (fileName, mode);
  if (gzfp == Z_NULL)
    {
      fprintf(stderr, "Error: Cannot Open the file %s\n", fileName);
      fflush(stderr);
      exit(0);
    }
  return gzfp;
}
/**********************************************/
double getTime(void)
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec+t.tv_usec/1000000.0;
}

/**********************************************/
char reverseComplementChar(char c)
{
  char ret;
  switch (c)
    {
    case 'A': 
      ret = 'T';
      break;
    case 'T':
      ret = 'A';
      break;
    case 'C':	
      ret = 'G';
      break;
    case 'G':
      ret = 'C';
      break;
    default:
      ret = 'N';
      break;
    }
  return ret;
}
/**********************************************/
void reverseComplement (char *seq, char *rcSeq , int length)
{
  int i;
  for (i=0; i<length; i++)
    {
      rcSeq[i]=reverseComplementChar (seq[length-1-i]) ;
    }
}
/**********************************************/
void * getMem(size_t size)
{
  void *ret;
  ret = malloc(size);
  if (ret == NULL){
    fprintf(stderr, "Cannot allocate memory. Currently addressed memory = %0.2f MB, requested memory = %0.2f MB.\nCheck the available main memory, and if you have user limits (ulimit -v).\n", getMemUsage(), (float)(size/1048576.0));
    exit(0);
  }
  memUsage+=size;
  return ret;
}
/**********************************************/
void reMem(void *ptr, size_t oldsize, size_t newsize)
{
  ptr = realloc(ptr, newsize);
  if (ptr == NULL){
    fprintf(stderr, "Cannot reallocate memory. Currently addressed memory = %0.2f MB, requested memory = %0.2f MB.\nCheck the available main memory, and if you have user limits (ulimit -v).\n", getMemUsage(), (double)newsize);
    exit(0);
  }
  memUsage+=newsize-oldsize;
}
/**********************************************/
void freeMem(void *ptr, size_t size)
{
  memUsage-=size;
  free(ptr);
}
/**********************************************/
double getMemUsage()
{
  return memUsage/1048576.0;
}
/**********************************************/
void reverse (char *seq, char *rcSeq , int length)
{
  int i;
  int l = length;
  if(l != strlen(seq))
    l = strlen(seq);
  for (i=0; i<l; i++)
    {
      rcSeq[i]=seq[l-1-i] ;
    }
  rcSeq[l] = '\0';

}
/**********************************************/
void stripPath(char *full, char **path, char **fileName)
{
  int i;
  int pos = -1;

  for (i=strlen(full)-1; i>=0; i--)
    {
      if (full[i]=='/')
	{
	  pos = i;
	  break;
	}

    }

  if (pos != -1)
    {
      sprintf(*fileName, "%s%c", (full+pos+1), '\0');
      full[pos+1]='\0';
      sprintf(*path,"%s%c", full, '\0');
    }
  else
    {
      sprintf(*fileName, "%s%c", full, '\0');
      sprintf(*path,"%c", '\0');
    }
}
