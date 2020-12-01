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
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include "Common.h"
#include "Reads.h"
#include "HashTable.h"

#define CHARCODE(a) (a=='A' ? 0 : (a=='C' ? 1 : (a=='G' ? 2 : (a=='T' ? 3 : 4))))

FILE *_r_fp1;
FILE *_r_fp2;
gzFile _r_gzfp1;
gzFile _r_gzfp2;
Read *_r_seq;
int _r_seqCnt;
int *_r_samplingLocs;

/**********************************************/
char *(*readFirstSeq)(char *);
char *(*readSecondSeq)(char *);
/**********************************************/
char *readFirstSeqTXT( char *seq )
{
  return fgets(seq, SEQ_MAX_LENGTH, _r_fp1);
}

/**********************************************/
char *readSecondSeqTXT( char *seq )
{
  return fgets(seq, SEQ_MAX_LENGTH, _r_fp2);
}
/**********************************************/
char *readFirstSeqGZ( char *seq )
{
  return gzgets(_r_gzfp1, seq, SEQ_MAX_LENGTH);
}

/**********************************************/
char *readSecondSeqGZ( char *seq )
{
  return gzgets(_r_gzfp2, seq, SEQ_MAX_LENGTH);
}
/**********************************************/
int toCompareRead(const void * elem1, const void * elem2)
{
  return strcmp(((Read *)elem1)->seq, ((Read *)elem2)->seq);	
}
/**********************************************/
int readAllReads(char *fileName1,
		 char *fileName2,
		 int compressed,
		 unsigned char *fastq,
		 unsigned char pairedEnd,
		 Read **seqList,
		 unsigned int *seqListSize)
{
  double startTime=getTime();
  
  char seq1[SEQ_MAX_LENGTH];
  char rseq1[SEQ_MAX_LENGTH];
  char name1[SEQ_MAX_LENGTH];
  char qual1[SEQ_MAX_LENGTH];
  char seq2[SEQ_MAX_LENGTH];
  char rseq2[SEQ_MAX_LENGTH];
  char name2[SEQ_MAX_LENGTH];
  char qual2[SEQ_MAX_LENGTH];

  char dummy[SEQ_MAX_LENGTH];
  char ch;
  int err1, err2;
  int nCnt;
  int discarded = 0;
  int seqCnt = 0;
  int maxCnt = 0;
  int i;
  Read *list = NULL;

  int clipped = 0;

  if (!compressed)
    {
      _r_fp1 = fileOpen( fileName1, "r");

      if (_r_fp1 == NULL)
	{
	  return 0;
	}

      ch = fgetc(_r_fp1);

      if ( pairedEnd && fileName2 != NULL )
	{
	  _r_fp2 = fileOpen ( fileName2, "r" );
	  if (_r_fp2 == NULL)
	    {
	      return 0;
	    }
	}
      else
	{
	  _r_fp2 = _r_fp1;
	}

      readFirstSeq = &readFirstSeqTXT;
      readSecondSeq = &readSecondSeqTXT;
    }
  else
    {

      _r_gzfp1 = fileOpenGZ (fileName1, "r");

      if (_r_gzfp1 == NULL)
	{
	  return 0;
	}

      ch = gzgetc(_r_gzfp1);

      if ( pairedEnd && fileName2 != NULL )
	{
	  _r_gzfp2 = fileOpenGZ ( fileName2, "r" );
	  if (_r_gzfp2 == NULL)
	    {
	      return 0;
	    }
	}
      else
	{
	  _r_gzfp2 = _r_gzfp1;
	}

      readFirstSeq = &readFirstSeqGZ;
      readSecondSeq = &readSecondSeqGZ;
    }

  if (ch == '>')
    *fastq = 0;
  else
    *fastq = 1;

  // Counting the number of lines in the file
  while (readFirstSeq(dummy)) maxCnt++;

  if (!compressed)
    {
      rewind(_r_fp1);
    }
  else
    {
      gzrewind(_r_gzfp1);
    }

  // Calculating the Maximum # of sequences
  if (*fastq)
    {
      if (maxCnt % 4 != 0){
        fprintf(stderr, "Input FASTQ file seems to be truncated. Number of lines: %d. Exiting.\n", maxCnt);
	exit (1);
      }
      maxCnt /= 4;
    }
  else
    {
      if (maxCnt % 2 != 0){
        fprintf(stderr, "Input FASTQ file seems to be truncated. Number of lines: %d. Exiting.\n", maxCnt);
	exit (1);
      }
      maxCnt /= 2;
    }

  if (pairedEnd && fileName2 != NULL )
    maxCnt *= 2;

  list = getMem(sizeof(Read)*maxCnt);

  while( readFirstSeq(name1) )
    {
      err1 = 0;
      err2 = 0;
      readFirstSeq(seq1);
      
      
      name1[strlen(name1)-1] = '\0';
      for (i=0; i<strlen(name1);i++)
	{
	  if (name1[i] == ' ')
	    {
	      name1[i] = '\0';
	      break;
	    }
	  
	}
      
      if (errThreshold == 255){
	if (cropSize > 0) 
	  {
	    errThreshold = (int) ceil(cropSize * 0.04);
	    fprintf(stderr, "Sequence length: %d bp. Error threshold is set to %d bp.\n", cropSize, errThreshold);
	  }
	else
	  {
	    errThreshold = (int) ceil((strlen(seq1)-1) * 0.04);
	    fprintf(stderr, "Sequence length: %d bp. Error threshold is set to %d bp.\n", ((int)strlen(seq1)-1), errThreshold);
	  }
	fprintf(stderr, "You can override this value using the -e parameter.\n");
      }
	
      
      if ( *fastq )
	{
	  readFirstSeq(dummy);
	  readFirstSeq(qual1);
	  qual1[strlen(qual1)-1] = '\0';
	}
      else
	{
	  sprintf(qual1, "*");
	}
	
      // Cropping
      if (cropSize > 0)
	{
	  seq1[cropSize] = '\0';
	  if ( *fastq )
	    qual1[cropSize] = '\0';
	}


      nCnt = 0;
      for (i=0; i<strlen(seq1); i++)
	{
	  seq1[i] = toupper (seq1[i]);
	  if (seq1[i] == 'N')
	    {
	      nCnt++;
	    }
	  else if (isspace(seq1[i]))
	    {

	      seq1[i] = '\0';
	      break;
	    }
	}

      if (nCnt > errThreshold)
	{
	  err1 = 1;
	}

      // Reading the second seq of pair-ends
      if (pairedEnd)
	{
	  readSecondSeq(name2);
	  readSecondSeq(seq2);
	  name2[strlen(name2)-1] = '\0';
	  for (i=0; i<strlen(name2);i++)
	    {
	      if (name2[i] == ' ')
		{
		  name2[i] = '\0';
		  break;
		}
	      
	    }
	  
	  if ( *fastq )
	    {
	      readSecondSeq(dummy);
	      readSecondSeq(qual2);

	      qual2[strlen(qual2)-1] = '\0';
	    }
	  else
	    {
	      sprintf(qual2, "*");
	    }
			
			
	  // Cropping
	  if (cropSize > 0)
	    {
	      seq2[cropSize] = '\0';
	      if ( *fastq )
		qual2[cropSize] = '\0';
	    }


	  nCnt = 0;
	  for (i=0; i<strlen(seq2); i++)
	    {
	      seq2[i] = toupper (seq2[i]);
	      if (seq2[i] == 'N')
		{
		  nCnt++;

		}
	      else if (isspace(seq2[i]))
		{
		  seq2[i] = '\0';
		}
	    }
	  if (nCnt > errThreshold)
	    {
	      err2 = 1;
	    }


	  if (strlen(seq1) < strlen(seq2)) { 
	    seq2[strlen(seq1)] = '\0'; 
	    if ( *fastq )
	      qual2[strlen(seq1)] = '\0'; 
	    if (!clipped) clipped = 2; 
	  }
	  else if (strlen(seq1) > strlen(seq2)){
	    seq1[strlen(seq2)] = '\0';
	    if ( *fastq )
	      qual1[strlen(seq2)] = '\0';
	    if (!clipped) clipped = 1;
	  }

	  if (clipped == 1 || clipped == 2){
	    fprintf(stderr, "[PE mode Warning] Sequence lengths are different,  read #%d is clipped to match.\n", clipped);
	    clipped = 3;
	  }
			


	}

      if (!pairedEnd && !err1)
	{
	  int _mtmp = strlen(seq1);
	  list[seqCnt].hits = getMem (1+3*_mtmp+3+strlen(name1)+1);
	  list[seqCnt].seq  = list[seqCnt].hits + 1;
	  list[seqCnt].rseq = list[seqCnt].seq + _mtmp+1;
	  list[seqCnt].qual = list[seqCnt].rseq + _mtmp+1;
	  list[seqCnt].name = list[seqCnt].qual + _mtmp+1;
                        

	  list[seqCnt].readNumber = seqCnt;			

	  reverseComplement(seq1, rseq1, _mtmp);

	  int i;

	  list[seqCnt].hits[0] = 0;

	  for (i=0; i<=_mtmp; i++)
	    {
	      list[seqCnt].seq[i] = seq1[i];
	      list[seqCnt].rseq[i] = rseq1[i] ;
	      list[seqCnt].qual[i] = qual1[i];
	    }
	  
	  list[seqCnt].rseq[_mtmp]=list[seqCnt].qual[_mtmp]='\0';	
	
	  sprintf(list[seqCnt].name,"%s%c", ((char*)name1)+1,'\0');

	  seqCnt++;

	}
      else if (pairedEnd && !err1 && !err2)
	{

	  // Naming Conventions X/1, X/2 OR X
	  int tmplen = strlen(name1);
	  if (strcmp(name1, name2) != 0)
	    {
	      tmplen = strlen(name1)-2;
	    }
		
	  //first seq
	  int _mtmp = strlen(seq1);
	  list[seqCnt].hits = getMem (1+3*_mtmp+3+tmplen+1);
	  list[seqCnt].seq = list[seqCnt].hits + 1;
	  list[seqCnt].rseq = list[seqCnt].seq + _mtmp+1;
	  list[seqCnt].qual = list[seqCnt].rseq + _mtmp+1;
	  list[seqCnt].name = list[seqCnt].qual + _mtmp+1;

	  list[seqCnt].readNumber = seqCnt;

	  reverseComplement(seq1, rseq1, _mtmp);

	  int i;

	  list[seqCnt].hits[0] = 0;

	  for (i=0; i<=_mtmp; i++)
	    {
	      list[seqCnt].seq[i] = seq1[i];
	      list[seqCnt].rseq[i] = rseq1[i] ;
	      list[seqCnt].qual[i] = qual1[i];
	    }


	  name1[tmplen]='\0';
	  list[seqCnt].rseq[_mtmp]=list[seqCnt].qual[_mtmp]='\0';

	  sprintf(list[seqCnt].name,"%s%c", ((char*)name1)+1,'\0');

	  seqCnt++;

	  //second seq
	  list[seqCnt].hits = getMem (1+3*_mtmp+3+tmplen+1);
	  list[seqCnt].seq = list[seqCnt].hits + 1;
	  list[seqCnt].rseq = list[seqCnt].seq + _mtmp+1;
	  list[seqCnt].qual = list[seqCnt].rseq + _mtmp+1;
	  list[seqCnt].name = list[seqCnt].qual + _mtmp+1;

	  list[seqCnt].readNumber = seqCnt;

	  reverseComplement(seq2, rseq2, _mtmp);
	  //rseq2[_mtmp] =  '\0';

	  list[seqCnt].hits[0] = 0;

	  for (i=0; i<=_mtmp; i++)
	    {
	      list[seqCnt].seq[i] = seq2[i];
	      list[seqCnt].rseq[i] = rseq2[i] ;
	      list[seqCnt].qual[i] = qual2[i];
	    }


	  name2[tmplen]='\0';
	  list[seqCnt].rseq[_mtmp]=list[seqCnt].qual[_mtmp]='\0';

	  sprintf(list[seqCnt].name,"%s%c", ((char*)name2)+1,'\0');

	  seqCnt++;

	}
      else
	{
	  discarded++;
	}
    }

  if (seqCnt > 0)
    {
      SEQ_LENGTH = strlen(list[0].seq);
    }
  else
    {
      fprintf(stderr, "ERROR: No reads can be found for mapping\n");
      return 0;
    }


  // Closing Files
  if (!compressed)
    {
      fclose(_r_fp1);
      if ( pairedEnd && fileName2 != NULL )
	{
	  fclose(_r_fp2);
	}
    }
  else
    {
      gzclose(_r_gzfp1);
      if ( pairedEnd && fileName2 != NULL)
	{
	  gzclose(_r_gzfp2);
	}
    }

  adjustQual(list, seqCnt);

  *seqList = list;
  *seqListSize = seqCnt;
	

  _r_seq = list;
  _r_seqCnt = seqCnt;

  if ( pairedEnd ) discarded *= 2;

  if (seqCnt>1)
    fprintf(stderr, "%d sequences are read in %0.2f. (%d discarded) [Mem:%0.2f M]\n", seqCnt, (getTime()-startTime), discarded, getMemUsage());
  else
    fprintf(stderr, "%d sequence is read in %0.2f. (%d discarded) [Mem:%0.2f M]\n", seqCnt, (getTime()-startTime), discarded, getMemUsage());

  return 1;
}


/**********************************************/
void loadSamplingLocations(int **samplingLocs, int * samplingLocsSize)
{
  int i;
  int samLocsSize = errThreshold + 1;
  int *samLocs = getMem(sizeof(int)*samLocsSize);

  for (i=0; i<samLocsSize; i++)
    {
      samLocs[i] = (SEQ_LENGTH / samLocsSize) *i;
      if ( samLocs[i] + WINDOW_SIZE > SEQ_LENGTH)
	samLocs[i] = SEQ_LENGTH - WINDOW_SIZE;
    }
  *samplingLocs = samLocs;
  *samplingLocsSize = samLocsSize;
  _r_samplingLocs = samLocs;
}

void finalizeReads(char *fileName)
{
  FILE *fp1=NULL;

  if (fileName != NULL)
    {
      fp1 = fileOpen(fileName, "w");
    }
  if (pairedEndMode)
    _r_seqCnt /=2;

  int i=0;
  for (i = 0; i < _r_seqCnt; i++)
    {
      if (pairedEndMode && _r_seq[2*i].hits[0] == 0 && _r_seq[2*i+1].hits[0] == 0  &&  strcmp(_r_seq[2*i].qual,"*")!=0)
	{
	  fprintf(fp1,"@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].qual, _r_seq[i*2].name, _r_seq[i*2+1].seq, _r_seq[i*2+1].qual);
	}
      else if (pairedEndMode && _r_seq[2*i].hits[0] == 0 && _r_seq[2*i+1].hits[0] == 0)
	{
	  fprintf(fp1, ">%s/1\n%s\n>%s/2\n%shits=%d\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].name, _r_seq[i*2+1].seq, _r_seq[2*i+1].hits[0]);
	}
      else if (!pairedEndMode && _r_seq[i].hits[0] == 0 && strcmp(_r_seq[i].qual, "*")!=0)
	{
	  fprintf(fp1,"@%s\n%s\n+\n%s\n", _r_seq[i].name, _r_seq[i].seq, _r_seq[i].qual);
	}
      else if (!pairedEndMode && _r_seq[i].hits[0] == 0)
	{
	  fprintf(fp1,">%s\n%s\n", _r_seq[i].name, _r_seq[i].seq);
	}
    }

  fclose(fp1);
  if (pairedEndMode)
    _r_seqCnt *= 2;

  for (i = 0; i < _r_seqCnt; i++)
    {
      freeMem(_r_seq[i].hits,0);
    }


  freeMem(_r_seq,0);
  freeMem(_r_samplingLocs,0);
}

void adjustQual(Read *list, int seqCnt){
  /* This function will automatically determine the phred_offset and readjust quality values if needed */
  int i,j,q, offset=64;
  int len = strlen(list[0].qual);
  
  for (i=0; i<10000 && i<seqCnt; i++){
    for (j=0;j<len;j++){
      q = (int) list[i].qual[j] - offset;
      if (q < 0){
	offset = 33;
	break;
      }
    }
    if (offset == 33)
      break;
  }
  
  if (offset == 64){
    fprintf(stderr, "[Quality Warning] Phred offset is 64. Readjusting to 33.\n");
    fflush(stderr);
    for (i=0;i<seqCnt;i++){
      for (j=0;j<len;j++){
	list[i].qual[j] -= 31;
      }
    }
  }
}


