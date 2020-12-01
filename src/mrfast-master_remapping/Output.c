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
#include <zlib.h>
#include <string.h>
#include "Common.h"
#include "Output.h"

FILE			*_out_fp;
gzFile			_out_gzfp;

char buffer[300000];
int bufferSize = 0;


void finalizeGZOutput()
{
  gzclose(_out_gzfp);
}

void finalizeTXOutput()
{
  fclose(_out_fp);
}


void gzOutputQ(SAM map)
{
  gzprintf(_out_gzfp, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s", 
	   map.QNAME, 
	   map.FLAG,
	   map.RNAME, 
	   map.POS,
	   map.MAPQ,
	   map.CIGAR,
	   map.MRNAME,
	   map.MPOS,
	   map.ISIZE,
	   map.SEQ,
	   map.QUAL);
	
  int i;

  for ( i = 0; i < map.optSize; i++)
    {
      switch (map.optFields[i].type)
	{
	case 'A':
	  gzprintf(_out_gzfp, "\t%s:%c:%c", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].cVal);
	  break;
	case 'i':
	  gzprintf(_out_gzfp, "\t%s:%c:%d", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].iVal);
	  break;
	case 'f':
	  gzprintf(_out_gzfp, "\t%s:%c:%f", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].fVal);
	  break;
	case 'Z':
	case 'H':
	  gzprintf(_out_gzfp, "\t%s:%c:%s", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].sVal);
	  break;
	}
    }

  if (readGroup[0] != 0)
    gzprintf(_out_gzfp, "\t%s:%c:%s", "RG", 'Z', readGroup);

  gzprintf(_out_gzfp, "\n");
}

void outputSAM(FILE *fp, SAM map)
{
  fprintf(fp, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
	  map.QNAME,
	  map.FLAG,
	  map.RNAME,
	  map.POS,
	  map.MAPQ,
	  map.CIGAR,
	  map.MRNAME,
	  map.MPOS,
	  map.ISIZE,
	  map.SEQ,
	  map.QUAL);


  int i;

  for ( i = 0; i < map.optSize; i++)
    {
      switch (map.optFields[i].type)
	{
	case 'A':
	  fprintf(fp, "\t%s:%c:%c", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].cVal);
	  break;
	case 'i':
	  fprintf(fp, "\t%s:%c:%d", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].iVal);
	  break;
	case 'f':
	  fprintf(fp, "\t%s:%c:%f", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].fVal);
	  break;
	case 'Z':
	case 'H':
	  fprintf(fp, "\t%s:%c:%s", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].sVal);
	  break;
	}
    }

  if (readGroup[0] != 0)
    fprintf(fp, "\t%s:%c:%s", "RG", 'Z', readGroup);

  fprintf(fp, "\n");

}

void outputQ(SAM map)
{

  fprintf(_out_fp, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s", 
	  map.QNAME, 
	  map.FLAG,
	  map.RNAME, 
	  map.POS,
	  map.MAPQ,
	  map.CIGAR,
	  map.MRNAME,
	  map.MPOS,
	  map.ISIZE,
	  map.SEQ,
	  map.QUAL);

	
  int i;

  for ( i = 0; i < map.optSize; i++)
    {
      switch (map.optFields[i].type)
	{
	case 'A':
	  fprintf(_out_fp, "\t%s:%c:%c", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].cVal);
	  break;
	case 'i':
	  fprintf(_out_fp, "\t%s:%c:%d", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].iVal);
	  break;
	case 'f':
	  fprintf(_out_fp, "\t%s:%c:%f", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].fVal);
	  break;
	case 'Z':
	case 'H':
	  fprintf(_out_fp, "\t%s:%c:%s", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].sVal);
	  break;
	}
    }

  if (readGroup[0] != 0)
    fprintf(_out_fp, "\t%s:%c:%s", "RG", 'Z', readGroup);
	
  fprintf(_out_fp, "\n");
}

int initOutput ( char *fileName, int compressed)
{
  if (compressed)
    {
      char newFileName[strlen(fileName)+4];
      sprintf(newFileName, "%s.gz", fileName);

      _out_gzfp = fileOpenGZ(newFileName, "w1f");
      if (_out_gzfp == Z_NULL)
	{
	  return 0;
	}
      finalizeOutput = &finalizeGZOutput;

      output = &gzOutputQ;
      SAMheaderGZ(_out_gzfp);
    }
  else
    {
      _out_fp = fileOpen(fileName, "w");
      if (_out_fp == NULL)
	{
	  return 0;
	}
      finalizeOutput = &finalizeTXOutput;
      output = &outputQ;
      SAMheaderTX(_out_fp, 1);
    }
  buffer[0] = '\0';
  return 1;
}

FILE* getOutputFILE() 
{
  if(_out_fp != NULL)
    return _out_fp;
  else if(_out_gzfp != NULL)
    return (FILE *)_out_gzfp;
  else
    return NULL;
}
void SAMheaderTX(FILE *outfp, int check)
{
  FILE *fp;
  char fainame[FILE_NAME_LENGTH];
  char chrom[FILE_NAME_LENGTH];
  int chromlen;
  char rest[FILE_NAME_LENGTH];
  char *ret = NULL;

  sprintf(fainame, "%s.fai",fileName[0]);
  fp = fopen(fainame, "r");

  if (fp != NULL){
    fprintf(outfp, "@HD\tVN:1.4\tSO:unknown\n");
    
    while (fscanf(fp, "%s\t%d\t", chrom, &chromlen) > 0){
      ret = fgets(rest, FILE_NAME_LENGTH, fp);
      fprintf(outfp, "@SQ\tSN:%s\tLN:%d\n", chrom, chromlen);
    }
    fclose(fp);

    if (readGroup[0] != 0 && sampleName[0] != 0)
      fprintf(outfp, "@RG\tID:%s\tSM:%s\tLB:%s\tPL:illumina\n", readGroup, sampleName, libName);
    
    fprintf(outfp, "@PG\tID:mrFAST\tPN:mrFAST\tVN:%s.%s\n",  versionNumber, versionNumberF);
  }
  else if (check){
    fprintf(stderr, "WARNING: %s.fai not found, the SAM file(s) will not have a header.\n", fileName[0]);
    fprintf(stderr, "You can generate the .fai file using samtools. Please place it in the same directory with the index to enable SAM headers.\n");
  }

  if (ret == NULL) 
    fprintf(stderr, "Reference genome index file read error.\n");
}

void SAMheaderGZ(gzFile outgzfp)
{
  FILE *fp;
  char fainame[FILE_NAME_LENGTH];
  char chrom[FILE_NAME_LENGTH];
  int chromlen;
  char rest[FILE_NAME_LENGTH];
  char *ret = NULL;

  sprintf(fainame, "%s.fai",fileName[0]);
  fp = fopen(fainame, "r");

  if (fp != NULL){
    gzprintf(outgzfp, "@HD\tVN:1.4\tSO:unknown\n");

    while (fscanf(fp, "%s\t%d\t", chrom, &chromlen) > 0){
      ret = fgets(rest, FILE_NAME_LENGTH, fp); 
      gzprintf(outgzfp, "@SQ\tSN:%s\tLN:%d\n", chrom, chromlen);
    }
    fclose(fp);

    if (readGroup[0] != 0 && sampleName[0] != 0)
      gzprintf(outgzfp, "@RG\tID:%s\tSM:%s\tLB:%s\tPL:illumina\n", readGroup, sampleName, libName);

    gzprintf(outgzfp, "@PG\tID:mrFAST\tPN:mrFAST\tVN:%s.%s\n",  versionNumber, versionNumberF);
  }
  else{
    fprintf(stderr, "WARNING: %s.fai not found, the SAM file(s) will not have a header.\n", fileName[0]);
    fprintf(stderr, "You can generate the .fai file using samtools. Please place it in the same directory with the index to enable SAM headers.\n");
  }

  if (ret == NULL) 
    fprintf(stderr, "Reference genome index file read error.\n");

}



