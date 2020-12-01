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


#ifndef __OUTPUT__
#define __OUTPUT__

#define	FORWARD 0
#define REVERSE 1

typedef struct
{
  char		*tag;
  char		type;
  char		cVal;
  int		iVal;
  float		fVal;
  char		*sVal;
} OPT_FIELDS;

typedef struct  
{
  char			*QNAME;
  short			FLAG;
  char			*RNAME;
  int			POS;
  unsigned char  	MAPQ;
  char			*CIGAR;
  char			*MRNAME;
  int			MPOS;
  int			ISIZE;
  char			*SEQ;
  char			*QUAL;
  int			optSize;
  OPT_FIELDS	        *optFields;
} SAM;

int initOutput(char *fileName, int compressed);
void (*finalizeOutput)();
void (*output)(SAM map);
void outputSAM(FILE *fp, SAM map);

FILE * getOutputFILE();
void SAMheaderTX(FILE *, int);
void SAMheaderGZ(gzFile);

#endif
