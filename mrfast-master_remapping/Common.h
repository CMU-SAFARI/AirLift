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


#ifndef __COMMON__
#define __COMMON__

#include <zlib.h>

#define SEQ_MAX_LENGTH		1000			// Seq Max Length
#define CONTIG_OVERLAP		2000 		// No. of characters overlapped between contings
#define CONTIG_NAME_SIZE	200			// Contig name max size
#define FILE_NAME_LENGTH	400			// Filename Max Length
#define DISCORDANT_CUT_OFF	800
#define MAX_OPEN_FILE		600		
#define MAX_TRANS_CHROMOSAL_OUTPUT 50
#define MAX_OEA_OUT		500

extern unsigned int		CONTIG_SIZE;
extern unsigned int		CONTIG_MAX_SIZE;


extern unsigned char	WINDOW_SIZE;		// WINDOW SIZE for indexing/searching
extern unsigned short	SEQ_LENGTH;		// Sequence(read) length

extern char				*versionNumber;
extern char				*versionNumberF;
extern unsigned char	mrFAST;


extern int				maxOEAOutput;
extern int				maxDiscordantOutput;
extern int				uniqueMode;
extern int				indexingMode;
extern int				searchingMode;
extern int				pairedEndModeMP;
extern int				pairedEndModePE;
extern int				pairedEndMode;
extern int				pairedEndDiscordantMode;
extern int 				transChromosomal;
extern int				pairedEndProfilingMode;
extern int				bestMode;
extern int				nosamMode;
extern int              debugMode;
extern int				seqCompressed;
extern int				outCompressed;
extern int				cropSize;
extern int				progressRep;
extern char 			*seqFile1;
extern char				*seqFile2;
extern char				*seqUnmapped;
extern char				*mappingOutput;
extern char 			*mappingOutputPath;
extern char				*unmappedOutput;
extern char             readGroup[FILE_NAME_LENGTH];
extern char             sampleName[FILE_NAME_LENGTH];
extern char             libName[FILE_NAME_LENGTH];
extern unsigned char	seqFastq;
extern unsigned char	errThreshold;
extern unsigned char	maxHits;	
extern int				minPairEndedDiscordantDistance;
extern int				maxPairEndedDiscordantDistance;
extern int				minPairEndedDistance;
extern int				maxPairEndedDistance;
extern char				fileName[2][FILE_NAME_LENGTH];
extern int				fileCnt;
extern long long		memUsage;

FILE	* fileOpen(char *fileName, char *mode);
gzFile	fileOpenGZ(char *fileName, char *mode);
double	getTime(void);
void	reverseComplement (char *seq, char *rcSeq , int length);
void	* getMem(size_t size);
void    reMem(void *, size_t, size_t);
void	freeMem(void * ptr, size_t size);
double	getMemUsage();
void 	reverse (char *seq, char *rcSeq , int length);
void 	stripPath(char *full, char **path, char **fileName);
#endif
