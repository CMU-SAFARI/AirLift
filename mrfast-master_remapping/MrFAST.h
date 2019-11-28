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

#ifndef __MR_FAST__
#define __MR_FAST__

#include "Reads.h"
#include <ctype.h>

#define MAP_CHUNKS 15
#define MAX_CIGAR_SIZE 100


// Pair is used to pre-processing and making the read index table
typedef struct
{
  int hv;
  int readNumber;
} Pair;

typedef struct 
{
  int loc;
  char dir;
  int err;
  float score;
  char md[MAX_CIGAR_SIZE];
  char cigar[MAX_CIGAR_SIZE];
  int cigarSize;
  int mdSize;
} FullMappingInfo;

typedef struct
{
  int loc;
  char dir;
  int err;
  float score;
  char md[MAX_CIGAR_SIZE];
  char chr[MAX_CIGAR_SIZE];
  char cigar[MAX_CIGAR_SIZE];
  int cigarSize;
  int mdSize;
  double tprob;
} BestFullMappingInfo;

typedef struct lc
{
  char md[MAP_CHUNKS][MAX_CIGAR_SIZE];
  int mdSize[MAP_CHUNKS];

  char cigar[MAP_CHUNKS][MAX_CIGAR_SIZE];
  int cigarSize[MAP_CHUNKS];

  int err[MAP_CHUNKS];
  int loc[MAP_CHUNKS];
  struct lc *next;
} MappingLocations;

typedef struct inf
{
  int size;
  MappingLocations *next;
} MappingInfo;


typedef struct 
{
  FullMappingInfo *mi;
  int size;
} FullMappingInfoLink;


extern long long			verificationCnt;
extern long long			mappingCnt;
extern long long			mappingCnt_BeforeAlignment;
extern long long			mappedSeqCnt;
extern long long			completedSeqCnt;

void initFAST(Read *, int, int *, int, char *);

void initVerifiedLocs();
void initLookUpTable();
void initBestMapping();


void finalizeFAST();
void finalizeBestSingleMapping();
void finalizeBestConcordantDiscordant();
void finalizeOEAReads(char *);


int mapAllSingleEndSeq();

void generateCigarFromMD(char *, int, char *);

/*
int backwardEditDistance2SSE2(char *a, int lena, char *b,int lenb);
int forwardEditDistance2SSE2(char *a, int lena, char *b,int lenb);
*/

double mapProb(int, char *, int, int);
int mapQ(int);

/*
int forwardEditDistanceSSE2G(char *a, int lena, char *b,int lenb);
int backwardEditDistanceSSE2G(char *a, int lena, char *b,int lenb);

int forwardEditDistance4SSE2(char *a, int lena, char *b,int lenb);
int backwardEditDistance4SSE2(char *a, int lena, char *b,int lenb);
*/

int forwardEditDistanceSSE2Extension(char *a, int lena, char *b,int lenb);
int backwardEditDistanceSSE2Extension(char *a, int lena, char *b,int lenb);


/***********************************/

/*
int verifySingleEndEditDistance(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength, 
				char *matrix, int *map_location);

int verifySingleEndEditDistance2(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength, 
				 char *matrix, int *map_location);

int verifySingleEndEditDistance4(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength, 
				 char *matrix, int *map_location);

*/

int verifySingleEndEditDistanceExtension(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength,
					 char *matrix, int *map_location);

// for fastHASH 
int compareEntrySize (const void *a, const void *b);											// fastHASH()
void mapSingleEndSeq(unsigned int *l1, int s1, int readNumber, int readSegment, int direction,	// fastHASH()
                     int index, key_struct* keys_input, int potential_key_number); 				// fastHASH()
void mapPairEndSeqList(unsigned int *l1, int s1, int readNumber, int readSegment, int direction,// fastHASH()
                       int index, key_struct* keys_input, int potential_key_number); 			// fastHASH(
void mapPairedEndSeq();
void outputPairedEnd();
void setFullMappingInfo(int readNumber, int loc, int dir, int err, int score,
			char *md, char * refName, char *cigar);

void outputAllTransChromosomal(int flag);
void outputPairedEndDiscPP();

#endif
