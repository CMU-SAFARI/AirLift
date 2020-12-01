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
#include <math.h>
#include <dirent.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>

#include "Common.h"
#include "Reads.h"
#include "HashTable.h"
#include "Output.h"
#include "MrFAST.h"
#include "RefGenome.h"

#define min(a,b) ((a)>(b)?(b):(a))
#define min3(a,b,c) ((a)>(b)?(b>c?c:b):(a>c?c:a))
#define CHARCODE(a) (a=='A' ? 0 : (a=='C' ? 1 : (a=='G' ? 2 : (a=='T' ? 3 : 4))))


double binomial_coefficient(int n, int k);
float calculateScore(int index, char *seq, char *qual, char *md);
unsigned char mrFAST = 1;
char *versionNumberF = "1.0";

long long verificationCnt = 0;
long long mappingCnt = 0;
long long mappedSeqCnt = 0;
long long mappingCnt_BeforeAlignment = 0;
long long completedSeqCnt = 0;
char *mappingOutput;
/**********************************************/
char *_msf_refGen = NULL;
int _msf_refGenLength = 0;
int _msf_refGenOffset = 0;
char *_msf_refGenName = NULL;

int _msf_refGenBeg;
int _msf_refGenEnd;

IHashTable *_msf_hashTable = NULL;

int *_msf_samplingLocs;
int *_msf_samplingLocsEnds;
int _msf_samplingLocsSize;

Read *_msf_seqList;
int _msf_seqListSize;

Pair *_msf_sort_seqList = NULL;

SAM _msf_output;

OPT_FIELDS *_msf_optionalFields;

char *_msf_op;

int *_msf_verifiedLocs = NULL;

char _msf_numbers[200][3];
char _msf_cigar[5];

MappingInfo *_msf_mappingInfo;

int *_msf_seqHits;
int _msf_openFiles = 0;
int _msf_maxLSize = 0;
int _msf_maxRSize = 0;

BestFullMappingInfo *bestHitMappingInfo;

/*************************/
int _msf_maxFile = 0;
char _msf_fileName[4000][200][2][FILE_NAME_LENGTH];
int _msf_fileCount[4000];

char *_msf_readHasConcordantMapping; //boolean if a read has concordant mapping

int *_msf_oeaMapping;
int *_msf_discordantMapping;



int scoreF[SEQ_MAX_LENGTH][SEQ_MAX_LENGTH];
int scoreB[SEQ_MAX_LENGTH][SEQ_MAX_LENGTH];
int score[SEQ_MAX_LENGTH][SEQ_MAX_LENGTH];
int direction1[SEQ_MAX_LENGTH][SEQ_MAX_LENGTH];
int direction2[SEQ_MAX_LENGTH][SEQ_MAX_LENGTH];

__m128i MASK;


/**************************************************Methods***************************************************/
int smallEditDistanceF(char *a, int lena, char *b, int lenb)
{
  int matrix[20][20];
  int i = 0;
  int j = 0;

  for(i = 0; i <= lena; i++)
    {
      matrix[0][i] = i;
    }

  for(i = 0; i <= lenb; i++)
    {
      matrix[i][0] = i;
    }

  for(i = 1; i <= lenb; i++)
    {
      for(j = 1; j <= lena; j++)
	{
	  matrix[i][j] = min3(matrix[i-1][j-1]+ (a[j-1] != b[i-1]),matrix[i][j-1]+1 ,matrix[i-1][j]+1);
	}
    }
  return (matrix[lenb][lena] > errThreshold ? -1 : matrix[lenb][lena]);
}

int smallEditDistanceB(char *a, int lena, char *b, int lenb)
{
  int matrix[20][20];
  int i = 0;
  int j = 0;

  for(i = 0; i <= lena; i++)
    {
      matrix[0][i] = i;
    }

  for(i = 0; i <= lenb; i++)
    {
      matrix[i][0] = i;
    }

  for(i = 1; i <= lenb; i++)
    {
      for(j = 1; j <= lena; j++)
	{
	  matrix[i][j] = min3(matrix[i-1][j-1]+ (*(a-j+1) != *(b-i+1)),matrix[i][j-1]+1 ,matrix[i-1][j]+1);
	}
    }

  return (matrix[lenb][lena] > errThreshold ? -1 : matrix[lenb][lena]);
}


void initLookUpTable()
{
  int i = 0;

  MASK = _mm_insert_epi16(MASK,1,0);
  MASK = _mm_insert_epi16(MASK,1,1);
  MASK = _mm_insert_epi16(MASK,1,2);
  MASK = _mm_insert_epi16(MASK,1,3);
  MASK = _mm_insert_epi16(MASK,1,4);
  MASK = _mm_insert_epi16(MASK,0,5);
  MASK = _mm_insert_epi16(MASK,0,6);
  MASK = _mm_insert_epi16(MASK,0,7);

  for(i = 0; i < errThreshold + 1; i++)
    {
      scoreF[0][i] = i;
      scoreF[i][0] = i;
    }

  for(i = 0 ; i < errThreshold + 1; i++)
    {
      scoreB[0][i] = i;
      scoreB[i][0] = i;
    }

}


inline int backwardEditDistanceSSE2Extension(char *a, int lena, char *b,
					     int lenb) {
  if (lenb == 0 || lena == 0)
    return 0;

  int i = 0;
  int j = 0;
  int k = 0;

  int i0;
  int i1;
  int i2;
  int i4;
  int i5;

  int e = 4;
  int mismatch = errThreshold;

  int minError = 2 * errThreshold;
  int index = 0;
  int tmpValue = 0;

  if (lenb <= e) {
    return smallEditDistanceB(a, lena, b, lenb);
  }

  __m128i R0, R1;
  __m128i Diag;
  __m128i Side1, Side2;
  __m128i Down1, Down2;
  __m128i tmp;
  __m128i SeqA, SeqB;
  __m128i Result;

  /* initialize */
  R0 = _mm_setzero_si128();
  R1 = _mm_setzero_si128();
  Diag = _mm_setzero_si128();
  Side1 = _mm_setzero_si128();
  Side2 = _mm_setzero_si128();
  Down1 = _mm_setzero_si128();
  Down2 = _mm_setzero_si128();
  SeqA = _mm_setzero_si128();
  SeqB = _mm_setzero_si128();
  Result = _mm_setzero_si128();
  /* end initialize */

  R1 = _mm_xor_si128(R1, R1);
  R0 = _mm_xor_si128(R0, R0);

  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag,minError,0);

  i0 = (a[0] != b[0]);
  i1 = min(i0, ( *(a-1)!=*b) ) + 1;
  i2 = min(i0,( a[0] != *(b-1) ) ) + 1;

  i0 = min3( i0+ ( *(a-1)!=*(b-1) ),i1+1,i2+1);
  i4 = min(i1, ( *(a-2)!=b[0] )+1) + 1;
  i5 = min(i2, (a[0] != *(b-2))+1) + 1;

  R1 = _mm_insert_epi16(R1, 3, 0);
  R1 = _mm_insert_epi16(R1, i1, 1);
  R1 = _mm_insert_epi16(R1, i2, 2);
  R1 = _mm_insert_epi16(R1, 3, 3);

  R0 = _mm_insert_epi16(R0, 4, 0);
  R0 = _mm_insert_epi16(R0, i4, 1);
  R0 = _mm_insert_epi16(R0, i0, 2);
  R0 = _mm_insert_epi16(R0, i5, 3);
  R0 = _mm_insert_epi16(R0, 4, 4);

  Side2 = _mm_xor_si128(Side2, Side2);
  Down2 = _mm_xor_si128(Down2, Down2);
  Down1 = _mm_xor_si128(Down1, Down1);
  Side1 = _mm_xor_si128(Side1, Side1);

  Side2 = _mm_insert_epi16(Side2,minError,0);
  Down1 = _mm_insert_epi16(Down1,minError,0);

  Side1 = _mm_insert_epi16(Side1,1,0);

  index = 0;
  for (j = 0; j < e; j++) {
    Side2 = _mm_slli_si128(Side2, 2);
    Side2 = _mm_insert_epi16(Side2,1,0);

    Down1 = _mm_slli_si128(Down1, 2);
    Down1 = _mm_insert_epi16(Down1,1,0);

    Down2 = _mm_slli_si128(Down2, 2);
    Down2 = _mm_insert_epi16(Down2,1,0);

    Side1 = _mm_slli_si128(Side1, 2);
    Side1 = _mm_insert_epi16(Side1,1,0);

    SeqA = _mm_slli_si128(SeqA, 2);
    SeqB = _mm_slli_si128(SeqB, 2);
    SeqA = _mm_insert_epi16(SeqA,*(a-index),0);
    SeqB = _mm_insert_epi16(SeqB,*(b-index),0);
    index++;
  }

  Down2 = _mm_slli_si128(Down2, 2);
  Down2 = _mm_insert_epi16(Down2,minError,0);

  index = 4;
  i = 5;

  int loopEnd = 2 * lenb - (e - 1);
  for (; i <= loopEnd; i++) {

    Diag = _mm_xor_si128(Diag, Diag);
    if (i % 2 == 0) {
      SeqA = _mm_slli_si128(SeqA, 2);
      SeqB = _mm_slli_si128(SeqB, 2);
      SeqA = _mm_insert_epi16(SeqA,*(a-(index)),0);
      SeqB = _mm_insert_epi16(SeqB,*(b-(index)),0);

      index++;

      tmp = _mm_shufflelo_epi16(SeqB,27);
      tmp = _mm_slli_si128(tmp, 2);
      tmpValue = _mm_extract_epi16(tmp, 5);
      tmp = _mm_insert_epi16(tmp, tmpValue, 0);

      Result = _mm_cmpeq_epi16(SeqA, tmp);
      Diag = _mm_andnot_si128(Result, MASK);

      R0 = _mm_min_epi16(_mm_add_epi16(R1,Side2), _mm_add_epi16(R0,Diag));
      R0 = _mm_min_epi16(R0, _mm_add_epi16(_mm_slli_si128(R1,2) ,Down2));

      if (_mm_extract_epi16(R0, 0) > errThreshold
	  && _mm_extract_epi16(R0, 1) > errThreshold
	  && _mm_extract_epi16(R0, 2) > errThreshold
	  && _mm_extract_epi16(R0, 3) > errThreshold
	  && _mm_extract_epi16(R0, 4) > errThreshold
	  && _mm_extract_epi16(R1, 0) > errThreshold
	  && _mm_extract_epi16(R1, 1) > errThreshold
	  && _mm_extract_epi16(R1, 2) > errThreshold
	  && _mm_extract_epi16(R1, 3) > errThreshold)
	return -1;

      if (i == 2 * lenb - e) {
	tmp = _mm_srli_si128(R0,2);
	for (k = 0; k < e - 1; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = _mm_extract_epi16(tmp,0);
      }

    }

    else {
      Result = _mm_cmpeq_epi16(SeqA, _mm_shufflelo_epi16(SeqB,27));
      Diag = _mm_andnot_si128(Result, MASK);

      R1 = _mm_min_epi16(_mm_add_epi16(_mm_srli_si128(R0,2),Side1), _mm_add_epi16(R1,Diag));
      R1 = _mm_min_epi16(R1, _mm_add_epi16(R0 ,Down1));

      if (i >= 2 * lenb - e) {
	tmp = _mm_srli_si128(R1,2);
	for (k = 0; k < e - 2; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = min(minError, _mm_extract_epi16(tmp,0));
      }
    }

  }

  j = 0;
  int tmpE = e;
  for (; j < 2 * (e - 2) + 1; j++) {

    Diag = _mm_xor_si128(Diag, Diag);
    //set the first element
    if (j == 0) {
      for (k = 0; k <= e - 1; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(lenb-1-k)) != *(a-((i-lenb)-1+k)),0);
      }

      R0 = _mm_min_epi16(_mm_add_epi16(R1,Side2), _mm_add_epi16(R0,Diag));
      R0 = _mm_min_epi16(R0, _mm_add_epi16(_mm_slli_si128(R1,2) ,Down2));

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < e - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    } else if (j % 2 == 0) {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(lenb-1-k)) != *(a-((i-lenb)-1+k)),0);
      }

      R0 = _mm_min_epi16(_mm_add_epi16(R1,Side2), _mm_add_epi16(R0,Diag));
      R0 = _mm_min_epi16(R0, _mm_add_epi16(_mm_slli_si128(R1,2) ,Down2));

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < tmpE - 1; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }

    else {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(lenb-1-k)) != *(a-((i-lenb)-1+k)),0);
      }

      R1 = _mm_min_epi16(_mm_add_epi16(_mm_srli_si128(R0,2),Side1), _mm_add_epi16(R1,Diag));
      R1 = _mm_min_epi16(R1, _mm_add_epi16(R0 ,Down1));

      tmp = _mm_srli_si128(R1,2);
      for (k = 0; k < tmpE - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }
    i++;
  }
  //Diag

  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag, 2*errThreshold, 0);
  Diag = _mm_insert_epi16(Diag, *(a-(lenb+e-2)) != *(b-(lenb-1)), 1);

  Side1 = _mm_insert_epi16(Side1,1,0);
  Side1 = _mm_insert_epi16(Side1,1,1);

  Down1 = _mm_insert_epi16(Down1, 2*errThreshold, 0);
  Down1 = _mm_insert_epi16(Down1, 1, 1);

  R1 = _mm_min_epi16(_mm_add_epi16(R0,Side1), _mm_add_epi16(_mm_slli_si128(R1,2),Diag));
  R1 = _mm_min_epi16(R1, _mm_add_epi16(_mm_slli_si128(R0,2),Down1));

  minError = min(minError, _mm_extract_epi16(R1,1));

  Diag = _mm_insert_epi16(Diag, *(a-(lenb+e-1)) != *(b-(lenb-1)), 0);
  Down1 = _mm_insert_epi16(Down1, 1, 0);

  R0 = _mm_min_epi16(_mm_add_epi16(R1,Down1), _mm_add_epi16(R0,Diag));
  R0 = _mm_min_epi16(R0, _mm_add_epi16(_mm_srli_si128(R1,2) ,Side1));

  minError = min(minError, _mm_extract_epi16(R0,0));

  if (minError > mismatch)
    return -1;
  return minError;
}



inline int forwardEditDistanceSSE2Extension(char *a, int lena, char *b,
					    int lenb) {
  if (lenb == 0 || lena == 0)
    return 0;

  int i = 0;
  int j = 0;
  int k = 0;

  int i0 = 0;
  int i1 = 0;
  int i2 = 0;
  int i4 = 0;
  int i5 = 0;

  int mismatch = errThreshold;
  int e = 4;

  int minError = 4 * mismatch + 1;
  int index = 0;
  int tmpValue = 0;

  if (lenb <= e) {
    return smallEditDistanceF(a, lena, b, lenb);
  }

  register __m128i R0, R1;
  __m128i Diag;
  __m128i Side1, Side2;
  __m128i Down1, Down2;
  __m128i tmp;
  register __m128i SeqA, SeqB;
  __m128i Result;

  __m128i tmpSeqA;
  __m128i tmpSeqB;

  /* initialize */
  R0 = _mm_setzero_si128();
  R1 = _mm_setzero_si128();
  Diag = _mm_setzero_si128();
  Side1 = _mm_setzero_si128();
  Side2 = _mm_setzero_si128();
  Down1 = _mm_setzero_si128();
  Down2 = _mm_setzero_si128();
  SeqA = _mm_setzero_si128();
  SeqB = _mm_setzero_si128();
  Result = _mm_setzero_si128();
  /* end initialize */

  R1 = _mm_xor_si128(R1, R1);
  R0 = _mm_xor_si128(R0, R0);

  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag,minError,0);

  i0 = (a[0] != b[0]);
  i1 = min(i0, (a[1]!=b[0])) + 1;
  i2 = min(i0,(a[0]!=b[1])) + 1;

  i0 = min3(i0+(a[1]!=b[1]),i1+1,i2+1);
  i4 = min(i1, (a[2]!=b[0])+1) + 1;
  i5 = min(i2, (a[0]!=b[2])+1) + 1;

  R1 = _mm_insert_epi16(R1, 3, 0);
  R1 = _mm_insert_epi16(R1, i1, 1);
  R1 = _mm_insert_epi16(R1, i2, 2);
  R1 = _mm_insert_epi16(R1, 3, 3);

  R0 = _mm_insert_epi16(R0, 4, 0);
  R0 = _mm_insert_epi16(R0, i4, 1);
  R0 = _mm_insert_epi16(R0, i0, 2);
  R0 = _mm_insert_epi16(R0, i5, 3);
  R0 = _mm_insert_epi16(R0, 4, 4);

  Side2 = _mm_xor_si128(Side2, Side2);
  Down2 = _mm_xor_si128(Down2, Down2);
  Down1 = _mm_xor_si128(Down1, Down1);
  Side1 = _mm_xor_si128(Side1, Side1);

  Side2 = _mm_insert_epi16(Side2,minError,0);
  Down1 = _mm_insert_epi16(Down1,minError,0);

  Side1 = _mm_insert_epi16(Side1,1,0);

  index = 0;
  for (j = 0; j < e; j++) {
    Side2 = _mm_slli_si128(Side2, 2);
    Side2 = _mm_insert_epi16(Side2,1,0);

    Down1 = _mm_slli_si128(Down1, 2);
    Down1 = _mm_insert_epi16(Down1,1,0);

    Down2 = _mm_slli_si128(Down2, 2);
    Down2 = _mm_insert_epi16(Down2,1,0);

    Side1 = _mm_slli_si128(Side1, 2);
    Side1 = _mm_insert_epi16(Side1,1,0);

    SeqA = _mm_slli_si128(SeqA, 2);
    SeqB = _mm_slli_si128(SeqB, 2);
    SeqA = _mm_insert_epi16(SeqA,a[index],0);
    SeqB = _mm_insert_epi16(SeqB,b[index],0);
    index++;
  }

  Down2 = _mm_slli_si128(Down2, 2);
  Down2 = _mm_insert_epi16(Down2,minError,0);

  index = 4;
  i = 5;

  int loopEnd = 2 * lenb - (e - 1);
  for (; i <= loopEnd; i++) {
    if (i % 2 == 0) {
      tmpSeqA = _mm_slli_si128(SeqA, 2);
      tmpSeqB = _mm_slli_si128(SeqB, 2);
      SeqA = _mm_insert_epi16(tmpSeqA,a[index],0);
      SeqB = _mm_insert_epi16(tmpSeqB,b[index],0);

      index++;

      tmp = _mm_shufflelo_epi16(SeqB,27);
      tmp = _mm_slli_si128(tmp, 2);
      tmpValue = _mm_extract_epi16(tmp, 5);
      tmp = _mm_insert_epi16(tmp, tmpValue, 0);

      Result = _mm_cmpeq_epi16(SeqA, tmp);
      Diag = _mm_andnot_si128(Result, MASK);

      R0 = _mm_min_epi16(_mm_add_epi16(R1,Side2), _mm_add_epi16(R0,Diag));
      R0 = _mm_min_epi16(R0, _mm_add_epi16(_mm_slli_si128(R1,2) ,Down2));

      if (_mm_extract_epi16(R0, 0) > errThreshold
	  && _mm_extract_epi16(R0, 1) > errThreshold
	  && _mm_extract_epi16(R0, 2) > errThreshold
	  && _mm_extract_epi16(R0, 3) > errThreshold
	  && _mm_extract_epi16(R0, 4) > errThreshold
	  && _mm_extract_epi16(R1, 0) > errThreshold
	  && _mm_extract_epi16(R1, 1) > errThreshold
	  && _mm_extract_epi16(R1, 2) > errThreshold
	  && _mm_extract_epi16(R1, 3) > errThreshold)
	return -1;

      if (i == 2 * lenb - e) {
	tmp = _mm_srli_si128(R0,2);
	for (k = 0; k < e - 1; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = _mm_extract_epi16(tmp,0);
      }

    }

    else {
      Result = _mm_cmpeq_epi16(SeqA, _mm_shufflelo_epi16(SeqB,27));
      Diag = _mm_andnot_si128(Result, MASK);

      R1 = _mm_min_epi16(_mm_add_epi16(_mm_srli_si128(R0,2),Side1), _mm_add_epi16(R1,Diag));
      R1 = _mm_min_epi16(R1, _mm_add_epi16(R0 ,Down1));

      if (i >= 2 * lenb - e) {
	tmp = _mm_srli_si128(R1,2);
	for (k = 0; k < e - 2; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = min(minError, _mm_extract_epi16(tmp,0));
      }
    }
  }

  j = 0;
  int tmpE = e;
  for (; j < 2 * (e - 2) + 1; j++) {

    Diag = _mm_xor_si128(Diag, Diag);
    //set the first element
    if (j == 0) {
      for (k = 0; k <= e - 1; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag = _mm_insert_epi16(Diag, b[lenb-1-k] != a[(i-lenb)-1+k],0);
      }

      R0 = _mm_min_epi16(_mm_add_epi16(R1,Side2), _mm_add_epi16(R0,Diag));
      R0 = _mm_min_epi16(R0, _mm_add_epi16(_mm_slli_si128(R1,2) ,Down2));

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < e - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    } else if (j % 2 == 0) {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag = _mm_insert_epi16(Diag, b[lenb-1-k] != a[(i-lenb)-1+k],0);
      }

      R0 = _mm_min_epi16(_mm_add_epi16(R1,Side2), _mm_add_epi16(R0,Diag));
      R0 = _mm_min_epi16(R0, _mm_add_epi16(_mm_slli_si128(R1,2) ,Down2));

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < tmpE - 1; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }

    else {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag = _mm_insert_epi16(Diag, b[lenb-1-k] != a[(i-lenb)-1+k],0);
      }

      R1 = _mm_min_epi16(_mm_add_epi16(_mm_srli_si128(R0,2),Side1), _mm_add_epi16(R1,Diag));
      R1 = _mm_min_epi16(R1, _mm_add_epi16(R0 ,Down1));

      tmp = _mm_srli_si128(R1,2);
      for (k = 0; k < tmpE - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }
    i++;
  }
  //Diag

  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag, minError, 0);
  Diag = _mm_insert_epi16(Diag, a[lenb+e-2] != b[lenb-1], 1);

  Side1 = _mm_insert_epi16(Side1,1,0);
  Side1 = _mm_insert_epi16(Side1,1,1);

  Down1 = _mm_insert_epi16(Down1, minError, 0);
  Down1 = _mm_insert_epi16(Down1, 1, 1);

  R1 = _mm_min_epi16(_mm_add_epi16(R0,Side1), _mm_add_epi16(_mm_slli_si128(R1,2),Diag));
  R1 = _mm_min_epi16(R1, _mm_add_epi16(_mm_slli_si128(R0,2),Down1));

  minError = min(minError, _mm_extract_epi16(R1,1));

  Diag = _mm_insert_epi16(Diag, a[lenb+e-1] != b[lenb-1], 0);
  Down1 = _mm_insert_epi16(Down1, 1, 0);

  R0 = _mm_min_epi16(_mm_add_epi16(R1,Down1), _mm_add_epi16(R0,Diag));
  R0 = _mm_min_epi16(R0, _mm_add_epi16(_mm_srli_si128(R1,2) ,Side1));

  minError = min(minError, _mm_extract_epi16(R0,0));

  if (minError > mismatch)
    return -1;
  return minError;
}



void initBestMapping(int totalReadNumber)
{
  int i = 0;
  bestHitMappingInfo = getMem(totalReadNumber * sizeof(BestFullMappingInfo));
  for (i = 0; i < totalReadNumber; i++) {
    bestHitMappingInfo[i].loc = -1;
    bestHitMappingInfo[i].tprob = 0.0; 
  }
}

void finalizeBestSingleMapping() 
{
  int i = 0;
  char *_tmpQual, *_tmpSeq;
  char rqual[SEQ_LENGTH + 1];
  rqual[SEQ_LENGTH] = '\0';

  for(i = 0; i < _msf_seqListSize; i++)
    {
      if(_msf_seqList[i].hits[0] != 0)
	{		
	  if (bestHitMappingInfo[i].dir)
	    {
	      reverse(_msf_seqList[i].qual, rqual, SEQ_LENGTH);
	      _tmpQual = rqual;
	      _tmpSeq = _msf_seqList[i].rseq;
	    }
	  else
	    {
	      _tmpQual = _msf_seqList[i].qual;
	      _tmpSeq = _msf_seqList[i].seq;
	    }

	  _msf_output.QNAME = _msf_seqList[i].name;
	  _msf_output.FLAG = 16 * bestHitMappingInfo[i].dir;
	  _msf_output.RNAME = bestHitMappingInfo[i].chr;

	  _msf_output.POS = bestHitMappingInfo[i].loc;

	  if (seqFastq)
	    _msf_output.MAPQ = mapQ(i);	  
	  else
	    _msf_output.MAPQ = 255;

	  _msf_output.CIGAR = bestHitMappingInfo[i].cigar;
	  _msf_output.MRNAME = "*";
	  _msf_output.MPOS = 0;
	  _msf_output.ISIZE = 0;

	  _msf_output.SEQ = _tmpSeq;
	  _msf_output.QUAL = _tmpQual;

	  _msf_output.optSize = 2;
	  _msf_output.optFields = _msf_optionalFields;

	  _msf_optionalFields[0].tag = "NM";
	  _msf_optionalFields[0].type = 'i';
	  _msf_optionalFields[0].iVal = bestHitMappingInfo[i].err;

	  _msf_optionalFields[1].tag = "MD";
	  _msf_optionalFields[1].type = 'Z';
	  _msf_optionalFields[1].sVal = bestHitMappingInfo[i].md;

	  output(_msf_output);
	}
    }
  freeMem(bestHitMappingInfo, _msf_seqListSize * sizeof(FullMappingInfo));
}
/**********************************************/
int compare(const void *a, const void *b) {
  return ((Pair *) a)->hv - ((Pair *) b)->hv;
}
/**********************************************/
void preProcessReads() {
  int i = 0;

  _msf_sort_seqList = getMem(_msf_seqListSize * sizeof(Pair));
  for (i = 0; i < _msf_seqListSize; i++) {
    _msf_sort_seqList[i].hv = hashVal(_msf_seqList[i].seq);
    _msf_sort_seqList[i].readNumber = i;    
  }
  
  qsort(_msf_sort_seqList, _msf_seqListSize, sizeof(Pair), compare);

}
/**********************************************/

int verifySingleEnd(int index, char* seq, int offset) {
  int curOff = 0;
  int i;

  char *ref;

  int err;
  int errCnt = 0;
  int errCntOff = 0;
  int NCntOff = 0;

  ref = _msf_refGen + index - 1;

  verificationCnt++;

  for (i = 0; i < SEQ_LENGTH; i++) {
    err = *ref != *seq;
    errCnt += err;
    if (errCnt > errThreshold) {

      return -1;
    }

    if (i >= _msf_samplingLocs[curOff]
	&& i <= _msf_samplingLocsEnds[curOff]) {
      errCntOff += err;
      NCntOff += (*seq == 'N');
    } else if (curOff < _msf_samplingLocsSize
	       && i >= _msf_samplingLocs[curOff + 1]) {

      if (errCntOff == 0 && NCntOff == 0 && offset > curOff) {
	return -1;
      }

      errCntOff = 0;
      NCntOff = 0;
      curOff++;

      if (i >= _msf_samplingLocs[curOff]) {
	errCntOff += err;
	NCntOff += (*seq == 'N');
      }
    }

    ref++;
    seq++;
  }
  return errCnt;
}

/*********************************************/
void initFAST(Read *seqList, int seqListSize, int *samplingLocs,
	      int samplingLocsSize, char *genFileName) {
  int i;

  if (_msf_optionalFields == NULL) {
    _msf_op = getMem(SEQ_LENGTH);
    if (pairedEndMode) {
      _msf_optionalFields = getMem(8 * sizeof(OPT_FIELDS));
    } else {
      _msf_optionalFields = getMem(2 * sizeof(OPT_FIELDS));
    }

    for (i = 0; i < 200; i++) {
      sprintf(_msf_numbers[i], "%d%c", i, '\0');
    }
    sprintf(_msf_cigar, "%dM", SEQ_LENGTH);
  }

  if (_msf_samplingLocsEnds == NULL) {
    _msf_samplingLocs = samplingLocs;
    _msf_samplingLocsSize = samplingLocsSize;

    _msf_samplingLocsEnds = getMem(sizeof(int) * _msf_samplingLocsSize);
    for (i = 0; i < _msf_samplingLocsSize; i++) {
      _msf_samplingLocsEnds[i] = _msf_samplingLocs[i] + WINDOW_SIZE - 1;
    }

    _msf_seqList = seqList;
    _msf_seqListSize = seqListSize;

    preProcessReads();

    _msf_oeaMapping = getMem(_msf_seqListSize * sizeof(int));
    for (i = 0; i < _msf_seqListSize; i++) {
      _msf_oeaMapping[i] = 0;
    }

    _msf_discordantMapping = getMem(_msf_seqListSize * sizeof(int));
    for (i = 0; i < _msf_seqListSize; i++) {
      _msf_discordantMapping[i] = 0;
    }

  }

  if (_msf_refGenName == NULL) {
    _msf_refGenName = getMem(4 * SEQ_LENGTH);
  }
  _msf_refGen = getRefGenome();
  _msf_refGenLength = strlen(_msf_refGen);

  _msf_refGenOffset = getRefGenomeOffset();
  snprintf(_msf_refGenName, 4 * SEQ_LENGTH, "%s%c", getRefGenomeName(), '\0');
  _msf_refGenName[strlen(getRefGenomeName())] = '\0';

  if (_msf_verifiedLocs != NULL) {
    freeMem(_msf_verifiedLocs, sizeof(int) * (_msf_refGenLength + 1));
  }

  _msf_verifiedLocs = (int *) getMem(sizeof(int) * (_msf_refGenLength + 1));

  for (i = 0; i <= _msf_refGenLength; i++)
    _msf_verifiedLocs[i] = _msf_seqListSize * 10 + 1;

  if (pairedEndMode && _msf_seqHits == NULL) {

    _msf_mappingInfo = getMem(seqListSize * sizeof(MappingInfo));

    for (i = 0; i < seqListSize; i++) {
      _msf_mappingInfo[i].next = NULL;
      _msf_mappingInfo[i].size = 0;
    }

    _msf_seqHits = getMem((_msf_seqListSize) * sizeof(int));

    for (i = 0; i < _msf_seqListSize; i++) {
      _msf_seqHits[i] = 0;
    }

    _msf_readHasConcordantMapping = getMem(
					   _msf_seqListSize / 2 * sizeof(char));
    for (i = 0; i < _msf_seqListSize / 2; i++) {
      _msf_readHasConcordantMapping[i] = 0;
    }

    initLoadingRefGenome(genFileName);
  }

  if (_msf_refGenOffset == 0) {
    _msf_refGenBeg = 1;
  } else {
    _msf_refGenBeg = CONTIG_OVERLAP - SEQ_LENGTH + 2 + errThreshold;
  }
  _msf_refGenEnd = _msf_refGenLength - SEQ_LENGTH + 1;

}
/**********************************************/
void finalizeFAST() {
  freeMem(_msf_seqHits, (_msf_seqListSize) * sizeof(int));
  freeMem(_msf_refGenName, 4 * SEQ_LENGTH);

  freeMem(_msf_sort_seqList, sizeof(int) * _msf_seqListSize);

}


int verifySingleEndEditDistanceExtension(int refIndex, char *lSeq,
					 int lSeqLength, char *rSeq, int rSeqLength, int segLength, char *matrix,
					 int *map_location) {
  int i = 0;

  char * ref;
  char * tempref;

  int rIndex = 0; //reference Index

  int error = 0;
  int error1 = 0;

  int error2 = 0;
  int error3 = 0;
  int totalError = 0;
  //int errorSegment = 0;

  int ERROR_BOUND = min(4, errThreshold);

  /*
    1: Up
    2: Side
    3: Diagonal Match
    4: Diagonal Mismatch
  */

  int min = 0;
  int minIndex1 = 0;
  int minIndex2 = 0;

  int directionIndex = 0;

  int size = 0;

  ref = _msf_refGen + refIndex - 1;
  tempref = _msf_refGen + refIndex - 1;

  if (lSeqLength != 0) {
    error3 = backwardEditDistanceSSE2Extension(ref - 1, lSeqLength,
					       lSeq + lSeqLength - 1, lSeqLength);
    if (error3 == -1) {
      return -1;
    }
  }

  if (rSeqLength != 0) {
    error2 = forwardEditDistanceSSE2Extension(ref + segLength, rSeqLength,
					      rSeq, rSeqLength);
    if (error2 == -1)
      return -1;
  }

  if (error2 + error3 > errThreshold)
    return -1;

  rIndex = 1;

  //int prevError = 0;

  int tempUp = 0;
  int tempDown = 0;

  int errorString = 0;

  int upValue;
  int diagValue;
  int sideValue;
  if (lSeqLength > ERROR_BOUND) {
    while (rIndex <= lSeqLength + ERROR_BOUND && lSeqLength != 0) {
      tempUp = (
		(rIndex - ERROR_BOUND) > 0 ?
		((rIndex > lSeqLength) ?
		 lSeqLength - ERROR_BOUND :
		 rIndex - ERROR_BOUND) :
		1);
      tempDown = (
		  (rIndex >= lSeqLength - ERROR_BOUND) ?
		  lSeqLength + 1 : rIndex + ERROR_BOUND + 1);
      for (i = tempUp; i < tempDown; i++) {
	errorString = (*(ref - rIndex) == *(lSeq + lSeqLength - i));

	upValue = scoreB[i - 1][rIndex] + 1;
	diagValue = scoreB[i - 1][rIndex - 1] + !errorString;
	sideValue = scoreB[i][rIndex - 1] + 1;

	if (i != tempUp && i != tempDown - 1)
	  scoreB[i][rIndex] = min3(sideValue, diagValue , upValue);

	else if ((i
		  == ((rIndex - ERROR_BOUND) > 0 ?
		      rIndex - ERROR_BOUND : 1))
		 && rIndex <= lSeqLength)
	  scoreB[i][rIndex] = min(sideValue, diagValue);
	else if (rIndex > lSeqLength && (i == lSeqLength - ERROR_BOUND))
	  scoreB[i][rIndex] = sideValue;
	else
	  scoreB[i][rIndex] = min(diagValue , upValue);

	if (i == tempUp)
	  error = scoreB[i][rIndex];
	else if (error > scoreB[i][rIndex])
	  error = scoreB[i][rIndex];
      }
      if (rIndex <= lSeqLength) {
	//errorSegment = error-prevError;
      }
      rIndex++;
    }

    if (lSeqLength != 0) {
      min = scoreB[lSeqLength][lSeqLength + ERROR_BOUND];
      minIndex1 = lSeqLength + ERROR_BOUND;

      // Find the Best error for all the possible ways.
      for (i = 1; i <= 2 * ERROR_BOUND; i++) {
	if (min >= scoreB[lSeqLength][lSeqLength + ERROR_BOUND - i]
	    && lSeqLength + ERROR_BOUND - i > 0) {
	  min = scoreB[lSeqLength][lSeqLength + ERROR_BOUND - i];
	  minIndex1 = lSeqLength + ERROR_BOUND - i;
	}
      }
      error = scoreB[lSeqLength][minIndex1];
    }
  } else {
    int j = 0;
    for (i = 1; i <= lSeqLength; i++) {
      for (j = 1; j <= lSeqLength; j++) {
	scoreB[i][j] =
	  min3(scoreB[i-1][j-1]+ (*(ref-j) != *(lSeq+lSeqLength-i) ),scoreB[i][j-1]+1 ,scoreB[i-1][j]+1);
      }
    }
    error = scoreB[lSeqLength][lSeqLength];
    minIndex1 = lSeqLength;

  }
  error1 = error;

  error = 0;
  //errorSegment = 0;

  directionIndex = lSeqLength;
  rIndex = minIndex1;

  *map_location = ((lSeqLength == 0) ? refIndex : refIndex - rIndex);

  ref = ref + segLength;

  if (rSeqLength != 0 && rSeqLength > ERROR_BOUND) {
    ERROR_BOUND = min(ERROR_BOUND, rSeqLength);

    if (rSeqLength == ERROR_BOUND) {
      for (i = 0; i < 2 * ERROR_BOUND; i++)
	scoreF[0][i] = i;
    }

    rIndex = 1;
    while (rIndex <= rSeqLength + ERROR_BOUND) {
      tempUp =
	(rIndex - ERROR_BOUND) > 0 ?
	((rIndex > rSeqLength) ?
	 rSeqLength - ERROR_BOUND :
	 rIndex - ERROR_BOUND) :
	1;
      tempDown = (
		  (rIndex >= rSeqLength - ERROR_BOUND) ?
		  rSeqLength + 1 : rIndex + ERROR_BOUND + 1);
      for (i = tempUp; i < tempDown; i++) {
	errorString = (*(ref + rIndex - 1) == *(rSeq + i - 1));
	upValue = scoreF[i - 1][rIndex] + 1;
	diagValue = scoreF[i - 1][rIndex - 1] + !errorString;
	sideValue = scoreF[i][rIndex - 1] + 1;

	if (i != tempUp && i != tempDown - 1)
	  scoreF[i][rIndex] = min3(sideValue, diagValue , upValue);
	else if ((i
		  == ((rIndex - ERROR_BOUND) > 0 ?
		      rIndex - ERROR_BOUND : 1))
		 && rIndex <= rSeqLength)
	  scoreF[i][rIndex] = min(sideValue, diagValue);
	else if (rIndex > rSeqLength && (i == rSeqLength - ERROR_BOUND))
	  scoreF[i][rIndex] = sideValue;
	else
	  scoreF[i][rIndex] = min(diagValue , upValue);

	if (i == tempUp)
	  error = scoreF[i][rIndex];
	if (error > scoreF[i][rIndex])
	  error = scoreF[i][rIndex];
      }
      if (rIndex <= rSeqLength) {
	//errorSegment = error;
      }
      rIndex++;
    }
    min = scoreF[rSeqLength][rSeqLength + ERROR_BOUND];
    minIndex2 = rSeqLength + ERROR_BOUND;

    // Find the Best error for all the possible ways.
    for (i = 1; i <= 2 * ERROR_BOUND; i++) {
      if (min > scoreF[rSeqLength][rSeqLength + ERROR_BOUND - i]
	  && rSeqLength + ERROR_BOUND - i > 0) {
	min = scoreF[rSeqLength][rSeqLength + ERROR_BOUND - i];
	minIndex2 = rSeqLength + ERROR_BOUND - i;
      }
    }
    error = scoreF[rSeqLength][minIndex2];
  } else {
    int j = 0;
    for (i = 1; i <= rSeqLength; i++) {
      for (j = 1; j <= rSeqLength; j++) {
	scoreF[i][j] =
	  min3(scoreF[i-1][j-1]+ (*(ref+j-1) != *(rSeq+i-1) ),scoreF[i][j-1]+1 ,scoreF[i-1][j]+1);
      }
    }
    error = scoreF[rSeqLength][rSeqLength];
    minIndex2 = rSeqLength;
  }

  totalError = error + error1;

  /* Farhad 08/07/2012 */
  if(totalError > errThreshold)
    return -1;
  /* Farhad 08/07/2012 */

  if (debugMode && totalError != error2 + error3) {
    for (i = 0; i < lSeqLength; i++)
      fprintf(stderr, "%c", *(tempref - 1 - i));
    fprintf(stderr, "\n");
    for (i = 0; i < lSeqLength; i++)
      fprintf(stderr, "%c", *(lSeq + i));
    fprintf(stderr, "\n");

    for (i = 0; i < rSeqLength; i++)
      fprintf(stderr, "%c", *(tempref + segLength + i));
    fprintf(stderr, "\n");

    for (i = 0; i < rSeqLength; i++)
      fprintf(stderr, "%c", *(rSeq + i));
    fprintf(stderr, "\n");

    fprintf(stderr, "ERROR=%d\n", totalError);
    fprintf(stderr, "ERROR_SSE=%d\n", error3 + error2);

    fprintf(stderr, "ERROR_SSE_back=%d E_SSE_forw=%d\n", error3, error2);
    fprintf(stderr, "ERROR_back=%d E_forw=%d\n", error1, error);

  }

  char matrixR[SEQ_MAX_LENGTH];
  char matrixL[SEQ_MAX_LENGTH];

  matrixR[0] = '\0';
  matrixL[0] = '\0';

  size = 0;
  directionIndex = rSeqLength;
  rIndex = minIndex2;

  /* Farhad 29/11/2012 */
  while (directionIndex > 0 || rIndex > 0) {
    if (directionIndex - rIndex == errThreshold) {
      if (scoreF[directionIndex][rIndex]
	  - scoreF[directionIndex - 1][rIndex] == 1) {
	matrixR[size] = *(rSeq + directionIndex - 1);
	size++;
	matrixR[size] = 'I';
	directionIndex--;
      } else if (scoreF[directionIndex][rIndex]
		 - scoreF[directionIndex - 1][rIndex - 1] == 1) {
	matrixR[size] = *(ref + rIndex - 1);
	rIndex--;
	directionIndex--;
      } else {
	matrixR[size] = 'M';
	rIndex--;
	directionIndex--;
      }

    } else if (rIndex - directionIndex == errThreshold) {
      if (scoreF[directionIndex][rIndex]
	  - scoreF[directionIndex][rIndex - 1] == 1) {
	matrixR[size] = *(ref + rIndex - 1);
	size++;
	matrixR[size] = 'D';
	rIndex--;
      } else if (scoreF[directionIndex][rIndex]
		 - scoreF[directionIndex - 1][rIndex - 1] == 1) {
	matrixR[size] = *(ref + rIndex - 1);
	rIndex--;
	directionIndex--;
      } else {
	matrixR[size] = 'M';
	rIndex--;
	directionIndex--;
      }
    } else {
      if (scoreF[directionIndex][rIndex]
	  - scoreF[directionIndex - 1][rIndex] == 1
	  && directionIndex != 0) {
	matrixR[size] = *(rSeq + directionIndex - 1);
	size++;
	matrixR[size] = 'I';
	directionIndex--;
      } else if (scoreF[directionIndex][rIndex]
		 - scoreF[directionIndex][rIndex - 1] == 1 && rIndex != 0) {
	matrixR[size] = *(ref + rIndex - 1);
	size++;
	matrixR[size] = 'D';
	rIndex--;
      } else if (scoreF[directionIndex][rIndex]
		 - scoreF[directionIndex - 1][rIndex - 1] == 1) {
	matrixR[size] = *(ref + rIndex - 1);
	rIndex--;
	directionIndex--;
      } else {
	matrixR[size] = 'M';
	rIndex--;
	directionIndex--;
      }
    }
    size++;
  }
  matrixR[size] = '\0';

  size = 0;
  directionIndex = lSeqLength;
  rIndex = minIndex1;

  while (directionIndex != 0 || rIndex != 0) {
    if (directionIndex - rIndex == errThreshold) {
      if (scoreB[directionIndex][rIndex]
	  - scoreB[directionIndex - 1][rIndex] == 1) {
	matrixL[size] = 'I';
	size++;
	matrixL[size] = *(lSeq + lSeqLength - directionIndex);
	directionIndex--;
      } else if (scoreB[directionIndex][rIndex]
		 - scoreB[directionIndex - 1][rIndex - 1] == 1) {
	matrixL[size] = *(tempref - rIndex);
	rIndex--;
	directionIndex--;
      } else {
	matrixL[size] = 'M';
	rIndex--;
	directionIndex--;
      }

    } else if (rIndex - directionIndex == errThreshold) {
      if (scoreB[directionIndex][rIndex]
	  - scoreB[directionIndex][rIndex - 1] == 1) {
	matrixL[size] = 'D';
	size++;
	matrixL[size] = *(tempref - rIndex);
	rIndex--;
      } else if (scoreB[directionIndex][rIndex]
		 - scoreB[directionIndex - 1][rIndex - 1] == 1) {
	matrixL[size] = *(tempref - rIndex);
	rIndex--;
	directionIndex--;
      } else {
	matrixL[size] = 'M';
	rIndex--;
	directionIndex--;
      }
    } else {
      if (scoreB[directionIndex][rIndex]
	  - scoreB[directionIndex - 1][rIndex] == 1
	  && directionIndex != 0) {
	matrixL[size] = 'I';
	size++;
	matrixL[size] = *(lSeq + lSeqLength - directionIndex);
	directionIndex--;
      } else if (scoreB[directionIndex][rIndex]
		 - scoreB[directionIndex][rIndex - 1] == 1 && rIndex != 0) {
	matrixL[size] = 'D';
	size++;
	matrixL[size] = *(tempref - rIndex);
	rIndex--;
      } else if (scoreB[directionIndex][rIndex]
		 - scoreB[directionIndex - 1][rIndex - 1] == 1) {
	matrixL[size] = *(tempref - rIndex);
	rIndex--;
	directionIndex--;
      } else {
	matrixL[size] = 'M';
	rIndex--;
	directionIndex--;
      }
    }
    size++;
  }
  matrixL[size] = '\0';

  char middle[SEQ_MAX_LENGTH];
  middle[0] = '\0';
  for (i = 0; i < segLength; i++)
    middle[i] = 'M';
  middle[segLength] = '\0';

  char rmatrixR[SEQ_MAX_LENGTH];

  reverse(matrixR, rmatrixR, strlen(matrixR));

  sprintf(matrix, "%s%s%s", matrixL, middle, rmatrixR);

  return totalError;

}



int addCigarSize(int cnt) {
  if (cnt < 10)
    return 1;
  else if (cnt < 100)
    return 2;
  return 3;
}

/*
  Generate Cigar from the back tracking matrix
*/
void generateCigar(char *matrix, int matrixLength, char *cigar) {
  int i = 0;

  int counterM = 0;
  int counterI = 0;
  int counterD = 0;

  int cigarSize = 0;

  cigar[0] = '\0';

  while (i < matrixLength) {
    if (matrix[i] == 'M') {
      counterM++;
      if (counterI != 0) {
	sprintf(cigar, "%s%dI", cigar, counterI);
	cigarSize += addCigarSize(counterI) + 1;
	cigar[cigarSize] = '\0';
	counterI = 0;
      } else if (counterD != 0) {
	sprintf(cigar, "%s%dD", cigar, counterD);
	cigarSize += addCigarSize(counterD) + 1;
	cigar[cigarSize] = '\0';
	counterD = 0;
      }
    } else if (matrix[i] == 'I') {
      if (counterM != 0) {
	sprintf(cigar, "%s%dM", cigar, counterM);
	cigarSize += addCigarSize(counterM) + 1;
	cigar[cigarSize] = '\0';
	counterM = 0;
      } else if (counterD != 0) {
	sprintf(cigar, "%s%dD", cigar, counterD);
	cigarSize += addCigarSize(counterD) + 1;
	cigar[cigarSize] = '\0';
	counterD = 0;
      }
      counterI++;
      i++;

    } else if (matrix[i] == 'D') {
      if (counterM != 0) {
	sprintf(cigar, "%s%dM", cigar, counterM);
	cigarSize += addCigarSize(counterM) + 1;
	cigar[cigarSize] = '\0';
	counterM = 0;
      } else if (counterI != 0) {
	sprintf(cigar, "%s%dI", cigar, counterI);
	cigarSize += addCigarSize(counterI) + 1;
	cigar[cigarSize] = '\0';
	counterI = 0;
      }

      counterD++;
      i++;

    } else {
      counterM++;
      if (counterI != 0) {
	sprintf(cigar, "%s%dI", cigar, counterI);
	cigarSize += addCigarSize(counterI) + 1;
	cigar[cigarSize] = '\0';
	counterI = 0;
      } else if (counterD != 0) {
	sprintf(cigar, "%s%dD", cigar, counterD);
	cigarSize += addCigarSize(counterD) + 1;
	cigar[cigarSize] = '\0';
	counterD = 0;
      }
    }
    i++;
  }

  if (counterM != 0) {
    sprintf(cigar, "%s%dM", cigar, counterM);
    cigarSize += addCigarSize(counterM) + 1;
    cigar[cigarSize] = '\0';
    counterM = 0;
  } else if (counterI != 0) {
    sprintf(cigar, "%s%dI", cigar, counterI);
    cigarSize += addCigarSize(counterI) + 1;
    cigar[cigarSize] = '\0';
    counterI = 0;
  } else if (counterD != 0) {
    sprintf(cigar, "%s%dD", cigar, counterD);
    cigarSize += addCigarSize(counterD) + 1;
    cigar[cigarSize] = '\0';
    counterD = 0;
  }

  cigar[cigarSize] = '\0';
}

/*
  Creates the Cigar output from the mismatching positions format  [0-9]+(([ACTGN]|\^[ACTGN]+)[0-9]+)*
*/
void generateCigarFromMD(char *mismatch, int mismatchLength, char *cigar) {
  int i = 0;
  int j = 0;

  int start = 0;
  int cigarSize = 0;

  cigar[0] = '\0';

  while (i < mismatchLength) {
    if (mismatch[i] >= '0' && mismatch[i] <= '9') {
      start = i;

      while (mismatch[i] >= '0' && mismatch[i] <= '9'
	     && i < mismatchLength)
	i++;

      int value = atoi(mismatch + start);
      for (j = 0; j < value - 1; j++) {
	cigar[cigarSize] = 'M';
	cigarSize++;
      }
      cigar[cigarSize] = 'M';
    } else if (mismatch[i] == '^') {
      cigar[cigarSize] = 'I';
      i++;
    } else if (mismatch[i] == '\'') {
      cigar[cigarSize] = 'D';
      i++;
    } else {
      cigar[cigarSize] = 'M';
      cigarSize++;
    }
    cigarSize++;
    i++;
  }
  cigar[cigarSize] = '\0';
}

void generateSNPSAM(char *matrix, int matrixLength, char *outputSNP) {

  int i = 0;

  int counterM = 0;
  int counterD = 0;

  char delete[100];

  int snpSize = 0;

  outputSNP[0] = '\0';
  delete[0] = '\0';

  while (i < matrixLength) {
    if (matrix[i] == 'M') {
      counterM++;
      if (counterD != 0) {
	delete[counterD] = '\0';
	counterD = 0;
	sprintf(outputSNP, "%s^%s", outputSNP, delete);
	snpSize += strlen(delete) + 1;
	outputSNP[snpSize] = '\0';
	delete[0] = '\0';
      }
    } else if (matrix[i] == 'D') {
      if (counterM != 0) {
	sprintf(outputSNP, "%s%d", outputSNP, counterM);
	snpSize += addCigarSize(counterM);
	outputSNP[snpSize] = '\0';
	counterM = 0;
	delete[counterD] = matrix[i + 1];
	i++;
	counterD++;
      } else if (counterD != 0) {
	delete[counterD] = matrix[i + 1];
	counterD++;
	i++;
      } else {
	delete[counterD] = matrix[i + 1];
	counterD++;
	i++;
      }
    } else if (matrix[i] == 'I') {
      if (counterM != 0) {
	// sprintf(outputSNP, "%s%d\0", outputSNP, counterM);
	//counterM++;
      } else if (counterD != 0) {
	delete[counterD] = '\0';
	sprintf(outputSNP, "%s^%s", outputSNP, delete);
	snpSize += strlen(delete) + 1;
	outputSNP[snpSize] = '\0';
	counterD = 0;
	delete[0] = '\0';
      }
      i++;

    } else {
      if (counterM != 0) {
	sprintf(outputSNP, "%s%d", outputSNP, counterM);
	snpSize += addCigarSize(counterM);
	outputSNP[snpSize] = '\0';
	counterM = 0;
      }
      if (counterD != 0) {
	delete[counterD] = '\0';
	counterD = 0;
	sprintf(outputSNP, "%s^%s", outputSNP, delete);
	snpSize += strlen(delete) + 1;
	outputSNP[snpSize] = '\0';
	delete[0] = '\0';
      }
      sprintf(outputSNP, "%s%c", outputSNP, matrix[i]);
      snpSize += 1;
      outputSNP[snpSize] = '\0';
    }
    i++;
  }

  if (counterM != 0) {
    sprintf(outputSNP, "%s%d", outputSNP, counterM);
    snpSize += addCigarSize(counterM);
    outputSNP[snpSize] = '\0';
    counterM = 0;
  } else if (counterD != 0) {
    delete[counterD] = '\0';
    sprintf(outputSNP, "%s^%s", outputSNP, delete);
    snpSize += strlen(delete) + 1;
    outputSNP[snpSize] = '\0';
    counterD = 0;
  }

  outputSNP[snpSize] = '\0';
}

/************************************************/
/* MrFAST with fastHASH: searchKey()			*/
/************************************************/
int searchKey(int target_coor, unsigned int* entry_coor, int entry_size) {
  if (entry_size <= 0)
    return 0;
  int lower_bound = 1;
  int upper_bound = entry_size;
  int mid = lower_bound + entry_size / 2;

  while (lower_bound < upper_bound) {
    if (entry_coor[mid] <= target_coor + errThreshold
	&& entry_coor[mid] >= target_coor - errThreshold)
      break;
    else if (entry_coor[mid] < target_coor)
      lower_bound = mid + 1;
    else
      upper_bound = mid - 1;
    mid = lower_bound + (upper_bound - lower_bound) / 2;
  }

  if (entry_coor[mid] <= target_coor + errThreshold
      && entry_coor[mid] >= target_coor - errThreshold) {
    return 1;
  } else
    return 0;
}

/************************************************/
/* direction = 0 forward						*/
/*  		   1 backward						*/
/************************************************/
/************************************************/
/* MrFAST with fastHASH: compareEntrySize()		*/
/************************************************/
void mapSingleEndSeq(unsigned int *l1, int s1, int readNumber, int readSegment,
		     int direction, int index, key_struct* keys_input,
		     int potential_key_number) {
  int j = 0;
  int z = 0;
  int *locs = (int *) l1;
  char *_tmpSeq, *_tmpQual;
  char rqual[SEQ_LENGTH + 1];

  int genLoc = 0;
  int leftSeqLength = 0;
  int rightSeqLength = 0;
  int middleSeqLength = 0;

  char matrix[SEQ_MAX_LENGTH];
  char editString[2 * SEQ_MAX_LENGTH];
  char cigar[MAX_CIGAR_SIZE];

  int key_number = SEQ_LENGTH / WINDOW_SIZE;
  int realLoc, readId;

  rqual[SEQ_LENGTH] = '\0';

  if (direction) {
    reverse(_msf_seqList[readNumber].qual, rqual, SEQ_LENGTH);
    _tmpQual = rqual;
    _tmpSeq = _msf_seqList[readNumber].rseq;
  } else {
    _tmpQual = _msf_seqList[readNumber].qual;
    _tmpSeq = _msf_seqList[readNumber].seq;
  }

   readId = 2 * readNumber + direction;

  for (z = 0; z < s1; z++) {
    int map_location = 0;
    int a = 0;
    int o = index;

    genLoc = locs[z];

    //hxin: If the read is at the beginning of the contig, and there are insertions
    //hxin: then genLoc - _msf_sampleingLoc[0] might be smaller than _refGenBeg
    realLoc = genLoc - _msf_samplingLocs[o];
    
    if (genLoc < _msf_samplingLocs[o]
	|| genLoc - _msf_samplingLocs[o] < _msf_refGenBeg
	|| genLoc - _msf_samplingLocs[o] > _msf_refGenEnd) {
      if (genLoc - _msf_samplingLocs[o] > _msf_refGenBeg - errThreshold)
	realLoc =  _msf_refGenBeg;
      else
	continue;
    }

    
    if (_msf_verifiedLocs[realLoc] == readId)
      continue;
    
    //Begin of long-K
    int mergeIdx = 2 * errThreshold + 1 - index;
    
    if (mergeIdx < potential_key_number) {
      if (!searchKey(
		     genLoc
		     + (keys_input[mergeIdx].key_number
			- keys_input[o].key_number) * WINDOW_SIZE,
		     keys_input[mergeIdx].key_entry,
		     keys_input[mergeIdx].key_entry_size)) {
	continue;
      }
    }
    //End of long-K
    

    // Adjacency Filtering Start ---------------------------------
    int skip_edit_distance = 0;
    int diff_num = 0;
    int ix = 0;
    for (ix = 0; ix < potential_key_number; ix++) {
      if (ix >= key_number - errThreshold) {
	break;
      }
      if (ix != o && ix != mergeIdx) { // Changed with long-K
	if (!searchKey(
		       genLoc
		       + (keys_input[ix].key_number
			  - keys_input[o].key_number)
		       * WINDOW_SIZE, keys_input[ix].key_entry,
		       keys_input[ix].key_entry_size)) {
	  diff_num++;
	  if (diff_num > errThreshold) {
	    skip_edit_distance = 1;
	    break;
	  }
	}
      }
    }
    // Adjacency Filtering End -----------------------------------

    int err = -1;
    map_location = 0;

        
    leftSeqLength = _msf_samplingLocs[o];
    middleSeqLength = WINDOW_SIZE;
    a = leftSeqLength + middleSeqLength;
    rightSeqLength = SEQ_LENGTH - a;
    

    /* CALKAN: skip alignment if it is a perfect match 
       this has to be re-addressed by Donghyuk

       if (diff_num == 0) {
       skip_edit_distance = 1;
       err = 0;
       sprintf(cigar, "%dM", SEQ_LENGTH);
       sprintf(editString, "%d", SEQ_LENGTH);
       }*/

    if (skip_edit_distance == 0) {
		mappingCnt_BeforeAlignment++;
      	err = verifySingleEndEditDistanceExtension(genLoc, _tmpSeq,
						   leftSeqLength, _tmpSeq + a, rightSeqLength,
						   middleSeqLength, matrix, &map_location);
    } 
    

    for (j = -errThreshold+1; j < errThreshold; j++) {
      if(genLoc-(readSegment*WINDOW_SIZE)+j >= _msf_refGenBeg &&
	 genLoc-(readSegment*WINDOW_SIZE)+j <= _msf_refGenEnd){
	_msf_verifiedLocs[genLoc-(readSegment*WINDOW_SIZE)+j] = readId;
      }
    }
      

    if (err != -1) {
      generateSNPSAM(matrix, strlen(matrix), editString);
      generateCigar(matrix, strlen(matrix), cigar);


      if (!bestMode) {
	mappingCnt++;
	
	_msf_seqList[readNumber].hits[0]++;
	_msf_output.QNAME = _msf_seqList[readNumber].name;
	_msf_output.FLAG = 16 * direction;
	_msf_output.RNAME = _msf_refGenName;
	_msf_output.POS = map_location + _msf_refGenOffset;
	_msf_output.MAPQ = 255;
	_msf_output.CIGAR = cigar;
	_msf_output.MRNAME = "*";
	_msf_output.MPOS = 0;
	_msf_output.ISIZE = 0;
	_msf_output.SEQ = _tmpSeq;
	_msf_output.QUAL = _tmpQual;
	  
	_msf_output.optSize = 2;
	_msf_output.optFields = _msf_optionalFields;
	  
	_msf_optionalFields[0].tag = "NM";
	_msf_optionalFields[0].type = 'i';
	_msf_optionalFields[0].iVal = err;
	  
	_msf_optionalFields[1].tag = "MD";
	_msf_optionalFields[1].type = 'Z';
	_msf_optionalFields[1].sVal = editString;
	  
	//output(_msf_output);
	  
	if (_msf_seqList[readNumber].hits[0] == 1) {
	  mappedSeqCnt++;
	}
	
	if (maxHits == 0) {
	  _msf_seqList[readNumber].hits[0] = 2;
	}
	
	if (maxHits != 0 && _msf_seqList[readNumber].hits[0] == maxHits) {
	  completedSeqCnt++;
	  break;
	}
      } 
      
      else  {  /* if mapped (err!=-1) and if it is best mode */
	mappingCnt++;
	_msf_seqList[readNumber].hits[0]++;
	
	if (_msf_seqList[readNumber].hits[0] == 1) {
	  mappedSeqCnt++;
	}
	
	if (maxHits == 0) {
	  _msf_seqList[readNumber].hits[0] = 2;
	}
	
	if (seqFastq)
	  bestHitMappingInfo[readNumber].tprob += mapProb(readNumber, editString, direction, err);
	
	if(err  < bestHitMappingInfo[readNumber].err || bestHitMappingInfo[readNumber].loc == -1)
	  {
	    setFullMappingInfo(readNumber, map_location + _msf_refGenOffset, direction, err, 0, editString, _msf_refGenName, cigar );
	  }
      }
    } 
  }
}


/************************************************/
/* MrFAST with fastHASH: compareEntrySize()		*/
/************************************************/
int compareEntrySize(const void *a, const void *b) {
  return ((*(key_struct*) a).key_entry_size
	  - (*(key_struct*) b).key_entry_size);
}

/************************************************/
/* MrFAST with fastHASH: mapAllSingleEndSeq()	*/
/************************************************/
int mapAllSingleEndSeq() {
  int i = 0;
  int j = 0;
  int k = 0;
  int it = 0;
  unsigned int *locs = NULL;
  int key_number = SEQ_LENGTH / WINDOW_SIZE;
  key_struct* sort_input = getMem(key_number * sizeof(key_struct));

  // Forward Mode
  for (i = 0; i < _msf_seqListSize; i++) {
    k = _msf_sort_seqList[i].readNumber;
    int available_key_num = 0;
    for (it = 0; it < key_number; it++) {
      locs = getCandidates(hashVal(_msf_seqList[k].seq + it * WINDOW_SIZE));

      if (locs != NULL) {
	sort_input[available_key_num].key_number = it;
	sort_input[available_key_num].key_entry = locs;
	sort_input[available_key_num].key_entry_size = locs[0];
	available_key_num++;
      }
    }

    int operating_key_num = _msf_samplingLocsSize;
    if (available_key_num < operating_key_num) {
      operating_key_num = available_key_num;
    }

    qsort(sort_input, available_key_num, sizeof(key_struct),
	  compareEntrySize);

    for (j = 0; j < operating_key_num; j++) {
      _msf_samplingLocs[j] = sort_input[j].key_number * WINDOW_SIZE;
      mapSingleEndSeq(sort_input[j].key_entry + 1,
		      sort_input[j].key_entry_size, k, sort_input[j].key_number,
		      0, j, sort_input, available_key_num);
    }
  }

  // Reverse Mode
  for (i = 0; i < _msf_seqListSize; i++) {
    k = _msf_sort_seqList[i].readNumber;
    int available_key_num = 0;
    for (it = 0; it < key_number; it++) {
      
      locs = getCandidates(hashVal(_msf_seqList[k].rseq + it * WINDOW_SIZE));

      if (locs != NULL) {
	sort_input[available_key_num].key_number = it;
	sort_input[available_key_num].key_entry = locs;
	sort_input[available_key_num].key_entry_size = locs[0];
	available_key_num++;
      }
    }

    qsort(sort_input, available_key_num, sizeof(key_struct),
	  compareEntrySize);

    int operating_key_num = _msf_samplingLocsSize;
    if (available_key_num < operating_key_num) {
      operating_key_num = available_key_num;
    }

    for (j = 0; j < operating_key_num; j++) {
      _msf_samplingLocs[j] = sort_input[j].key_number * WINDOW_SIZE;
      mapSingleEndSeq(sort_input[j].key_entry + 1,
		      sort_input[j].key_entry_size, k, sort_input[j].key_number,
		      1, j, sort_input, available_key_num);
    }
  }
  freeMem(sort_input, key_number * sizeof(key_struct));

  return 1;
}


/************************************************/
int compareOut(const void *a, const void *b) {
  FullMappingInfo *aInfo = (FullMappingInfo *) a;
  FullMappingInfo *bInfo = (FullMappingInfo *) b;
  return aInfo->loc - bInfo->loc;
}

/************************************************/

/************************************************/
/* direction = 0 forward						*/
/*  		   1 backward						*/
/************************************************/
/************************************************/
/* MrFAST with fastHASH: mapPairEndSeqList()	*/
/************************************************/
void mapPairEndSeqList(unsigned int *l1, int s1, int readNumber,
		       int readSegment, int direction, int index, key_struct* keys_input,
		       int potential_key_number) {
  int z = 0;
  int *locs = (int *) l1;
  char *_tmpSeq;

  char matrix[SEQ_MAX_LENGTH];
  char editString[2 * SEQ_MAX_LENGTH];
  char cigar[MAX_CIGAR_SIZE];


  int leftSeqLength = 0;
  int rightSeqLength = 0;
  int middleSeqLength = 0;

  int genLoc = 0;
  int r = readNumber;

  char d = (direction == 1) ? -1 : 1;

  int readId = 2 * readNumber + direction;
  int key_number = SEQ_LENGTH / WINDOW_SIZE;

  if (d == -1) {
    _tmpSeq = _msf_seqList[readNumber].rseq;
  } 
  else {
    _tmpSeq = _msf_seqList[readNumber].seq;
  }

  for (z = 0; z < s1; z++) {
    int map_location = 0;
    int a = 0;
    int o = index;
    genLoc = locs[z];

    leftSeqLength = _msf_samplingLocs[o];
    middleSeqLength = WINDOW_SIZE;
    a = leftSeqLength + middleSeqLength;
    rightSeqLength = SEQ_LENGTH - a;


    /*    
    //hxinPE: If the read is at the beginning of the contig, and there are insertions
    //hxinPE: then genLoc - _msf_sampleingLoc[0] might be smaller than _refGenBeg
    int realLoc = genLoc - _msf_samplingLocs[o];
    if (genLoc < _msf_samplingLocs[o]
	|| genLoc - _msf_samplingLocs[o] < _msf_refGenBeg
	|| genLoc - _msf_samplingLocs[o] > _msf_refGenEnd) {
      if (genLoc - _msf_samplingLocs[o] > _msf_refGenBeg - errThreshold)
	realLoc =  _msf_refGenBeg;
      else
	continue;
    }

    if (_msf_verifiedLocs[realLoc] == readId || _msf_verifiedLocs[realLoc] == -readId )
      continue;
    */

    
    if (genLoc - leftSeqLength < _msf_refGenBeg
	|| genLoc + rightSeqLength + middleSeqLength > _msf_refGenEnd
	|| _msf_verifiedLocs[genLoc - _msf_samplingLocs[o]] == readId
	|| _msf_verifiedLocs[genLoc - _msf_samplingLocs[o]] == -readId)
      continue;
   
    //Begin of long-K
    int mergeIdx = 2 * errThreshold + 1 - index;
    
    if (mergeIdx < potential_key_number) {
      if (!searchKey(
		     genLoc
		     + (keys_input[mergeIdx].key_number
			- keys_input[o].key_number) * WINDOW_SIZE,
		     keys_input[mergeIdx].key_entry,
		     keys_input[mergeIdx].key_entry_size)) {
	continue;
      }
    }
    //End of long-K



    // Adjacency Filtering Start ---------------------------------
    int skip_edit_distance = 0;
    int diff_num = 0;
    int ix = 0;
    for (ix = 0; ix < potential_key_number; ix++) {
      if (ix >= key_number - errThreshold) {
	break;
      }
      if (ix != o && ix != mergeIdx) { // Changed with long-K
	if (!searchKey(
		       genLoc
		       + (keys_input[ix].key_number
			  - keys_input[o].key_number)
		       * WINDOW_SIZE, keys_input[ix].key_entry,
		       keys_input[ix].key_entry_size)) {
	  diff_num++;
	  if (diff_num > errThreshold) {
	    skip_edit_distance = 1;
	    break;
	  }
	}
      }
    }
    // Adjacency Filtering End -----------------------------------

    int err = -1;
    map_location = 0;

    if (skip_edit_distance == 0) {
     
      err = verifySingleEndEditDistanceExtension(genLoc, _tmpSeq,
						 leftSeqLength, _tmpSeq + a, rightSeqLength,
						 middleSeqLength, matrix, &map_location);
	
    } else {
      err = -1; 
    }

    int j = 0;

    for (j = -errThreshold+1; j < errThreshold; j++) {
      if(genLoc-(readSegment*WINDOW_SIZE)+j >= _msf_refGenBeg &&
	 genLoc-(readSegment*WINDOW_SIZE)+j <= _msf_refGenEnd)
	_msf_verifiedLocs[genLoc-(readSegment*WINDOW_SIZE)+j] = readId;
    }


    if (err != -1) {
      int i = 0;
      /* calkan counter */
      mappingCnt++;
   
      generateSNPSAM(matrix, strlen(matrix), editString);
      generateCigar(matrix, strlen(matrix), cigar);
      MappingLocations *parent = NULL;
      MappingLocations *child = _msf_mappingInfo[r].next;
      genLoc = map_location + _msf_refGenOffset;

      for (i = 0; i < (_msf_mappingInfo[r].size / MAP_CHUNKS); i++) {
	parent = child;
	child = child->next;
      }

      if (child == NULL) {
	MappingLocations *tmp = getMem(sizeof(MappingLocations));
	tmp->next = NULL;
	tmp->loc[0] = genLoc * d; // d is required: DHL
	tmp->err[0] = err;
	tmp->cigarSize[0] = strlen(cigar);
	sprintf(tmp->cigar[0], "%s", cigar);
	tmp->mdSize[0] = strlen(editString);
	sprintf(tmp->md[0], "%s", editString);

	if (parent == NULL)
	  _msf_mappingInfo[r].next = tmp;
	else
	  parent->next = tmp;
      } else {
	if (strlen(cigar) > SEQ_LENGTH
	    || strlen(editString) > SEQ_LENGTH) {
	  fprintf(stderr, 
		 "ERROR in %d read size(After mapping) exceeds cigar=%d md =%d cigar=%s md =%s\n",
		 r, (int) strlen(cigar), (int) strlen(editString),
		 cigar, editString);
	}
	child->loc[_msf_mappingInfo[r].size % MAP_CHUNKS] = genLoc * d;
	child->err[_msf_mappingInfo[r].size % MAP_CHUNKS] = err;
	child->cigarSize[_msf_mappingInfo[r].size % MAP_CHUNKS] = strlen(cigar);
	sprintf(child->cigar[_msf_mappingInfo[r].size % MAP_CHUNKS], "%s", cigar);
	child->mdSize[_msf_mappingInfo[r].size % MAP_CHUNKS] = strlen(editString);
	sprintf(child->md[_msf_mappingInfo[r].size % MAP_CHUNKS], "%s", editString);
      }
      _msf_mappingInfo[r].size++;
    } 
  }
}


/************************************************/
/* MrFAST with fastHASH: mapPairedEndSeq()		*/
/************************************************/
void mapPairedEndSeq() {
  // DHL: Changed Start
  int i = 0;
  int j = 0;
  int k = 0;
  unsigned int *locs = NULL;
  int key_number = SEQ_LENGTH / WINDOW_SIZE;
  key_struct* sort_input = getMem(key_number * sizeof(key_struct));

  // Forward Mode
  for (i = 0; i < _msf_seqListSize; i++) {
    k = _msf_sort_seqList[i].readNumber;
    int available_key_num = 0;
    int it = 0;
    for (it = 0; it < key_number; it++) {
      //int key_hash = 
      locs = getCandidates(hashVal(_msf_seqList[k].seq + it * WINDOW_SIZE));
      if (locs != NULL) {
	sort_input[available_key_num].key_number = it;
	sort_input[available_key_num].key_entry = locs;
	sort_input[available_key_num].key_entry_size = locs[0];
	available_key_num++;
      }
    }

    int operating_key_num = _msf_samplingLocsSize;
    if (available_key_num < operating_key_num) {
      operating_key_num = available_key_num;
    }

    qsort(sort_input, available_key_num, sizeof(key_struct),
	  compareEntrySize);

    for (j = 0; j < operating_key_num; j++) {
      _msf_samplingLocs[j] = sort_input[j].key_number * WINDOW_SIZE;
      mapPairEndSeqList(sort_input[j].key_entry + 1,
			sort_input[j].key_entry_size, k, sort_input[j].key_number,
			0, j, sort_input, available_key_num);
    }
  }

  // Reverse Mode
  for (i = 0; i < _msf_seqListSize; i++) {
    k = _msf_sort_seqList[i].readNumber;
    int key_number = SEQ_LENGTH / WINDOW_SIZE;
    int available_key_num = 0;
    int it = 0;
    for (it = 0; it < key_number; it++) {
      //int key_hash = 
      locs = getCandidates(hashVal(_msf_seqList[k].rseq + it * WINDOW_SIZE));

      if (locs != NULL) {
	sort_input[available_key_num].key_number = it;
	sort_input[available_key_num].key_entry = locs;
	sort_input[available_key_num].key_entry_size = locs[0];
	available_key_num++;
      }
    }

    int operating_key_num = _msf_samplingLocsSize;
    if (available_key_num < operating_key_num) {
      operating_key_num = available_key_num;
    }

    qsort(sort_input, available_key_num, sizeof(key_struct),
	  compareEntrySize);

    for (j = 0; j < operating_key_num; j++) {
      _msf_samplingLocs[j] = sort_input[j].key_number * WINDOW_SIZE;
      mapPairEndSeqList(sort_input[j].key_entry + 1,
			sort_input[j].key_entry_size, k, sort_input[j].key_number,
			1, j, sort_input, available_key_num);
    }
  }
  freeMem(sort_input, key_number * sizeof(key_struct));
  // DHL: Changed End

  char fname1[FILE_NAME_LENGTH];
  char fname2[FILE_NAME_LENGTH];
  MappingLocations *cur;
  int tmpOut;
  int lmax = 0, rmax = 0;

  sprintf(fname1, "%s__%s__%s__%d__1.tmp", mappingOutputPath, _msf_refGenName,
	  mappingOutput, _msf_openFiles);
  sprintf(fname2, "%s__%s__%s__%d__2.tmp", mappingOutputPath, _msf_refGenName,
	  mappingOutput, _msf_openFiles);

  FILE* out;
  FILE* out1 = fileOpen(fname1, "w");
  FILE* out2 = fileOpen(fname2, "w");

  _msf_openFiles++;

  for (i = 0; i < _msf_seqListSize; i++) {
    if (i % 2 == 0) {
      out = out1;
      if (lmax < _msf_mappingInfo[i].size) {
	lmax = _msf_mappingInfo[i].size;
      }
    } else {
      out = out2;
      if (rmax < _msf_mappingInfo[i].size) {
	rmax = _msf_mappingInfo[i].size;
      }
    }
    tmpOut = fwrite(&(_msf_mappingInfo[i].size), sizeof(int), 1, out);
    if (_msf_mappingInfo[i].size > 0) {
      cur = _msf_mappingInfo[i].next;
      for (j = 0; j < _msf_mappingInfo[i].size; j++) {
	if (j > 0 && j % MAP_CHUNKS == 0) {
	  cur = cur->next;
	}
	if(debugMode && (cur->cigarSize[j % MAP_CHUNKS] > SEQ_LENGTH || cur->mdSize[j % MAP_CHUNKS] > SEQ_LENGTH))
	  {
	    fprintf(stderr, "ERROR in %d read size exceeds cigar=%d md =%d cigar=%s md =%s\n", i,  cur->cigarSize[j % MAP_CHUNKS], cur->mdSize[j % MAP_CHUNKS], cur->cigar[j % MAP_CHUNKS], cur->md[j % MAP_CHUNKS]);	
	  }

	tmpOut = fwrite(&(cur->loc[j % MAP_CHUNKS]), sizeof(int), 1,
			out);
	tmpOut = fwrite(&(cur->err[j % MAP_CHUNKS]), sizeof(int), 1,
			out);
	tmpOut = fwrite(&(cur->cigarSize[j % MAP_CHUNKS]), sizeof(int),
			1, out);
	tmpOut = fwrite((cur->cigar[j % MAP_CHUNKS]), sizeof(char),
			(cur->cigarSize[j % MAP_CHUNKS]), out);
	tmpOut = fwrite(&(cur->mdSize[j % MAP_CHUNKS]), sizeof(int), 1,
			out);
	tmpOut = fwrite((cur->md[j % MAP_CHUNKS]), sizeof(char),
			(cur->mdSize[j % MAP_CHUNKS]), out);
      }
      _msf_mappingInfo[i].size = 0;
    }
  }
  _msf_maxLSize += lmax;
  _msf_maxRSize += rmax;
  tmpOut++;


  fclose(out1);
  fclose(out2);
}

void outputPairFullMappingInfo(FILE *fp, int readNumber) {

  char *seq1, *seq2, *rseq1, *rseq2, *qual1, *qual2;
  char rqual1[SEQ_LENGTH + 1], rqual2[SEQ_LENGTH + 1];

  rqual1[SEQ_LENGTH] = rqual2[SEQ_LENGTH] = '\0';

  seq1 = _msf_seqList[readNumber * 2].seq;
  rseq1 = _msf_seqList[readNumber * 2].rseq;
  qual1 = _msf_seqList[readNumber * 2].qual;

  reverse(_msf_seqList[readNumber * 2].qual, rqual1, SEQ_LENGTH);

  seq2 = _msf_seqList[readNumber * 2 + 1].seq;
  rseq2 = _msf_seqList[readNumber * 2 + 1].rseq;
  qual2 = _msf_seqList[readNumber * 2 + 1].qual;

  reverse(_msf_seqList[readNumber * 2 + 1].qual, rqual2, SEQ_LENGTH);

  if (bestHitMappingInfo[readNumber * 2].loc == -1
      && bestHitMappingInfo[readNumber * 2 + 1].loc == -1)
    return;
  else {

    char *seq;
    char *qual;
    char d1;
    char d2;
    int isize;
    int proper = 0;
    // ISIZE CALCULATION
    // The distance between outer edges
    isize = abs(
		bestHitMappingInfo[readNumber * 2].loc
		- bestHitMappingInfo[readNumber * 2 + 1].loc)
      + SEQ_LENGTH - 2;

    if (bestHitMappingInfo[readNumber * 2].loc
	- bestHitMappingInfo[readNumber * 2 + 1].loc > 0) {
      isize *= -1;
    }
    d1 = (bestHitMappingInfo[readNumber * 2].dir == -1) ? 1 : 0;
    d2 = (bestHitMappingInfo[readNumber * 2 + 1].dir == -1) ? 1 : 0;

    if (d1) {
      seq = rseq1;
      qual = rqual1;
    } else {
      seq = seq1;
      qual = qual1;
    }
    if ((bestHitMappingInfo[readNumber * 2].loc
	 < bestHitMappingInfo[readNumber * 2 + 1].loc && !d1 && d2)
	|| (bestHitMappingInfo[readNumber * 2].loc
	    > bestHitMappingInfo[readNumber * 2 + 1].loc && d1
	    && !d2)) {
      proper = 2;
    } else {
      proper = 0;
    }

    _msf_output.POS = bestHitMappingInfo[readNumber * 2].loc;
    _msf_output.MPOS = bestHitMappingInfo[readNumber * 2 + 1].loc;
    _msf_output.FLAG = 1 + proper + 16 * d1 + 32 * d2 + 64;
    _msf_output.ISIZE = isize;
    _msf_output.SEQ = seq;
    _msf_output.QUAL = qual;
    _msf_output.QNAME = _msf_seqList[readNumber * 2].name;
    _msf_output.RNAME = bestHitMappingInfo[readNumber * 2].chr;
    if (seqFastq)
      _msf_output.MAPQ = mapQ(readNumber * 2) + mapQ(readNumber * 2 + 1);
    else
      _msf_output.MAPQ = 255;
    _msf_output.CIGAR = bestHitMappingInfo[readNumber * 2].cigar;
    _msf_output.MRNAME = "=";

    _msf_output.optSize = 2;
    _msf_output.optFields = _msf_optionalFields;

    _msf_optionalFields[0].tag = "NM";
    _msf_optionalFields[0].type = 'i';
    _msf_optionalFields[0].iVal = bestHitMappingInfo[readNumber * 2].err;

    _msf_optionalFields[1].tag = "MD";
    _msf_optionalFields[1].type = 'Z';
    _msf_optionalFields[1].sVal = bestHitMappingInfo[readNumber * 2].md;

    
    output(_msf_output);

    if (d2) {
      seq = rseq2;
      qual = rqual2;
    } else {
      seq = seq2;
      qual = qual2;
    }

    _msf_output.POS = bestHitMappingInfo[readNumber * 2 + 1].loc;
    _msf_output.MPOS = bestHitMappingInfo[readNumber * 2].loc;
    _msf_output.FLAG = 1 + proper + 16 * d2 + 32 * d1 + 128;
    _msf_output.ISIZE = -isize;
    _msf_output.SEQ = seq;
    _msf_output.QUAL = qual;
    _msf_output.QNAME = _msf_seqList[readNumber * 2].name;
    _msf_output.RNAME = bestHitMappingInfo[readNumber * 2].chr;
    if (seqFastq)
      _msf_output.MAPQ = mapQ(readNumber * 2) + mapQ(readNumber * 2 + 1);
    else
      _msf_output.MAPQ = 255;
    _msf_output.CIGAR = bestHitMappingInfo[readNumber * 2 + 1].cigar;
    _msf_output.MRNAME = "=";

    _msf_output.optSize = 2;
    _msf_output.optFields = _msf_optionalFields;

    _msf_optionalFields[0].tag = "NM";
    _msf_optionalFields[0].type = 'i';
    _msf_optionalFields[0].iVal =
      bestHitMappingInfo[readNumber * 2 + 1].err;

    _msf_optionalFields[1].tag = "MD";
    _msf_optionalFields[1].type = 'Z';
    _msf_optionalFields[1].sVal = bestHitMappingInfo[readNumber * 2 + 1].md;

    output(_msf_output);
  }
}

/*
  Find the closet one to the c
  @return 0: if the x1 is closer to c
  1: if the x2 is closer to c
  2: if both distance are equal
  -1: if error
*/
int findNearest(int x1, int x2, int c) {

  if (abs(x1 - c) > abs(x2 - c))
    return 0;
  else if (abs(x1 - c) < abs(x2 - c))
    return 1;
  else if (abs(x1 - c) == abs(x2 - c))
    return 2;
  else
    return -1;
}



void finalizeBestConcordantDiscordant() {
  int i = 0;

  for (i = 0; i < _msf_seqListSize / 2; i++) {
    outputPairFullMappingInfo(NULL, i);
  }
  freeMem(bestHitMappingInfo, _msf_seqListSize * sizeof(FullMappingInfo));
}

double mapProb(int readNumber, char *md, int dir, int err){
  int i = 0;
  int mdlen = strlen(md);
  char buf[MAX_CIGAR_SIZE];
  int j = 0;

  double phred = 0.0;
  int errloc = 0;
  int errcnt = 0; //since I cannot calculate deletion base quality
 

  buf[0] = 0;

  if (err == 0) 
    return 1.0;

  while (i<mdlen){
    if (isdigit(md[i])) buf[j++]=md[i++];
    else if (isalpha(md[i])){
      /* mismatch */
      errcnt++;
      buf[j] = '\0'; 
      if (j != 0)
	errloc += atoi(buf);
      else if (i!=0)
	errloc++;

      j=0; buf[0]=0;

      if (dir)
	phred += (double) (_msf_seqList[readNumber].qual[SEQ_LENGTH-errloc-1] - 33);
      else
	phred += (double) (_msf_seqList[readNumber].qual[errloc] - 33);

      i++;
    }
    
    else if (md[i]=='^'){
      /* insertion to the read / deletion from reference  */
      if (j!=0){
	buf[j]=0;
	errloc += atoi(buf);
	buf[0] = 0;
      }
      j=0; 
      i++; /* pass ^ */
      while (isalpha(md[i++])) j++;
      errloc += j;
      j = 0;
    }
  }

  double indel_prob = 1; 
  if (errcnt != err)
    indel_prob = 0.0002 * (err - errcnt);

  return pow(10, -1 * (phred / 10)) * indel_prob;

}

int mapQ(int readNumber)
{
  int mapqual;
  double mapprob;

  mapprob = mapProb(readNumber, bestHitMappingInfo[readNumber].md, 
		    bestHitMappingInfo[readNumber].dir, bestHitMappingInfo[readNumber].err); 

  if (mapprob == bestHitMappingInfo[readNumber].tprob)
    mapqual = 40;

  else 
    mapqual =  (int) (round(-10.0 * log10(1 - (mapprob / bestHitMappingInfo[readNumber].tprob))));
  
  if (mapqual > 40) mapqual = 40;

  return mapqual;

}

void setFullMappingInfo(int readNumber, int loc, int dir, int err, int score,
			char *md, char * refName, char *cigar) {
  bestHitMappingInfo[readNumber].loc = loc;
  bestHitMappingInfo[readNumber].dir = dir;
  bestHitMappingInfo[readNumber].err = err;
  bestHitMappingInfo[readNumber].score = score;

  strncpy(bestHitMappingInfo[readNumber].md, md, strlen(md) + 1);

  /*
  if (bestHitMappingInfo[readNumber].chr == NULL)
    bestHitMappingInfo[readNumber].chr = (char *) getMem(sizeof(char) * (strlen(refName)+1));
  else if (strlen(bestHitMappingInfo[readNumber].chr) < strlen(refName)){
    freeMem(bestHitMappingInfo[readNumber].chr, (strlen(bestHitMappingInfo[readNumber].chr)+1));
    bestHitMappingInfo[readNumber].chr = (char *) getMem(sizeof(char) * (strlen(refName)+1));
  }
  */

  strncpy(bestHitMappingInfo[readNumber].chr, refName, strlen(refName) + 1);
  strncpy(bestHitMappingInfo[readNumber].cigar, cigar, strlen(cigar) + 1);
}

void setPairFullMappingInfo(int readNumber, FullMappingInfo mi1,
			    FullMappingInfo mi2) {

  bestHitMappingInfo[readNumber * 2].loc = mi1.loc;
  bestHitMappingInfo[readNumber * 2].dir = mi1.dir;
  bestHitMappingInfo[readNumber * 2].err = mi1.err;
  bestHitMappingInfo[readNumber * 2].score = mi1.score;

  /*
  if (bestHitMappingInfo[readNumber * 2].chr == NULL)
    bestHitMappingInfo[readNumber * 2].chr = (char *) getMem(sizeof(char) * (strlen(_msf_refGenName)+1));
  else if (strlen(bestHitMappingInfo[readNumber * 2].chr) < strlen(_msf_refGenName)){
    freeMem(bestHitMappingInfo[readNumber * 2].chr, (strlen(bestHitMappingInfo[readNumber * 2].chr)+1));
    bestHitMappingInfo[readNumber * 2].chr = (char *) getMem(sizeof(char) * (strlen(_msf_refGenName)+1));
    }
  */

  snprintf(bestHitMappingInfo[readNumber * 2].chr, strlen(_msf_refGenName)+1, "%s",
	   _msf_refGenName);

  strncpy(bestHitMappingInfo[readNumber * 2].md, mi1.md, strlen(mi1.md) + 1);
  strncpy(bestHitMappingInfo[readNumber * 2].cigar, mi1.cigar,
	  strlen(mi1.cigar) + 1);

  bestHitMappingInfo[readNumber * 2 + 1].loc = mi2.loc;
  bestHitMappingInfo[readNumber * 2 + 1].dir = mi2.dir;
  bestHitMappingInfo[readNumber * 2 + 1].err = mi2.err;
  bestHitMappingInfo[readNumber * 2 + 1].score = mi2.score;

  /*
  if (bestHitMappingInfo[readNumber * 2 + 1].chr == NULL)
    bestHitMappingInfo[readNumber * 2 + 1].chr = (char *) getMem(sizeof(char) * (strlen(_msf_refGenName)+1));
  else if (strlen(bestHitMappingInfo[readNumber * 2 + 1].chr) < strlen(_msf_refGenName)){
    freeMem(bestHitMappingInfo[readNumber * 2 + 1].chr, (strlen(bestHitMappingInfo[readNumber * 2 + 1].chr)+1));
    bestHitMappingInfo[readNumber * 2 + 1].chr = (char *) getMem(sizeof(char) * (strlen(_msf_refGenName)+1));
  }
  */

  snprintf(bestHitMappingInfo[readNumber * 2 + 1].chr, strlen(_msf_refGenName)+1, "%s",
	   _msf_refGenName);

  strncpy(bestHitMappingInfo[readNumber * 2 + 1].md, mi2.md,
	  strlen(mi2.md) + 1);
  strncpy(bestHitMappingInfo[readNumber * 2 + 1].cigar, mi2.cigar,
	  strlen(mi2.cigar) + 1);

}

/**********************************************/
void outputPairedEnd() {
  int i = 0;

  char cigar[MAX_CIGAR_SIZE];

  int tmpOut;

  FILE* in1[_msf_openFiles];
  FILE* in2[_msf_openFiles];

  char fname1[_msf_openFiles][FILE_NAME_LENGTH];
  char fname2[_msf_openFiles][FILE_NAME_LENGTH];

  // discordant
  FILE *out = NULL, *out1 = NULL;

  char fname3[FILE_NAME_LENGTH];
  char fname4[FILE_NAME_LENGTH];

  int meanDistanceMapping = 0;

  char rqual1[SEQ_LENGTH + 1];
  char rqual2[SEQ_LENGTH + 1];
  int tmp = 0;

  loadRefGenome(&_msf_refGen, &_msf_refGenName, &tmpOut);

  if (pairedEndDiscordantMode) {
    sprintf(fname3, "%s__%s__disc", mappingOutputPath, mappingOutput);
    sprintf(fname4, "%s__%s__oea", mappingOutputPath, mappingOutput);
    out = fileOpen(fname3, "a");
    out1 = fileOpen(fname4, "a");
  }

  FullMappingInfo *mi1 = getMem(sizeof(FullMappingInfo) * _msf_maxLSize);
  FullMappingInfo *mi2 = getMem(sizeof(FullMappingInfo) * _msf_maxRSize);

  _msf_fileCount[_msf_maxFile] = 0;
  for (i = 0; i < _msf_openFiles; i++) {
    sprintf(fname1[i], "%s__%s__%s__%d__1.tmp", mappingOutputPath,
	    _msf_refGenName, mappingOutput, i);
    sprintf(_msf_fileName[_msf_maxFile][_msf_fileCount[_msf_maxFile]][0],
	    "%s", fname1[i]);

    sprintf(fname2[i], "%s__%s__%s__%d__2.tmp", mappingOutputPath,
	    _msf_refGenName, mappingOutput, i);
    sprintf(_msf_fileName[_msf_maxFile][_msf_fileCount[_msf_maxFile]][1],
	    "%s", fname2[i]);

    in1[i] = fileOpen(fname1[i], "r");
    in2[i] = fileOpen(fname2[i], "r");
    _msf_fileCount[_msf_maxFile]++;
  }
  _msf_maxFile++;

  int size;
  int j, k;
  int size1, size2;

  meanDistanceMapping =
    (pairedEndDiscordantMode == 1) ?
    (minPairEndedDiscordantDistance
     + maxPairEndedDiscordantDistance) / 2 + SEQ_LENGTH :
    (minPairEndedDistance + maxPairEndedDistance) / 2
    + SEQ_LENGTH;

  for (i = 0; i < _msf_seqListSize / 2; i++) {
    size1 = size2 = 0;
    for (j = 0; j < _msf_openFiles; j++) {
      tmpOut = fread(&size, sizeof(int), 1, in1[j]);
      if (size > 0) {
	for (k = 0; k < size; k++) {
	  mi1[size1 + k].dir = 1;
	  tmpOut = fread(&(mi1[size1 + k].loc), sizeof(int), 1,
			 in1[j]);
	  tmpOut = fread(&(mi1[size1 + k].err), sizeof(int), 1,
			 in1[j]);

	  tmpOut = fread(&(mi1[size1 + k].cigarSize), sizeof(int), 1,
			 in1[j]);
	  tmpOut = fread((mi1[size1 + k].cigar), sizeof(char),
			 mi1[size1 + k].cigarSize, in1[j]);
	  mi1[size1 + k].cigar[mi1[size1 + k].cigarSize] = '\0';

	  tmpOut = fread(&(mi1[size1 + k].mdSize), sizeof(int), 1,
			 in1[j]);
	  tmpOut = fread((mi1[size1 + k].md), sizeof(char),
			 (mi1[size1 + k].mdSize), in1[j]);
	  mi1[size1 + k].md[mi1[size1 + k].mdSize] = '\0';

	  if (mi1[size1 + k].loc < 1) {
	    mi1[size1 + k].loc *= -1;
	    mi1[size1 + k].dir = -1;
	  }
	}
	qsort(mi1 + size1, size, sizeof(FullMappingInfo), compareOut);
	size1 += size;
      }
    }

    for (j = 0; j < _msf_openFiles; j++) {
      tmpOut = fread(&size, sizeof(int), 1, in2[j]);
      if (size > 0) {
	for (k = 0; k < size; k++) {
	  mi2[size2 + k].dir = 1;
	  tmpOut = fread(&(mi2[size2 + k].loc), sizeof(int), 1,
			 in2[j]);
	  tmpOut = fread(&(mi2[size2 + k].err), sizeof(int), 1,
			 in2[j]);

	  tmpOut = fread(&(mi2[size2 + k].cigarSize), sizeof(int), 1,
			 in2[j]);
	  tmpOut = fread((mi2[size2 + k].cigar), sizeof(char),
			 mi2[size2 + k].cigarSize, in2[j]);
	  mi2[size2 + k].cigar[mi2[size2 + k].cigarSize] = '\0';

	  tmpOut = fread(&(mi2[size2 + k].mdSize), sizeof(int), 1,
			 in2[j]);
	  tmpOut = fread((mi2[size2 + k].md), sizeof(char),
			 mi2[size2 + k].mdSize, in2[j]);
	  mi2[size2 + k].md[mi2[size2 + k].mdSize] = '\0';

	  if (mi2[size2 + k].loc < 1) {
	    mi2[size2 + k].loc *= -1;
	    mi2[size2 + k].dir = -1;
	  }
	}
	qsort(mi2 + size2, size, sizeof(FullMappingInfo), compareOut);
	size2 += size;
      }
    }

    int lm, ll, rl, rm;
    int pos = 0;

    if (pairedEndDiscordantMode) {

      for (j = 0; j < size1; j++) {
	lm = mi1[j].loc - maxPairEndedDiscordantDistance + 1;
	ll = mi1[j].loc - minPairEndedDiscordantDistance + 1;
	rl = mi1[j].loc + minPairEndedDiscordantDistance - 1;
	rm = mi1[j].loc + maxPairEndedDiscordantDistance - 1;

	while (pos < size2 && mi2[pos].loc < lm) {
	  pos++;
	}

	k = pos;
	while (k < size2 && mi2[k].loc <= rm) {
	  if (mi2[k].loc <= ll || mi2[k].loc >= rl) {
	    if ((mi1[j].loc < mi2[k].loc && mi1[j].dir == 1
		 && mi2[k].dir == -1)
		|| (mi1[j].loc > mi2[k].loc && mi1[j].dir == -1
		    && mi2[k].dir == 1)) {
	      _msf_seqList[i * 2].hits[0] = 1;
	      _msf_seqList[i * 2 + 1].hits[0] = 1;

	      if (nosamMode != 0) {
		size1 = 0;
		size2 = 0;
	      }

	      break;
	    }
	  }
	  k++;
	}
      }

      _msf_seqHits[i * 2] += size1;
      _msf_seqHits[i * 2 + 1] += size2;

      if (_msf_seqHits[i * 2 + 1] * _msf_seqHits[i * 2]
	  > DISCORDANT_CUT_OFF && nosamMode != 0) {
	_msf_seqList[i * 2].hits[0] = 1;
	_msf_seqList[i * 2 + 1].hits[0] = 1;
	size1 = 0;
	size2 = 0;
      }

      int rNo = 0;
      int loc = 0;
      int err = 0;
      float sc = 0;
      char l = 0;

      //write the OEA data
      if (_msf_seqHits[i * 2] == 0){
	for (k = 0;
	     k < size2 && _msf_oeaMapping[i * 2 + 1] < maxOEAOutput;
	     k++) {
	  rNo = i * 2 + 1;
	  loc = mi2[k].loc * mi2[k].dir;
	  err = mi2[k].err;
	  sc = mi2[k].score;

	  l = strlen(_msf_refGenName);

	  tmp = fwrite(&rNo, sizeof(int), 1, out1);

	  tmp = fwrite(&l, sizeof(char), 1, out1);
	  tmp = fwrite(_msf_refGenName, sizeof(char), l, out1);

	  tmp = fwrite(&loc, sizeof(int), 1, out1);
	  tmp = fwrite(&err, sizeof(int), 1, out1);
	  tmp = fwrite(&sc, sizeof(float), 1, out1);

	  if (mi2[k].cigarSize > SEQ_LENGTH || mi2[k].cigarSize <= 0)
	    fprintf(stderr, "ERROR  CIGAR size=%d %s\n", mi2[k].cigarSize,
		   _msf_seqList[i * 2 + 1].seq);

	  tmp = fwrite(&(mi2[k].cigarSize), sizeof(int), 1, out1);
	  tmp = fwrite((mi2[k].cigar), sizeof(char), mi2[k].cigarSize,
		       out1);

	  tmp = fwrite(&(mi2[k].mdSize), sizeof(int), 1, out1);
	  tmp = fwrite((mi2[k].md), sizeof(char), mi2[k].mdSize,
		       out1);

	  _msf_oeaMapping[i * 2 + 1]++;
	}
      }
      if (_msf_seqHits[i * 2 + 1] == 0){ 
	for (j = 0; j < size1 && _msf_oeaMapping[i * 2] < maxOEAOutput;
	     j++) {
	  rNo = i * 2;
	  loc = mi1[j].loc * mi1[j].dir;
	  err = mi1[j].err;
	  sc = mi1[j].score;

	  l = strlen(_msf_refGenName);

	  tmp = fwrite(&rNo, sizeof(int), 1, out1);

	  tmp = fwrite(&l, sizeof(char), 1, out1);
	  tmp = fwrite(_msf_refGenName, sizeof(char), l, out1);

	  tmp = fwrite(&loc, sizeof(int), 1, out1);
	  tmp = fwrite(&err, sizeof(int), 1, out1);
	  tmp = fwrite(&sc, sizeof(float), 1, out1);

	  if (mi1[j].cigarSize > SEQ_LENGTH || mi1[j].cigarSize <= 0)
	    fprintf(stderr, "ERROR %d %s\n", mi1[j].cigarSize,
		   _msf_seqList[i * 2 + 1].seq);

	  tmp = fwrite(&(mi1[j].cigarSize), sizeof(int), 1, out1);
	  tmp = fwrite((mi1[j].cigar), sizeof(char), mi1[j].cigarSize,
		       out1);

	  tmp = fwrite(&(mi1[j].mdSize), sizeof(int), 1, out1);
	  tmp = fwrite((mi1[j].md), sizeof(char), mi1[j].mdSize,
		       out1);

	  _msf_oeaMapping[i * 2]++;
	}
      }
    }

    char *seq1, *seq2, *rseq1, *rseq2, *qual1, *qual2;

    rqual1[SEQ_LENGTH] = '\0';
    rqual2[SEQ_LENGTH] = '\0';
    rqual1[0] = '\0';
    rqual2[0] = '\0';

    seq1 = _msf_seqList[i * 2].seq;
    rseq1 = _msf_seqList[i * 2].rseq;
    qual1 = _msf_seqList[i * 2].qual;

    strncpy(rqual1, _msf_seqList[i * 2].qual, SEQ_LENGTH);

    seq2 = _msf_seqList[i * 2 + 1].seq;
    rseq2 = _msf_seqList[i * 2 + 1].rseq;
    qual2 = _msf_seqList[i * 2 + 1].qual;

    strncpy(rqual2, _msf_seqList[i * 2 + 1].qual, SEQ_LENGTH);

    if (pairedEndDiscordantMode) {
      for (k = 0; k < size1; k++) {
	mi1[k].score = calculateScore(mi1[k].loc,
				      (mi1[k].dir == -1) ? rseq1 : seq1,
				      (mi1[k].dir == -1) ? rqual1 : qual1, mi1[k].cigar);
      }

      for (k = 0; k < size2; k++) {
	mi2[k].score = calculateScore(mi2[k].loc,
				      (mi2[k].dir == -1) ? rseq2 : seq2,
				      (mi2[k].dir == -1) ? rqual2 : qual2, mi2[k].cigar);
      }

    }

    
    /* CALKAN MAPQ FOR PE */
    if (seqFastq){
      for (j = 0; j < size1; j++) {
	if (mi1[j].err != 0){
	  bestHitMappingInfo[i*2].tprob += mapProb(i*2, mi1[j].md, mi1[j].dir, mi1[j].err);
	}
      }
      for (k = 0; k < size2; k++) {
	if (mi2[k].err != 0){
	  bestHitMappingInfo[i*2+1].tprob += mapProb((i*2+1), mi2[k].md, mi2[k].dir, mi2[k].err);
	}
      }
    }
    
    if (pairedEndDiscordantMode) {
      for (j = 0; j < size1; j++) {
	for (k = 0; k < size2; k++) {
	  if (
	      (pairedEndModePE && 
	       (
		(mi2[k].loc - mi1[j].loc
		 >= minPairEndedDiscordantDistance
		 && mi2[k].loc - mi1[j].loc
		 <= maxPairEndedDiscordantDistance
		 && mi1[j].dir > 0 && mi2[k].dir < 0)
		
		||
		
		(mi1[j].loc - mi2[k].loc
		 >= minPairEndedDiscordantDistance
		 && mi1[j].loc - mi2[k].loc
		 <= maxPairEndedDiscordantDistance
		 && mi1[j].dir < 0 && mi2[k].dir > 0)
		) ) 

	      ||  // CALKAN MPPE
	      
	      (pairedEndModeMP &&
	       (
		(mi2[k].loc - mi1[j].loc
		 >= minPairEndedDiscordantDistance
		 && mi2[k].loc - mi1[j].loc
		 <= maxPairEndedDiscordantDistance
		 && mi1[j].dir < 0 && mi2[k].dir > 0)
		
		||
		
		(mi1[j].loc - mi2[k].loc
		 >= minPairEndedDiscordantDistance
		 && mi1[j].loc - mi2[k].loc
		 <= maxPairEndedDiscordantDistance
		 && mi1[j].dir > 0 && mi2[k].dir < 0)
		) ) 
	      
	      )
	    {
	    

	      //POSSIBLE CONCORDANT
	      if(_msf_readHasConcordantMapping[i] == 0)
		{
		  setPairFullMappingInfo(i, mi1[j], mi2[k]);
		  _msf_readHasConcordantMapping[i] = 1;
		  _msf_seqList[i * 2].hits[0] = 1;
		  _msf_seqList[i * 2 + 1].hits[0] = 1;
		} else {
		if (bestHitMappingInfo[i * 2].err
		    + bestHitMappingInfo[i * 2 + 1].err
		    >= mi1[j].err + mi2[k].err) {

		  if (bestHitMappingInfo[i * 2].err
		      + bestHitMappingInfo[i * 2 + 1].err
		      == mi1[j].err + mi2[k].err
		      && findNearest(
				     abs(
					 bestHitMappingInfo[i * 2 + 1].loc
					 - bestHitMappingInfo[i * 2].loc),
				     abs(mi2[k].loc - mi1[j].loc),
				     meanDistanceMapping) == 0) {
		    continue;
		  }
		  setPairFullMappingInfo(i, mi1[j], mi2[k]);
		}
	      }
	    }
	  //DISCORDANT TO TEMP FILE FOR POST PROCESSING
	  else if (_msf_readHasConcordantMapping[i] == 0
		   && _msf_seqHits[i * 2] != 0
		   && _msf_seqHits[i * 2 + 1] != 0) {

	    int rNo = i;
	    int loc = mi1[j].loc * mi1[j].dir;
	    int err = mi1[j].err;
	    float sc = mi1[j].score;

	    char l = strlen(_msf_refGenName);

	    if (_msf_discordantMapping[i * 2]
		< maxDiscordantOutput) {

	      tmp = fwrite(&rNo, sizeof(int), 1, out);

	      tmp = fwrite(&l, sizeof(char), 1, out);
	      tmp = fwrite(_msf_refGenName, sizeof(char), l, out);

	      tmp = fwrite(&loc, sizeof(int), 1, out);
	      tmp = fwrite(&err, sizeof(int), 1, out);
	      tmp = fwrite(&sc, sizeof(float), 1, out);

	      tmp = fwrite(&(mi1[j].cigarSize), sizeof(int), 1,
			   out);
	      tmp = fwrite((mi1[j].cigar), sizeof(char),
			   mi1[j].cigarSize, out);

	      tmp = fwrite(&(mi1[j].mdSize), sizeof(int), 1, out);
	      tmp = fwrite((mi1[j].md), sizeof(char),
			   mi1[j].mdSize, out);

	      loc = mi2[k].loc * mi2[k].dir;
	      err = mi2[k].err;
	      sc = mi2[k].score;

	      tmp = fwrite(&loc, sizeof(int), 1, out);
	      tmp = fwrite(&err, sizeof(int), 1, out);
	      tmp = fwrite(&sc, sizeof(float), 1, out);

	      tmp = fwrite(&(mi2[k].cigarSize), sizeof(int), 1,
			   out);
	      tmp = fwrite((mi2[k].cigar), sizeof(char),
			   mi2[k].cigarSize, out);

	      tmp = fwrite(&(mi2[k].mdSize), sizeof(int), 1, out);
	      tmp = fwrite((mi2[k].md), sizeof(char),
			   mi2[k].mdSize, out);

	      _msf_discordantMapping[i * 2]++;
	    }
	    //SET THE BEST DISCORDANT
	    //BEGIN {Farhad Hormozdiari}
	    if (bestHitMappingInfo[i * 2].loc == -1
		&& bestHitMappingInfo[i * 2 + 1].loc == -1
		&& _msf_readHasConcordantMapping[i] == 0) {
	      setPairFullMappingInfo(i, mi1[j], mi2[k]);
	      _msf_seqList[i * 2].hits[0] = 1;
	      _msf_seqList[i * 2 + 1].hits[0] = 1;
	    } else if (bestHitMappingInfo[i * 2].err
		       + bestHitMappingInfo[i * 2 + 1].err
		       >= mi1[j].err + mi2[k].err
		       && _msf_readHasConcordantMapping[i] == 0) {
	      if (bestHitMappingInfo[i * 2].err
		  + bestHitMappingInfo[i * 2 + 1].err
		  == mi1[j].err + mi2[k].err
		  && findNearest(
				 abs(
				     bestHitMappingInfo[i * 2 + 1].loc
				     - bestHitMappingInfo[i * 2].loc),
				 abs(mi1[j].loc - mi2[k].loc),
				 meanDistanceMapping) == 0) {
		continue;
	      }
	      setPairFullMappingInfo(i, mi1[j], mi2[k]);
	    }
	    //END {Farhad Hormozdiari}
	  }
	}
      }
    } else {
      for (j = 0; j < size1; j++) {
	for (k = 0; k < size2; k++) {
	  if ((mi2[k].loc - mi1[j].loc >= minPairEndedDistance
	       && mi2[k].loc - mi1[j].loc <= maxPairEndedDistance
	       && mi1[j].dir > 0 && mi2[k].dir < 0)
	      || (mi1[j].loc - mi2[k].loc >= minPairEndedDistance
		  && mi1[j].loc - mi2[k].loc
		  <= maxPairEndedDistance
		  && mi1[j].dir < 0 && mi2[k].dir > 0)) {
	    char *seq;
	    char *qual;
	    char d1;
	    char d2;
	    int isize;
	    int proper = 0;
	    // ISIZE CALCULATION
	    // The distance between outer edges
	    isize = abs(mi1[j].loc - mi2[k].loc) + SEQ_LENGTH - 2;
	    if (mi1[j].loc - mi2[k].loc > 0) {
	      isize *= -1;
	    }

	    d1 = (mi1[j].dir == -1) ? 1 : 0;
	    d2 = (mi2[k].dir == -1) ? 1 : 0;

	    //SET THE READ HAS CONCORDANT MAPPING
	    _msf_readHasConcordantMapping[i] = 1;

	    if (d1) {
	      seq = rseq1;
	      qual = rqual1;
	    } else {
	      seq = seq1;
	      qual = qual1;
	    }

	    if ((mi1[j].loc < mi2[k].loc && !d1 && d2)
		|| (mi1[j].loc > mi2[k].loc && d1 && !d2)) {
	      proper = 2;
	    } else {
	      proper = 0;
	    }

	    _msf_output.POS = mi1[j].loc;
	    _msf_output.MPOS = mi2[k].loc;
	    _msf_output.FLAG = 1 + proper + 16 * d1 + 32 * d2 + 64;
	    _msf_output.ISIZE = isize;
	    _msf_output.SEQ			= seq;
	    _msf_output.QUAL		= qual;
	    _msf_output.QNAME = _msf_seqList[i * 2].name;
	    _msf_output.RNAME = _msf_refGenName;
	    _msf_output.MAPQ = 255;
	    _msf_output.CIGAR = cigar;
	    _msf_output.MRNAME = "=";

	    _msf_output.optSize = 2;
	    _msf_output.optFields = _msf_optionalFields;

	    _msf_optionalFields[0].tag = "NM";
	    _msf_optionalFields[0].type = 'i';
	    _msf_optionalFields[0].iVal = mi1[j].err;

	    _msf_optionalFields[1].tag = "MD";
	    _msf_optionalFields[1].type = 'Z';
	    _msf_optionalFields[1].sVal = mi1[j].md;

	    if (!bestMode)
	      output(_msf_output);

	    if (d2) {
	      seq = rseq2;
	      qual = rqual2;
	    } else {
	      seq = seq2;
	      qual = qual2;
	    }

	    _msf_output.POS = mi2[k].loc;
	    _msf_output.MPOS = mi1[j].loc;
	    _msf_output.FLAG = 1 + proper + 16 * d2 + 32 * d1 + 128;
	    _msf_output.ISIZE = -isize;
	    _msf_output.SEQ		= seq;
	    _msf_output.QUAL		= qual;
	    _msf_output.QNAME = _msf_seqList[i * 2].name;
	    _msf_output.RNAME = _msf_refGenName;
	    _msf_output.MAPQ = 255;
	    _msf_output.CIGAR = cigar;
	    _msf_output.MRNAME = "=";

	    _msf_output.optSize = 2;
	    _msf_output.optFields = _msf_optionalFields;

	    _msf_optionalFields[0].tag = "NM";
	    _msf_optionalFields[0].type = 'i';
	    _msf_optionalFields[0].iVal = mi2[k].err;
	    

	    _msf_optionalFields[1].tag = "MD";
	    _msf_optionalFields[1].type = 'Z';
	    _msf_optionalFields[1].sVal = mi2[k].md;

	    if (!bestMode)
	      output(_msf_output);
	    //SET THE BEST CONCORDANT
	    //BEGIN {Farhad Hormozdiari}
	    if (bestHitMappingInfo[i * 2].loc == -1
		&& bestHitMappingInfo[i * 2 + 1].loc == -1) {
	      setPairFullMappingInfo(i, mi1[j], mi2[k]);
	    } else {
	      if (bestHitMappingInfo[i * 2].err
		  + bestHitMappingInfo[i * 2 + 1].err
		  >= mi1[j].err + mi2[k].err) {

		if (bestHitMappingInfo[i * 2].err
		    + bestHitMappingInfo[i * 2 + 1].err
		    == mi1[j].err + mi2[k].err
		    && findNearest(
				   abs(
				       bestHitMappingInfo[i * 2
							  + 1].loc
				       - bestHitMappingInfo[i
							    * 2].loc),
				   abs(mi2[k].loc - mi1[j].loc),
				   meanDistanceMapping) == 0) {
		  continue;
		}
		setPairFullMappingInfo(i, mi1[j], mi2[k]);
	      }
	    }
	    //END   {Farhad Hormozdiari}
	  }
	}
      }

    }
  }

  if (pairedEndDiscordantMode) {
    fclose(out);
    fclose(out1);
  }

  for (i = 0; i < _msf_openFiles; i++) {
    fclose(in1[i]);
    fclose(in2[i]);

    unlink(fname1[i]);
    unlink(fname2[i]);
  }

  tmp++;

  freeMem(mi1, sizeof(FullMappingInfo) * _msf_maxLSize);
  freeMem(mi2, sizeof(FullMappingInfo) * _msf_maxRSize);

  _msf_openFiles = 0;

  /* calkan counter */
  int unmappedCnt = 0;
  for (i = 0; i < _msf_seqListSize; i++) {
    if (_msf_seqHits[i] == 0) unmappedCnt++;
  }

  mappedSeqCnt = _msf_seqListSize - unmappedCnt; 

}

/**********************************************/
/**********************************************/
/**********************************************/
/**********************************************/
float str2int(char *str, int index1, int index2) {
  char tmp[SEQ_MAX_LENGTH];
  strncpy(tmp, &str[index1], index2 - index1);
  tmp[index2 - index1] = '\0';
  return atol(tmp);
}

double binomial_coefficient(int n, int k){
  double ret;
  int i;
  ret = 1.0;
  
  for (i=0; i<k; i++){
    ret *= (n - i);
    ret /= (k - i);
  }
 
  return ret; 

}

float calculateScore(int index, char *seq, char *qual, char *md) {
  int i;
  int j;
  char *ref;
  char *ver;

  ref = _msf_refGen + index - 1;
  ver = seq;
  float score = 1;

  char tmp[2 * SEQ_MAX_LENGTH];
  int value = 0;
  int end = 0;
  int index1 = 0;
  int index2 = 0;

  i = 0;
  while (1) {

    if (i >= strlen(md))
      break;

    index1 = i;

    while (md[i] >= '0' && md[i] <= '9') {
      i++;
    }

    index2 = i;

    value = str2int(md, index1, index2);

    if (md[i] == 'M') {
      for (j = 0; j < value; j++) {
	tmp[end] = 'M';
	end++;
      }
    } else if (md[i] == 'I') {
      for (j = 0; j < value; j++) {
	tmp[end] = 'I';
	end++;
      }

    } else if (md[i] == 'D') {
      for (j = 0; j < value; j++) {
	tmp[end] = 'D';
	end++;
      }
    }
    i++;
  }

  tmp[end] = '\0';

  j = 0;

  for (i = 0; i < end; i++) {
    if (tmp[i] == 'M') {
      if (*ref != *ver) {
	score *= 0.001 + 1 / pow(10, ((qual[j] - 33) / 10.0));
      }

      ref++;
      ver++;
      j++;
    } else if (tmp[i] == 'I') {
      ver++;
      j++;
      score *= 0.0003;  // 0.0001 + 0.0002;  0.0001: indel rate in normal human, 0.0002: indel error rate in Illumina
    } else if (tmp[i] == 'D') {
      ref++;
      score *= 0.0003; // 0.0001 + 0.0002
    }
  }

  return score;
}

int matoi(char *str, int start, int end) {
  int i = 0;
  char tmp[SEQ_MAX_LENGTH];

  for (i = 0; i < end - start; i++)
    tmp[i] = str[start + i];
  tmp[i] = '\0';

  return atoi(tmp);
}

void convertCigarToMatrix(char *cigar, int cigar_size, char * matrix) {
  int i = 0;
  int j = 0;

  int start = 0;
  int size = 0;

  matrix[0] = '\0';

  while (i < cigar_size) {
    if (cigar[i] >= '0' && cigar[i] <= '9') {
      start = i;

      while (cigar[i] >= '0' && cigar[i] <= '9' && i < cigar_size)
	i++;

      int value = matoi(cigar, start, i);
      for (j = 0; j < value; j++) {
	if (cigar[i] == 'M')
	  matrix[size] = 'M';
	else if (cigar[i] == 'D')
	  matrix[size] = 'D';
	else if (cigar[i] == 'I')
	  matrix[size] = 'I';
	size++;
      }
    }
    i++;
  }
  matrix[size] = '\0';
}

void convertMDToMatrix(char *md, int md_size, char * matrix) {
  int i = 0;
  int j = 0;

  int start = 0;
  int size = 0;

  matrix[0] = '\0';

  while (i < md_size) {
    if (md[i] >= '0' && md[i] <= '9') {
      start = i;

      while (md[i] >= '0' && md[i] <= '9' && i < md_size)
	i++;

      int value = matoi(md, start, i);
      for (j = 0; j < value; j++) {
	matrix[size] = 'M';
	size++;
      }
      i--;
    } else if (md[i] == '^') {
      matrix[size] = 'D';
      size++;
    } else {
      matrix[size] = md[i];
      size++;
    }
    i++;
  }
  matrix[size] = '\0';
}

void convertMDCigarToMatrix(char *cigar, int cigar_size, char *md, int md_size,
			    char *matrix) {
  int i = 0;
  int j = 0;

  int size = 0;

  char tmp1[SEQ_MAX_LENGTH];
  char tmp2[SEQ_MAX_LENGTH];
  convertMDToMatrix(md, md_size, tmp2);

  convertCigarToMatrix(cigar, cigar_size, tmp1);

  while (i < strlen(tmp1)) {
    if (tmp1[i] == 'M') {
      if (j < strlen(tmp2)) {
	if (tmp2[j] == 'M') {
	  matrix[size] = 'M';
	  size++;
	}
	if (tmp2[j] != 'M') {
	  matrix[size] = tmp2[j];
	  size++;
	}
      } else {
	matrix[size] = 'M';
	size++;
      }
    } else if (tmp1[i] == 'D') {
      matrix[size] = 'D';
      size++;
      j++;
      matrix[size] = tmp2[j];
      size++;

    } else if (tmp1[i] == 'I') {
      matrix[size] = 'I';
      size++;
    }

    i++;
    if (j < strlen(tmp2))
      j++;
  }

  if (strlen(tmp1))

    matrix[size] = '\0';

}

void convertInsertion(char * in_matrix, char * seq, char *out_matrix) {
  int i = 0;
  int j = 0;
  int size = 0;

  while (i < strlen(in_matrix)) {
    if (in_matrix[i] == 'M') {
      out_matrix[size] = 'M';
      size++;
      j++;
    } else if (in_matrix[i] == 'D') {
      out_matrix[size] = 'D';
      size++;

      i++;
      j++;

      out_matrix[size] = seq[j];
      j++;
      size++;
    } else if (in_matrix[i] == 'I') {
      out_matrix[size] = 'I';
      size++;
      out_matrix[size] = seq[j];
      size++;
      j++;
    } else {
      out_matrix[size] = in_matrix[i];
      size++;
      j++;
    }
    i++;
  }
  out_matrix[size] = '\0';
}

/**********************************************/
void outputPairedEndDiscPP() {
  char tmp_matrix1[SEQ_MAX_LENGTH];
  char tmp_matrix2[SEQ_MAX_LENGTH];

  char matrix1[SEQ_MAX_LENGTH];
  char matrix2[SEQ_MAX_LENGTH];

  char cigar1[MAX_CIGAR_SIZE];
  char editString1[2 * SEQ_MAX_LENGTH];

  char cigar2[MAX_CIGAR_SIZE];
  char editString2[2 * SEQ_MAX_LENGTH];
  char seq1[SEQ_LENGTH + 1];

  char seq2[SEQ_LENGTH + 1];

  char genName[SEQ_LENGTH];
  char fname1[FILE_NAME_LENGTH];
  char fname2[FILE_NAME_LENGTH];
  char l;
  int l_size;
  int loc1, loc2;
  int err1, err2;
  char dir1, dir2;
  float sc1, sc2, lsc = 0;
  int flag = 0;
  int rNo, lrNo = -1;
  int tmp;
  FILE *in, *out;

  sprintf(fname1, "%s__%s__disc", mappingOutputPath, mappingOutput);
  sprintf(fname2, "%s%s_DIVET.vh", mappingOutputPath, mappingOutput);

  in = fileOpen(fname1, "r");
  out = fileOpen(fname2, "w");

  if (in != NULL) {
    flag = fread(&rNo, sizeof(int), 1, in);
  } else {
    flag = 0;
  }

  seq1[SEQ_LENGTH] = '\0';
  seq2[SEQ_LENGTH] = '\0';

  while (flag) {
    tmp = fread(&l, sizeof(char), 1, in);
    tmp = fread(genName, sizeof(char), l, in);
    genName[(int) l] = '\0';
    tmp = fread(&loc1, sizeof(int), 1, in);
    tmp = fread(&err1, sizeof(int), 1, in);
    tmp = fread(&sc1, sizeof(float), 1, in);

    tmp = fread(&l_size, sizeof(int), 1, in);
    tmp = fread(cigar1, sizeof(char), l_size, in);
    cigar1[(int) l_size] = '\0';

    tmp = fread(&l_size, sizeof(int), 1, in);
    tmp = fread(editString1, sizeof(char), l_size, in);
    editString1[(int) l_size] = '\0';

    tmp = fread(&loc2, sizeof(int), 1, in);
    tmp = fread(&err2, sizeof(int), 1, in);
    tmp = fread(&sc2, sizeof(float), 1, in);

    tmp = fread(&l_size, sizeof(int), 1, in);
    tmp = fread(cigar2, sizeof(char), l_size, in);
    cigar2[(int) l_size] = '\0';

    tmp = fread(&l_size, sizeof(int), 1, in);
    tmp = fread(editString2, sizeof(char), l_size, in);
    editString2[(int) l_size] = '\0';

    convertMDCigarToMatrix(cigar1, strlen(cigar1), editString1,
			   strlen(editString1), tmp_matrix1);
    convertMDCigarToMatrix(cigar2, strlen(cigar2), editString2,
			   strlen(editString2), tmp_matrix2);

    /* CALKAN: GO OVER THIS VERY CAREFULLY FOR PE vs MP */

    if (_msf_readHasConcordantMapping[rNo] == 0 && _msf_discordantMapping[rNo * 2] < maxDiscordantOutput ) {

      dir1 = dir2 = 'F';

      strncpy(seq1, _msf_seqList[rNo * 2].seq, SEQ_LENGTH);
      strncpy(seq2, _msf_seqList[rNo * 2 + 1].seq, SEQ_LENGTH);

      if (loc1 < 0) {
	dir1 = 'R';
	loc1 = -loc1;

	strncpy(seq1, _msf_seqList[rNo * 2].rseq, SEQ_LENGTH);
      }

      if (loc2 < 0) {
	dir2 = 'R';
	loc2 = -loc2;

	strncpy(seq2, _msf_seqList[rNo * 2 + 1].rseq, SEQ_LENGTH);
      }

      convertInsertion(tmp_matrix1, seq1, matrix1);
      convertInsertion(tmp_matrix2, seq2, matrix2);

      if (rNo != lrNo) {
	int j;
	for (j = 0; j < SEQ_LENGTH; j++) {
	  lsc += _msf_seqList[rNo * 2].qual[j]
	    + _msf_seqList[rNo * 2 + 1].qual[j];
	}
	lsc /= 2 * SEQ_LENGTH;
	lsc -= 33;
	lrNo = rNo;
      }

      char event = '\0';

      if (dir1 == dir2) {
	event = 'V';
      } 
      else {
	if (pairedEndModePE && loc1 < loc2 && dir1 == 'R' && dir2 == 'F') 
	  event = 'E';
	else if (pairedEndModeMP && loc1 < loc2 && dir1 == 'F' && dir2 == 'R') 
	  event = 'E';
	else if (pairedEndModePE && loc2 < loc1 && dir1 == 'F' && dir2 == 'R') 
	  event = 'E';
	else if (pairedEndModeMP && loc2 < loc1 && dir1 == 'R' && dir2 == 'F') 
	  event = 'E';
	else if (abs(loc2 - loc1) >= maxPairEndedDiscordantDistance) 
	  event = 'D';
	else 
	  event = 'I';	    	 
      }

      _msf_seqList[rNo * 2].hits[0] = 2;
      fprintf(out,
	      "%s\t%s\t%d\t%d\t%c\t=\t%d\t%d\t%c\t%c\t%d\t%0.0f\t%e\n",
	      _msf_seqList[rNo * 2].name, genName, loc1,
	      (loc1 + SEQ_LENGTH - 1), dir1, loc2,
	      (loc2 + SEQ_LENGTH - 1), dir2, event, (err1 + err2),
	      lsc, sc1 * sc2);
      
      //	      lsc, sc1 * sc2 * binomial_coefficient(2 * SEQ_LENGTH, (err1 + err2)));

    }
    flag = fread(&rNo, sizeof(int), 1, in);
  }

  tmp++;

  fclose(in);
  fclose(out);

  unlink(fname1);
}

void finalizeOEAReads(char *fileName) {
  FILE *fp_out1;
  FILE * in;

  char genName[SEQ_LENGTH];

  char fname1[FILE_NAME_LENGTH];
  char fname2[FILE_NAME_LENGTH];

  char l = 0;
  int loc1 = 0;

  int err1;

  char d;

  float sc1 = 0;
  int flag = 0;
  int rNo = -1;
  int tmp = 0;

  int cigarSize = 0;
  int mdSize = 0;

  char cigar[SEQ_LENGTH + 1];
  char md[SEQ_LENGTH + 1];

  char *seq1, *seq2, *qual1, *qual2;
  char rqual1[SEQ_LENGTH + 1];

  seq1 = NULL;
  seq2 = NULL;
  qual1 = NULL;
  qual2 = NULL;

  rqual1[0] = '\0';

  sprintf(fname1, "%s%s_OEA.sam", mappingOutputPath, mappingOutput);

  fp_out1 = fileOpen(fname1, "w");

  SAMheaderTX(fp_out1, 0);
  in = NULL;
  if (pairedEndDiscordantMode) {
    sprintf(fname2, "%s__%s__oea", mappingOutputPath, mappingOutput);

    in = fileOpen(fname2, "r");
  }

  if (in != NULL) {
    flag = fread(&rNo, sizeof(int), 1, in);
  } else {
    flag = 0;
  }

  while (flag) {
    cigar[0] = '\0';
    md[0] = '\0';

    tmp = fread(&l, sizeof(char), 1, in);
    tmp = fread(genName, sizeof(char), l, in);

    genName[(int) l] = '\0';

    tmp = fread(&loc1, sizeof(int), 1, in);
    tmp = fread(&err1, sizeof(int), 1, in);
    tmp = fread(&sc1, sizeof(float), 1, in);

    tmp = fread(&cigarSize, sizeof(int), 1, in);
    tmp = fread(cigar, sizeof(char), cigarSize, in);

    cigar[cigarSize] = '\0';

    tmp = fread(&mdSize, sizeof(int), 1, in);
    tmp = fread(md, sizeof(char), mdSize, in);
    md[mdSize] = '\0';

    d = 1;

    if (loc1 < 0) {
      d = -1;
      loc1 *= -1;

      seq1 = _msf_seqList[rNo].rseq;
      reverse(_msf_seqList[rNo].qual, rqual1, SEQ_LENGTH);
      rqual1[SEQ_LENGTH] = '\0';
      qual1 = rqual1;
    } else {
      seq1 = _msf_seqList[rNo].seq;
      qual1 = _msf_seqList[rNo].qual;
      qual1[SEQ_LENGTH] = '\0';
    }

    if (rNo % 2 == 0) {
      seq2 = _msf_seqList[rNo + 1].seq;
      qual2 = _msf_seqList[rNo + 1].qual;
      qual2[SEQ_LENGTH] = '\0';
    } else {
      seq2 = _msf_seqList[rNo - 1].seq;
      qual2 = _msf_seqList[rNo - 1].qual;
      qual2[SEQ_LENGTH] = '\0';
    }

    
    if (_msf_seqHits[rNo] != 0 && _msf_seqHits[rNo] < maxOEAOutput
	&& _msf_seqHits[(rNo % 2 == 0) ? rNo + 1 : rNo - 1] == 0) {
      _msf_output.POS = loc1;
      _msf_output.MPOS = 0;
      _msf_output.FLAG =
	(rNo % 2 == 0) ? 1 + 4 + 32 * d + 128 : 1 + 8 + 16 * d + 64;
      _msf_output.ISIZE = 0;
      _msf_output.SEQ = seq1;
      _msf_output.QUAL = qual1;
      _msf_output.QNAME = _msf_seqList[rNo].name;
      _msf_output.RNAME = genName;
      _msf_output.MAPQ = 255;
      _msf_output.CIGAR = cigar;
      _msf_output.MRNAME = "=";

      _msf_output.optSize = 4;
      _msf_output.optFields = _msf_optionalFields;

      _msf_optionalFields[0].tag = "NM";
      _msf_optionalFields[0].type = 'i';
      _msf_optionalFields[0].iVal = err1;

      _msf_optionalFields[1].tag = "MD";
      _msf_optionalFields[1].type = 'Z';
      _msf_optionalFields[1].sVal = md;

      //for the OEA reads
      _msf_optionalFields[2].tag = "NS";
      _msf_optionalFields[2].type = 'Z';
      _msf_optionalFields[2].sVal = seq2;

      _msf_optionalFields[3].tag = "NQ";
      _msf_optionalFields[3].type = 'Z';
      _msf_optionalFields[3].sVal = qual2;

      outputSAM(fp_out1, _msf_output);

      _msf_seqList[rNo].hits[0] = -1;
      _msf_seqList[(rNo % 2 == 0) ? rNo + 1 : rNo - 1].hits[0] = -1;
    }
    else if(_msf_seqHits[rNo] != 0 && _msf_seqHits[(rNo % 2 == 0) ? rNo + 1 : rNo - 1] == 0)
      {
	_msf_seqList[rNo].hits[0] = -1;
	_msf_seqList[(rNo % 2 == 0) ? rNo + 1 : rNo - 1].hits[0] = -1;
      }
    flag = fread(&rNo, sizeof(int), 1, in);
  }

  tmp++;

  fclose(in);
  unlink(fname2);

  fclose(fp_out1);
}

/* disabled until completed 


   void outputTransChromosomal(char *fileName1, char *fileName2, FILE * fp_out)
   {
   int i = 0;
   int j = 0;
   int k = 0;

   char *index;

   int size1 = 0;
   int size2 = 0;

   FILE *fp1 = NULL;
   FILE *fp2 = NULL;

   char geneFileName1[FILE_NAME_LENGTH];
   char geneFileName2[FILE_NAME_LENGTH];

   FullMappingInfoLink *miL = getMem(_msf_seqListSize * sizeof(FullMappingInfoLink));
   FullMappingInfoLink *miR = getMem(_msf_seqListSize * sizeof(FullMappingInfoLink));


   if(fileName1 != NULL && fileName2 != NULL)
   {

   fp1 = fileOpen(fileName1, "r");
   fp2 = fileOpen(fileName2, "r");

   index = strstr(fileName1, "__");
   strncpy(geneFileName1, index + 2 * sizeof(char), strstr(index + 2, "__") - index - 2);
   geneFileName1[strstr(index + 2, "__") - index - 2] = '\0';

   index = strstr(fileName2, "__");
   strncpy(geneFileName2, index + 2 * sizeof(char), strstr(index + 2, "__") - index - 2);
   geneFileName2[strstr(index + 2, "__") - index - 2] = '\0';


   for(i = 0; i < _msf_seqListSize / 2; i++)
   {
   fread(&size1, sizeof(int), 1, fp1);
   fread(&size2, sizeof(int), 1, fp2);

   miL[i].mi = getMem(size1 * sizeof(FullMappingInfo) );
   miR[i].mi = getMem(size2 * sizeof(FullMappingInfo) );

   miL[i].size = size1;
   miR[i].size = size2;

   for(j = 0; j < size1; j++)
   {
   fread(&(miL[i].mi[j].loc), sizeof(int), 1, fp1);

   fread (&(miL[i].mi[j].err), sizeof(int), 1, fp1);

   fread (&(miL[i].mi[j].cigarSize), sizeof(int), 1, fp1);
   fread ((miL[i].mi[j].cigar), sizeof(char), miL[i].mi[j].cigarSize+1, fp1);

   fread (&(miL[i].mi[j].mdSize), sizeof(int), 1, fp1);
   fread ((miL[i].mi[j].md), sizeof(char), miL[i].mi[j].mdSize+1, fp1);

   miL[i].mi[j].dir = 1;
   if(miL[i].mi[j].loc < 1)
   {
   miL[i].mi[j].loc *= -1;
   miL[i].mi[j].dir = -1;
   }
   }
   for(k = 0; k < size2; k++)
   {
   fread(&(miR[i].mi[k].loc), sizeof(int), 1, fp2);

   fread (&(miR[i].mi[k].err), sizeof(int), 1, fp2);

   fread (&(miR[i].mi[k].cigarSize), sizeof(int), 1, fp2);
   fread ((miR[i].mi[k].cigar), sizeof(char), miR[i].mi[k].cigarSize+1, fp2);

   fread (&(miR[i].mi[k].mdSize), sizeof(int), 1, fp2);
   fread ((miR[i].mi[k].md), sizeof(char), miR[i].mi[k].mdSize+1, fp2);

   miR[i].mi[k].dir = 1;
   if(miR[i].mi[k].loc < 1)
   {
   miR[i].mi[k].loc *= -1;
   miR[i].mi[k].dir = -1;
   }
   }
   if(_msf_readHasConcordantMapping[i] == 0 && size1 != 0 && size2 != 0 && (size1 * size2 < MAX_TRANS_CHROMOSAL_OUTPUT))
   {
   int d1 = 0;
   int d2 = 0;
   char *seq, *qual;
   char *seq1, *seq2, *rseq1, *rseq2, *qual1, *qual2;
   char rqual1[SEQ_LENGTH+1], rqual2[SEQ_LENGTH+1];
   rqual1[SEQ_LENGTH] = rqual2[SEQ_LENGTH] = '\0';
   seq1 = _msf_seqList[i*2].seq;
   rseq1 = _msf_seqList[i*2].rseq;
   qual1 = _msf_seqList[i*2].qual;
   reverse(_msf_seqList[i*2].qual, rqual1, SEQ_LENGTH);

   seq2 = _msf_seqList[i*2+1].seq;
   rseq2 = _msf_seqList[i*2+1].rseq;
   qual2 = _msf_seqList[i*2+1].qual;
   reverse(_msf_seqList[i*2+1].qual, rqual2, SEQ_LENGTH);

   for(j = 0; j < size1; j++)
   {
   d1 = (miL[i].mi[j].dir == -1)?1:0;

   if ( d1 )
   {
   seq = rseq1;
   qual = rqual1;
   }
   else
   {
   seq = seq1;
   qual = qual1;
   }

   for(k = 0; k < size2; k++)
   {

   d2 = (miR[i].mi[k].dir == -1)?1:0;

   _msf_output.POS                 = miL[i].mi[j].loc;
   _msf_output.MPOS                = miR[i].mi[k].loc;
   _msf_output.FLAG                = 0;
   _msf_output.ISIZE               = 0;
   _msf_output.SEQ                 = seq;
   _msf_output.QUAL                = qual;
   _msf_output.QNAME               = _msf_seqList[i*2].name;
   _msf_output.RNAME               = geneFileName1;
   _msf_output.MAPQ                = 255;
   _msf_output.CIGAR               = miL[i].mi[j].cigar;
   _msf_output.MRNAME              = "=";

   _msf_output.optSize     = 2;
   _msf_output.optFields   = _msf_optionalFields;

   _msf_optionalFields[0].tag = "NM";
   _msf_optionalFields[0].type = 'i';
   _msf_optionalFields[0].iVal = miL[i].mi[j].err;

   _msf_optionalFields[1].tag = "MD";
   _msf_optionalFields[1].type = 'Z';
   _msf_optionalFields[1].sVal = miL[i].mi[j].md;


   if ( d2 )
   {
   seq = rseq2;
   qual = rqual2;
   }
   else
   {
   seq = seq2;
   qual = qual2;
   }

   outputSAM(fp_out, _msf_output);


   _msf_output.POS                 = miR[i].mi[k].loc;
   _msf_output.MPOS                = miL[i].mi[j].loc;
   _msf_output.FLAG                = 0;
   _msf_output.ISIZE               = 0;
   _msf_output.SEQ                 = seq;
   _msf_output.QUAL                = qual;
   _msf_output.QNAME               = _msf_seqList[i*2+1].name;
   _msf_output.RNAME               = geneFileName2;
   _msf_output.MAPQ                = 255;
   _msf_output.CIGAR               = miR[i].mi[k].cigar;
   _msf_output.MRNAME              = "=";

   _msf_output.optSize     = 2;
   _msf_output.optFields   = _msf_optionalFields;

   _msf_optionalFields[0].tag = "NM";
   _msf_optionalFields[0].type = 'i';
   _msf_optionalFields[0].iVal = miR[i].mi[k].err;

   _msf_optionalFields[1].tag = "MD";
   _msf_optionalFields[1].type = 'Z';
   _msf_optionalFields[1].sVal = miR[i].mi[k].md;

   outputSAM(fp_out, _msf_output);

   }
   }
   }
   }

   }

   for(i = 0; i < _msf_seqListSize / 2; i++)
   {
   freeMem(miL[i].mi, miL[i].size * sizeof(FullMappingInfo));
   freeMem(miR[i].mi, miR[i].size * sizeof(FullMappingInfo));
   }

   freeMem(miL, _msf_seqListSize * sizeof(FullMappingInfoLink));
   freeMem(miR, _msf_seqListSize * sizeof(FullMappingInfoLink));

   fclose(fp1);
   fclose(fp2);
   }

*/

/*
  if flag is 1 it will output all the possible trans chromsal mapping
  otherwise only tmp file will be delete

*/

void outputAllTransChromosomal(int flag) {
  return;
  /* disabled until completed

     int i = 0;
     int j = 0;
     int k = 0;
     int l = 0;

     FILE *fp_out = NULL;
     char fname1[FILE_NAME_LENGTH];

     if(flag)
     {
     fp_out = fileOpen(fname1, "w");

     sprintf(fname1, "%s%s_TRANSCHROMOSOMAL", mappingOutputPath, mappingOutput);

     i  = 0;
     for(j = i+1; j < _msf_maxFile; j++)
     {
     if(i != j)
     {
     for(k = 0; k < _msf_fileCount[i]; k++)
     {
     for(l = 0; l < _msf_fileCount[j]; l++)
     {
     outputTransChromosomal(_msf_fileName[i][k][0], _msf_fileName[j][l][1], fp_out);
     }// for l
     }// for k
     }// if
     }// for j
     }

     for(i = 0; i < _msf_maxFile; i++)
     {
     for(j = 0; j < _msf_fileCount[i]; j++)
     {
     unlink(_msf_fileName[i][j][0]);
     unlink(_msf_fileName[i][j][1]);
     }
     }
     if(flag)
     fclose(fp_out);
  */
}

