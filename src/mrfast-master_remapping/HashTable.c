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
#include "Common.h"
#include "RefGenome.h"
#include "HashTable.h"
/**********************************************/
FILE		*_ih_fp			= NULL;
IHashTable	*_ih_hashTable		= NULL;
int 		_ih_maxHashTableSize	= 0;
char		*_ih_refGen		= NULL;
char		*_ih_refGenName		= NULL;
long long	_ih_memUsage		= 0;
int		_ih_refGenOff		= 0;
/**********************************************/

int hashVal(char *seq)
{
  int i=0;
  int val=0, numericVal=0;

  while(i<WINDOW_SIZE)
    {
      switch (seq[i])
	{
	case 'A':
	  numericVal = 0;
	  break;
	case 'C':
	  numericVal = 1;
	  break;
	case 'G' :
	  numericVal = 2;
	  break;
	case 'T':
	  numericVal = 3;
	  break;
	default:
	  return -1;
	  break;
	}
      val = (val << 2)|numericVal;
      i++;
    }
  return val;
}


/**********************************************/
void freeIHashTableContent(IHashTable *hashTable, unsigned int maxSize)
{
  int i=0;
  for (i=0; i<maxSize; i++)
    {
      if (hashTable[i].locs != NULL)
	{
	  freeMem(hashTable[i].locs, (hashTable[i].locs[0]+1)*(sizeof(unsigned int)));
	  hashTable[i].locs = NULL;
	}
    }
}
/**********************************************/
void initSavingIHashTable(char *fileName)
{
  int tmp;
  _ih_fp = fileOpen(fileName, "w");
  unsigned char bIndex = 0;				// Bisulfite Index
	
  // First Two bytes are indicating the type of the index & window size
  tmp = fwrite(&bIndex, sizeof(bIndex), 1, _ih_fp);
  tmp = fwrite(&WINDOW_SIZE, sizeof(WINDOW_SIZE), 1, _ih_fp);
	
  if (tmp == 0)
    fprintf(stderr, "Write error while initializing hash table.\n");

}
/**********************************************/
void finalizeSavingIHashTable()
{
  fclose(_ih_fp);
}
/**********************************************/
void saveIHashTable(IHashTable *hashTable,  unsigned int size, unsigned int maxSize, char *refGen, char *refGenName, int refGenOffset)
{
  int tmp;
	
  // Every Chunk starts with a byte indicating whether it has extra info;
  unsigned char extraInfo = 0;
  tmp = fwrite (&extraInfo, sizeof(extraInfo), 1, _ih_fp);

  short len = strlen(refGenName);
  tmp = fwrite(&len, sizeof(len), 1, _ih_fp);
  tmp = fwrite(refGenName, sizeof(char), len, _ih_fp);

  tmp = fwrite(&refGenOffset, sizeof(refGenOffset), 1, _ih_fp);
	
  unsigned int refGenLength = strlen(refGen);
  tmp = fwrite(&refGenLength, sizeof(refGenLength), 1, _ih_fp);
	
  tmp = fwrite(refGen, sizeof(char), refGenLength, _ih_fp);
  tmp = fwrite(&size, sizeof(size), 1, _ih_fp);

  int i=0,j=0;
  unsigned char cnt=0;
  for (i=0; i<maxSize; i++)
    {
      if (hashTable[i].locs != NULL)
	{
	  tmp = fwrite(&i, sizeof(i), 1, _ih_fp);

	  if (hashTable[i].locs[0] < 250)
	    {
	      cnt = hashTable[i].locs[0];
	      tmp = fwrite(&cnt, sizeof(cnt), 1, _ih_fp);
	    }
	  else
	    {
	      cnt =0;
	      tmp = fwrite (&cnt, sizeof(cnt), 1, _ih_fp);
	      tmp = fwrite (&(hashTable[i].locs[0]), sizeof(hashTable[i].locs[0]), 1, _ih_fp);
	    }

	  for (j=1; j<=hashTable[i].locs[0]; j++)
	    tmp = fwrite(&(hashTable[i].locs[j]), sizeof(hashTable[i].locs[j]), 1, _ih_fp);
	}
    }

  if (tmp == 0)
    fprintf(stderr, "Write error while saving hash table.\n");
}
/**********************************************/
unsigned int addIHashTableLocation(IHashTable *hashTable, int hv, int location)
{
  unsigned int	sizeInc				= 0;

  if (hashTable[hv].locs == NULL)
    {
      sizeInc = 1;
      hashTable[hv].locs = getMem (sizeof(unsigned int)*2);
      hashTable[hv].locs[0]=1;
      hashTable[hv].locs[1]=location;
    }
  else
    {
      int size = hashTable[hv].locs[0];
      int i;
      unsigned int *tmp = getMem( (size + 2) * sizeof(unsigned int) );

      for (i = 0; i <= size; i++)
	{
	  tmp[i] = hashTable[hv].locs[i];
	}
      size++;
      tmp[0] = size;
      tmp[size] = location;
		
      freeMem(hashTable[hv].locs, (hashTable[hv].locs[0]*(sizeof(unsigned int))));
      hashTable[hv].locs = tmp;
    }
  return sizeInc;
}
/**********************************************/
void generateIHashTable(char *fileName, char *indexName)
{
  double          startTime           = getTime();
  unsigned int	hashTableSize		= 0;
  unsigned int 	hashTableMaxSize	= pow(4, WINDOW_SIZE);
  IHashTable		*hashTable			= getMem(sizeof(IHashTable)*hashTableMaxSize);
  char 			*refGenName;
  char			*refGen;
  int				refGenOff			= 0;
  int i, hv, l, flag;


  for ( i = 0; i < hashTableMaxSize; i++)
    {
      hashTable[i].locs = NULL;
    }

  //Loading Fasta File
  if (!initLoadingRefGenome(fileName))
    return;		
  initSavingIHashTable(indexName);
	
  fprintf(stderr, "Generating Index from %s", fileName);
  fflush(stderr);

  char *prev = getMem (CONTIG_NAME_SIZE);
  prev[0]='\0';

  do
    {
      flag = 	 loadRefGenome (&refGen, &refGenName, &refGenOff);	

      if ( strcmp(prev, refGenName) != 0)
	{
	  fprintf(stderr, "\n - %s ", refGenName);
	  fflush(stderr);
	  sprintf(prev, "%s", refGenName);
	}
      else
	{
	  fprintf(stderr, ".");
	  fflush(stderr);
	}
		
      l = strlen(refGen) - WINDOW_SIZE;

      for (i=0; i < l; i++)
	{
	  hv = hashVal(refGen + i);
	  if (hv != -1)
	    {
	      hashTableSize += addIHashTableLocation (hashTable, hv, i+1);
	    }
	}

      saveIHashTable(hashTable, hashTableSize, hashTableMaxSize, refGen, refGenName, refGenOff);
      freeIHashTableContent(hashTable, hashTableMaxSize);
      hashTableSize = 0;
    } while (flag);

  freeMem(prev, CONTIG_NAME_SIZE);
  freeMem(hashTable, sizeof(IHashTable)*hashTableMaxSize);

  finalizeLoadingRefGenome();
  finalizeSavingIHashTable();

  fprintf(stderr, "\nDONE in %0.2fs!\n", (getTime()-startTime));
}

/**********************************************/
void finalizeLoadingIHashTable()
{
  freeIHashTableContent(_ih_hashTable, _ih_maxHashTableSize);
  freeMem(_ih_hashTable, sizeof(IHashTable)* _ih_maxHashTableSize);
  freeMem(_ih_refGen, strlen(_ih_refGen)+1) ;
  freeMem(_ih_refGenName, strlen(_ih_refGenName)+1);
  fclose(_ih_fp);
}

/**********************************************/
int  loadIHashTable(double *loadTime)
{
  double startTime = getTime();
  unsigned char extraInfo = 0;
  short len;
  unsigned int refGenLength;
  unsigned int hashTableSize;
  unsigned int tmpSize;
  int tmp;
  int i=0,j=0;

  if ( fread(&extraInfo, sizeof(extraInfo), 1, _ih_fp) != sizeof(extraInfo) )
    return 0;

  freeIHashTableContent(_ih_hashTable, _ih_maxHashTableSize);
  freeMem(_ih_refGen, strlen(_ih_refGen)+1) ;
  freeMem(_ih_refGenName, strlen(_ih_refGenName)+1);

  // Reading Chr Name
  tmp = fread(&len, sizeof(len), 1, _ih_fp);
  if (tmp == 0){
    fprintf(stderr, "Read error while loading hash table.\n");
    exit(0);
  }

  _ih_refGenName = getMem(sizeof(char)* (len+1));
  tmp = fread(_ih_refGenName, sizeof(char), len, _ih_fp);
  _ih_refGenName [len] ='\0';

  tmp = fread(&_ih_refGenOff, sizeof (_ih_refGenOff), 1, _ih_fp);

  // Reading Size and Content of Ref Genome
  tmp = fread(&refGenLength, sizeof(refGenLength), 1, _ih_fp);
  _ih_refGen = getMem(sizeof(char)*(refGenLength+1));
  tmp = fread(_ih_refGen, sizeof(char), refGenLength, _ih_fp);
  _ih_refGen[refGenLength]='\0';

  //Reading Hashtable Size and Content
  tmp = fread(&hashTableSize, sizeof(hashTableSize), 1, _ih_fp);
	
  unsigned int hv;
  unsigned char cnt=0;
  //unsigned int  index = 1;
  //unsigned int  tmpValue = 0;

  for (i=0; i<hashTableSize; i++)
    {
      //index = 1;
      tmp = fread(&hv, sizeof(hv), 1, _ih_fp);
      tmp = fread(&cnt, sizeof(cnt), 1, _ih_fp);

      if (cnt>0)
	{
	  tmpSize = cnt;
	}
      else
	{
	  tmp = fread(&tmpSize, sizeof(tmpSize), 1, _ih_fp);
	}
		
      _ih_hashTable[hv].locs = getMem( sizeof(unsigned int)* (tmpSize+1) );
      _ih_hashTable[hv].locs[0] = tmpSize;

		
      for (j=1; j<=tmpSize; j++) 
	tmp = fread(&(_ih_hashTable[hv].locs[j]), sizeof(_ih_hashTable[hv].locs[j]), 1, _ih_fp);
		
      /*

	if(tmpSize != 0)
	{
	tmp = fread(&(_ih_hashTable[hv].locs[index]), sizeof(_ih_hashTable[hv].locs[index]), 1, _ih_fp);
	index++;
	}

	for (j=2; j<=tmpSize; j++)	
	{
	tmp = fread(&tmpValue, sizeof(_ih_hashTable[hv].locs[index-1]), 1, _ih_fp);
	if(abs(tmpValue-_ih_hashTable[hv].locs[index-1]) >  errThreshould)
	{
	_ih_hashTable[hv].locs[index] = tmpValue;
	index++;
	}
	}
	_ih_hashTable[hv].locs[0] = index-1;
		
      */
    }

  *loadTime = getTime()-startTime;
  return 1;
}
/**********************************************/
unsigned int *getIHashTableCandidates(int hv)
{
  if ( hv != -1 )
    return _ih_hashTable[hv].locs;
  else 
    return NULL;
}
/**********************************************/
/**********************************************/
/**********************************************/
void configHashTable()
{
  if (WINDOW_SIZE <= 15)
    {
      generateHashTable = &generateIHashTable;
      loadHashTable = &loadIHashTable;
      finalizeLoadingHashTable = &finalizeLoadingIHashTable;
      getCandidates = &getIHashTableCandidates;
    }
}
/**********************************************/
int initLoadingHashTable(char *fileName)
{
  int i;	
  unsigned char bsIndex;
  int tmp; 

  _ih_fp = fileOpen(fileName, "r");	

  if (_ih_fp == NULL)
    return 0;

  tmp = fread(&bsIndex, sizeof(bsIndex), 1, _ih_fp);

  if (tmp == 0){
    fprintf(stderr, "Read error while initializing hash table loading.\n");
    exit(0);
  }

  if (bsIndex)
    {
      fprintf(stderr, "Error: Wrong Type of Index indicated");
      return 0;
    }
	
  tmp = fread(&WINDOW_SIZE, sizeof(WINDOW_SIZE), 1, _ih_fp);

  configHashTable();

  if (_ih_maxHashTableSize != pow(4, WINDOW_SIZE))
    {

      if (_ih_hashTable != NULL)
	{
	  freeIHashTableContent(_ih_hashTable, _ih_maxHashTableSize);
	  freeMem(_ih_hashTable, sizeof(IHashTable)* _ih_maxHashTableSize);
	  freeMem(_ih_refGen, strlen(_ih_refGen)+1) ;
	  freeMem(_ih_refGenName, strlen(_ih_refGenName)+1);
	}

      _ih_maxHashTableSize = pow(4, WINDOW_SIZE);

      _ih_hashTable = getMem (sizeof(IHashTable) * _ih_maxHashTableSize);
      for (i=0; i<_ih_maxHashTableSize; i++)
	_ih_hashTable[i].locs = NULL;
      _ih_refGen = getMem(1);
      _ih_refGen[0]='\0';
      _ih_refGenName = getMem(1);
      _ih_refGenName[0] = '\0';
    }

  return 1;
}
/**********************************************/
char *getRefGenome()
{
  return _ih_refGen;

}
/**********************************************/
char *getRefGenomeName()
{
  return _ih_refGenName;

}
/**********************************************/
int getRefGenomeOffset()
{
  return _ih_refGenOff;

}
/**********************************************/
HashTable *getHashTable()
{
  return NULL;
}
