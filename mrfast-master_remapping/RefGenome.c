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
#include "Common.h"
#include "RefGenome.h"

FILE *_rg_fp;
char *_rg_gen;
char *_rg_name;
int _rg_offset;
int _rg_contGen; 

/**********************************************/
int initLoadingRefGenome(char *fileName)
{
  char ch;
  _rg_fp = fileOpen (fileName, "r");
  if (fscanf(_rg_fp, "%c", &ch))
    {
      if (ch == '>')
	{
	  _rg_contGen = 0;
	  _rg_offset = 0;
	  _rg_gen = getMem(CONTIG_MAX_SIZE);
	  _rg_name = getMem(CONTIG_NAME_SIZE);
	  return 1;
	}
    }
  return 0;
}
/**********************************************/
void finalizeLoadingRefGenome()
{
  freeMem(_rg_gen, CONTIG_MAX_SIZE);
  freeMem(_rg_name, CONTIG_NAME_SIZE); 
  fclose(_rg_fp);
}
/**********************************************/
int loadRefGenome(char **refGen, char **refGenName, int *refGenOff)
{
  char ch;
  int i;
  int returnVal = 0;
  int actualSize=0;
  int size;
  char *tmp;
	
  // New Contig 
  if (!_rg_contGen)
    {
      size = 0;
      tmp = fgets(_rg_name, SEQ_MAX_LENGTH, _rg_fp);

      if (tmp==NULL)
	fprintf(stderr, "Error reading the contig.\n");
      
      int k;
      for (k=0; k<strlen(_rg_name);k++)
	{
	  if (_rg_name[k] == ' ')
	    {
	      _rg_name[k]='\0';
	      break;
	    }
	}
    }
  else
    {
      size=strlen(_rg_gen);
      for( i = 0 ; i < CONTIG_OVERLAP ; i++ )
	{
			
	  _rg_gen[i] = _rg_gen[size-CONTIG_OVERLAP+i];
	  if (_rg_gen[i] != 'N')
	    actualSize++;
	}
      size = CONTIG_OVERLAP;
    }
  while( fscanf(_rg_fp, "%c", &ch) > 0 )
    {
      if (ch == '>')
	{
	  _rg_contGen = 0;
	  returnVal = 1;
	  break;
	}
      else if (!isspace(ch))
	{
	  ch = toupper(ch);
	  _rg_gen[size++] = ch;
	  if (ch != 'N')
	    {
	      actualSize++;
	    }
	  if (actualSize == CONTIG_SIZE || size == CONTIG_MAX_SIZE)
	    {
	      _rg_contGen = 1;
	      returnVal=1;
	      break;
	    }
	}

    }

  _rg_gen[size] = '\0';

  for (i=strlen(_rg_name)-1; i >= 0; i--)
    if (!isspace(_rg_name[i]))
      break;
  _rg_name[i+1] = '\0';

  *refGenOff = _rg_offset;
  *refGenName = _rg_name;
  *refGen = _rg_gen;
  
  if (_rg_contGen == 1)
    {
      _rg_offset += size-CONTIG_OVERLAP;
    }
  else
    {
      _rg_offset = 0;
    }
	
	
  return returnVal;
}
/**********************************************/
