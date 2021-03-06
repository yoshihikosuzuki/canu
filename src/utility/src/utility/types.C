
/******************************************************************************
 *
 *  This file is part of meryl-utility, a collection of miscellaneous code
 *  used by Meryl, Canu and others.
 *
 *  This software is based on:
 *    'Canu' v2.0              (https://github.com/marbl/canu)
 *  which is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "types.H"

//  In hex, a 128-bit integer needs 32 digits.
//  In dec, a 128-bit integer needs 39 digits.  (it is 340,282,366,920,938,463,463,374,607,431,768,211,456)
//
//  We'll just allocate 64 digits and be done with it.  Until we want to
//  print that 128-bit integer as binary.  (that we overallocate space makes
//  conversion to decimal a little bit easier)

char   dec[10]     = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' };
char   hex[16]     = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f' };
char   str[16][64] = { 0 };
uint32 pos         =   0;



template<typename uintType>
char const *
toHex(uintType v) {
  char   *ret = str[pos++];
  uint32  w = sizeof(uintType) * 2;
  uint32  p = w;
  uint32  s = 0;

  if (pos >= 16)
    pos = 0;

  while (p > 0) {
    p -= 1;
    ret[p] = hex[ (v >> s) & 0xf ];
    s += 4;
  }

  ret[w] = 0;

  return(ret);
}

template char const *toHex<uint128>(uint128 v);
template char const *toHex<uint64> (uint64  v);
template char const *toHex<uint32> (uint32  v);
template char const *toHex<uint16> (uint16  v);



template<typename uintType>
char const *
toDec(uintType v) {
  char   *ret = str[pos++];
  uint32  p   = 64;
  uint32  x   = 0;

  if (pos >= 16)
    pos = 0;

  if (v == 0) {
    ret[0] = '0';
    ret[1] =  0;
  }

  else {
    while (v > 0) {             //  Write the number, low-order
      p -= 1;                   //  digits first, to the end
      ret[p] = dec[ v % 10 ];   //  of the string.
      v /= 10;
    }

    for (x=0; p<64; x++, p++)   //  Shift the string so it
      ret[x] = ret[p];          //  starts at the origin.

    ret[x] = 0;
  }

  return(ret);
}

template char const *toDec<uint128>(uint128 v);
template char const *toDec<uint64> (uint64  v);
template char const *toDec<uint32> (uint32  v);
template char const *toDec<uint16> (uint16  v);
