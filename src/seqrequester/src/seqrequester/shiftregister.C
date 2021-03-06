
/******************************************************************************
 *
 *  This file is part of seqrequester, a tool for summarizing, extracting,
 *  generating and modifying DNA sequences.
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

#include "seqrequester.H"
#include "bits.H"



void
shiftRegisterParameters::initialize(void) {
  uint32 srLen = strlen(sr);
  uint32 snLen = strlen(svmin);
  uint32 sxLen = strlen(svmax);
  bool   fail  = false;

  fail |=  (order == 0);

  fail |= ((srLen > 0) && (snLen > 0) && (srLen != snLen));
  fail |= ((srLen > 0) && (sxLen > 0) && (srLen != sxLen));
  fail |= ((snLen > 0) && (sxLen > 0) && (snLen != sxLen));

  fail |= ((srLen > 0) && (srLen != order));
  fail |= ((snLen > 0) && (snLen != order));
  fail |= ((sxLen > 0) && (snLen != order));

  if (fail == true) {
    fprintf(stderr, "ERROR: order %u\n", order);
    fprintf(stderr, "ERROR: sr    %s len %u\n", sr,    srLen);
    fprintf(stderr, "ERROR: svmin %s len %u\n", svmin, snLen);
    fprintf(stderr, "ERROR: svmax %s len %u\n", svmax, sxLen);
  }

  assert(fail == false);

  if (srLen > 0)   std::reverse(sr,    sr    + srLen);
  if (snLen > 0)   std::reverse(svmin, svmin + snLen);
  if (sxLen > 0)   std::reverse(svmax, svmax + sxLen);
}



uint64
shiftRegisterParameters::getEncodedSR(void) {
  uint64  r = 0llu;

  if (sr[0] == 0) {
    r = 1llu;
  }

  else {
    for (uint32 ii=0; ii<order; ii++) {
      r <<= 2;
      r  |= sr[order-1-ii] - '0';
    }
  }

  return(r);
}



uint64
shiftRegisterParameters::getCycleLen(void) {
  return(0llu);
}



uint64
shiftRegisterParameters::getCycleMax(void) {
  return(1llu << (2 * order));
}



uint64
shiftRegisterParameters::getEncodedSVmin(void) {
  uint64  r = 0llu;

  if (svmin[0] == 0) {
    r   = 1llu;
    r <<= 2 * order - 2;
  }

  else {
    for (uint32 ii=0; ii<order; ii++) {
      r <<= 2;
      r  |= svmin[order-1-ii] - '0';
    }
  }

  return(r);
}



uint64
shiftRegisterParameters::getEncodedSVmax(void) {
  uint64  r = 0llu;

  if (svmax[0] == 0) {
    r   = 1llu;
    r <<= 2 * order;
    r  -= 1;
  }

  else {
    for (uint32 ii=0; ii<order; ii++) {
      r <<= 2;
      r  |= svmax[order-1-ii] - '0';
    }
  }

  return(r);
}



uint64
shiftRegisterParameters::getEncodedSVmask(void) {
  return((1llu << (2 * order)) - 1);
}




void  searchShiftRegisterFast(shiftRegisterParameters &srPar);
void  searchShiftRegisterSlow(shiftRegisterParameters &srPar);
void  emitShiftRegisterFast(shiftRegisterParameters &srPar);

void  emitShiftRegisterSlow(shiftRegisterParameters &srPar) {
}



void
doShiftRegister(shiftRegisterParameters &srPar) {

  srPar.initialize();

  fprintf(stderr, "VERSION 7\n");

  if      ((srPar.search == true) && (srPar.fast == true))
    searchShiftRegisterFast(srPar);

  else if ((srPar.search == true) && (srPar.fast == false))
    searchShiftRegisterSlow(srPar);

  else if ((srPar.search == false) && (srPar.fast == true))
    emitShiftRegisterFast(srPar);

  else if ((srPar.search == true) && (srPar.fast == false))
    emitShiftRegisterSlow(srPar);
}

