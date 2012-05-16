#include "HiggsAnalysis/HiggsTo2photons/interface/Tools.h"

int fexist(char *filename) {
  struct stat buffer ;
  if (stat( filename, &buffer ))
    return 1;
  return 0;
}
