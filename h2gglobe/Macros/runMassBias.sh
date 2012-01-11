#!/bin/csh

rm BkgBias.root

setenv M 115
while ( $M <= 150 )
  if ( $M != 145 ) root -b -l -q 'massBias.C('${M}')'
  setenv M `expr $M + 5`
end
