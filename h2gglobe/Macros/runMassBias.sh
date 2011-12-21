#!/bin/csh

rm BkgBias.root

root -b -l -q 'massBias.C('115',true)'
setenv M 120
while ( $M <= 150 )
  if ( $M != 145 ) root -b -l -q 'massBias.C('${M}')'
  setenv M `expr $M + 5`
end
