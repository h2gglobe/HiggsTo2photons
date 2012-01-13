#!/bin/csh

if ($1 == "") then
  echo "output directory name required as argument"
  exit
endif

#cmsenv
setenv M 115
while ( $M <= 150 )
  if ( $M != 145 ) root -b -l -q 'backgroundModelPlots.C('${M}',1,"'${1}'")'
  setenv M `expr $M + 5`
end
