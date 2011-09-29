#!/bin/csh

if ($1 == "") then
  echo "output directory name required as argument"
  exit
endif

cmsenv
setenv M 115
while ( $M <= 150 )
  setenv N `expr $M - 105`
  setenv N `expr $N / 5`
  setenv N `expr $N - 1`
  if ( $M == 150 ) setenv N `expr $N - 1`
  sed s/backgroundModelPlots/backgroundModelPlots\_"${M}"/g backgroundModelPlots.C > temp1
  sed s/mlow\_cat2/mlow\_cat"${N}"/g temp1 > temp2
  sed s/msig\_cat2/msig\_cat"${N}"/g temp2 > temp3
  sed s/mhigh\_cat2/mhigh\_cat"${N}"/g temp3 > temp4
  sed s/gg\_120/gg\_"${M}"/g temp4 > temp5
  sed s/120\_cat0/"${M}"\_cat0/g temp5 > temp6
  sed s/in\=120/in\="${M}"/g temp6 > backgroundModelPlots_${M}.C
  rm temp1 temp2 temp3 temp4 temp5 temp6
  if ( $M != 145 ) root -b -l -q 'backgroundModelPlots_${M}.C(1,"'${1}'")'
  rm backgroundModelPlots_${M}.C
  setenv M `expr $M + 5`
end

setenv OUTDIR $1
setenv M 115
while ( $M <= 150 )
  setenv N `expr $M - 105`
  setenv N `expr $N / 5`
  setenv N `expr $N - 1`
  if ( $M == 150 ) setenv N `expr $N - 1`
  if ( $M != 120 ) cp -p ${OUTDIR}/120/gifs/sigmaM_cat.gif ${OUTDIR}/${M}/gifs/
  setenv M `expr $M + 5`
end
