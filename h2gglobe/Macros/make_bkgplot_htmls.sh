#!/bin/csh

if ($1 == "") then
  echo "output directory name required as argument"
  exit
endif

./make_bkgplot_html.sh $1 mass
./make_bkgplot_html.sh $1 bdtOut_ada
./make_bkgplot_html.sh $1 bdtOut_grad
./make_bkgplot_html.sh $1 bdtOutBin_ada
./make_bkgplot_html.sh $1 bdtOutBin_grad
./make_bkgplot_html.sh $1 sigmaM_cat
./make_bkgplot_html.sh $1 sigmaM
./make_bkgplot_html.sh $1 sigmaMOverM
./make_bkgplot_html.sh $1 deltaMOverMH
./make_bkgplot_html.sh $1 deltaMOverSigmaM
./make_bkgplot_html.sh $1 deltaEta
./make_bkgplot_html.sh $1 deltaPhi
./make_bkgplot_html.sh $1 eta
./make_bkgplot_html.sh $1 helicityAngle
./make_bkgplot_html.sh $1 ptOverM
./make_bkgplot_html.sh $1 maxeta
./make_bkgplot_html.sh $1 pho1_eta
./make_bkgplot_html.sh $1 pho2_eta
./make_bkgplot_html.sh $1 pho1_ptOverMH
./make_bkgplot_html.sh $1 pho2_ptOverMH
./make_bkgplot_html.sh $1 pho_minr9
