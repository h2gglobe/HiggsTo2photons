#!/bin/csh

echo $2

setenv GIFSIZE 1400
if (${2} == "mass") setenv GIFSIZE 1000

setenv OUTDIR ${1}
if (! -d ${OUTDIR} ) mkdir ${OUTDIR}
cd ${OUTDIR}

setenv M 115
while ( $M <= 150 )

    if (! -d $M ) mkdir $M
    cd $M
    if (! -d gifs ) mkdir gifs

    if (-e ${2}.html) rm -f ${2}.html

cat > ${2}.html <<@EOF
<html>
<body>
<big>
Sideband configuration:&nbsp;&nbsp;
<table border="0">
<tr>
<td><a href="../../7percent_1sideband/${M}/${2}.html">7percent_1sideband</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="../../7percent_1sideband_loose/${M}/${2}.html">7percent_1sideband_loose</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="../../7percent_1sideband_veryloose/${M}/${2}.html">7percent_1sideband_veryloose</a>&nbsp;&nbsp;&nbsp;</td>
</tr>
<tr>
<td><a href="../../2percent_3sidebands/${M}/${2}.html">2percent_3sidebands</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="../../2percent_3sidebands_loose/${M}/${2}.html">2percent_3sidebands_loose</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="../../2percent_3sidebands_veryloose/${M}/${2}.html">2percent_3sidebands_veryloose</a>&nbsp;&nbsp;&nbsp;</td>
</tr>
</table> 
<br>
Mass hypothesis:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="../115/${2}.html">115</a>&nbsp;&nbsp;&nbsp;
<a href="../120/${2}.html">120</a>&nbsp;&nbsp;&nbsp;
<a href="../125/${2}.html">125</a>&nbsp;&nbsp;&nbsp;
<a href="../130/${2}.html">130</a>&nbsp;&nbsp;&nbsp;
<a href="../135/${2}.html">135</a>&nbsp;&nbsp;&nbsp;
<a href="../140/${2}.html">140</a>&nbsp;&nbsp;&nbsp;
<a href="../150/${2}.html">150</a>&nbsp;&nbsp;&nbsp;
<br><br>
<table border="0">
<tr>
<td>BDT Output:</td>
<td><a href="bdtOutBin_grad_nominalbins.html">Gradient boost</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="bdtOutBin_ada_nominalbins.html">Adaptive boost</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="bdtOutBin_grad_sob_nominalbins.html">Gradient boost (S/B binning)</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="bdtOutBin_ada_sob_nominalbins.html">Adaptive boost (S/B binning)</a>&nbsp;&nbsp;&nbsp;</td>
</tr>
<tr>
<td>BDT Input variables:</td>
<td><a href="deltaMOverMH.html">deltaMOverMH</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="sigmaMOverM.html">sigmaMOverM</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="deltaPhi.html">deltaPhi</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho1_ptOverMH.html">pho1_ptOverMH</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho2_ptOverMH.html">pho2_ptOverMH</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="ptOverMH.html">ptOverMH</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="eta.html">eta</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="maxeta.html">maxeta</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho1_eta.html">pho1_eta</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho2_eta.html">pho2_eta</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho_minr9.html">pho_minr9</a>&nbsp;&nbsp;&nbsp;</td>
</tr>
<tr>
<td>BDT inputs NOT used:</td>
<td><a href="pho1_pt.html">pho1_pt</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho2_pt.html">pho2_pt</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="sigmaMOverMH.html">sigmaMOverMH</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="deltaMOverSigmaM.html">deltaMOverSigmaM</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="deltaEta.html">deltaEta</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="cosDeltaPhi.html">cosDeltaPhi</a>&nbsp;&nbsp;&nbsp;</td>
</tr>
<tr>
<td>Mass distribution:</td>
<td><a href="mass.html">Di-Photon mass</a>&nbsp;&nbsp;&nbsp;</td>
</tr>
</table> 
<a href="gifs/${2}.gif"><img src="gifs/${2}.gif" width="${GIFSIZE}"></a>
</body>
</html>
@EOF

    cd ..
    setenv M `expr $M + 5`
end

#<a href="sigmaM_cat.html">sigmaM in categories</a> (for mass hypothesis = 120)

#<td>BDT Output:</td>
#<td><a href="bdtOutBin_grad.html">Gradient boost</a>&nbsp;&nbsp;&nbsp;</td>
#<td><a href="bdtOutBin_ada.html">Adaptive boost</a>&nbsp;&nbsp;&nbsp;</td>
#<td><a href="bdtOutBin_grad_sob.html">Gradient boost (S/B binning)</a>&nbsp;&nbsp;&nbsp;</td>
#<td><a href="bdtOutBin_ada_sob.html">Adaptive boost (S/B binning)</a>&nbsp;&nbsp;&nbsp;</td>
#</tr>
#<tr>
