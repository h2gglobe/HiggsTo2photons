#!/bin/csh

echo $2

setenv GIFSIZE 1400
if (${2} == "mass") setenv GIFSIZE 1000
if (${2} == "bdtOutBin_grad_biascorrected") setenv GIFSIZE 1000
if (${2} == "bdtOutBin_ada_biascorrected") setenv GIFSIZE 1000

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
Mass hypothesis:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
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
<td>Final BDT with mass output:</td>
<td><a href="bdtOutBin_grad.html">Gradient boost</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="bdtOutBin_ada.html">Adaptive boost</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="bdtOutBin_grad_biascorrected.html">Gradient boost, bias corrected</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="bdtOutBin_ada_biascorrected.html">Adaptive boost, bias corrected</a>&nbsp;&nbsp;&nbsp;</td>
</tr>
<tr>
<td>Final BDT with mass input:</td>
<td><a href="bdtoutput.html">MIT BDT output</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="deltaMOverMH.html">deltaM/MH</a>&nbsp;&nbsp;&nbsp;</td>
</tr>
<tr>
<td>BDT input variables:</td>
<td><a href="pho1_phoidMva.html">lead photon ID MVA output</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho2_phoidMva.html">sublead photon ID MVA output</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="sigmaMOverM.html">sigmaM/Mgg (right vertex)</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="sigmaMOverM_wrongVtx.html">sigmaM/Mgg (wrong vertex)</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="vtxProb.html">vertex probability</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="deltaPhi.html">deltaPhi</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="cosDeltaPhi.html">cosDeltaPhi</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="ptOverM.html">diphoton pt/Mgg</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho1_ptOverM.html">lead pt/Mgg</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho2_ptOverM.html">sublead pt/Mgg</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho1_eta.html">pho1_eta</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho2_eta.html">pho2_eta</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho_minr9.html">pho_minr9</a>&nbsp;&nbsp;&nbsp;</td>
</tr>
<tr>
<td>Other variables:</td>
<td><a href="maxeta.html">maxeta</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="deltaEta.html">deltaEta</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="eta.html">diphoton eta</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="ptOverMH.html">diphoton pt/MH</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho1_ptOverMH.html">lead pt/MH</a>&nbsp;&nbsp;&nbsp;</td>
<td><a href="pho2_ptOverMH.html">sublead pt/MH</a>&nbsp;&nbsp;&nbsp;</td>
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
