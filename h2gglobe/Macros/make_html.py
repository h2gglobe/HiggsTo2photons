import os
import fnmatch
import sys

def printLinks(outFile,webpath,folder,syst):
  outFile.write('<center>\n <p> \n <font size="5">Signal Interpolation Diagnostics</font>\n </p> \n')
  outFile.write('<script language="Javascript"> \n document.write("Results from interpolation run at: " + document.lastModified + ""); \n </SCRIPT> <br> \n ')
  outFile.write('Run on file: '+folder+'<br>\n')
  outFile.write('Output file with interpolated histograms: '+folder+'_interpolated.root <br>\n')
  outFile.write('<a href=\"'+webpath+'plots/BDTInterpolationDiagnostics.txt\">Diagnostics File</a><br>\n')
  outFile.write('<a href=\"'+webpath+'plots/grad/david/plots.html\">Signal, background and data plots</a><br>\n')
  outFile.write('<a href=\"'+webpath+'plots/ada/fracs/plots.html\">Up, down and interpolated signal plots</a><br>\n')
  outFile.write('<a href=\"'+webpath+'plots/ada/systs/plots.html\">Fractional systematic up and down templates from mass histogram </a><br>\n')
  outFile.write('<a href=\"'+webpath+'plots/syst.html\">Comparison of interpolated to true</a><br>\n')
  if syst==0:
    outFile.write('BDT type: \n')
    outFile.write('<a href=\"../../ada/david/plots.html\">ada</a>\n')
    outFile.write('<a href=\"../../grad/david/plots.html\">grad</a><br> </center>\n')

user = sys.argv[1]
inPath = sys.argv[2]
folder = sys.argv[3]
webpath='https://www.hep.ph.ic.ac.uk/~'+user+'/h2g/MVA/SigInt/Diagnostics/'+inPath+'/'
genFile=open('plots/plots.html','w')
genFile.write('<center>\n <p> \n <font size="5">Signal Interpolation Diagnostics</font>\n </p> \n')
genFile.write('<script language="Javascript"> \n document.write("Run at: " + document.lastModified + ""); \n </SCRIPT> <br> \n')
genFile.write('Run on file '+folder+'<br>')
genFile.write('Output file with interpolated histograms: '+folder+'_interpolated.root <br>\n')
genFile.write('<a href=\"BDTInterpolationDiagnostics.txt\">Diagnostics File</a><br>\n')
genFile.write('<a href=\"ada/david/plots.html\">Signal, background and data plots</a><br>\n')
genFile.write('<a href=\"ada/fracs/plots.html\">Up, down and interpolated signal plots</a><br>\n')
genFile.write('<a href=\"ada/systs/plots.html\">Fractional systematic up and down templates from mass histogram </a><br>\n')
genFile.write('<a href=\"syst.html\">Comparison of interpolated to true</a> </center>\n')

htmlFiles = []

for root, dirs, files in os.walk('plots'):
  for filename in fnmatch.filter(files,'*.png'):
    if root not in htmlFiles:
      htmlFiles.append(root)
    #print root
    #print os.path.join(root,filename)

systFile = open('plots/syst.html','w')
printLinks(systFile,webpath,folder,1)

for htmlName in htmlFiles:
  htmlFile = open(htmlName+'/plots.html','w')
  printLinks(htmlFile,webpath,folder,0)
  for root, dirnames, filenames in os.walk(htmlName):
    filenames.sort()
    for files in filenames:
      link=webpath+str(htmlName)+'/'+str(files)
      if 'syst' in files and 'Fracs' not in files:
        systFile.write('<a href='+link+'><img height=\"400\" src=\"'+link+'\"></a>\n')
      else:
        htmlFile.write('<a href='+link+'><img height=\"400\" src=\"'+link+'\"></a>\n')
