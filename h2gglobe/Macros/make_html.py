import os
import fnmatch

def printLinks(outFile,webpath):
  outFile.write('<p> \n <center>\n <font size="5">Signal Interpolation Diagnostics</font>\n </center> \n </p> \n')
  outFile.write('<a href=\"../../BDTInterpolationDiagnostics.txt\">Diagnostics File</a><br>\n')
  outFile.write('<a href=\"../david/plots.html\">Signal, background and data plots</a><br>\n')
  outFile.write('<a href=\"../fracs/plots.html\">Up, down and interpolated signal plots</a><br>\n')
  outFile.write('BDT type: \n')
  outFile.write('<a href=\"../../ada/david/plots.html\">ada</a>\n')
  outFile.write('<a href=\"../../grad/david/plots.html\">grad</a><br>\n')

webpath='https://www.hep.ph.ic.ac.uk/~mk1009/h2g/MVA/SigInt/Diagnostics/'
genFile=open('plots/plots.html','w')
genFile.write('<p> \n <center>\n <font size="5">Signal Interpolation Diagnostics</font>\n </center> \n </p> \n')
genFile.write('<a href=\"BDTInterpolationDiagnostics.txt\">Diagnostics File</a><br>\n')
genFile.write('<a href=\"ada/david/plots.html\">Signal, background and data plots</a><br>\n')
genFile.write('<a href=\"ada/fracs/plots.html\">Up, down and interpolated signal plots</a><br>\n')

htmlFiles = []

for root, dirs, files in os.walk('plots'):
  for filename in fnmatch.filter(files,'*.png'):
    if root not in htmlFiles:
      htmlFiles.append(root)
    #print root
    #print os.path.join(root,filename)

for htmlName in htmlFiles:
  htmlFile = open(htmlName+'/plots.html','w')
  printLinks(htmlFile,webpath)
  for root, dirnames, filenames in os.walk(htmlName):
    filenames.sort()
    for files in filenames:
      link=webpath+str(htmlName)+'/'+str(files)
      #print link
      htmlFile.write('<a href='+link+'><img height=\"400\" src=\"'+link+'\"></a>\n')
