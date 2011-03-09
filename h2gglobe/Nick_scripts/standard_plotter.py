import ROOT
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetLineWidth(2)
import sys,os

r = ROOT.TRandom3()
current_lumi    = 36.1
data_name 	= 'hist/Data.root'

print "Higgs Analysis - Standard Plotter"
print "Plotting all histograms available in ./hist"

file_names	= {'hist/DYEEM20_hist.root'	:[0,1636.621	   	]
		  ,'hist/DYEEM1020_hist.root'	:[0,726.965	   	]
		  ,'hist/QCD40_hist.root'	:[1,544.862	   	]
		  ,'hist/GJet20_hist.root'	:[1,2395.580	   	]
		  ,'hist/Born_10_hist.root'	:[2,2211.781	   	]
		  ,'hist/Born_25_hist.root'	:[2,24025.257	   	]
		  ,'hist/Born_250_hist.root'	:[2,67685208.127   	]
		  ,'hist/Box_10_hist.root'	:[3,2227.736	   	]
		  ,'hist/Box_25_hist.root'	:[3,62871.867	   	]
		  ,'hist/Box_250_hist.root'	:[3,3795528846.154 	]
		  ,'hist/WZTTH115_hist.root'	:[4,0.001*40475893.478   ]
		  ,'hist/VBF115_hist.root'	:[4,0.001*38717590.830   ]
		  ,'hist/GGH115_hist.root'	:[4,0.001*2848260.736    ]
		  ,'hist/WZTTH120_hist.root'	:[5,0.001*43907180.221   ]
		  ,'hist/VBF120_hist.root'	:[5,0.001*38470186.499   ]
		  ,'hist/GGH120_hist.root'	:[5,0.001*2939587.092    ]
		  ,'hist/WZTTH130_hist.root'	:[6,0.001*56914996.108   ]
		  ,'hist/VBF130_hist.root'	:[6,0.001*41722136.164   ]
		  ,'hist/GGH130_hist.root'	:[6,0.001*3427365.075    ]
		  }

print "Using files/type/weight :"
for k in file_names.keys(): print k, file_names[k] 

# should have one color and title per Channel
titles_colors = {'0':['Z/#gamma* #rightarrow ee',ROOT.kRed+1	]
		,'1':['QCD & #gamma + Jet',ROOT.kOrange-3	]
	 	,'2':['Born',ROOT.kGreen+3			]
		,'3':['Box',ROOT.kMagenta+3			]
		,'4':['Higgs M-115',ROOT.kBlue-2		]
		,'5':['Higgs M-120',ROOT.kTeal-2		]
		,'6':['Higgs M-130',ROOT.kCyan-2		]
	 	}

signal_index = [4,5,6]
#Build linked files lists:
for f in file_names: 
  if not (os.path.isfile(f)):
   sys.exit("No file found named %s"%f)

linked_files = [[filter(lambda x: (file_names[x])[0]==i,file_names),
		 i in signal_index] \
		for i in range(len(file_names))]

file_list = [[[ROOT.TFile(f) for f in F[0]],F[1]] for F in linked_files]
file_list = filter(lambda x: len(x[0]) > 0, file_list)

if os.path.isfile(data_name):  	
  data_file = ROOT.TFile(data_name)
else:
  sys.exit('No file found names %s'%data_name)

if len(titles_colors) != len(file_list): 
  sys.exit("Incorrect number of colors to channels!") 

keys 	  = [[f.GetListOfKeys() for f in F[0]] 
	      for F in file_list]
data_k	  = data_file.GetListOfKeys()

tmp=[0 for kk in file_list]

for i in range(len(keys[0][0])):


	leg=ROOT.TLegend(0.55,0.6,0.85,0.85)
	leg.SetFillColor(0)
	leg.SetBorderSize(0)
	
	stack = ROOT.THStack()	
	c = ROOT.TCanvas()
	c.SetLogy()
	
	to_be_stacked = []

	for l in range(len(file_list)):
	  for j in range(len(file_list[l])):

		f = file_list[l][0][j]
		h = keys[l][j][i].ReadObj()
		w = (file_names[f.GetName()])[1]

		h.Scale(current_lumi/w)
		if j == 0: 
			tmp[l] = h.Clone()
		else:  tmp[l].Add(h)

	  tmp[l].SetFillColor(titles_colors[str(l)][1])
	  leg.AddEntry(tmp[l],titles_colors[str(l)][0],'F')
          #tmp[l].SetLineWidth(1.5)

	  to_be_stacked.append(tmp[l])
	

	data_hist = data_k[i].ReadObj()
	n_data	= data_hist.Integral()

	stacked_sum = 0
	for l,h in enumerate(to_be_stacked):
	  if not file_list[l][1]:
	    stacked_sum+=h.Integral()

	for h in to_be_stacked:
	  h.Scale(n_data/stacked_sum)
	  stack.Add(h)  

	data_hist.SetMarkerStyle(21)
	data_hist.SetMarkerSize(0.7)
	data_hist.Sumw2()

	leg.AddEntry(data_hist, 'Data','PLE')
	
	data_hist.Draw()
	stack.Draw('same')
	data_hist.Draw('same')
	leg.Draw()

	c.SetGrid(True)	
	c.SaveAs('plots/'+h.GetName()+'.pdf')

#EOF



  

