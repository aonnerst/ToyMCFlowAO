import sys
sys.path.append("JPyPlotRatio");

import numpy as np
import ROOT
import pickle
import JPyPlotRatio
import matplotlib.ticker as plticker

#-----Read in files------------------
data = "graphs":ROOT.TFile("out_GraphsForRatio.root","read");

#-----Plot---------------------------
NC       = 3;
NMethod  = 3;
gr_Names =np.array(["SP","TP","EP"]);
plots=np.empty(NC+1,dtype=int);
centBins = np.array([0,5,20,30,50,60]);
rlimits  = [(-1.1,5.1),(-1.1,5.1)];

plot = JPyPlotRatio.JPyPlotRatio(panels=(1,3),disableRatio=[],
	legendLoc=(0.3,0.27),
	tickLabelSize=11,
	ratioBounds=rlimits,# for nrow
	panelLabel={i:"{}-{}%".format(centBins[i],centBins[i+1]) for i in range(0,6)},panelLabelLoc=(0.37,0.15),panelLabelAlign="left",ylabel="$v_n$");
#plot.EnableLatex(True);

#-----Get the tGraphs from files-----
for ic in NC:
	for im in NMethod:
		x,y,xerr,yerr =JPyPlotRatio.TGraphErrorsToNumpy(data.Get("gr_pv{}_{}".format(ic,gr_Names[im])));

for ic in NC:
	x1,y1,xerr1,yerr1 =JPyPlotRatio.TGraphErrorsToNumpy(data.Get("gr_vnin_{}".format(ic)));

# one centrality per pad
for ic in NC:
	plots[ic]=plot.Add(ic,(x1,y1,xerr1,yerr1),nlabel="hej");
	plots[NC+1]=plot.Add(ic,(x1,y1,xerr1,yerr1), nlabel="input");
	plot.Ratio()
