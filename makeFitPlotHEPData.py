import sys, string, copy
from hepdata_lib import Submission
from hepdata_lib import Uncertainty as U
from hepdata_lib import Variable as V
from hepdata_lib import Table as T
from hepdata_lib import Submission as S
from hepdata_lib import RootFileReader as R

# Given a mystery array of uncertainties from hepdata
# Make sure the format is correct for symmetric or not
def makeUncArray(u):

    length = len(u)
    newUnc = []
    isSymm = None

    # Empty array, get the hell out of here
    if length == 0:
        return u, True

    else:
        for unc in u:
            # If the array is filled with tuples
            # still need to determine if symm or not
            if type(unc) == tuple:
                # If tuples are symm, then make simple array
                if (abs(unc[0]) - abs(unc[1])) / abs(unc[0]) < 1e-6:
                    newUnc.append(abs(unc[0]))
                    isSymm = True
                else: 
                    newUnc.append((unc[0], unc[1]))
                    isSymm = False
            else:
                newUnc.append(unc)
                isSymm = True

    return newUnc, isSymm

def makeFitPlotHEPData(sub):
    years = ["Y16", "Y17", "Y18A", "Y18B"]
    procs = ["Fit", "sigRefHist1", "sigRefHist2", "Nobs"]
    nnBins = ["D1", "D2", "D3", "D4"]
    tmap = {"Y16"  : None,
            "Y17"  : None,
            "Y18A" : None,
            "Y18B" : None
    }
    
    # Map histo name to name for presentation
    names = {"Fit" : "Bkg Fit", "sigRefHist1" : "RPV $m_{\\tilde{t}}$ = 450 GeV", "sigRefHist2" : "SYY $m_{\\tilde{t}}$ = 850 GeV", "Nobs" : "N observed"}
    lumiMap = {"Y16" : 35.9, "Y17" : 41.5, "Y18A" : 21.1, "Y18B" : 38.7}
    files = {"Y16" : R("inputs/Fit_RPV450Combo16b.root"), "Y18A" : R("inputs/Fit_RPV450Combo18preb.root"),
             "Y17" : R("inputs/Fit_RPV450Combo17b.root"), "Y18B" : R("inputs/Fit_RPV450Combo18postb.root")
    }
    
    # Define one instance of the independent variable, here Njets
    xvar = V("$N_{jets}$-$S_{\\textrm{NN}}$ bin", is_independent=True, is_binned=False, units="")
    xvar.values = list(range(1,25)) 

    # Create a unique table for each SNN, Njets bin
    letters = list(string.ascii_lowercase); count = 0
    for year, t in sorted(tmap.items()):
        tab = T("Figure 4%s"%(letters[count]))
        tab.description = "Fitted background prediction and observed data counts for %s as functions of $N_{\\text{jets}}$ in each of the four $S_{\\textrm{NN}}$ bins. The signal" \
                          " distributions normalized to the predicted cross section for the RPV model with $m_{\\tilde{t}}$ = 450 GeV" \
                          " and the stealth SYY model with $m_{\\tilde{t}}$ = 850 GeV are shown for comparison."%(year.replace("Y", "20"))
        tab.location = "Data from Figure 4, located on page 11"
        tmap[year] = copy.deepcopy(tab)
        tmap[year].add_variable(xvar)
        count += 1
    
    for year in sorted(years):
        for proc in sorted(procs):
            varArr = []; uncArr = []

            # Make a variable object to hold one of the histos or graphs
            # and give it the y vals from the object
            v = V("%s"%(names[proc]), is_independent=False, is_binned=False, units="")
            v.add_qualifier("SQRT(S)", 13, "TeV")
            v.add_qualifier("LUMINOSITY", lumiMap[year], "fb$^{-1}$")

            for nnBin in sorted(nnBins): 
                h = None
    
                # For the fit and data TGraphs read them as graphs
                # Otherwise, read everything as a simple 1D histo
                if "Fit" in proc or "Nobs" in proc: 
                    h = files[year].read_graph("%s_%s_%s"%(proc, year.replace("A","pre").replace("B","post"), nnBin))
                else:
                    h = files[year].read_hist_1d("%s_%s_%s"%(proc, year.replace("A","pre").replace("B","post"), nnBin))

                varArr += h["y"]
    
                # For the fit and data graphs, there is an associated uncertainty
                # So add that uncertainty to the Fit and Nobs variables
                if proc == "Fit" or proc == "Nobs":
                    uncArr += h["dy"]

            v.values = varArr

            if proc == "Fit" or proc == "Nobs":
                uncVals, isSymm = makeUncArray(uncArr)
                unc = U("unc.", is_symmetric=isSymm)
                unc.values = uncVals
                v.add_uncertainty(unc)

            # Add each sig1, sig2, fit, data variable to the corresponding table
            tmap[year].add_variable(v)

    for year, t in sorted(tmap.items()):
        sub.add_table(t)
    
if __name__ == "__main__":
    sub = Submission()
    makeFitPlotHEPData(sub)
    sub.create_files("output")
