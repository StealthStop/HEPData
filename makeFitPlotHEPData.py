import sys, string, copy
sys.path.insert(0, "/uscms_data/d3/jhiltb/susy/CMSSW_11_2_0/src/sandbox/HEPData/hepdata_lib")
from hepdata_lib import Submission
from hepdata_lib import Uncertainty as U
from hepdata_lib import Variable as V
from hepdata_lib import Table as T
from hepdata_lib import Submission as S
from hepdata_lib import RootFileReader as R

# If the Variable object is present in the Table
# just do nothing, don't add again!
def addUniqueVar(t, v):

    alreadyAdded = False
    for var in t.variables:

        if var.name == v.name:
            alreadyAdded = True
            break

    if not alreadyAdded:
        t.add_variable(v)

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

def makeFitPlotHEPData():
    # Dictionary to keep track of all the Variable objects
    vmap = {"Y16"  : {"Fit"         : {"D1" : None, "D2" : None, "D3" : None, "D4" : None},
                      "sigRefHist1" : {"D1" : None, "D2" : None, "D3" : None, "D4" : None},
                      "sigRefHist2" : {"D1" : None, "D2" : None, "D3" : None, "D4" : None},
                      "Nobs"        : {"D1" : None, "D2" : None, "D3" : None, "D4" : None}},
            "Y17"  : {"Fit"         : {"D1" : None, "D2" : None, "D3" : None, "D4" : None},
                      "sigRefHist1" : {"D1" : None, "D2" : None, "D3" : None, "D4" : None},
                      "sigRefHist2" : {"D1" : None, "D2" : None, "D3" : None, "D4" : None},
                      "Nobs"        : {"D1" : None, "D2" : None, "D3" : None, "D4" : None}},
            "Y18A" : {"Fit"         : {"D1" : None, "D2" : None, "D3" : None, "D4" : None},
                      "sigRefHist1" : {"D1" : None, "D2" : None, "D3" : None, "D4" : None},
                      "sigRefHist2" : {"D1" : None, "D2" : None, "D3" : None, "D4" : None},
                      "Nobs"        : {"D1" : None, "D2" : None, "D3" : None, "D4" : None}},
            "Y18B" : {"Fit"         : {"D1" : None, "D2" : None, "D3" : None, "D4" : None},
                      "sigRefHist1" : {"D1" : None, "D2" : None, "D3" : None, "D4" : None},
                      "sigRefHist2" : {"D1" : None, "D2" : None, "D3" : None, "D4" : None},
                      "Nobs"        : {"D1" : None, "D2" : None, "D3" : None, "D4" : None}}
    }
    
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
    
    sub = Submission()
    
    # Create a unique table for each SNN, Njets bin
    letters = list(string.ascii_lowercase); count = 0
    for year, t in sorted(tmap.items()):
        tab = T("Figure 4%s"%(letters[count]))
        tab.description = "Fitted background prediction and observed data counts for %s as functions of $N_{\\text{jets}}$ in each of the four $S_{\\textrm{NN}}$ bins. The signal" \
                          " distributions normalized to the predicted cross section for the RPV model with $m_{\\tilde{t}}$ = 450 GeV" \
                          " and the stealth SYY model with $m_{\\tilde{t}}$ = 850 GeV are shown for comparison."%(year.replace("Y", "20"))
        tab.location = "Data from Figure 4, located on page 11"
        tmap[year] = copy.deepcopy(tab)
        count += 1
    
    xvar = None
    
    for year, procd in vmap.items():
        for proc, nnBind in procd.items():
            for nnBin, v in nnBind.items():
    
                h = None
    
                # For the fit and data TGraphs read them as graphs
                # Otherwise, read everything as a simple 1D histo
                if "Fit" in proc or "Nobs" in proc: 
                    h = files[year].read_graph("%s_%s_%s"%(proc, year.replace("A","pre").replace("B","post"), nnBin))
                else:
                    h = files[year].read_hist_1d("%s_%s_%s"%(proc, year.replace("A","pre").replace("B","post"), nnBin))
    
                # Define one instance of the independent variable, here Njets
                if "sig" in proc:
                    xvar = V("$N_{jets}$", is_independent=True, is_binned=False, units="")
                    temp = h["x"]
                    for j in range(0,len(temp)):
                        temp[j] += 6.5
                        temp[j] = int(temp[j])
                    xvar.values = temp 
    
                # Make a variable object to hold one of the histos or graphs
                # and give it the y vals from the object
                v = V("%s %s"%(nnBin,names[proc]), is_independent=False, is_binned=False, units="")
                v.values = h["y"]
                v.add_qualifier("SQRT(S)", 13, "TeV")
                v.add_qualifier("LUMINOSITY", lumiMap[year], "fb$^{-1}$")
    
                # For the fit and data graphs, there is an associated uncertainty
                # So add that uncertainty to the Fit and Nobs variables
                if proc == "Fit":
                    temp = h["dy"]
                    uncArr, isSymm = makeUncArray(temp)
                    unc = U("unc.", is_symmetric=isSymm)
                    unc.values = uncArr 
                    v.add_uncertainty(unc)
    
                elif proc == "Nobs":
                    temp = h["dy"]
                    uncArr, isSymm = makeUncArray(temp)
                    unc = U("unc.", is_symmetric=isSymm)
                    unc.values = uncArr 
                    v.add_uncertainty(unc)
    
                # Add each sig1, sig2, fit, data variable to the corresponding table
                tmap[year].add_variable(v)
    
                # Add the xvar to each year, SNN table just once!
                if xvar is not None:
                    addUniqueVar(tmap[year], xvar)
    
    for year, t in tmap.items():
        sub.add_table(t)
    
    sub.create_files("output")

if __name__ == "__main__":
    makeFitPlotHEPData()
