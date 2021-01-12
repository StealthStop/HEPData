import numpy as np
from hepdata_lib import Submission
from hepdata_lib import RootFileReader
from hepdata_lib import Table
from hepdata_lib import Variable, Uncertainty

def addLimitPlot(submission, config):
    table = Table(config["name"])
    table.description = config["description"]
    table.location = config["location"]
    table.keywords["observables"] = ["SIG"]
    table.keywords["reactions"] = ["P P --> TOP --> tt + 6j"]
    table.add_image(config["image"])
    
    reader = RootFileReader(config["inputData"])
    data = reader.read_limit_tree()
    stop_pair_Br           = np.array([10.00, 4.43, 2.15, 1.11,  0.609, 0.347, 0.205, 0.125, 0.0783, 0.0500, 0.0326, 0.0216, 0.0145, 0.00991, 0.00683, 0.00476, 0.00335, 0.00238, 0.00170, 0.00122, 0.000887, 0.000646, 0.000473])
    stop_pair_Br1SPpercent = np.array([6.65, 6.79, 6.99, 7.25,  7.530, 7.810, 8.120, 8.450, 8.8000, 9.1600, 9.5300, 9.9300, 10.3300, 10.76, 11.2, 11.65, 12.12, 12.62, 13.13, 13.66, 14.21, 14.78, 15.37])
    stop_pair_unc = stop_pair_Br*stop_pair_Br1SPpercent/100.0
    stop_pair_up = stop_pair_Br+stop_pair_unc
    stop_pair_down = stop_pair_Br-stop_pair_unc

    nData = len(data)
    for mass_id in range(0,nData):
        data[mass_id][1:] = stop_pair_Br[mass_id]*data[mass_id][1:] 
    
    #####################################################################################    
    d = Variable("Top squark mass", is_independent=True, is_binned=False, units="GeV")
    d.values = data[:,0]
    
    sig = Variable("Top squark cross section", is_independent=False, is_binned=False, units="pb")
    sig.values = np.array(stop_pair_Br[:nData])
    sig.add_qualifier("Limit", "")
    sig.add_qualifier("SQRT(S)", 13, "TeV")
    sig.add_qualifier("LUMINOSITY", 137, "fb$^{-1}$")
    
    obs = Variable("Observed cross section upper limit at 95% CL", is_independent=False, is_binned=False, units="pb")
    obs.values = data[:,6]
    obs.add_qualifier("Limit", "Observed")
    obs.add_qualifier("SQRT(S)", 13, "TeV")
    obs.add_qualifier("LUMINOSITY", 137, "fb$^{-1}$")
    
    exp = Variable("Expected cross section upper limit at 95% CL", is_independent=False, is_binned=False, units="pb")
    exp.values = data[:,3]
    exp.add_qualifier("Limit", "Expected")
    exp.add_qualifier("SQRT(S)", 13, "TeV")
    exp.add_qualifier("LUMINOSITY", 137, "fb$^{-1}$")
    
    unc_sig = Uncertainty("1 s.d.", is_symmetric=False)
    unc_sig.set_values_from_intervals(zip(stop_pair_up[:nData], stop_pair_down[:nData]), nominal=sig.values)
    sig.add_uncertainty(unc_sig)
    
    # +/- 1 sigma
    unc_1s = Uncertainty("1 s.d.", is_symmetric=False)
    unc_1s.set_values_from_intervals(zip(data[:,2], data[:,4]), nominal=exp.values)
    exp.add_uncertainty(unc_1s)
    
    # +/- 2 sigma
    unc_2s = Uncertainty("2 s.d.", is_symmetric=False)
    unc_2s.set_values_from_intervals(zip(data[:,1], data[:,5]), nominal=exp.values)
    exp.add_uncertainty(unc_2s)
    
    table.add_variable(d)
    table.add_variable(sig)
    table.add_variable(obs)
    table.add_variable(exp)
    submission.add_table(table)

if __name__ == "__main__":
    submission = Submission()
    config = {}
    
    #Add RPV Combo limit plot
    config["name"] = "Figure 6a"
    config["description"] = "Exclusion limits on the product of the production cross section and the branching fraction for a new spin-2 resonance decaying to WW, as a function of the resonance mass hypothesis."    
    config["location"] = "Data from Figure 6a, located on page 14."
    config["image"] = "limit_inputs/sigBrLim_RPV_Combo_Jun15_2020_CLs_Observed.pdf"
    config["inputData"] = "limit_inputs/higgsCombineCombo.AsymptoticLimits.merged.MODELRPV.root"
    addLimitPlot(submission, config)

    submission.create_files("output")
    
