import numpy as np
from hepdata_lib import Submission
from hepdata_lib import RootFileReader
from hepdata_lib import Table
from hepdata_lib import Variable, Uncertainty
from makeFitPlotHEPData import makeFitPlotHEPData

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

def makeCutFlow(submission, config):
    table = Table(config["name"])
    table.description = config["description"]
    table.location = config["location"]
    #table.keywords["observables"] = ["SIG"]
    #table.keywords["reactions"] = ["P P --> TOP --> tt + 6j"]
    
    data1 = config["data1"]
    error1 = config["error1"]
    data2 = config["data2"]
    error2 = config["error2"]
    
    #####################################################################################    
    d = Variable("Step", is_independent=True, is_binned=False, units="")
    d.values = np.array(list(i for i in range(0,data1.size)))

    cuts = Variable("Selection requirement", is_independent=False, is_binned=False, units="")
    cuts.values = config["cutnames"]
    cuts.add_qualifier("SQRT(S)", 13, "TeV")
    cuts.add_qualifier("LUMINOSITY", config["lumi"], "fb$^{-1}$")
    
    obs1 = Variable("RPV $m_{\\tilde{t}}$ = 450 GeV", is_independent=False, is_binned=False, units="")
    obs1.values = data1
    obs1.add_qualifier("SQRT(S)", 13, "TeV")
    obs1.add_qualifier("LUMINOSITY", config["lumi"], "fb$^{-1}$")
    
    unc_obs1 = Uncertainty("1 s.d.", is_symmetric=True)
    unc_obs1.values = error1
    obs1.add_uncertainty(unc_obs1)

    obs2 = Variable("SYY $m_{\\tilde{t}}$ = 850 GeV", is_independent=False, is_binned=False, units="")
    obs2.values = data2
    obs2.add_qualifier("SQRT(S)", 13, "TeV")
    obs2.add_qualifier("LUMINOSITY", config["lumi"], "fb$^{-1}$")
    
    unc_obs2 = Uncertainty("1 s.d.", is_symmetric=True)
    unc_obs2.values = error2
    obs2.add_uncertainty(unc_obs2)

    table.add_variable(d)
    table.add_variable(cuts)
    table.add_variable(obs1)
    table.add_variable(obs2)
    submission.add_table(table)
    
    
if __name__ == "__main__":
    submission = Submission()
    config = {}

    #Add fit Yeld tables
    makeFitPlotHEPData(submission)
    
    #Add RPV Combo limit plot
    config["name"] = "Figure 6a"    
    config["description"] = "Expected and observed 95% CL upper limit on the top squark pair production cross section as a function of the top squark mass for the RPV SUSY models."
    config["location"] = "Data from Figure 6a, located on page 13."
    config["image"] = "inputs/sigBrLim_RPV_Combo_Jun15_2020_CLs_Observed.pdf"
    config["inputData"] = "inputs/higgsCombineCombo.AsymptoticLimits.merged.MODELRPV.root"
    addLimitPlot(submission, config)
    
    config["name"] = "Figure 6b"
    config["description"] = "Expected and observed 95% CL upper limit on the top squark pair production cross section as a function of the top squark mass for the stealth SYY SUSY models."
    config["location"] = "Data from Figure 6b, located on page 13."
    config["image"] = "inputs/sigBrLim_SYY_Combo_Jun15_2020_CLs_Observed.pdf"
    config["inputData"] = "inputs/higgsCombineCombo.AsymptoticLimits.merged.MODELSYY.root"
    addLimitPlot(submission, config)

    #Add cutflows for signal
    config["name"] = "Signal cutflow 2016"
    config["description"] = "Cut flow of the signal region selection decribed in section II."
    config["location"] = "Cut flow of the two signals shown in Figure 4a. $M(l,b)$ refers to the invariant mass of the system formed by the b-tagged jet and the lepton. At each stage of the cut flow all event weights used for the signal region are applied."
    config["cutnames"] = ['#No cuts', '#Trigger and data quality', '#Exactly 1 lepton', '#$H_{T} > 300$ GeV', '#At least 1 b-tagged jet', '#$50 < M(l,b) < 250$ GeV', '#$N_{jets} \geq 7$']
    config["data1"] = np.array([39587.9, 10994, 8102.54, 8033.69, 6414.19, 5679.5, 3955])    
    config["error1"] = [51.3583, 26.7699, 22.9344, 22.8371, 20.256, 19.0251, 15.8792]
    config["data2"] = np.array([770.303, 232.683, 159.505, 159.505, 133.883, 112.255, 105.48])
    config["error2"] = [3.75411, 2.04207, 1.68552, 1.68552, 1.53699, 1.40427, 1.3614]
    config["lumi"] = 35.9
    makeCutFlow(submission, config)

    config["name"] = "Signal cutflow 2017"
    config["description"] = "Cut flow of the signal region selection decribed in section II."
    config["location"] = "Cut flow of the two signals shown in Figure 4b. $M(l,b)$ refers to the invariant mass of the system formed by the b-tagged jet and the lepton. At each stage of the cut flow all event weights used for the signal region are applied."
    config["cutnames"] = ['#No cuts', '#Trigger and data quality', '#Exactly 1 lepton', '#$H_{T} > 300$ GeV', '#At least 1 b-tagged jet', '#$50 < M(l,b) < 250$ GeV', '#$N_{jets} \geq 7$']
    config["data1"] = np.array([45076.9, 12001.1, 8947.06, 8872.36, 6816.34, 6053.54, 4267.69])
    config["error1"] = [59.4338, 29.9482, 25.7483, 25.6416, 21.9583, 20.6329,17.3188]
    config["data2"] = np.array([878.772, 254.262, 178.874, 178.874, 142.104, 120.355, 113.85])
    config["error2"] = [4.34625, 2.2916, 1.91071, 1.91071, 1.66191, 1.52037, 1.48061]
    config["lumi"] = 41.5
    makeCutFlow(submission, config)

    config["name"] = "Signal cutflow 2018A"
    config["description"] = "Cut flow of the signal region selection decribed in section II."
    config["location"] = "Cut flow of the two signals shown in Figure 4c. $M(l,b)$ refers to the invariant mass of the system formed by the b-tagged jet and the lepton. At each stage of the cut flow all event weights used for the signal region are applied."
    config["cutnames"] = ['#No cuts', '#Trigger and data quality', '#Exactly 1 lepton', '#$H_{T} > 300$ GeV', '#At least 1 b-tagged jet', '#$50 < M(l,b) < 250$ GeV', '#$N_{jets} \geq 7$']
    config["data1"] =  np.array([23126.8, 6224.31, 4679.45, 4641.48, 3638.17, 3239.1, 2267.07])
    config["error1"] = [25.7698, 13.1026, 11.3365, 11.2918, 9.93458, 9.35931, 7.82141]
    config["data2"] =  np.array([450.339, 134.459, 94.6625, 94.6625, 76.4676, 65.412, 61.3689])
    config["error2"] = [1.87485, 1.01652, 0.854426, 0.854426, 0.766858, 0.69681, 0.67013]
    config["lumi"] = 21.1
    makeCutFlow(submission, config)
    
    config["name"] = "Signal cutflow 2018B"
    config["description"] = "Cut flow of the signal region selection decribed in section II."
    config["location"] = "Cut flow of the two signals shown in Figure 4d. $M(l,b)$ refers to the invariant mass of the system formed by the b-tagged jet and the lepton. At each stage of the cut flow all event weights used for the signal region are applied."
    config["cutnames"] = ['#No cuts', '#Trigger and data quality', '#Exactly 1 lepton', '#$H_{T} > 300$ GeV', '#At least 1 b-tagged jet', '#$50 < M(l,b) < 250$ GeV', '#$N_{jets} \geq 7$']
    config["data1"] =  np.array([42544.4, 8039.05, 6045.75, 6000.3, 4726.61, 4210.44, 2898.18])
    config["error1"] = [47.3966, 20.309, 17.5823, 17.519, 15.4541, 14.5699, 12.0853]
    config["data2"] =  np.array([827.605, 172.116, 121.12,  121.12, 97.4285, 83.1049, 77.4309])
    config["error2"] =  [3.44402, 1.57136, 1.32318, 1.32318, 1.18439, 1.06967, 1.02218]
    config["lumi"] = 38.7
    makeCutFlow(submission, config)
    
    submission.create_files("output")
    
