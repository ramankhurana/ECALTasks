from ROOT import * 
import os
import sys 
from sets import Set


def maskedTT(filename, histname, txtfile):
    fout = open(txtfile, "w")
    f = TFile(filename,"READ")
    eta_phi_map_tp = f.Get(histname) 
    masked_TT_list = []
    for ibinx in range(1,eta_phi_map_tp.GetNbinsX()+1):
        if abs(eta_phi_map_tp.GetXaxis().GetBinCenter(ibinx)) < 18: continue 
        for ibiny in range(1,eta_phi_map_tp.GetNbinsY()+1):
            if ( str(eta_phi_map_tp.GetBinContent(ibinx,ibiny)) == str(0.0) ) :
                TT_= str(eta_phi_map_tp.GetXaxis().GetBinCenter(ibinx)) + " " + str(eta_phi_map_tp.GetYaxis().GetBinCenter(ibiny)) + "\n"
                TT_l = str(eta_phi_map_tp.GetXaxis().GetBinCenter(ibinx)) + " " + str(eta_phi_map_tp.GetYaxis().GetBinCenter(ibiny))
                fout.write(TT_)
                masked_TT_list.append(TT_l)
    fout.close()
    masked_TT_set = Set(masked_TT_list)
    return masked_TT_set





masked_tp = maskedTT("ttplot_debug_306456_tp32_oldped.root", "ieta_vs_iphi_TP", "masked_TT.txt")

masked_etp = maskedTT("ttplot_debug_306456_tp32_oldped.root", "ieta_vs_iphi_ETP", "masked_ETT.txt")

#print masked_tp
#print masked_etp

masked_all = masked_tp | masked_etp 

#print masked_all

f_all = open("masted_TT_all.txt","w")

for iele in masked_all:
    f_all.write(iele)
    f_all.write("\n")
    
f_all.close()



