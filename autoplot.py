#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pyrootplots.ROOTDataLoader import ROOTDataLoader
from pyrootplots.Cut import Cut
from pyrootplots.ROCCurve import ROCCurve
from pyrootplots.Criterion import Criterion
from pyrootplots.Condition import Condition
import seaborn
from matplotlib import pyplot as plt
import numpy as np

def plotROCCurve(pdf_ref_sig, pdf_ref_bkg, pdfs_sig, pdfs_bkg):
    ref_sig = len(pdf_ref_sig["event"])
    ref_bkg = len(pdf_ref_bkg["event"])
    ref_bkg = 10829
    sig = {i: len(pdf_sig["event"]) for (i, pdf_sig) in pdfs_sig.items()}
    bkg = {i: len(pdf_bkg["event"]) for (i, pdf_bkg) in pdfs_bkg.items()}
    criteria = {"0": Criterion("Sum_all_muons_pt", Condition.GREATER_THAN_OR_EQUAL_TO, [160.0]),
                "3": Criterion("DR_dimuon", Condition.LESS_THAN_OR_EQUAL_TO, [1.8]),
                "4": Criterion("DR_dimuon", Condition.LESS_THAN_OR_EQUAL_TO, [1.2]),
                "5": Criterion("DR_dimuon", Condition.LESS_THAN_OR_EQUAL_TO, [0.8]),
                "6": Criterion("DPhi_dimuon", Condition.LESS_THAN_OR_EQUAL_TO, [1.57]),
                "7": Criterion("DPhi_dimuon", Condition.LESS_THAN_OR_EQUAL_TO, [1.2]),
                "8": Criterion("DPhi_dimuon", Condition.LESS_THAN_OR_EQUAL_TO, [0.8])}
    print(*[f"{key}: {val}" for (key, val) in criteria.items()], sep="\n")
    print(sig)
    print(bkg)
    sig_scale = 0.06855
    bkg_scale = 0.17946
    print(ref_sig)
    print(ref_bkg)
    cuts = [Cut(criterion, ref_sig, sigEvtAfter, {"ttbar": ref_bkg}, {"ttbar": bkgEvtAfter}, sig_scale, bkg_scale, scale=True)
            for (criterion, sigEvtAfter, bkgEvtAfter) in zip(criteria.values(), sig.values(), bkg.values())]
    print(*cuts, sep="\n")
    print(*cuts[:1], sep="\n")
    print(*cuts[1:4], sep="\n")
    print(*cuts[4:], sep="\n")
    ROCCurve(cuts[:1],  scatter={"marker": "x"}, txt=map(lambda cut: cut.criterion, cuts[:1])).plot(whichBkg="ttbar")
    ROCCurve(cuts[1:4], scatter={"marker": "x"}, txt=map(lambda cut: cut.criterion, cuts[1:4])).plot(whichBkg="ttbar")
    ROCCurve(cuts[4:],  scatter={"marker": "x"}, txt=map(lambda cut: cut.criterion, cuts[4:])).plot(whichBkg="ttbar")
    plt.show()

if __name__ == "__main__":
    signal_fname = {}
    ttbar_fname = {}
    for i in (0, 3, 4, 5, 6, 7, 8):
        signal_fname[str(i)] = f"analysed/analysed_TT2L2Nu_2_Muons_signal.root_{i}"
        ttbar_fname[str(i)]  = f"analysed/analysed_TT2L2Nu_2_Muons_ttbar.root_{i}"
    variables = ["Selected_muon_deltaR_min_neutral",
                 "dimuon_deltaR",
                 "St",
                 "DeltaPhi_dimuon",
                 "Vectorial_sum_dimuon_mass",
                 "Vectorial_sum_bl_mass",
                 "DeltaPhi_dimuon",
                 "Selected_muon_sum_all_muons_pt",
                 "event"]
    signal = {}
    ttbar = {}
    for i, fname in signal_fname.items():
        print(i, fname)
        rdl = ROOTDataLoader(fname, "outputTree2", variables)
        pdf = rdl.getPandasDataFrame()
        signal[str(i)] = pdf
    for i, fname in ttbar_fname.items():
        print(i, fname)
        rdl = ROOTDataLoader(fname, "outputTree2", variables)
        pdf = rdl.getPandasDataFrame()
        ttbar[str(i)] = pdf
    rdl_ref_signal = ROOTDataLoader("analysed/analysed_TT2L2Nu_2_Muons_signal.root", "outputTree2", variables)
    pdf_ref_signal = rdl_ref_signal.getPandasDataFrame()
    rdl_ref_ttbar  = ROOTDataLoader("analysed/analysed_TT2L2Nu_2_Muons_ttbar.root", "outputTree2", variables)
    pdf_ref_ttbar  = rdl_ref_ttbar.getPandasDataFrame()
    plotROCCurve(pdf_ref_signal, pdf_ref_ttbar, signal, ttbar)

    """
    print(pdf.keys())
    fig, axs = plt.subplots(3, 3)
    for i, var in enumerate(variables):
        j = i%3
        k = i//3
        seaborn.histplot(pdf[var], bins=30, ax=axs[j, k])
    plt.show()
    """
