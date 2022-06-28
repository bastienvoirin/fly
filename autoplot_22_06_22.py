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
    #ref_bkg = 10829
    #print(ref_sig)
    #print(ref_bkg)
    sig_scale = 0.06855
    bkg_scale = 0.17946
    sig = {i: (len(pdf_sig["event"]), ref_sig) for (i, pdf_sig) in pdfs_sig.items()}
    bkg = {i: (len(pdf_bkg["event"]), ref_bkg) for (i, pdf_bkg) in pdfs_bkg.items()}
    #print("sig:", sig)
    #print("bkg:", bkg)
    sig = {i: len(pdf_sig["event"]) / ref_sig for (i, pdf_sig) in pdfs_sig.items()}
    bkg = {i: len(pdf_bkg["event"]) / ref_bkg for (i, pdf_bkg) in pdfs_bkg.items()}
    #print("sig:", sig)
    #print("bkg:", bkg)

    fig, ax = plt.subplots()
    ax.set_title("ROC curve")
    ax.set_xlabel("B(after) / B(before)")
    ax.set_ylabel("S(after) / S(before)")
    ax.set_xticks(np.linspace(0.0, 1.0, 11))
    ax.set_yticks(np.linspace(0.0, 1.0, 11))
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.0])
    ax.scatter(bkg.values(), sig.values())
    for i, (xi, yi, txt) in enumerate(zip(bkg.values(), sig.values(), sig.keys())):
        S = yi * ref_sig * sig_scale
        B = xi * ref_bkg * bkg_scale
        SBR = 100 * S/B
        print(f"{i}: S = {S:.2f}, B = {B:.2f}, S/B = {SBR:.4f}%")
        ax.annotate(f"{txt}",#\nS = {S:.2f}, B = {B:.2f}\nS/B = {SBR:.4f}%",
                    xy=(xi, yi), textcoords="offset points", va="center",
                    xytext=(4, 0))
    plt.show()

if __name__ == "__main__":
    signal_fname = {}
    ttbar_fname = {}
    wz_fname = {}
    for i in range(10):
        signal_fname[str(i)] = f"analysed/analysed_TT2L2Nu_2_Muons_signal.root_{i}"
        ttbar_fname[str(i)]  = f"analysed/analysed_TT2L2Nu_2_Muons_ttbar.root_{i}"
        wz_fname[str(i)]     = f"analysed/analysed_TT2L2Nu_2_Muons_WZ.root_{i}"
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
    wz = {}
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
    for i, fname in wz_fname.items():
        print(i, fname)
        rdl = ROOTDataLoader(fname, "outputTree2", variables)
        pdf = rdl.getPandasDataFrame()
        wz[str(i)] = pdf
    rdl_ref_signal = ROOTDataLoader("analysed/analysed_TT2L2Nu_2_Muons_signal.root_9", "outputTree2", variables)
    pdf_ref_signal = rdl_ref_signal.getPandasDataFrame()
    rdl_ref_ttbar  = ROOTDataLoader("analysed/analysed_TT2L2Nu_2_Muons_ttbar.root_9", "outputTree2", variables)
    pdf_ref_ttbar  = rdl_ref_ttbar.getPandasDataFrame()
    rdl_ref_wz  = ROOTDataLoader("analysed/analysed_TT2L2Nu_2_Muons_WZ.root_9", "outputTree2", variables)
    pdf_ref_wz  = rdl_ref_wz.getPandasDataFrame()
    print("\nttbar")
    plotROCCurve(pdf_ref_signal, pdf_ref_ttbar, signal, ttbar)
    print("\nWZ")
    plotROCCurve(pdf_ref_signal, pdf_ref_wz, signal, wz)

    """
    print(pdf.keys())
    fig, axs = plt.subplots(3, 3)
    for i, var in enumerate(variables):
        j = i%3
        k = i//3
        seaborn.histplot(pdf[var], bins=30, ax=axs[j, k])
    plt.show()
    """
