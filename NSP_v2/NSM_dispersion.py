# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Nonlinear Static Procedure with dispersion

# <markdowncell>

# The user should change only section 1 "Define Option"

# <codecell>

"""
Created on Tue May 13 18:04:57 2014

@author: chiaracasotto
"""
# Clear existing variables
def clearall():
    all = [var for var in globals() if var[0] != "_"]
    for var in all:
        del globals()[var]
clearall()

# Import functions
import matplotlib.pyplot as plt
import numpy as np
import os
from common.get_data import read_data
from common.fragility_process import fragility_process
from common.export_outputs import export
from common.vulnerability_process import vulnerability_process

pi = 3.141592653589793
plt.close("all")
cd = os.getcwd()

# <headingcell level=2>

# 1. Define Options

# <markdowncell>

# >####an_type (analysis type) <br />
# 0: elastoplastic Pushover, inelastic displacement obtained from Ruiz-Garcia Miranda [1] <br />
# 1: any pushover shape (bilinear, trilinear, quadrilinear) inelastic displacement obtained from spo2ida [2] <br />
# 2: any pushover shape (bilinear, trilinear, quadrilinear) inelastic displacement obtained from Dolsek and Fajfar [3] <br />
# 
# >####in_type (input_type)
# 0: input is idealised pushover curve <br />
# 1: input is pushover curve <br />

# <codecell>

an_type = 2
in_type = 1

# <markdowncell>

# >####vuln
# 0: derive fragility curves <br />
# 1: derive vulnerability curves from fragility curves <br />
# ####g
# Units of g compatible with T and displacement
# ####iml
# Array of intensity measure level used to discretise the fragility curves and get discrete vulnerability curves

# <codecell>

vuln = 0
g = 9.81
iml = np.linspace(0.1,2,50)

# <markdowncell>

# >####plotflag 
# 4 integers for each plot: idealised pushover curve, 16%-50%-84% ida curves, fragility curves, vulenarbility curve <br />
# 1: plot <br />
# 0: don't plot <br />
# >####linew
# Line width for plots
# >####fontsize
# Fontsize used for labels, graphs etc.
# >####units
# List of 3 strings defining the displacements, forces and Sa units as ['[kN]' '[m]' ['m/s^2]]

# <codecell>

plotflag = [1, 1, 1, 1]
linew = 2
fontsize = 10
units = ['[m]', '[kN]', '[g]']

# <markdowncell>

# >####N
# Number of points per segment of IDA curve derived with spo2ida (only spo2ida procedure)
# >####MC 
# Number of Monte Carlo simulations to account for uncertainty in damage thresholds (only spo2ida procedure)
# >####Tc, Td
# Constant accel-constant velocity  and constant velocity-constant displacement corner periods of a Newmark-Hall type spectrum. Default
# values <br />
# Tc = 0.5 <br />
# Td = 1.8 <br />

# <codecell>

N = 10
MC = 10
Tc = 0.5
Td = 1.8

# <headingcell level=2>

# 2. Read data from csv input file and Process data

# <markdowncell>

# Obtain: First period, first modal participation factor, weights associated to each building, roof disp. at each limit state, idealised pushover parameters, dispersion in the limit states

# <codecell>

plot_feature = [plotflag, linew, fontsize, units, iml]
[T, Gamma, w, dcroof, SPO, bUthd, noBlg, Tav, Sa_ratios] = read_data(in_type,an_type,linew,fontsize,units,plotflag[0])

# <headingcell level=2>

# 3.a Derive Fragility

# <codecell>

[log_meanSa, log_stSa] = fragility_process(vuln, an_type, T, Gamma, w, dcroof, SPO, bUthd, noBlg, g, MC, Sa_ratios, plot_feature, N, Tc, Td)

# <headingcell level=3>

# 3.a.2 Plot and Export fragility

# <codecell>

export(vuln, plot_feature, log_meanSa, log_stSa)

# <headingcell level=2>

# 3.b Derive Vulnerability curves from damage-to-loss ratios

# <codecell>

[LR50, bLR] = vulnerability_process(vuln, an_type, T, Gamma, w, dcroof, SPO, bUthd, noBlg, g, MC, Sa_ratios, plot_feature, N, Tc, Td)

# <headingcell level=3>

# 3.b.2 Derive and export Vulnerability

# <codecell>

export(vuln, plot_feature, LR50, bLR)

# <codecell>


