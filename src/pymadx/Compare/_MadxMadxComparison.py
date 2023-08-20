"""
Functions for comparing the optical distributions of two
BDSIM models.

Functions for plotting the individual optical functions, and an
eighth, helper function ''compare_all_optics``, to plot display all
seven in one go.
"""
import datetime as _datetime
from matplotlib.backends.backend_pdf import PdfPages as _PdfPages

from ._CompareCommon import _MakePlotter

# Predefined lists of tuples for making the standard plots,
# format = (optical_var_name, optical_var_error_name, legend_name)

_BETA = [("BETX", r'$\beta_{x}$'),
         ("BETY", r'$\beta_{y}$')]

_ALPHA = [("ALFX", r"$\alpha_{x}$"),
          ("ALFY", r"$\alpha_{y}$")]

_DISP = [("DXBETA", r"$D_{x}$ / $\beta$"),
         ("DYBETA", r"$D_{y}$ / $\beta$")]

_DISP_P = [("DPXBETA", r"$D_{p_{x}}$"),
           ("DPYBETA", r"$D_{p_{y}}$")]

_SIGMA = [("SIGMAX", r"$\sigma_{x}$"),
          ("SIGMAY", r"$\sigma_{y}$")]

_SIGMA_P = [("SIGMAXP", r"$\sigma_{xp}$"),
            ("SIGMAYP", r"$\sigma_{yp}$")]

_MEAN = [("X", r"$\bar{x}$"),
         ("Y", r"$\bar{y}$")]

PlotBeta   = _MakePlotter(_BETA, "S / m", r"$\beta_{x,y}$ / m", "Beta")
PlotAlpha  = _MakePlotter(_ALPHA, "S / m", r"$\alpha_{x,y}$ / m", "Alpha")
PlotDisp   = _MakePlotter(_DISP, "S / m", r"$D_{x,y} / m$", "Dispersion")
PlotDispP  = _MakePlotter(_DISP_P, "S / m", r"$D_{p_{x},p_{y}}$ / m", "Momentum_Dispersion")
PlotSigma  = _MakePlotter(_SIGMA, "S / m", r"$\sigma_{x,y}$ / m", "Sigma")
PlotSigmaP = _MakePlotter(_SIGMA_P, "S / m", r"$\sigma_{xp,yp}$ / rad", "SigmaP")
PlotMean   = _MakePlotter(_MEAN, "S / m", r"$\bar{x}, \bar{y}$ / m", "Mean")

def MADXVsMADX(first, second, first_name=None,
               second_name=None, saveAll=True, 
               outputFileName=None, **kwargs):
    """
    Display all the optical function plots for the two input optics files.
    """
    figures = [
    PlotBeta(first, second, first_name=first_name,
             second_name=second_name, **kwargs),
    PlotAlpha(first, second, first_name=first_name,
              second_name=second_name, **kwargs),
    PlotDisp(first, second, first_name=first_name,
             second_name=second_name, **kwargs),
    PlotDispP(first, second, first_name=first_name,
              second_name=second_name, **kwargs),
    PlotSigma(first, second, first_name=first_name,
              second_name=second_name, **kwargs),
    PlotSigmaP(first, second, first_name=first_name,
               second_name=second_name, **kwargs),
    PlotMean(first, second, first_name=first_name,
             second_name=second_name, **kwargs),
    ]

    if saveAll:
        if outputFileName is not None:
            output_filename = outputFileName
        else:
            output_filename = "optics-report.pdf"

        with _PdfPages(output_filename) as pdf:
            for figure in figures:
                pdf.savefig(figure)
            d = pdf.infodict()
            d['Title'] = "{} VS {} Optical Comparison".format(first_name, second_name)
            d['CreationDate'] = _datetime.datetime.today()
        print("Written ", output_filename)
