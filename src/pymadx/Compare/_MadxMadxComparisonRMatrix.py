"""
Functions for comparing the optical distributions of two
BDSIM models.

Functions for plotting the individual optical functions, and an
eighth, helper function ''compare_all_optics``, to plot display all
seven in one go.
"""
import datetime as _datetime
from matplotlib.backends.backend_pdf import PdfPages as _PdfPages

from ._CompareCommon import _MakePlotterWithScale

_HORIZONTAL = [('RE11', 1.0,  'RE11 / mm'),
               ('RE12', 1.0,  'RE12 / mm'),
               ('RE16', 10.0, 'RE16 / %')]

_VERTICAL = [('RE33', 1.0,  'RE33 / mm'),
             ('RE34', 1.0,  'RE34 / mm'),
             ('RE36', 10.0, 'RE36 / %')]

PlotHorizontal = _MakePlotterWithScale(_HORIZONTAL, "S / m", "x / mm", "Horizontal", zeroLine=True)
PlotVertical   = _MakePlotterWithScale(_VERTICAL,   "S / m", "x / mm", "Vertical", zeroLine=True)

def MadxVsMadxRMatrix(first, second, first_name=None,
                      second_name=None, sOffsetSecond=0,
                      output_filename=None, **kwargs):
    """
    Display vertical and horizontal RMatrix components for two files.
    """
    figures = [
        PlotHorizontal(first, second, first_name, second_name, sOffsetSecond, **kwargs),
        PlotVertical(first, second, first_name, second_name, sOffsetSecond, **kwargs)
    ]

    if output_filename is not None:
        if 'pdf' in output_filename:
            with _PdfPages(output_filename) as pdf:
                for figure in figures:
                    pdf.savefig(figure)
                d = pdf.infodict()
                d['Title'] = "{} VS {} RMatrix Optical Comparison".format(first_name, second_name)
                d['CreationDate'] = _datetime.datetime.today()
            print("Written ", output_filename)
        else:
            for figure, label in zip(figures, ['x', 'y']):
                figure.savefig(output_filename+label+".png", dpi=500)
