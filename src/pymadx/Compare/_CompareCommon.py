from .. import Data as _Data
from .. import Plot as _Plot

from os import path as _path

import matplotlib.pyplot as _plt

def _LoadTfsInfput(tfs_in, name):
    """
    Return tfs_in as a Tfs instance, which should either be a path
    to a TFS file or a Tfs instance, and in either case, generate a
    name if None is provided, and return that as well.
    """
    if isinstance(tfs_in, str):
        if not _path.isfile(tfs_in):
            raise IOError("file \"{}\" not found!".format(tfs_in))
        name = (_path.splitext(_path.basename(tfs_in))[0]
                if name is None else name)
        return _Data.Tfs(tfs_in), name
    try:
        name = tfs_in.filename if name is None else name
        return tfs_in, name
    except AttributeError:
        raise TypeError("Expected Tfs input is neither a "
                        "file path nor a Tfs instance: {}".format(tfs_in))


# use closure to avoid tonnes of boilerplate code as happened with
# MadxBdsimComparison.py
def _MakePlotter(plot_info_tuples, x_label, y_label, title):
    def f_out(first, second, first_name=None, second_name=None, **kwargs):
        """first and second should be tfs files."""
        first, first_name = _LoadTfsInfput(first, first_name)
        second, second_name = _LoadTfsInfput(second, second_name)

        plot = _plt.figure(title, figsize=(9,5), **kwargs)
        # Loop over the variables in plot_info_tuples and draw the plots.
        colours = ["tab:blue", "tab:orange"]
        for colour, (var, legend_name) in zip(colours, plot_info_tuples):
            _plt.plot(first.GetColumn('S'), first.GetColumn(var), c=colour,
                      label="{}: {}".format(legend_name, first_name), **kwargs)
            _plt.plot(second.GetColumn('S'), second.GetColumn(var), '--', c=colour,
                      label="{}: {}".format(legend_name, second_name), **kwargs)

        # Set axis labels and draw legend
        axes = _plt.gcf().gca()
        axes.set_ylabel(y_label)
        axes.set_xlabel(x_label)
        axes.legend(loc='best')

        _Plot.AddMachineLatticeToFigure(plot, first)
        _plt.show(block=False)
        return plot
    return f_out


def _MakePlotterWithScale(plot_info_tuples, x_label, y_label, title, zeroLine=False):
    def f_out(first, second, first_name=None, second_name=None, dSSecond=0, **kwargs):
        """first and second should be tfs files."""
        first, first_name = _LoadTfsInfput(first, first_name)
        second, second_name = _LoadTfsInfput(second, second_name)
        cmap = _plt.get_cmap('tab20')
        plot = _plt.figure(title, figsize=(9,5), **kwargs)
        # Loop over the variables in plot_info_tuples and draw the plots.
        for i, (var, scale, legend_name) in enumerate(plot_info_tuples):
            _plt.plot(first.GetColumn('S'),
                      first.GetColumn(var)*scale,
                      c=cmap(2*i),
                      label="{}: {}".format(legend_name, first_name),
                      **kwargs)
            _plt.plot(second.GetColumn('S')+dSSecond,
                      second.GetColumn(var)*scale,
                      c=cmap(2*i + 1), ls='--',
                      label="{}: {}".format(legend_name, second_name),
                      **kwargs)
        if zeroLine:
            s = first.GetColumn('S')
            _plt.plot([s[0], s[-1]], [0, 0], c='grey', alpha=0.3)

        # Set axis labels and draw legend
        axes = _plt.gcf().gca()
        axes.set_ylabel(y_label)
        axes.set_xlabel(x_label)
        axes.legend(loc='best')

        _Plot.AddMachineLatticeToFigure(plot, first)
        _plt.show(block=False)
        return plot
    return f_out