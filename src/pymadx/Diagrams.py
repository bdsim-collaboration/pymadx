"""
Various plots for madx TFS files using the pymadx Tfs class
"""
from collections import defaultdict as _defaultdict
import numpy as _np
import matplotlib.patches as _patches
import matplotlib.pyplot as _plt
from matplotlib.collections import PatchCollection as _PatchCollection

import pymadx.Data as _Data


defaultElementColours = {'DRIFT': u'#c0c0c0',
                         'QUADRUPOLE': u'#d10000',
                         'RBEND': u'#0066cc',
                         'SBEND': u'#0066cc',
                         'HKICKER': u'#4c33b2',
                         'VKICKER': u'#ba55d3',
                         'SOLENOID': u'#ff8800',
                         'RCOLLIMATOR': 'k',
                         'ECOLLIMATOR': 'k',
                         'COLLIMATOR': 'k',
                         'SEXTUPOLE': u'#ffcc00',
                         'OCTUPOLE': u'#00994c'
                         }

def Rotate(points, angle, origin=None):
    """
    Apply a rotation to an array of [x,y] pairs about (0,0) or an optional origin point.

    :param points: numpy array of [x,y] points (any length)
    :type points: numpy.array([x1,y1], [x2,y2],...])
    :param angle: angle to rotate by in radians. Positive is anti-clockwise in the x-y plane.
    :type angle: float, int
    :param origin: optional [x,y] to rotate the data about - default [0,0]
    :type origin: np.array([x,y])
    """
    if origin is None:
        origin = _np.array([[0, 0]])
    c, s = _np.cos(-angle), _np.sin(-angle)
    R = _np.array(((c, -s), (s, c)))
    return (points - origin) @ R + origin

def RotateTranslate(points, angle, translation, rotationOrigin=None):
    """
    Apply a rotation and then add a translation to all points.

    See Rotate()

    :param translation: [x,y] to add to all points after rotation
    :type translation: np.array([x,y])
    """
    rotated = Rotate(points, angle, rotationOrigin)
    rotoTranslated = rotated + translation
    return rotoTranslated

def Polygon(edges, colour, alpha):
    """
    Return a polygon patch from a list of points as a numpy array.
    """
    return _patches.Polygon(edges, colour=colour, fill=True, alpha=alpha)

def _CoilPolygonsDipoleH(xend, yend, length, rotation, roto_translation, dx=0, coil_dict=None, colour="#b87333",
                         alpha=1.0, greyOut=False, greyColour='#c0c0c0'):
    if greyOut:
        colour = greyColour
    cl = coil_dict['coil_length']
    if cl == 0:
        return []
    cw = coil_dict['coil_width']
    cdx = coil_dict['coil_dx']
    arc1 = _np.array([[_np.sin(theta), _np.cos(theta)] for theta in _np.linspace(0, _np.pi/2, 5)])
    r = 0.5*cl
    arc_top = arc1*r + _np.array([cl-r, (0.5*cw) - r ])
    arc_bot = _np.array(arc_top)[::-1]
    arc_bot[:,1] *= -1 # flip y
    coil_offset = _np.array([0, cdx + dx])
    edges_out = _np.array([[0, 0.5*cw], *arc_top, *arc_bot, [0, -0.5*cw]]) + coil_offset
    edges_in = _np.array(edges_out)
    edges_in[:,0] *= -1 # flip x
    edges_in[:,0] -= length
    global_offset = _np.array([xend, yend])
    edges_out = Rotate(edges_out, rotation) + global_offset
    edges_in = Rotate(edges_in, rotation) + global_offset
    edges_out = roto_translation * edges_out
    edges_in = roto_translation * edges_in
    return [_patches.Polygon(edges_out, color=colour, fill=True, alpha=alpha),
            _patches.Polygon(edges_in, color=colour, fill=True, alpha=alpha)]

def _CoilPolygonsQuad(xend, yend, length, rotation, roto_translation, alpha, coil_dict, colour="#b87333",
                      greyOut=False, greyColour='#c0c0c0'):
    if greyOut:
        colour = greyColour
    cl = coil_dict['coil_length']
    if cl == 0:
        return []
    cw = coil_dict['coil_width']
    arc1 = _np.array([[_np.sin(theta), _np.cos(theta)] for theta in _np.linspace(0, _np.pi/2, 5)])
    r = 0.5*cl
    arc_bot = arc1*r + _np.array([cl-r, -r])
    arc_top = _np.array(arc_bot)
    arc_top[:,1] *= -1
    edges_out = _np.array([*arc_bot, [cl, -0.5*cw], [0, -0.5*cw], [0, 0.5*cw], [cl, 0.5*cw], *(arc_top[::-1])])
    edges_in = _np.array(edges_out)
    edges_in[:,0] *= -1 # flip x
    edges_in[:,0] -= length
    global_offset = _np.array([xend, yend])
    edges_out = Rotate(edges_out, rotation) + global_offset
    edges_in = Rotate(edges_in, rotation) + global_offset
    edges_out = roto_translation * edges_out
    edges_in = roto_translation * edges_in
    return [_patches.Polygon(edges_out, color=colour, fill=True, alpha=alpha),
            _patches.Polygon(edges_in, color=colour, fill=True, alpha=alpha)]

def _CurvedLine(x0, y0, xydir, angle, arcLength, stepSize=0.1):
    return [[x0, y0]]

def _SBend():
    return None # polygon

def _Rectangle(xend, yend, half_width, length, rotation, roto_translation, colour, alpha, dx=0, dy=0):
    """
    dx, dy are in curvilinear x,y so are applied to plot y,x respectively for an ZX plot.
    """
    edges = _np.array([[0, 0.5 * half_width + dx],
                       [0, -0.5 * half_width + dx],
                       [-length, -0.5 * half_width + dx],
                       [-length, 0.5 * half_width + dx]])
    edges = RotateTranslate(edges, rotation, _np.array([xend, yend]))
    edges = roto_translation * edges
    return _patches.Polygon(edges, color=colour, fill=True, alpha=alpha)

def _Collimator(xend, yend, width, length, rotation, roto_translation, colour, alpha):
    edges = _np.array([[0,            0.3*width], [0,           -0.3*width], [-0.3*length, -0.3*width],
                       [-0.3*length, -0.5*width], [-0.7*length, -0.5*width], [-0.7*length, -0.3*width],
                       [-length,     -0.3*width], [-length,      0.3*width], [-0.7*length,  0.3*width],
                       [-0.7*length,  0.3*width], [-0.7*length,  0.5*width], [-0.3*length,  0.5*width],
                       [-0.3*length,  0.3*width]])
    edges = RotateTranslate(edges, rotation, _np.array([xend, yend]))
    edges = roto_translation * edges
    return _patches.Polygon(edges, color=colour, fill=True, alpha=alpha)

def _UpdateParams(elementDict, typeDict, element, params, insideFactor):
    n = element['NAME']
    allKeys = set(params.keys())
    allowedTypeKeys = set(allKeys)
    for k, prms in elementDict.items():
        if k in n:
            params.update(prms)
            eleKeys = set(prms.keys())
            allowedTypeKeys -= eleKeys
            break
    for k, prms in typeDict.items():
        if k in n:
            keysToTake = allowedTypeKeys.intersection(set(prms.keys()))
            subDict = {l:prms[l] for l in keysToTake}
            params.update(subDict)
            break
    insideFactor2 = 1 if params['inside'] else -1
    params['dx'] *= insideFactor * insideFactor2
    params['coil_dx'] *= insideFactor * insideFactor2
    return params

def Survey2DZX(survey_tfsfile, ax=None, greyOut=False, offsetRotoTranslation=None, title='', outputfilename=None,
               zOffset=0, defaultWidth=0.5, elementDict=None, typeDict=None, funcDict=None, maskNames=None,
               ignoreNames=None, resolution=0.1, defaultCoilLength=0.15, pipeRadius=None, pipeMaskRanges=None,
               invisibleAxes=False, arrowsDy=0, arrowsDx=0, defaultAlpha=None):
    """
    Draw a schematic of a beamline in the Z-X plane. Allow styling of different elements by type or by name.

    :param survey_tfsfile: list of tfs files as strings or already loaded pymadx.Data.Tfs objects.
    :type survey_tfsfile: str, pymadx.Data.Tfs
    :param ax: matplotlib axes object
    :type ax: matplotlib.axes.Axes, None
    :param greyOut: whether to plot the whole line in greyscale
    :type greyOut: bool
    :param offsetRotoTranslation: Optional roto-translation to apply before plotting to offset the whole diagram. Can
    be the name of the element in the survey file also.
    :type offsetRotoTranslation: None, str, pymadx.Data.RotoTranslation2D
    :param title: title of the plot
    :type title: str
    :param outputfilename: output filename to save the plot to
    :type outputfilename: str
    :param zOffset: z-offset when drawing in matplotlib. Higher is more to the front. Use in single integer steps.
    :type zOffset: int
    :param defaultWidth: default width of each magnet in metres.
    :type defaultWidth: float
    :param elementDict: dictionary of element_name to dict of parameters
    :type elementDict: dict(str, dict(str, str|float|int))
    :param typeDict: dictionary of element_name to dict of parameters for a given type, e.g. QUADRUPOLE
    :type typeDict: dict(str, dict(str, str|float|int))
    :param funcDict: dictionary of element_name to function to draw it
    :type funcDict: dict(str, function)
    :param maskNames: list of elements to grey out as grey with alpha = 0.1
    :type maskNames: list(str)
    :param ignoreNames: list of elements to ignore and not draw at all
    :type ignoreNames: list(str)
    :param resolution: sampling distance along curved arcs in metres
    :type resolution: float
    :param defaultCoilLength: default coil length along S in metres
    :type defaultCoilLength: float
    :param pipeRadius: pipe radius (i.e. half-width) in metres
    :type pipeRadius: None, float
    :param pipeMaskRanges: list of tuples of (start, end) to hide the beam pipe
    :type pipeMaskRanges: None, list(tuple(float, float))
    :param invisibleAxes: whether to hide the bounding box of the matplotlib axes
    :type invisibleAxes: bool
    :param arrowsDy: fraction of figure height to move the orientation arrows (positive is up in the plot)
    :type arrowsDy: float
    :param arrowsDx: fraction of figure width to move the orientation arrows (positive is down in the plot)
    :type arrowsDx: float
    :param defaultAlpha: default alpha value to use for the whole beamline if not otherwise specified per element
    :type defaultAlpha: None, float


    Order of precedence is: funcDict, elementDict, typeDict. i.e. funcDict replaces any information given
    in the other two. An elementDict always overrides any typeDict information.
    """
    survey = _Data.CheckItsTfs(survey_tfsfile)

    # need these list / dicts but can't have a mutable type as a default argument
    def _InitialiseList(var):
        return [] if var is None else var
    def _InitialiseDict(var):
        return {} if var is None else var
    elementDict = _InitialiseDict(elementDict)
    typeDict = _InitialiseDict(typeDict)
    funcDict = _InitialiseDict(funcDict)
    maskNames = _InitialiseList(maskNames)
    ignoreNames = _InitialiseList(ignoreNames)
    greyColour = u'#c0c0c0'

    if pipeMaskRanges is None:
        pipeMaskRanges = [[survey[0]['S']-100, survey[0]['S']-99]]
    buildPipes = pipeRadius is not None
    _pipeOn = True

    if offsetRotoTranslation is not None:
        if isinstance(offsetRotoTranslation, str):
            rt = survey.GetRotoTranslationFromElementZX(offsetRotoTranslation).Inverse()
        elif isinstance(offsetRotoTranslation, _Data.RotoTranslation2D):
            rt = offsetRotoTranslation.Inverse()
        else:
            raise ValueError("offsetRotoTranslation must be None, str or a RotoTranslation2D")
    else:
        rt = _Data.RotoTranslation2D()

    if not ax:
        f = _plt.figure()
        ax = f.add_subplot(111)
    else:
        f = ax.get_figure()

    X, Z, ang = 0, 0, 0
    XZdir = _np.array([0,1])
    axisLine = []

    # loop over elements and prepare patches
    # patches are turned into patch collection which is more efficient later
    quads, bends, kickers, collimators, sextupoles = [], [], [], [], []
    octupoles, multipoles, solenoids, other, coils = [], [], [], [], []
    pipes = []

    _defaults = _defaultdict(lambda: r"#cccccc", defaultElementColours)  # default is grey

    for i, e in enumerate(survey):
        if i == 0:
            X = e['X']
            Z = e['Z']
            ang = e['THETA']
            axisLine.append([Z,X])
        l = e['L']
        name = e['NAME']
        if "BEGIN.VAC" in name or "BEG.VAC" in name:
            _pipeOn = True
        elif "END.VAC" in name:
            _pipeOn = False

        if l == 0:
            continue

        angle = e['ANGLE']
        insideFactor = -1 * _np.sign(angle)
        # draw line from previous point until now
        Xend, Zend = e['X'], e['Z']
        pt = [[Zend, Xend]] if angle == 0 else _CurvedLine(Z, X, XZdir, angle, l, resolution)
        axisLine.extend(pt)

        if name in ignoreNames:
            continue

        alpha = 0.1 if name in maskNames else 1.0
        if defaultAlpha != None:
            alpha = defaultAlpha
        # tolerate very slight tilts, e.g. earth's curvature corrections
        if 'TILT' in e:
            vertical = abs(e['TILT'] - 0.5*_np.pi) < 0.02 * _np.pi
        else:
            vertical = False

        # draw element
        kw = e['KEYWORD']
        params = {'colour' : _defaults[kw],
                  'width' : defaultWidth,
                  'height' : defaultWidth,
                  'dx' : 0, # curvilinear x
                  'dy' : 0, # curvilinear y
                  'coil_length' : defaultCoilLength,
                  'coil_width' : 0.7*defaultWidth,
                  'coil_height' : 0,
                  'coil_dx' : 0, # curvilinear x
                  'coil_dy' : 0, # curvilinear x
                  'coil_edge' : 0,
                  'inside' : True,
                  'style' : 'normal',
                  }
        params = _UpdateParams(elementDict, typeDict, e, params, insideFactor)
        #params['coil_dx']# *= inside
        c = params['colour']
        w = params['width']
        h = params['height']
        th = e['THETA']
        dx = params['dx']
        dy = params['dy']
        coils_this_mag = []
        if vertical:
            dx, dy = dy, dx

        # for everything
        if greyOut:
            c = greyColour
            alpha = 1.0

        if name in funcDict:
            # delegate function gives back list of polygons as x,y coords in plot
            edges = funcDict[name](locals())
            edges = rt * edges
            return _patches.Polygon(edges, color=c, fill=True, alpha=alpha)
        elif kw == 'DRIFT':
            # don't deal with pole faces - just put behind everything
            if _pipeOn and buildPipes:
                ignore = any([a <= e['S'] <= b for (a,b) in pipeMaskRanges])
                if not ignore:
                    pipes.append(_Rectangle(Zend, Xend, pipeRadius, l, th, rt, c, alpha))
        elif kw == 'QUADRUPOLE':
            quads.append(_Rectangle(Zend, Xend, w, l, th, rt, c, alpha, dx, dy))
            coils_this_mag.extend(_CoilPolygonsQuad(Zend, Xend, l, th, rt, alpha, params))
        elif kw == 'RBEND':
            a = 0 if vertical else 0.5*angle
            ang = th+0.5*a
            bends.append(_Rectangle(Zend, Xend, w, l, ang, rt, c, alpha, dx, dy))
            coils_this_mag.extend(_CoilPolygonsDipoleH(Zend, Xend, l, ang, rt, dx, params, alpha=alpha, greyOut=greyOut))
        # elif kw == 'SBEND':
        #     bends.append(DrawBend(element, c, alpha, dx, dy)) #blue
        #     coils_this_mag.extend(_CoilPolygonsDipoleH(Zend, Xend, l, th, alpha, dx, params))
        elif kw in ['HKICKER', 'VKICKER', 'KICKER']:
            kickers.append(_Rectangle(Zend, Xend, w, l, th, rt, c, alpha, dx, dy))
            coils_this_mag.extend(_CoilPolygonsDipoleH(Zend, Xend, l, th, rt, dx, params, alpha=alpha, greyOut=greyOut))
        elif kw == 'SOLENOID':
            solenoids.append(_Rectangle(Zend, Xend, w, l, e['THETA'], rt, c, alpha, dx, dy))
        elif kw in ['RCOLLIMATOR', 'ECOLLIMATOR', 'COLLIMATOR']:
            if params['style'] == 'fancy':
                cc = c if greyOut else u'#606060'
                collimators.append(_Collimator(Zend, Xend, w, l, e['THETA'], rt, cc, alpha))
            else:
                collimators.append(_Rectangle(Zend, Xend, w, l, e['THETA'], rt, c, alpha))
        elif kw == 'SEXTUPOLE':
            sextupoles.append(_Rectangle(Zend, Xend, w, l, e['THETA'], rt, c, alpha, dx, dy)) #yellow
        elif kw == 'OCTUPOLE':
            octupoles.append(_Rectangle(Zend, Xend, w, l, e['THETA'], rt, c, alpha, dx, dy)) #green
        else:
            #unknown so make light in alpha
            if l > 0.1:
                other.append(_Rectangle(Zend, Xend, w, l, e['THETA'], rt, c, alpha)) #light grey

        if len(coils_this_mag) > 0:
            coils.extend(coils_this_mag)

        # update the coords at the incoming end for the next loop iteration
        X = Xend
        Z = Zend

    axisLine.extend([[Z, X]])

    zo = zOffset * 30
    ax.add_collection(_PatchCollection(bends, match_original=True, zorder=20+zo))
    ax.add_collection(_PatchCollection(quads, match_original=True, zorder=19+zo))
    ax.add_collection(_PatchCollection(kickers, match_original=True, zorder=18+zo))
    ax.add_collection(_PatchCollection(collimators, match_original=True, zorder=16+zo))
    ax.add_collection(_PatchCollection(sextupoles, match_original=True, zorder=15+zo))
    ax.add_collection(_PatchCollection(octupoles, match_original=True, zorder=14+zo, edgecolor=None))
    ax.add_collection(_PatchCollection(multipoles, match_original=True, zorder=13+zo, edgecolor=None))
    ax.add_collection(_PatchCollection(other, match_original=True, zorder=12+zo, edgecolor=None))
    ax.add_collection(_PatchCollection(solenoids, match_original=True, zorder=11+zo))
    ax.add_collection(_PatchCollection(coils, match_original=True, zorder=10+zo))
    ax.add_collection(_PatchCollection(pipes, match_original=True, zorder=9+zo))

    # axis of machine over the top always
    axisLine = rt * _np.array(axisLine)
    ax.plot(axisLine[:, 0], axisLine[:, 1], c='k', zorder=21+zo, alpha=0.5, lw=1)

    _plt.suptitle(title, size='x-large')
    _plt.xlabel('Z (m)')
    _plt.ylabel('X (m)')
    _plt.tight_layout()

    if invisibleAxes:
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # orientation arrows
        fratio = f.get_size_inches()
        xOverY = fratio[0]/fratio[1]
        w = 0.003
        ys = arrowsDy
        xs = arrowsDx
        _plt.arrow(0.02-0.5*w+xs, 0.02+ys, 0.07/xOverY, 0, width=w*xOverY, head_width=3*w*xOverY, head_length=1*w*xOverY, transform=ax.transAxes, color='k', capstyle='butt')
        _plt.text(0.02+(0.13/xOverY)+xs, 0.02+ys, 'Z (m)', transform=ax.transAxes)
        _plt.arrow(0.02+xs, 0.02-(0.5*w*xOverY)+ys, 0, 0.07, width=w, head_width=3*w, head_length=3*w*xOverY, transform=ax.transAxes, color='k', capstyle='butt')
        _plt.text(0.02+xs, 0.13+ys, 'X (m)', transform=ax.transAxes)

    if outputfilename is not None:
        if 'png' in outputfilename:
            _plt.savefig(outputfilename, dpi=500)
        else:
            _plt.savefig(outputfilename)

    return f,ax
