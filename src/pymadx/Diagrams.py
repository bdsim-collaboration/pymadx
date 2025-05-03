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
    Apply a rotation and then add a translation to all points. See Rotate()

    :param translation: [x,y] to add to all points after rotation
    :type translation: np.array([x,y])

    """
    rotated = Rotate(points, angle, rotationOrigin)
    rotoTranslated = rotated + translation
    return rotoTranslated

def Polygon(edges, colour, alpha):
    """
    Return a polygon patch from a list of points as a numpy array.

    :param edges: numpy array of [x,y] edges
    :type edges: numpy.array([x1,y1], [x2,y2],...)
    :param colour: colour of the polygon
    :type colour: str
    :param alpha: transparency of the polygon
    :type alpha: float

    """
    return _patches.Polygon(edges, colour=colour, fill=True, alpha=alpha)

def _CoilPolygonsDipoleZX(xend, yend, length, rotation, roto_translation, dx=0, coil_dict=None, colour="#b87333",
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

def _CoilPolygonsDipoleZY(xend, yend, length, rotation, roto_translation, dy=0, coil_dict=None, colour="#b87333",
                          alpha=1.0, greyOut=False, greyColour='#c0c0c0'):
    if greyOut:
        colour = greyColour
    cl = coil_dict['coil_length']
    if cl == 0:
        return []
    ch = coil_dict['coil_height']
    cdy = coil_dict['coil_dy']
    global_offset = _np.array([xend, yend])
    gap = 0.3*ch
    hg = 0.5*gap
    hh = 0.5*ch
    l = length

    edges_top_out = _np.array([[0, hg], [0, hh], [cl, hh], [cl, hg]])
    edges_bot_out = _np.array([[0, -hg], [0, -hh], [cl, -hh], [cl, -hg]])
    edges_top_in  = _np.array([[-l, hg], [-l, hh], [-l-cl, hh], [-l-cl, hg]])
    edges_bot_in  = _np.array([[-l, -hg], [-l, -hh], [-l-cl, -hh], [-l-cl, -hg]])

    def tf(points):
        points[:,1] += cdy
        p = Rotate(points, rotation) + global_offset
        p = roto_translation * p
        return p

    return [_patches.Polygon(tf(edges_top_out), color=colour, fill=True, alpha=alpha),
            _patches.Polygon(tf(edges_bot_out), color=colour, fill=True, alpha=alpha),
            _patches.Polygon(tf(edges_top_in), color=colour, fill=True, alpha=alpha),
            _patches.Polygon(tf(edges_bot_in), color=colour, fill=True, alpha=alpha)]

def _CoilPolygonsQuad(xend, yend, length, rotation, roto_translation, alpha, coil_length, y_height, colour="#b87333",
                      greyOut=False, greyColour='#c0c0c0'):
    if greyOut:
        colour = greyColour
    cl = coil_length
    if cl == 0:
        return []
    cw = y_height
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

def _CoilPolygonsQuadZX(xend, yend, length, rotation, roto_translation, coil_dict, colour="#b87333", alpha=1.0,
                        greyOut=False, greyColour='#c0c0c0'):
    cl = coil_dict['coil_length']
    cw = coil_dict['coil_width']
    return _CoilPolygonsQuad(xend, yend, length, rotation, roto_translation, alpha, cl, cw, colour, greyOut, greyColour)

def _CoilPolygonsQuadZY(xend, yend, length, rotation, roto_translation, coil_dict, colour="#b87333", alpha=1.0,
                        greyOut=False, greyColour='#c0c0c0'):
    cl = coil_dict['coil_length']
    cw = coil_dict['coil_height']
    return _CoilPolygonsQuad(xend, yend, length, rotation, roto_translation, alpha, cl, cw, colour, greyOut, greyColour)

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


class _SurveyDiagram:
    def __init__(self, survey_tfsfile, ax=None, greyOut=False, offsetRotoTranslation=None, title='', outputfilename=None,
                 zOffset=0, defaultWidth=0.5, elementDict=None, typeDict=None, funcDict=None, maskNames=None,
                 ignoreNames=None, resolution=0.1, defaultCoilLength=0.15, pipeRadius=None, pipeMaskRanges=None,
                 invisibleAxes=False, arrowsDy=0, arrowsDx=0, defaultAlpha=None):
        self.survey = _Data.CheckItsTfs(survey_tfsfile)

        if not ax:
            self.f = _plt.figure()
            self.ax = self.f.add_subplot(111)
        else:
            self.ax = ax
            self.f = ax.get_figure()

        self.greyOut = greyOut
        self.offsetRotoTranslation = offsetRotoTranslation
        self.title = title
        self.outputfilename = outputfilename
        self.zOffset = zOffset
        self.defaultWidth = defaultWidth

        # need these list / dicts but can't have a mutable type as a default argument
        def _InitialiseList(var):
            return [] if var is None else var
        def _InitialiseDict(var):
            return {} if var is None else var
        self.elementDict = _InitialiseDict(elementDict)
        self.typeDict = _InitialiseDict(typeDict)
        self.funcDict = _InitialiseDict(funcDict)
        self.maskNames = _InitialiseList(maskNames)
        self.ignoreNames = _InitialiseList(ignoreNames)

        self.resolution = resolution
        self.defaultCoilLength = defaultCoilLength
        self.pipeRadius = pipeRadius

        if pipeMaskRanges is None:
            self.pipeMaskRanges = [[self.survey[0]['S'] - 100, self.survey[0]['S'] - 99]]
        else:
            self.pipeMaskRanges = pipeMaskRanges
        self.buildPipes = pipeRadius is not None
        self._pipeOn = True

        self.invisibleAxes = invisibleAxes
        self.arrowsDy = arrowsDy
        self.arrowsDx = arrowsDx
        self.defaultAlpha = defaultAlpha

        self.greyColour = u'#c0c0c0'

        self._SetupRotoTranslation()
        self.vars = self._PlotVariables()

        self._plane_global_angle = 0
        self._flip_geometry_parameters = False

    def _GetRTFromElement(self, offrot):
        raise NotImplementedError()

    def _SetupRotoTranslation(self):
        if self.offsetRotoTranslation is not None:
            if isinstance(self.offsetRotoTranslation, str):
                self.rt = self._GetRTFromElement(self.offsetRotoTranslation).Inverse()
            elif isinstance(self.offsetRotoTranslation, _Data.RotoTranslation2D):
                self.rt = self.offsetRotoTranslation.Inverse()
            else:
                raise ValueError("offsetRotoTranslation must be None, str or a RotoTranslation2D")
        else:
            self.rt = _Data.RotoTranslation2D()

    def _GetCoilsDipole(self, bends_in_plot_plane, x_end, y_end, l, xy_ang, rt, params, alpha):
        raise NotImplementedError()

    def _GetCoilsQuadrupole(self, xend, yend, l, xy_ang, rt, params, alpha):
        raise NotImplementedError()

    def _PlotVariables(self):
        raise NotImplementedError()

    def ApplyPlaneFlip(self, params):
        if not self._flip_geometry_parameters:
            return params
        p_new = dict(params)
        p_new['width'] = params['height']
        p_new['height'] = params['width']
        p_new['dx'] = params['dy']
        p_new['dy'] = params['dx']
        p_new['coil_width'] = params['coil_height']
        p_new['coil_height'] = params['coil_width']
        p_new['coil_dx'] = params['coil_dy']
        p_new['coil_dy'] = params['coil_dx']
        return p_new

    def Plot(self):
        # coordinates in the plot
        x, y, ang = 0, 0, 0
        xy_dir = _np.array([0,1])
        axisLine = []

        # loop over elements and prepare patches
        # patches are turned into patch collection which is more efficient later
        quads, bends, kickers, collimators, sextupoles = [], [], [], [], []
        octupoles, multipoles, solenoids, other, coils = [], [], [], [], []
        pipes = []

        _defaults = _defaultdict(lambda: r"#cccccc", defaultElementColours)  # default is grey

        rt = self.rt  # shortcut
        for i, e in enumerate(self.survey):
            if i == 0:
                x = e[self.vars['x']]
                y = e[self.vars['y']]
                #xy_ang = e[self.vars['angle']]
                axisLine.append([x, y])
            l = e['L']
            name = e['NAME']
            if "BEGIN.VAC" in name or "BEG.VAC" in name:
                _pipeOn = True
            elif "END.VAC" in name:
                _pipeOn = False

            if l == 0:
                continue

            angle = e['ANGLE']

            # tolerate very slight tilts, e.g. earth's curvature corrections
            # tilt may not be in the survey file... then have to assume we don't know
            tilt = e['TILT'] if 'TILT' in e else 0
            bends_in_plot_plane = (abs(tilt - self._plane_global_angle) < 0.02 * _np.pi)
            # change in angle this element induces in this plot plane
            d_angle = angle if bends_in_plot_plane else 0

            insideFactor = -1 * _np.sign(d_angle)
            # draw line from previous point until now
            x_end, y_end,  = e[self.vars['x']], e[self.vars['y']]
            pt = [[x_end, y_end]] if d_angle == 0 else _CurvedLine(x, y, xy_dir, d_angle, l, self.resolution)
            axisLine.extend(pt)

            if name in self.ignoreNames:
                continue

            alpha = 0.1 if name in self.maskNames else 1.0
            if self.defaultAlpha != None:
                alpha = self.defaultAlpha

            # draw element
            kw = e['KEYWORD']
            params = {'colour' : _defaults[kw],
                      'width' : self.defaultWidth,
                      'height' : self.defaultWidth,
                      'dx' : 0, # curvilinear x
                      'dy' : 0, # curvilinear y
                      'coil_length' : self.defaultCoilLength,
                      'coil_width' : 0.7*self.defaultWidth,
                      'coil_height' : 0.7*self.defaultWidth,
                      'coil_dx' : 0, # curvilinear x
                      'coil_dy' : 0, # curvilinear x
                      'coil_edge' : 0,
                      'inside' : True,
                      'style' : 'normal',
                      }
            params = _UpdateParams(self.elementDict, self.typeDict, e, params, insideFactor)
            params = self.ApplyPlaneFlip(params)
            c = params['colour']
            w = params['width']
            h = params['height']
            xy_ang = e[self.vars['angle']]
            dx = params['dx']
            dy = params['dy']
            coils_this_mag = []

            # for everything
            if self.greyOut:
                c = self.greyColour
                alpha = 1.0

            if name in self.funcDict:
                # delegate function gives back list of polygons as x,y coords in plot
                edges = self.funcDict[name](locals())
                edges = self.rt * edges
                return _patches.Polygon(edges, color=c, fill=True, alpha=alpha)
            elif kw == 'DRIFT':
                # don't deal with pole faces - just put behind everything
                if self._pipeOn and self.buildPipes:
                    ignore = any([a <= e['S'] <= b for (a,b) in self.pipeMaskRanges])
                    if not ignore:
                        pipes.append(_Rectangle(x_end, y_end, self.pipeRadius, l, xy_ang, rt, c, alpha))
            elif kw == 'QUADRUPOLE':
                quads.append(_Rectangle(x_end, y_end, w, l, xy_ang, rt, c, alpha, dx, dy))
                coils_this_mag.extend(self._GetCoilsQuadrupole(x_end, y_end, l, xy_ang, rt, params, alpha))
            elif kw == 'RBEND':
                ang = xy_ang + 0.5*d_angle
                bends.append(_Rectangle(x_end, y_end, w, l, ang, rt, c, alpha, dx, dy))
                coils_this_mag.extend(self._GetCoilsDipole(bends_in_plot_plane, x_end, y_end, l, ang, rt, params, alpha))
            # elif kw == 'SBEND':
            #     bends.append(DrawBend(element, c, alpha, dx, dy)) #blue
            #     coils_this_mag.extend(_CoilPolygonsDipoleH(Zend, Xend, l, th, alpha, dx, params))
            elif kw in ['HKICKER', 'VKICKER', 'KICKER']:
                kickers.append(_Rectangle(x_end, y_end, w, l, xy_ang, rt, c, alpha, dx, dy))
                coils_this_mag.extend(self._GetCoilsDipole(bends_in_plot_plane, x_end, y_end, l, ang, rt, params, alpha))
            elif kw == 'SOLENOID':
                solenoids.append(_Rectangle(x_end, y_end, w, l, xy_ang, rt, c, alpha, dx, dy))
            elif kw in ['RCOLLIMATOR', 'ECOLLIMATOR', 'COLLIMATOR']:
                if params['style'] == 'fancy':
                    cc = c if self.greyOut else u'#606060'
                    collimators.append(_Collimator(x_end, y_end, w, l, xy_ang, rt, cc, alpha))
                else:
                    collimators.append(_Rectangle(x_end, y_end, w, l, xy_ang, rt, c, alpha))
            elif kw == 'SEXTUPOLE':
                sextupoles.append(_Rectangle(x_end, y_end, w, l, xy_ang, rt, c, alpha, dx, dy)) #yellow
            elif kw == 'OCTUPOLE':
                octupoles.append(_Rectangle(x_end, y_end, w, l, xy_ang, rt, c, alpha, dx, dy)) #green
            else:
                #unknown so make light in alpha
                if l > 0.1:
                    other.append(_Rectangle(x_end, y_end, w, l, xy_ang, rt, c, alpha)) #light grey

            if len(coils_this_mag) > 0:
                coils.extend(coils_this_mag)

            # update the coords at the incoming end for the next loop iteration
            x = x_end
            y = y_end

        axisLine.extend([[x, y]])

        zo = self.zOffset * 30
        self.ax.add_collection(_PatchCollection(bends, match_original=True, zorder=20+zo))
        self.ax.add_collection(_PatchCollection(quads, match_original=True, zorder=19+zo))
        self.ax.add_collection(_PatchCollection(kickers, match_original=True, zorder=18+zo))
        self.ax.add_collection(_PatchCollection(collimators, match_original=True, zorder=16+zo))
        self.ax.add_collection(_PatchCollection(sextupoles, match_original=True, zorder=15+zo))
        self.ax.add_collection(_PatchCollection(octupoles, match_original=True, zorder=14+zo, edgecolor=None))
        self.ax.add_collection(_PatchCollection(multipoles, match_original=True, zorder=13+zo, edgecolor=None))
        self.ax.add_collection(_PatchCollection(other, match_original=True, zorder=12+zo, edgecolor=None))
        self.ax.add_collection(_PatchCollection(solenoids, match_original=True, zorder=11+zo))
        self.ax.add_collection(_PatchCollection(coils, match_original=True, zorder=10+zo))
        self.ax.add_collection(_PatchCollection(pipes, match_original=True, zorder=9+zo))

        # axis of machine over the top always
        axisLine = rt * _np.array(axisLine)
        self.ax.plot(axisLine[:, 0], axisLine[:, 1], c='k', zorder=21+zo, alpha=0.5, lw=1)

        _plt.suptitle(self.title, size='x-large')
        _plt.xlabel('Z (m)')
        _plt.ylabel('X (m)')
        _plt.tight_layout()

        if self.invisibleAxes:
            self.ax.get_xaxis().set_visible(False)
            self.ax.get_yaxis().set_visible(False)
            self.ax.spines['top'].set_visible(False)
            self.ax.spines['bottom'].set_visible(False)
            self.ax.spines['left'].set_visible(False)
            self.ax.spines['right'].set_visible(False)
            # orientation arrows
            fratio = self.f.get_size_inches()
            xOverY = fratio[0]/fratio[1]
            w = 0.003
            ys = self.arrowsDy
            xs = self.arrowsDx
            _plt.arrow(0.02-0.5*w+xs, 0.02+ys, 0.07/xOverY, 0, width=w*xOverY, head_width=3*w*xOverY,
                       head_length=1*w*xOverY, transform=self.ax.transAxes, color='k', capstyle='butt')
            _plt.text(0.02+(0.13/xOverY)+xs, 0.02+ys, 'Z (m)', transform=self.ax.transAxes)
            _plt.arrow(0.02+xs, 0.02-(0.5*w*xOverY)+ys, 0, 0.07, width=w, head_width=3*w,
                       head_length=3*w*xOverY, transform=self.ax.transAxes, color='k', capstyle='butt')
            _plt.text(0.02+xs, 0.13+ys, 'X (m)', transform=self.ax.transAxes)

        if self.outputfilename is not None:
            if 'png' in self.outputfilename:
                _plt.savefig(self.outputfilename, dpi=500)
            else:
                _plt.savefig(self.outputfilename)

        return self.f, self.ax


class _SurveyDiagramZX(_SurveyDiagram):
    def __init__(self, *args, **kwargs):
        super(_SurveyDiagramZX, self).__init__(*args, **kwargs)
        self._plane_global_angle = 0
        self._flip_geometry_parameters = False

    def _GetRTFromElement(self, offrot):
        return self.survey.GetRotoTranslationFromElementZX(self.offsetRotoTranslation)

    def _PlotVariables(self):
        return {'x' : 'Z', 'y' : 'X', 'angle' : 'THETA'}

    def _GetCoilsDipole(self, bends_in_plot_plane,  x_end, y_end, l, xy_ang, rt, params, alpha):
        dx = params['dx']
        if bends_in_plot_plane:
            return _CoilPolygonsDipoleZX(x_end, y_end, l, xy_ang, rt, dx, params, alpha=alpha,
                                         greyOut=self.greyOut, greyColour=self.greyColour)
        else:
            return _CoilPolygonsDipoleZY(x_end, y_end, l, xy_ang, rt, dx, params, alpha=alpha,
                                         greyOut=self.greyOut, greyColour=self.greyColour)

    def _GetCoilsQuadrupole(self, xend, yend, l, xy_ang, rt, params, alpha):
        return _CoilPolygonsQuadZY(xend, yend, l, xy_ang, rt, params, alpha=alpha,
                                   greyOut=self.greyOut, greyColour=self.greyColour)


class _SurveyDiagramZY(_SurveyDiagram):
    def __init__(self, *args, **kwargs):
        super(_SurveyDiagramZY, self).__init__(*args, **kwargs)
        self._plane_global_angle = _np.pi * 0.5
        self._flip_geometry_parameters = True

    def _GetRTFromElement(self, offrot):
        return self.survey.GetRotoTranslationFromElementZY(self.offsetRotoTranslation)

    def _PlotVariables(self):
        return {'x' : 'Z', 'y' : 'Y', 'angle' : 'PHI'}

    def _GetCoilsDipole(self, bends_in_plot_plane, x_end, y_end, l, xy_ang, rt, params, alpha):
        dy = params['dy']
        if bends_in_plot_plane:
            return _CoilPolygonsDipoleZX(x_end, y_end, l, xy_ang, rt, dy, params, alpha=alpha,
                                         greyOut=self.greyOut, greyColour=self.greyColour)
        else:
            return _CoilPolygonsDipoleZY(x_end, y_end, l, xy_ang, rt, dy, params, alpha=alpha,
                                         greyOut=self.greyOut, greyColour=self.greyColour)

    def _GetCoilsQuadrupole(self, xend, yend, l, xy_ang, rt, params, alpha):
        return _CoilPolygonsQuadZX(xend, yend, l, xy_ang, rt, params, alpha=alpha,
                                   greyOut=self.greyOut, greyColour=self.greyColour)


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
    :param offsetRotoTranslation: Optional roto-translation to apply before plotting to offset the whole diagram. Can be the name of the element in the survey file also.
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
    a = _SurveyDiagramZX(survey_tfsfile, ax, greyOut, offsetRotoTranslation, title, outputfilename,
                         zOffset, defaultWidth, elementDict, typeDict, funcDict, maskNames,
                         ignoreNames, resolution, defaultCoilLength, pipeRadius, pipeMaskRanges,
                         invisibleAxes, arrowsDy, arrowsDx, defaultAlpha)
    return a.Plot()


def Survey2DZY(survey_tfsfile, ax=None, greyOut=False, offsetRotoTranslation=None, title='', outputfilename=None,
               zOffset=0, defaultWidth=0.5, elementDict=None, typeDict=None, funcDict=None, maskNames=None,
               ignoreNames=None, resolution=0.1, defaultCoilLength=0.15, pipeRadius=None, pipeMaskRanges=None,
               invisibleAxes=False, arrowsDy=0, arrowsDx=0, defaultAlpha=None):
    """
    Draw a schematic of a beamline in the Z-Y plane. Allow styling of different elements by type or by name.

    :param survey_tfsfile: list of tfs files as strings or already loaded pymadx.Data.Tfs objects.
    :type survey_tfsfile: str, pymadx.Data.Tfs
    :param ax: matplotlib axes object
    :type ax: matplotlib.axes.Axes, None
    :param greyOut: whether to plot the whole line in greyscale
    :type greyOut: bool
    :param offsetRotoTranslation: Optional roto-translation to apply before plotting to offset the whole diagram. Can be the name of the element in the survey file also.
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
    a = _SurveyDiagramZY(survey_tfsfile, ax, greyOut, offsetRotoTranslation, title, outputfilename,
                         zOffset, defaultWidth, elementDict, typeDict, funcDict, maskNames,
                         ignoreNames, resolution, defaultCoilLength, pipeRadius, pipeMaskRanges,
                         invisibleAxes, arrowsDy, arrowsDx, defaultAlpha)
    return a.Plot()