import pymadx
import os

import matplotlib.pyplot as plt

def _fn(filename):
    return os.path.join(os.path.dirname(__file__), "test_input", filename)

def test_survey2d_ZX():
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    pymadx.Diagrams.Survey2DZX(s)

def test_survey2d_ZX_offset():
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    offset = s.GetRotoTranslationFromElementZX("XWCA.X0410404")
    pymadx.Diagrams.Survey2DZX(s, offsetRotoTranslation=offset)

def test_survey2d_ZX_offset_by_name():
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    pymadx.Diagrams.Survey2DZX(s, offsetRotoTranslation="XWCA.X0410404")

def test_survey2d_ZX_grey_out():
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    pymadx.Diagrams.Survey2DZX(s, greyOut=True)

def test_survey2d_ZX_mask_names():
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    maskNames = ['MTN.X0400003', 'MTN.X0400007']
    pymadx.Diagrams.Survey2DZX(s, maskNames=maskNames)

def test_survey2d_ZX_fancy():
    #funcDict = {'MSN' : pymadx.Plot.MSNPatches}
    elemDict = {'MBXHC.X0410117' : {'inside': False},
                'MBXHC.X0410121' : {'inside': False},
                'MBXHC.X0410124' : {'inside': False},
                'MBXHC.X0410132' : {'inside': False},
                'MBXHC.X0410135' : {'inside': False},
                'MBXHC.X0410139' : {'inside': False},
                }
    typeDict = {'MCA' : {'width': 1.25, 'height': 1.25, 'dx': 0.14, 'colour': r'#a7d9b0', 'coil_length': 0.31,
                         'coil_width': 0.8, 'coil_dx': -0.14},
                'MSN' : {'width': 0.47, 'height': 0.3, 'dx': 0.11, 'coil_length': 0.1, 'coil_width': 0.25, 'coil_dx': -0.11},
                'MTN' : {'width': 1.2, 'height': 0.69, 'colour':r'#b2d6c0', 'coil_length': 0.1, 'coil_width': 0.6},
                'MBXHC' : {'width': 1.246, 'height': 1.25, 'colour': r'#a7d9b0', 'dx': 0.141,
                           'coil_length': 0.3, 'coil_width':0.9, 'coil_edge': 0.3, 'coil_dx': -0.141},
                'MBNV': {'width': 0.6, 'height': 1.1, 'colour': r'#f29010', 'coil_length': 0.22, 'coil_width': 0.6},
                'MBNH': {'width': 1.1, 'height': 0.6, 'colour': r'#f29010', 'coil_length': 0.22, 'coil_width': 0.64},
                'MBW' : {'width': 0.88, 'height': 0.47, 'colour': r'#d4340c', 'coil_length': 0.1, 'coil_width': 0.5},
                'MCWH' : {'width': 0.85, 'height': 1.0, 'colour': r'#a7d9b0', 'coil_length': 0.2, 'coil_width': 0.6},
                'MCWV': {'width': 1.0, 'height': 0.85, 'colour': r'#a7d9b0', 'coil_length': 0.2, 'coil_width': 0.6},
                'QNL' : {'width': 0.6, 'height': 0.8, 'colour': r'#fcd34c'},
                'QPL' : {'width': 1.1, 'height': 1.1, 'colour': r'#1e84eb', 'coil_length': 0.23, 'coil_width': 0.8},
                'QSL' : {'width': 0.32, 'height': 0.52},
                'QTS' : {'width': 0.6, 'height': 0.8, 'colour': r'#82bdb9'},
                'QWL' : {'width': 0.66, 'height': 0.84, 'colour': r'#d4340c', 'coil_length': 0.14, 'coil_width': 0.6},
                'LSX' : {'width':0.4, 'height':0.4, 'colour':r'#ddc8e9'},
                'XCHV' : {'style': 'fancy'}}
    typeDict['MCB'] = dict(typeDict['MCA']) # copy it
    typeDict['MBXGD'] = dict(typeDict['MCWV'])
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    maskRanges = [[9.2, 20.1]]

    f, ax = pymadx.Diagrams.Survey2DZX(s, typeDict=typeDict, elementDict=elemDict, pipeRadius=0.07, pipeMaskRanges=maskRanges)

    h8 = pymadx.Data.Tfs(_fn("h8-survey.tfs"))
    pymadx.Diagrams.Survey2DZX(h8, ax=ax, typeDict=typeDict)

    p42 = pymadx.Data.Tfs(_fn("p42-survey.tfs"))
    pymadx.Diagrams.Survey2DZX(p42, ax=ax, typeDict=typeDict)

def test_survey2d_ZY():
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    pymadx.Diagrams.Survey2DZY(s)

def test_survey2d_ZY_offset():
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    offset = s.GetRotoTranslationFromElementZY("XWCA.X0410404")
    pymadx.Diagrams.Survey2DZY(s, offsetRotoTranslation=offset)

def test_survey2d_ZY_offset_by_name():
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    pymadx.Diagrams.Survey2DZY(s, offsetRotoTranslation="XWCA.X0410404")

def test_survey2d_ZY_grey_out():
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    pymadx.Diagrams.Survey2DZY(s, greyOut=True)

def test_survey2d_ZY_mask_names():
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    maskNames = ['MTN.X0400003', 'MTN.X0400007']
    pymadx.Diagrams.Survey2DZY(s, maskNames=maskNames)

def test_survey2d_ZY_fancy():
    #funcDict = {'MSN' : pymadx.Plot.MSNPatches}
    elemDict = {'MBXHC.X0410117' : {'inside': False},
                'MBXHC.X0410121' : {'inside': False},
                'MBXHC.X0410124' : {'inside': False},
                'MBXHC.X0410132' : {'inside': False},
                'MBXHC.X0410135' : {'inside': False},
                'MBXHC.X0410139' : {'inside': False},
                }
    typeDict = {'MCA' : {'width': 1.25, 'height': 1.25, 'dx': 0.14, 'colour': r'#a7d9b0', 'coil_length': 0.31,
                         'coil_width': 0.8, 'coil_dx': -0.14},
                'MSN' : {'width': 0.47, 'height': 0.3, 'dx': 0.11, 'coil_length': 0.1, 'coil_width': 0.25, 'coil_dx': -0.11},
                'MTN' : {'width': 1.2, 'height': 0.69, 'colour':r'#b2d6c0', 'coil_length': 0.1, 'coil_width': 0.6},
                'MBXHC' : {'width': 1.246, 'height': 1.25, 'colour': r'#a7d9b0', 'dx': 0.141,
                           'coil_length': 0.3, 'coil_width':0.9, 'coil_edge': 0.3, 'coil_dx': -0.141},
                'MBNV': {'width': 0.6, 'height': 1.1, 'colour': r'#f29010', 'coil_length': 0.22, 'coil_width': 0.6, 'coil_height': 0.64},
                'MBNH': {'width': 1.1, 'height': 0.6, 'colour': r'#f29010', 'coil_length': 0.22, 'coil_width': 0.64}, 'coil_height': 0.6,
                'MBW' : {'width': 0.88, 'height': 0.47, 'colour': r'#d4340c', 'coil_length': 0.1, 'coil_width': 0.5},
                'MCWH' : {'width': 0.85, 'height': 1.0, 'colour': r'#a7d9b0', 'coil_length': 0.2, 'coil_width': 0.6},
                'MCWV': {'width': 1.0, 'height': 0.85, 'colour': r'#a7d9b0', 'coil_length': 0.2, 'coil_width': 0.6},
                'QNL' : {'width': 0.6, 'height': 0.8, 'colour': r'#fcd34c'},
                'QPL' : {'width': 1.1, 'height': 1.1, 'colour': r'#1e84eb', 'coil_length': 0.23, 'coil_width': 0.8, 'coil_height': 0.8},
                'QSL' : {'width': 0.32, 'height': 0.52},
                'QTS' : {'width': 0.6, 'height': 0.8, 'colour': r'#82bdb9'},
                'QWL' : {'width': 0.66, 'height': 0.84, 'colour': r'#d4340c', 'coil_length': 0.14, 'coil_width': 0.6},
                'LSX' : {'width':0.4, 'height':0.4, 'colour':r'#ddc8e9'},
                'XCHV' : {'style': 'fancy'}}
    typeDict['MCB'] = dict(typeDict['MCA']) # copy it
    typeDict['MBXGD'] = dict(typeDict['MCWV'])
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    maskRanges = [[9.2, 20.1]]

    f, ax = pymadx.Diagrams.Survey2DZY(s, typeDict=typeDict, elementDict=elemDict, pipeRadius=0.07, pipeMaskRanges=maskRanges)

    h8 = pymadx.Data.Tfs(_fn("h8-survey.tfs"))
    pymadx.Diagrams.Survey2DZY(h8, ax=ax, typeDict=typeDict)

    p42 = pymadx.Data.Tfs(_fn("p42-survey.tfs"))
    pymadx.Diagrams.Survey2DZY(p42, ax=ax, typeDict=typeDict)


# plots for manual
"""
pth = "../docs/source/figures/diagrams_"
test_survey2d_ZX()
plt.savefig(pth+"zx.png", dpi=400)
test_survey2d_ZX_offset()
plt.savefig(pth+"zx_offset.png", dpi=400)
test_survey2d_ZX_fancy()
plt.savefig(pth+"zx_fancy.png", dpi=400)
test_survey2d_ZX_grey_out()
plt.savefig(pth+"zx_grey.png", dpi=400)
test_survey2d_ZX_offset_by_name()

test_survey2d_ZY()
plt.savefig(pth+"zy.png", dpi=400)
test_survey2d_ZY_offset()
plt.savefig(pth+"zy_offset.png", dpi=400)
test_survey2d_ZY_fancy()
plt.savefig(pth+"zy_fancy.png", dpi=400)
test_survey2d_ZY_grey_out()
plt.savefig(pth+"zy_grey.png", dpi=400)
test_survey2d_ZY_offset_by_name()

h6 = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
h8 = pymadx.Data.Tfs(_fn("h8-survey.tfs"))
fig, ax = pymadx.Diagrams.Survey2DZX(h6)
pymadx.Diagrams.Survey2DZX(h8, ax=ax, zOffset=2)
plt.savefig(pth+"zx_multiple.png", dpi=400)
fig2, ax2 = pymadx.Diagrams.Survey2DZY(h6)
pymadx.Diagrams.Survey2DZY(h8, ax=ax2, zOffset=2)
plt.savefig(pth+"zy_multiple.png", dpi=400)

plt.show()
"""