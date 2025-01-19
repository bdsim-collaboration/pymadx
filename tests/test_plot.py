import os
import pymadx


def _fn(filename):
    return os.path.join(os.path.dirname(__file__), "test_input", filename)


def test_survey2d_ZX():
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    pymadx.Plot.Survey2DZX(s)


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
                'LSX' : {'width':0.4, 'height':0.4, 'colour':r'#ddc8e9'}}
    typeDict['MCB'] = dict(typeDict['MCA']) # copy it
    typeDict['MBXGD'] = dict(typeDict['MCWV'])
    oldmbn = r'#e39351'
    s = pymadx.Data.Tfs(_fn("h6-survey.tfs"))
    f, ax = pymadx.Plot.Survey2DZX(s, typeDict=typeDict, elementDict=elemDict)

    h8 = pymadx.Data.Tfs(_fn("h8-survey.tfs"))
    pymadx.Plot.Survey2DZX(h8, ax=ax, typeDict=typeDict)

    p42 = pymadx.Data.Tfs(_fn("p42-survey.tfs"))
    globalRotation = [4332.609, -638.9902819, 0.116276]
    pymadx.Plot.Survey2DZX(p42, ax=ax, typeDict=typeDict, globalRotation=globalRotation)

#test_survey2d_ZX()
test_survey2d_ZX_fancy()