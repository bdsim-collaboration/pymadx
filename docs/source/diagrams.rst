========
Diagrams
========

The :code:`pymadx.Diagrams` module provides machine diagrams. These use a MADX
survey output file (TFS format) and require all columns to be present. i.e. do
not select individual columns. For example, in MADX::

    select, flag=survey, clear;
    survey, file="survey-h6-from-zero.tfs";


Z-X Plane
---------

.. code-block::

    s = pymadx.Data.Tfs("survey_file.tfs")
    pymadx.Diagrams.Survey2DZX(s)


This will use the default colours and widths for the various components.

Offsetting
**********

An overall global offset and rotation can be specified in two ways. Firstly,
by naming an element. Secondly, by providing a roto-translation yourself, potentially
from another survey if combining machine diagrams.

.. code-block::

    s = pymadx.Data.Tfs("h6-survey.tfs")
    pymadx.Diagrams.Survey2DZX(s, offsetRotoTranslation="XWCA.X0410404")

.. code-block::

    s = pymadx.Data.Tfs("h6-survey.tfs")
    offset = s.GetRotoTranslationFromElementZX("XWCA.X0410404")
    pymadx.Diagrams.Survey2DZX(s, offsetRotoTranslation=offset)


Multiple Machines
-----------------

More than one machine diagram can be combined by providing a matplotlib axis
instance to draw into. The z-order can also be specified with a higher integer
number meaning in front.

.. code-block::

    h6 = pymadx.Data.Tfs("h6_survey_file.tfs")
    h8 = pymadx.Data.Tfs("h6_survey_file.tfs")
    fig, ax = pymadx.Diagrams.Survey2DZX(h6)
    pymadx.Diagrams.Survey2DZX(h8, ax=ax, zorder=2)

