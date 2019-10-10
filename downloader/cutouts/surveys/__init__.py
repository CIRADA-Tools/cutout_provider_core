# abstract base class components
from .survey_abc import processing_status

from .vlass import VLASS

def get_sexadecimal_string(position): 
    from .survey_abc import SurveyABC
    return SurveyABC.get_sexadecimal_string(position)
