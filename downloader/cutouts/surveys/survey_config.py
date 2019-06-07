import sys

import yaml as yml

from .survey_filters import wise_filters
from .survey_filters import grizy_filters

# supported suverys (nb: cf., SurveyConfig::self.supported_surveys)
from surveys.nvss      import NVSS
from surveys.first     import FIRST
from surveys.wise      import WISE
from surveys.sdss      import SDSS
from surveys.vlass     import VLASS
from surveys.panstarrs import PanSTARRS

class SurveyConfig:
    def __init__(self,yml_configuration_file):
        self.config = yml.load(open(yml_configuration_file,'r'))

        print(f"config: {self.config}")
        print()

        # define supported_surveys
        self.supported_surveys = (
            FIRST.__name__,
            NVSS.__name__,
            VLASS.__name__,
            WISE.__name__,
            PanSTARRS.__name__,
            SDSS.__name__,
            # TODO: Handle 2MASS case (i.e., number prefixed name -- python no like)
        )
        print(f"supported_surveys: {self.supported_surveys}")

        self.survey_filter_sets = list()
        for supported_survey in self.supported_surveys:
            survey_filters  = eval(f"{supported_survey}.get_filters()")
            if survey_filters:
                self.survey_filter_sets.append({supported_survey: [f for f in survey_filters]})
        print(f"survey_filters: {self.survey_filter_sets}")


        # set survey_names
        self.survey_block = self.config['cutouts']['surveys']
        self.survey_names = self.__extract_surveys_names(self.survey_block)

        print(f"get_survey_names: {self.get_survey_names()}")


    def __match_filters(self,survey,filters):
        def get_supported_filters(survey):
            filters = list()
            if self.has_survey(survey):
                for survey_filters in self.survey_filter_sets:
                    s = [k for k in survey_filters.keys()][0]
                    if s.lower() == survey.lower():
                      filters = survey_filters[s]
            return filters

        if isinstance(filters,str):
            filters = [filters]

        matched = list()
        for filter in filters:
            found = False
            for supported_filter in get_supported_filters(survey):
                if supported_filter.name.lower() == filter.lower():
                    found = True
                    break
            if found:
                matched.append(supported_filter)
                found = False
            else:
                self.__print(f"WARNING: '{survey}' filter '{filter}' is not supported!")
        
        return matched


    def __print(self,string,show_caller=False):
        prefix = type(self).__name__ + (f"({sys._getframe(1).f_code.co_name})" if show_caller else "")
        prefixed_string = "\n".join([f"{prefix}: {s}" for s in string.splitlines()])
        print(prefixed_string)


    def __extract_surveys_names(self,config_surveys_block):
        def extract_surveys_names_from_config_surveys_block(config_surveys_block):
            survey_names = list()
            for survey in config_surveys_block:
                print(f"survey: {survey}") # TODO: remove
                if isinstance(survey,str):
                    survey_names.append(survey)
                elif isinstance(survey,dict):
                    survey_names.append(list(survey.keys())[0])
                else:
                    # Whoops!
                    # TODO: Add output indicating a corrupt yml cfg file and 'raise'
                    pass
            return survey_names
        surveys = extract_surveys_names_from_config_surveys_block(config_surveys_block)
        supported = list()
        for survey in surveys:
            if self.__is_supported_survey(survey):
                for supported_survey in self.supported_surveys:
                    if supported_survey.lower() == survey.lower():
                        supported.append(supported_survey)
                        break
            else:
                self.__print(f"WARNING: Survey '{survey}' is not supported!")
        return supported


    def __is_supported_survey(self,survey):
        for supported_survey in self.supported_surveys:
            if survey.lower() == supported_survey.lower():
                return True
        return False


    def has_survey(self,survey):
        for survey_name in self.survey_names:
            if survey.lower() == survey_name.lower():
                return True
        self.__print(f"WARNING: Survey '{survey}' not found!")
        return False


    def has_filters(self,survey):
        if self.has_survey(survey):
            for s in self.survey_block:
                name = s if isinstance(s,str) else [k for k in s.keys()][0]
                if name.lower() == survey.lower() and isinstance(s,dict):
                   survey_parameters = s[name]
                   if isinstance(survey_parameters,dict) and ('filters' in survey_parameters.keys()):
                       return True
                   break
        return False


    def get_supported_survey(self):
        return self.supported_surveys


    def get_survey_names(self):
        return self.survey_names
    

    def get_filters(self,survey):
        filters = list()
        if self.has_filters(survey):
            for s in self.survey_block:
                name = s if isinstance(s,str) else [k for k in s.keys()][0]
                if name.lower() == survey.lower() and isinstance(s,dict):
                    filters = s[name]['filters']
                    if isinstance(filters,str):
                        filters = [filters.lower()]
                    elif isinstance(filters,list):
                        filters = [f.lower() for f in filters]
                    return self.__match_filters(survey,filters)
        return filters

