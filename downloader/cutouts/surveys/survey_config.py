#system
import sys

# configuration
import csv
import yaml as yml

# astropy libs
from astropy import units as u
from astropy.coordinates import SkyCoord

# filters used by various surveys
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

        # set the filters
        self.survey_filter_sets = list()
        for supported_survey in self.supported_surveys:
            survey_filters  = eval(f"{supported_survey}.get_filters()")
            if survey_filters:
                self.survey_filter_sets.append({supported_survey: [f for f in survey_filters]})

        # set survey_names
        self.survey_block = self.config['cutouts']['surveys']
        self.survey_names = self.__extract_surveys_names(self.survey_block)

        # set the cutout size
        self.size_arcmin = self.config['cutouts']['box_size_armin'] * u.arcmin

        # set targets
        self.targets = list()
        for coords_csv_file in self.config['cutouts']['ra_deg_deg_csv_files']:
            sources = self.__csv_to_dict(coords_csv_file)
    
            # make all keys lower case to effectively reference case-variants of RA and Dec.
            sources = [{k.lower(): v for k,v in s.items()} for s in sources]
    
            # extract position information
            self.targets.extend([{
                'coord': SkyCoord(x['ra'], x['dec'], unit=(u.deg, u.deg)),
                'size':  self.size_arcmin
            } for x in sources])


    def __csv_to_dict(self,filename):
        entries = []
    
        with open(filename, 'r') as infile:
            c = csv.DictReader(infile)
            for entry in c:
                entries.append(entry)
    
        return entries

    def get_targets(self):
        return self.targets


    def __get_target_list(self,config_file):
        config  = yml.load(open(config_file,'r'))['cutouts']
    
        targets = list()
        size = config['box_size_armin'] * u.arcmin
        for coord_csv_file in config['ra_deg_deg_csv_files']:
            sources = csv_to_dict(coord_csv_file)
    
            # make all keys lower case to effectively reference case-variants of RA and Dec.
            sources = [{k.lower(): v for k,v in s.items()} for s in sources]
    
            # extract position information
            targets.extend([
                {
                    'coord': SkyCoord(x['ra'], x['dec'], unit=(u.deg, u.deg)),
                    'size': size
                }
                for x in sources])
    
        return targets


    def __print(self,string,show_caller=False):
        prefix = type(self).__name__ + (f"({sys._getframe(1).f_code.co_name})" if show_caller else "")
        prefixed_string = "\n".join([f"{prefix}: {s}" for s in string.splitlines()])
        print(prefixed_string)


    def __extract_surveys_names(self,config_surveys_block):
        def extract_surveys_names_from_config_surveys_block(config_surveys_block):
            survey_names = list()
            for survey in config_surveys_block:
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
                   if isinstance(survey_parameters,dict) and \
                        ('filters' in survey_parameters.keys()) and \
                        (len(self.__match_filters(name,survey_parameters['filters'])) > 0):
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


    def get_processing_stack(self):
        # TODO: Fix the repeating message,
        #
        #         > In [354]: cfg = SurveyConfig("config_debug.yml") 
        #         > ...
        #         > In [355]: cfg.get_processing_stack()
        #         > 
        #         > SurveyConfig: INSTANTIATING: FIRST()
        #         > SurveyConfig: INSTANTIATING: NVSS()
        #         > SurveyConfig: WARNING: 'VLASS' filter 'foo' is not supported!
        #         > SurveyConfig: INSTANTIATING: VLASS()
        #         > VLASS: => Using CADC cutout server!
        #         > SurveyConfig: INSTANTIATING: WISE(filter=wise_filters.w1)
        #         > SurveyConfig: WARNING: 'PanSTARRS' filter 'k' is not supported!
        #         > SurveyConfig: WARNING: 'PanSTARRS' filter 'k' is not supported!
        #         > SurveyConfig: WARNING: 'PanSTARRS' filter 'k' is not supported!
        #         > SurveyConfig: INSTANTIATING: PanSTARRS(filter=grizy_filters.i)
        #         > SurveyConfig: INSTANTIATING: SDSS(filter=grizy_filters.g)
        #         > SurveyConfig: INSTANTIATING: SDSS(filter=grizy_filters.r)
        #         > SurveyConfig: INSTANTIATING: SDSS(filter=grizy_filters.i)
        #         > Out[355]:
        #         > [<surveys.first.FIRST at 0x11a76acc0>,
        #         >  <surveys.nvss.NVSS at 0x11a76af28>,
        #         >  <surveys.vlass.VLASS at 0x11a76a860>,
        #         >  <surveys.wise.WISE at 0x11a76ac50>,
        #         >  <surveys.panstarrs.PanSTARRS at 0x11a76af98>,
        #         >  <surveys.sdss.SDSS at 0x11a76ada0>,
        #         >  <surveys.sdss.SDSS at 0x11a76a7b8>,
        #         >  <surveys.sdss.SDSS at 0x11a76a9b0>]
        #         > 
        #         > In [356]:
        #
        #       problem.
        processing_stack = list()
        for survey_name in self.get_survey_names():
            if self.has_filters(survey_name):
                for filter in self.get_filters(survey_name):
                    self.__print(f"INSTANTIATING: {survey_name}(filter={filter})")
                    processing_stack.append(eval(f"{survey_name}(filter={filter})"))
            else:
                self.__print(f"INSTANTIATING: {survey_name}()")
                processing_stack.append(eval(f"{survey_name}()"))
        return processing_stack 



