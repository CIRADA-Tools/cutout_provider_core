#system
import os
import sys

# utilities
import re

# add module paths that are two levels up from here
this_source_file_dir = re.sub(r"(.*/).*$",r"\1",os.path.realpath(__file__))
sys.path.append(this_source_file_dir+"../..")

# import the vospace space module to get the data-subdir configuration
from vospace.hierarchy import LocalCutoutDirs

from random import shuffle

# configuration
import csv
import yaml as yml

# processing
from .survey_abc import processing_status as ProcStatus
from .survey_abc import get_sexadecimal_string

# astropy libs
from astropy import units as u
from astropy.coordinates import SkyCoord

# filters used by various surveys
from .survey_filters import wise_filters
from .survey_filters import grizy_filters

# supported suverys (nb: cf., SurveyConfig::self.supported_surveys)
from .nvss      import NVSS
from .first     import FIRST
from .wise      import WISE
from .sdss      import SDSS
from .vlass     import VLASS
from .panstarrs import PanSTARRS

def get_cutout_filename(position,size,survey,filter=None,extension=None):
    coords = get_sexadecimal_string(position)
    size   = re.sub(r"\.?0+$","","%f" % size.to(u.arcmin).value)
    filter = (lambda f: '' if f is None else f"-{f.name}")(filter)
    return f"J{coords}_s{size}arcmin_{survey}{filter}{'.%s' % extension if not extension is None else ''}"

class SurveyConfig:
    def __init__(self,yml_configuration_file):
        # get config file
        self.config = yml.load(open(yml_configuration_file,'r'))

        # get relative path to config file
        relative_path = re.sub(r"[^/]+$","",yml_configuration_file)

        # get cutout dir hierarchy class
        self.local_dirs = LocalCutoutDirs()

        # file flush utility flag
        self.is_force_flush = False


        #
        # S U V E R Y   S E T U P
        #

        # define supported_surveys
        self.supported_surveys = (
            FIRST.__name__,
            NVSS.__name__,
            VLASS.__name__,
            WISE.__name__,
            PanSTARRS.__name__,
            SDSS.__name__,
            # TODO (Issue #13): Handle 2MASS case (i.e., number prefixed name -- python no like)
        )

        # make sure supported_surveys are defined in the hierarchy class
        for survey in self.supported_surveys:
            if not self.local_dirs.has_survey(survey):
                self.__print(f"WARNING: '{survey}' not in {type(self.local_dirs).__name__}() class configuration: removing from list...")
                self.supported_surveys = tuple(s for s in self.supported_surveys if s != survey)

        # set the filters
        self.survey_filter_sets = list()
        for supported_survey in self.supported_surveys:
            survey_filters  = eval(f"{supported_survey}.get_supported_filters()")
            if survey_filters:
                self.survey_filter_sets.append({supported_survey: [f for f in survey_filters]})

        # set survey_names
        self.survey_block = self.config['cutouts']['surveys']
        self.survey_names = self.__extract_surveys_names(self.survey_block)

        # set the cutout size
        self.size_arcmin = self.config['cutouts']['box_size_arcmin'] * u.arcmin


        #
        # T A R G E T   C O N F I G U R A T I O N
        #

        # set targets
        self.targets = list()
        for coords_csv_file in self.config['cutouts']['ra_dec_deg_csv_files']:
            sources = self.__csv_to_dict(relative_path+coords_csv_file)
    
            # make all keys lower case to effectively reference case-variants of RA and Dec.
            sources = [{k.lower(): v for k,v in s.items()} for s in sources]
    
            # extract position information
            self.targets.extend([{
                'coord': SkyCoord(x['ra'], x['dec'], unit=(u.deg, u.deg)),
                'size':  self.size_arcmin
            } for x in sources])


        #
        # E N V I R O N M E N T   C O N F I G U R A T I O N
        #

        # get the configuration block
        configuration = self.config['configuration']

        # set the data output dir
        data_root = self.__sanitize_path(configuration['local_root'])
        if bool(re.match('/',data_root)): # absolute path case
           out_dir = data_root
        elif bool(re.match('~/',data_root)): # home path case
           out_dir = os.path.expanduser(data_root)
        else: # relative path case
           out_dir = relative_path+data_root
        self.local_dirs.set_local_root(out_dir)
        self.out_dirs = {s: self.local_dirs.get_survey_dir(s) for s in self.survey_names}
        for out_dir in self.out_dirs.values():
            try:
                os.makedirs(out_dir)
            except FileExistsError:
                self.__print(f"Using FITS output dir: {out_dir}")
            else:
                self.__print(f"Created FITS output dir: {out_dir}")

        # set the overwrite file parameter
        if 'overwrite' in configuration:
            self.overwrite = configuration['overwrite']
        else:
            self.overwrite = False

        # set the flush file parameter
        if 'flush' in configuration:
            self.flush = configuration['flush']
        else:
            self.flush = False


    def __sanitize_path(self,path):
        # clean up repeating '/'s with a trailing '/' convention
        return re.sub(r"(/+|/*$)",r"/",path)


    def __csv_to_dict(self,filename):
        entries = []
    
        with open(filename, 'r') as infile:
            c = csv.DictReader(infile)
            for entry in c:
                entries.append(entry)
    
        return entries


    def __get_target_list(self,config_file):
        config  = yml.load(open(config_file,'r'))['cutouts']
    
        targets = list()
        size = config['box_size_arcmin'] * u.arcmin
        for coord_csv_file in config['ra_dec_deg_csv_files']:
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


    def __print(self,string,show_caller=False,is_suspend_on_force_flush=False):
        if string is None or (self.is_force_flush and is_suspend_on_force_flush):
            return
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
                    # TODO (Issue #13): Add output indicating a corrupt yml cfg file and 'raise'
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
                # TODO (Issue #13): This should really be done in init to prevent the potential of repeated message...
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
            # TODO (Issue #13): we need to do this test in init, otherwise the message repeats
            #else:
            #    self.__print(f"WARNING: '{survey}' filter '{filter}' is not supported!")
        
        return set(matched)


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
    

    def get_supported_filters(self,survey):
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


    def get_survey_targets(self):
        return self.targets


    def get_survey_class_stack(self):
        class_stack = list()
        for survey_name in self.get_survey_names():
            if self.has_filters(survey_name):
                for filter in self.get_supported_filters(survey_name):
                    self.__print(f"USING_SUVERY_CLASS: {survey_name}(filter={filter})",is_suspend_on_force_flush=True)
                    class_stack.append(f"{survey_name}(filter={filter})")
            else:
                self.__print(f"USING_SUVERY_CLASS: {survey_name}()",is_suspend_on_force_flush=True)
                class_stack.append(f"{survey_name}()")
        return class_stack 


    def set_overwrite(self,overwrite=True):
        if isinstance(overwrite, bool):
            self.overwrite = overwrite
        return self.overwrite


    def get_overwrite(self):
        return self.overwrite


    def set_flush(self,flush=True):
        if isinstance(flush, bool):
            self.flush = flush
        return self.flush


    def get_flush(self):
        return self.flush


    def __set_force_flush(self):
        self.is_force_flush = True


    def __unset_force_flush(self):
        self.is_force_flush = False


    def force_flush(self):
        self.__print("Flushing...")
        self.__set_force_flush()
        self.get_procssing_stack()
        self.__unset_force_flush()
        self.__print("[done]")


    def get_cutout_filename(self, position,size,survey,filter=None,extension=None):
        return get_cutout_filename(position,size,survey,filter,extension)


    def get_procssing_stack(self):
        # ra-dec-size cutout targets
        survey_targets = self.get_survey_targets()

        # survey-class stack
        survey_classes = self.get_survey_class_stack()

        # ok, let's build the cutout-fetching processing stack
        pid = 0 # task tracking id
        procssing_stack = list()
        for survey_class in survey_classes:
            for survey_target in survey_targets:
                survey_instance = eval(survey_class)

                # ra-dec-size cutout target
                task = dict(survey_target)

                # add survey instance for processing stack
                task['survey'] = survey_instance

                # define the fits output filename
                #coords = survey_instance.get_sexadecimal_string(task['coord'])
                #size = re.sub(r"\.?0+$","","%f" % task['size'].value)
                #survey = type(task['survey']).__name__
                #filter = (lambda f: '' if f is None else f"-{f.name}")(survey_instance.get_filter_setting())
                #task['filename'] = f"{self.out_dirs[survey]}J{coords}_s{size}arcmin_{survey}{filter}.fits"
                coords = task['coord']
                size   = task['size']
                survey = type(task['survey']).__name__
                filter = survey_instance.get_filter_setting()
                path = self.out_dirs[survey]
                task['filename'] = f"{path}{self.get_cutout_filename(coords,size,survey,filter,'fits')}"

                # just flush the file/s if in flush mode

                # set task pid
                task['pid'] = pid

                if self.is_force_flush:
                    self.__print(ProcStatus.flush(task['filename'],is_all=True))
                elif self.overwrite or self.flush or (not ProcStatus.is_processed(task['filename'])):
                    # push the task onto the processing stack
                    procssing_stack.append(task)

                    ## set the task id ...
                    #survey_instance.set_pid(pid)
                    # increment task pid
                    pid += 1

                    if self.overwrite or self.flush:
                        # flush unreprocessable files
                        self.__print(ProcStatus.flush(task['filename'],is_all=self.flush))
                else:
                    files = ProcStatus.get_file_listing(task['filename'])
                    self.__print(f"File{'s' if len(files) > 1 else ''} {files} exist{'' if len(files) > 1 else 's'}; overwrite={self.overwrite}, skipping...")

        # randomize processing stack to minimize server hits...
        shuffle(procssing_stack)

        self.__print(f"CUTOUT PROCESSNING STACK SIZE: {pid}",is_suspend_on_force_flush=True)
        return procssing_stack

