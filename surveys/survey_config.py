#system
import os
import sys
import io

# utilities
import re

# add module paths that are two levels up from here
this_source_file_dir = re.sub(r"(.*/).*$",r"\1",os.path.realpath(__file__))
sys.path.append(this_source_file_dir+"../..")

# import the vospace space module to get the data-subdir configuration
# from .hierarchy import LocalCutoutDirs

from random import shuffle

# configuration
import csv
import yaml as yml

# processing
from .survey_abc import processing_status as ProcStatus
from .toolbox import extractCoordfromString, readCoordsFromFile, get_cutout_filename

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


class SurveyConfig:
    def __init__(self, yml_file_or_list):
        # set up outdirs
        self.relative_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/'
        #self.local_dirs = LocalCutoutDirs() #ONLY USED FOR HEIRARCHY

        # defaults
        self.overwrite = True
        self.flush = True
        self.survey_block = None
        # file flush utility flag
        self.is_force_flush = False

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

        # set the filters
        self.survey_filter_sets = list()
        for supported_survey in self.supported_surveys:
            survey_filters  = eval(f"{supported_survey}.get_supported_filters()")
            if survey_filters:
                self.survey_filter_sets.append({supported_survey: [f for f in survey_filters]})

        if type(yml_file_or_list)==list:
            out_dir = self.configure_single_source(yml_file_or_list)
        else:
            out_dir = self.configure_batch(yml_file_or_list)

        ############# old way ONLY USED FOR HEIRARCHY files###########
        #self.local_dirs.set_local_root(out_dir)
        #self.out_dirs = {s: self.local_dirs.get_survey_dir(s) for s in self.survey_names}
        ###############################################################
        # same config for sinlge or batch
        self.out_dirs = {s: out_dir + s for s in self.survey_names}
        for out_dir in self.out_dirs.values():
            try:
                os.makedirs(out_dir)
            except FileExistsError:
                self.__print(f"Using FITS output dir: {out_dir}")
            else:
                self.__print(f"Created FITS output dir: {out_dir}")

    def configure_batch(self, yml_file):
        file_obj = open(yml_file,'r')
        self.config = yml.load(file_obj, Loader=yml.SafeLoader)
        file_obj.close()
        # get relative path to config file
        #relative_path = re.sub(r"[^/]+$","",yml_file_or_list)
        # set survey_names
        self.survey_block = self.config['cutouts']['surveys']
        self.survey_names = self.__extract_surveys_names(self.survey_block)
        # set the cutout size
        self.size_arcmin = self.config['cutouts']['box_size_arcmin'] * u.arcmin

        # set all targets from csv of targets
        self.set_batch_targets(self.config['cutouts']['ra_dec_deg_csv_files'])

        # E N V I R O N M E N T   C O N F I G U R A T I O N
        # get the configuration block
        configuration = self.config['configuration']
        # set the data output dir
        data_root = self.__sanitize_path(configuration['local_root'])
        # set the overwrite file parameter
        if 'overwrite' in configuration:
            self.overwrite = configuration['overwrite']
        # set the flush file parameter
        if 'flush' in configuration:
            self.flush = configuration['flush']
        if bool(re.match('/',data_root)): # absolute path case
           out_dir = data_root
        elif bool(re.match('~/',data_root)): # home path case
           out_dir = os.path.expanduser(data_root)
        else: # relative path case
           out_dir = self.relative_path+data_root
        return out_dir

    def configure_single_source(self, survey_list):
        if not survey_list: #no surveys specified then use all
            self.survey_names = list(self.supported_surveys)
        else:
            self.survey_names = self.__check_supported(survey_list)
        self.overwrite = True #update single CUTOUTS
        return self.__sanitize_path(self.relative_path+'data')

    def set_single_target_params(self, single_target, size, is_name=False):
        self.size_arcmin = size * u.arcmin
        if not single_target:
            raise Exception("No Target provided!")
        self.targets = [{'position': extractCoordfromString(single_target, is_name), 'size': self.size_arcmin}]

    def set_batch_targets(self, csv_files):
        # set targets in list of dicts
        coords = list()
        for coords_csv_file in csv_files:
            with open(self.relative_path +coords_csv_file, newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                locations, errors = readCoordsFromFile(reader)
                coords.extend(locations)
        self.targets = [dict(item, size=self.size_arcmin) for item in coords]

    def __sanitize_path(self,path):
        # clean up repeating '/'s with a trailing '/' convention
        return re.sub(r"(/+|/*$)",r"/",path)

    # def __csv_to_dict(self,filename):
    #     entries = []
    #     with open(filename, 'r') as infile:
    #         c = csv.DictReader(infile)
    #         for entry in c:
    #             entries.append(entry)
    #     return entries

    def __print(self,string,show_caller=False,is_suspend_on_force_flush=False):
        if string is None or (self.is_force_flush and is_suspend_on_force_flush):
            return
        prefix = type(self).__name__ + (f"({sys._getframe(1).f_code.co_name})" if show_caller else "")
        prefixed_string = "\n".join([f"{prefix}: {s}" for s in string.splitlines()])
        print(prefixed_string)

    # check list of surveys and return list of supported surveys
    def __check_supported(self, surveys):
        supported = list()
        for survey in surveys:
            added = False
            for supported_survey in self.supported_surveys:
                if survey.lower() == supported_survey.lower():
                    supported.append(supported_survey)
                    added = True
                    break
            # TODO (Issue #13): This should really be done in init to prevent the potential of repeated message...
            if not added:
                self.__print(f"WARNING: Survey '{survey}' is not supported!")
        return supported

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
        return self.__check_supported(surveys)

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
        self.__print(f"WARNING: Survey '{survey}' nottargets found!")
        return False

    def has_filters(self,survey):
        if self.has_survey(survey) and self.survey_block:
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

    def get_supported_survey(self):
        return self.supported_surveys

    def get_survey_names(self):
        return self.survey_names

    def get_survey_targets(self):
        return self.targets

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

    def get_survey_class_stack(self):
        class_stack = list()
        for survey_name in self.get_survey_names():
            if self.has_filters(survey_name):
                for filter in self.get_supported_filters(survey_name):
                    self.__print(f"USING_SURVEY_CLASS: {survey_name}(filter={filter})",is_suspend_on_force_flush=True)
                    class_stack.append(f"{survey_name}(filter={filter})")
            else:
                self.__print(f"USING_SURVEY_CLASS: {survey_name}()",is_suspend_on_force_flush=True)
                class_stack.append(f"{survey_name}()")
        return class_stack

    def force_flush(self):
        self.__print("Flushing...")
        self.__set_force_flush()
        self.get_procssing_stack()
        self.__unset_force_flush()
        self.__print("[done]")

    # main initial processing step
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
                survey = type(task['survey']).__name__
                filter = survey_instance.get_filter_setting()
                path = self.out_dirs[survey]
                radius = task['size']/2
                task['filename'] = f"{path}/{get_cutout_filename(task['position'],radius,survey,filter,'fits')}"
                # set task pid
                task['pid'] = pid

                if self.is_force_flush:
                    self.__print(ProcStatus.flush(task['filename'],is_all=True))
                elif self.overwrite or self.flush or (not ProcStatus.is_processed(task['filename'])):
                    # push the task onto the processing stack
                    procssing_stack.append(task)
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
