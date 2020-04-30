#system
import os, sys, io, glob
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
from core.survey_abc import processing_status as ProcStatus
from core.toolbox import *
# astropy libs
from astropy import units as u
from astropy.coordinates import SkyCoord
# filters used by various surveys
from core.survey_filters import grizy_filters, wise_filters, ugriz_filters
# supported suverys (nb: cf., SurveyConfig::self.supported_surveys)
from core.nvss      import NVSS
from core.first     import FIRST
from core.wise      import WISE
from core.sdss      import SDSS
from core.vlass     import VLASS
from core.panstarrs import PANSTARRS


# THIS IS ONLY USED FOR THE COMMAND LINE INTERFACE CAN BE REMOVED FROM CORE
#
class CLIConfig:
    # define Class var

    def __init__(self, surveys, data_out='data_out', group_by=None):
        # set up outdirs
        # self.relative_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/'
        #self.local_dirs = LocalCutoutDirs() #ONLY USED FOR HEIRARCHY
        # defaults
        self.overwrite = False
        self.survey_filter_sets = {} #None # this is to keep track of requested survey filters
        self.supported_surveys = (
            FIRST.__name__,
            NVSS.__name__,
            VLASS.__name__,
            WISE.__name__,
            PANSTARRS.__name__,
            SDSS.__name__,
            GLEAM.__name__,
        )
        self.group_by = group_by
        self.survey_names = []
        self.set_survey_filter_sets(surveys)
        ############# old way ONLY USED FOR HEIRARCHY files###########
        #self.local_dirs.set_local_root(out_dir)
        #self.out_dirs = {s: self.local_dirs.get_survey_dir(s) for s in self.survey_names}
        ###############################################################
        # same config for sinlge or batch
        self.out_dirs = {s: os.path.join(data_out,s) for s in self.survey_names}
        if "PANSTARRS" in self.out_dirs.keys():
            self.out_dirs["PANSTARRS"] = self.out_dirs["PANSTARRS"].replace("PANSTARRS", "PanSTARRS")
        for out_dir in self.out_dirs.values():
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
                self.__print(f"Created FITS output dir: {out_dir}")

    def set_overwrite(self,overwrite=True):
        if isinstance(overwrite, bool):
            self.overwrite = overwrite
        return self.overwrite

    def get_overwrite(self):
        return self.overwrite

    def get_survey_targets(self):
        return self.targets

    def __print(self,string,show_caller=False):
        if string is None:
            return
        prefix = type(self).__name__ + (f"({sys._getframe(1).f_code.co_name})" if show_caller else "")
        prefixed_string = "\n".join([f"{prefix}: {s}" for s in string.splitlines()])

    def flush_old_survey_data(self):
        for survey, dir in self.out_dirs.items():
            if glob.glob(dir+'/*'):
                print(f"Flushing Old Data from {dir}")
                for file in glob.glob(dir+'/*'):
                    os.remove(file)
            else:
                print(f"NO DATA TO FLUSH AT {dir}")
        return

    def set_survey_filter_sets(self, survey_list):
        if not survey_list: #no surveys specified then use all
            self.survey_names = list(self.supported_surveys)
            self.survey_filter_sets = {s:set() for s in self.survey_names}
        else:
            self.survey_names = []
            for survey in survey_list:
                filters = set()
                # user can use square or round brackets
                rund = re.search('\(([^)]+)\)',survey)
                square = re.search('\[([^)]+)\]',survey)
                if rund:
                    filters = [str(f) for f in rund[0].replace('(','').replace(')','').split(',')]
                    survey = survey.replace(rund[0],'')
                elif square:
                    filters = [str(f) for f in square[0].replace('[','').replace(']','').split(',')]
                    survey = survey.replace(square[0],'')
                survey = survey.upper()
                if survey in self.supported_surveys:
                    safe_filters = self.match_filters(survey,filters)
                    self.survey_filter_sets[survey] = safe_filters
                    self.survey_names.append(survey)
                else:
                    print(f"Survey {survey} NOT Supported!")
        return

    def set_single_target_params(self, single_target, size, is_name=False):
        self.size_arcmin = size * u.arcmin
        if not single_target:
            raise Exception("No Target provided!")
        self.targets = [{'position': extractCoordfromString(single_target, is_name), 'size': self.size_arcmin}]

    def set_batch_targets(self, csv_files, relative_path, size):
        self.size_arcmin = size * u.arcmin
        # set targets in list of dicts
        coords = list()
        for coords_csv_file in csv_files:
            with open(relative_path+coords_csv_file, newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                locations, errors = readCoordsFromFile(reader)
                coords.extend(locations)
        self.targets = [dict(item, size=self.size_arcmin) for item in coords]

    def match_filters(self,survey,filters):
        # Get class from globals and create an instance
        survey_class = globals()[survey.upper()]()
        if isinstance(filters,str):
            filters = [filters]
        matched = list()
        for filter in filters:
            found = False
            # Get the function (from the instance) that we need to call
            for supported_filter in getattr(survey_class, 'get_supported_filters')():
            # for supported_filter in get_supported_filters(survey):
                if supported_filter.name.lower() == filter.lower():
                    found = True
                    break
            if found:
                matched.append(supported_filter)
                found = False
            else:
                self.__print(f"WARNING: '{survey}' filter '{filter}' is not supported!")
        return set(matched)

    def get_survey_class_stack(self):
        class_stack = list()
        for survey_name in self.survey_names:
            if len(self.survey_filter_sets[survey_name])>0: #list of filters is non-empty
            # if self.has_filters(survey_name):
                for filter in self.survey_filter_sets[survey_name]:
                    self.__print(f"USING_SURVEY_CLASS: {survey_name}(filter={filter})")
                    class_stack.append(f"{survey_name}(filter={filter})")
            else:
                self.__print(f"USING_SURVEY_CLASS: {survey_name}()")
                class_stack.append(f"{survey_name}()")
        return class_stack

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
                # ra-dec-size cutout target
                task = dict(survey_target)
                task['survey'] = eval(survey_class) # add survey instance to processing stack
                survey = type(task['survey']).__name__
                task['survey'].set_out_dir(self.out_dirs[survey]) #set where to store output
                task['survey'].overwrite = self.overwrite
                # filter = task['survey'].get_filter_setting()
                # radius = task['size']/2
                task['group_by'] = self.group_by
                # set task pid
                task['pid'] = pid
                # if self.overwrite or (not ProcStatus.is_processed(task['filename'])):
                # push the task onto the processing stack
                procssing_stack.append(task)
                # increment task pid
                pid += 1
                # else:
                #     files = ProcStatus.get_file_listing(task['filename'])
                #     self.__print(f"File{'s' if len(files) > 1 else ''} {files} exist{'' if len(files) > 1 else 's'}; overwrite={self.overwrite}, skipping...")
        # randomize processing stack to minimize server hits...
        shuffle(procssing_stack)

        self.__print(f"CUTOUT PROCESSNING STACK SIZE: {pid}")
        return procssing_stack
