import os
import re

import yaml as yml


from abc import ABC, abstractmethod
class HierarchyABC(ABC):
    def __init__(self):
        super().__init__()

        # load the config file, indepent of where the source code is run...
        this_source_file_dir = re.sub(r"(.*/).*$",r"\1",os.path.realpath(__file__))
        self.config = yml.load(open(this_source_file_dir+"hierarchy.yml"), Loader=yml.FullLoader)

        # set local (relative) and remote (absolute) paths
        self.local_root  = self.__sanitize_path(self.config['root']['local'])
        self.remote_root = self.__sanitize_path(self.config['root']['remote'])
        self.data_subdir = self.__sanitize_path(self.config['data_subdir'])

        # get the base data directory hierarchy
        self.hierarchy = self.config['hierarchy']

        # compile a list of support surveys
        self.surveys = [s.lower() for k in [list(self.hierarchy[b].keys()) for b in self.hierarchy.keys()] for s in k]

    def __sanitize_path(self,path):
        # clean up repeating '/'s with a trailing '/' convention
        return re.sub(r"(/+|/*$)",r"/",path)


    def __get_base_survey_path(self, survey):
        if self.has_survey(survey):
            for band in self.hierarchy.keys():
                if survey.lower() in self.hierarchy[band]:
                    return self.__sanitize_path(band+'/'+survey.lower()+'/'+self.data_subdir)
        return None


    def set_local_root(self,local_root):
        self.local_root = self.__sanitize_path(local_root)
        return self


    def get_local_root(self):
        return self.local_root


    def get_remote_root(self):
        return self.remote_root


    def has_survey(self, survey):
        return survey.lower() in self.surveys


    @abstractmethod
    def get_survey_dir(self, survey):
        pass


    def get_survey_dirs(self):
        return [self.get_survey_dir(s) for s in self.surveys]


class LocalDirs(HierarchyABC):
    def __init__(self):
        super().__init__()


    def get_survey_dir(self, survey):
        path = self._HierarchyABC__get_base_survey_path(survey)
        if not path is None:
            return self.get_local_root()+path
        return path


class RemoteDirs(HierarchyABC):
    def __init__(self):
        super().__init__()


    def get_survey_dir(self, survey):
        path = self._HierarchyABC__get_base_survey_path(survey)
        if not path is None:
            return self.get_remote_root()+path
        return path


class LocalCutoutDirs(LocalDirs):
    def __init__(self):
        super().__init__()
        self.sub_dir = self._HierarchyABC__sanitize_path(self.config['catagories']['cutouts']['sub_dir'])


    def get_survey_dir(self, survey):
        path = super().get_survey_dir(survey)
        if not path is None:
            return path+self.sub_dir
        return None


class RemoteCutouts(RemoteDirs):
    def __init__(self):
        super().__init__()
        self.sub_dir = self._HierarchyABC__sanitize_path(self.config['catagories']['cutouts']['sub_dir'])


    def get_survey_dir(self, survey):
        path = super().get_survey_dir(survey)
        if not path is None:
            return path+self.sub_dir
        return None


if __name__ == '__main__':
    #print()
    #print(" *** HierarchyABC() Class Test ***")
    #h = HierarchyABC()
    #print(f"DATA_ROOT: {h.get_local_root()}")
    #print(f"VOSPACE_ROOT: {h.get_remote_root()}")
    #print(f"Has FirSt: {h.has_survey('FirSt')}")
    #print(f"Has Bart: {h.has_survey('Bart')}")

    print()
    print(" *** LocalDirs(HierarchyABC) Class Test ***")
    local = LocalDirs()
    print(f"get_survey_dir('First'): {local.get_survey_dir('First')}")
    print("get_survey_dirs(): \n> "+"\n> ".join(local.get_survey_dirs()))

    print()
    print(" *** RemoteDirs(HierarchyABC) Class Test ***")
    remote = RemoteDirs()
    print(f"get_survey_dir('PanSTARRS'): {remote.get_survey_dir('PanSTARRS')}")
    print("get_survey_dirs(): \n> "+"\n> ".join(remote.get_survey_dirs()))

    print()
    print(" *** LocalCutoutDirs(LocalDirs(HierarchyABC)) Class Test ***")
    local_cutouts = LocalCutoutDirs()
    print(f"get_survey_dir('SDss'): {local_cutouts.get_survey_dir('SDss')}")
    print("get_survey_dirs(): \n> "+"\n> ".join(local_cutouts.get_survey_dirs()))

    print()
    print(" *** RemoteCutouts(RemoteDirs(HierarchyABC)) Class Test ***")
    remote_cutouts = RemoteCutouts()
    print("get_survey_dirs(): \n> "+"\n> ".join(remote_cutouts.get_survey_dirs()))




    print()
