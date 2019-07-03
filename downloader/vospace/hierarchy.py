import re

import yaml as yml


class Hierarchy():
    def __init__(self):
        self.config = yml.load(open("hierarchy.yml"))
        self.local_root  = self.__sanitize_path(self.config['root']['local'])
        self.remote_root = self.__sanitize_path(self.config['root']['remote'])
        self.hierarchy = self.config['hierarchy']
        self.surveys = [s.lower() for k in [list(self.hierarchy[b].keys()) for b in self.hierarchy.keys()] for s in k]

    def __sanitize_path(self,path):
        return re.sub(r"(/+|/*$)",r"/",path)

    def __get_survey_path(self, survey):
        if self.has_survey(survey):
            for band in self.hierarchy.keys():
                if survey.lower() in self.hierarchy[band]:
                    return self.__sanitize_path(band+'/'+survey.lower())
        return None


    def get_local_root(self):
        return self.local_root


    def get_remote_root(self):
        return self.remote_root


    def has_survey(self, survey):
        return survey.lower() in self.surveys


class LocalDirs(Hierarchy):
    def __init__(self):
        super().__init__()


    def get_survey_path(self, survey):
        path = self._Hierarchy__get_survey_path(survey)
        if not path is None:
            return self.get_local_root()+path
        return path


class RemoteDirs(Hierarchy):
    def __init__(self):
        super().__init__()


    def get_survey_path(self, survey):
        path = self._Hierarchy__get_survey_path(survey)
        if not path is None:
            return self.get_remote_root()+path
        return path


class LocalCutouts(LocalDirs):
    def __init__(self):
        super().__init__()
        self.sub_dir = self._Hierarchy__sanitize_path(self.config['catagories']['cutouts']['sub_dir'])


    def get_survey_path(self, survey):
        path = super().get_survey_path(survey)
        if not path is None:
            return path+self.sub_dir
        return None


class RemoteCutouts(RemoteDirs):
    def __init__(self):
        super().__init__()
        self.sub_dir = self._Hierarchy__sanitize_path(self.config['catagories']['cutouts']['sub_dir'])


    def get_survey_path(self, survey):
        path = super().get_survey_path(survey)
        if not path is None:
            return path+self.sub_dir
        return None


if __name__ == '__main__':
    print()
    print(" *** Hierarchy() Class Test ***")
    h = Hierarchy()
    print(f"DATA_ROOT: {h.get_local_root()}")
    print(f"VOSPACE_ROOT: {h.get_remote_root()}")
    print(f"Has FirSt: {h.has_survey('FirSt')}")
    print(f"Has Bart: {h.has_survey('Bart')}")

    print()
    print(" *** LocalDirs(Hierarchy) Class Test ***")
    local = LocalDirs()
    print(f"get_survey_path('First'): {local.get_survey_path('First')}")

    print()
    print(" *** RemoteDirs(Hierarchy) Class Test ***")
    local = RemoteDirs()
    print(f"get_survey_path('PanSTARRS'): {local.get_survey_path('PanSTARRS')}")

    print()
    print(" *** LocalCutouts(LocalDirs(Hierarchy)) Class Test ***")
    local_cutouts = LocalCutouts()
    print(f"get_survey_path('SDss'): {local_cutouts.get_survey_path('SDss')}")

    print()
    print(" *** RemoteCutouts(RemoteDirs(Hierarchy)) Class Test ***")
    remote_cutouts = RemoteCutouts()
    print(f"get_survey_path('SDss'): {remote_cutouts.get_survey_path('SDss')}")




    print()
