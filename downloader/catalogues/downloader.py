# notes: 
# [1] http://www.canfar.net/en/docs/storage/
# [2] https://python-future.org/compatible_idioms.html
from __future__ import print_function

import os
import shutil
import re
import yaml
import argparse
from vos import Client


class downloader(Client):
    def __init__(self):
        super(downloader,self).__init__()

        # extact yaml configuration information
        self.config = yaml.load(open("hierarchy.yml"))
        self.remote_root_dir = self.config['root']['remote']
        self.local_root_dir = self.config['root']['local']
        self.dps_root_dir = self.config['root']['dps']
        self.docker_compose_file = self.get_dps_root_dir()+"docker-compose.yml"

        #print(self.get_dps_root_dirs())
        #print 
        #print(self.get_dps_root_dir())


    def __get_full_repo_path(self,relative_repo_path):
        def get_this_source_file_directory():
            def path(fullpath):
                return re.sub(r"(.*/).*$",r"\1",fullpath)
            return path(os.path.realpath(__file__))
        return re.sub("%s" % self.local_root_dir,"",get_this_source_file_directory())+relative_repo_path

    def get_remote_root_dir(self):
        return self.remote_root_dir

    def get_local_root_dir(self):
        return self.__get_full_repo_path(self.local_root_dir)

    def get_dps_root_dir(self):
        return self.__get_full_repo_path(self.dps_root_dir)

    def get_remote_catalogue_dir(self):
        return self.config['catalogues']['remote']

    def get_local_catalogue_dir(self):
        return self.config['catalogues']['local']

    def get_dps_catalogue_dir(self):
        return self.config['catalogues']['dps']

    def __get_hierarchy(self):
        def get_h(r,s):
            if not isinstance(s,dict):
                return r+"/"
            h = list()
            for k in s.keys():
                f = get_h(k,s[k])
                if isinstance(f,list):
                    h.extend(["{0}/{1}".format(r,e) for e in f])
                else:
                    h.append("{0}/{1}".format(r,get_h(k,s[k])))
            return h
        return [re.sub(r"^/","",e) for e in get_h("",self.config['hierarchy'])]

    def get_remote_root_dirs(self):
        return [self.get_remote_root_dir()+e for e in self.__get_hierarchy()]

    def get_local_root_dirs(self):
        return [self.get_local_root_dir()+e for e in self.__get_hierarchy()]

    def get_dps_root_dirs(self):
        return [self.get_dps_root_dir()+e for e in self.__get_hierarchy()]

    def __create_dirs_if_not_exist(self,dirs,is_repo):
        if is_repo:
            for dir in dirs:
                if not os.path.isdir(dir):
                    print("> creating:",dir)
                    os.makedirs(dir)
        else:
            # TODO/OPTIONAL: Make it work for CADC/CANFAR's VOSpace (i.e., use vmkdir); may require work,
            # as it appears to have a limited library.
            pass
        return dirs

    def get_remote_catalogue_dirs(self):
        return self.__create_dirs_if_not_exist([e+self.get_remote_catalogue_dir() for e in self.get_remote_root_dirs()],is_repo=False)

    def get_local_catalogue_dirs(self):
        return self.__create_dirs_if_not_exist([e+self.get_local_catalogue_dir() for e in self.get_local_root_dirs()],is_repo=True)

    def get_dps_catalogue_dirs(self):
        return self.__create_dirs_if_not_exist([e+self.get_dps_catalogue_dir() for e in self.get_dps_root_dirs()],is_repo=True)

    def __get_file_configuration(self,root_dir=None):
        files = list()
        with open(self.docker_compose_file,'r') as dcf:
            if root_dir in [self.get_remote_root_dir(),self.get_local_root_dir()]:
                root = root_dir
                catalogue_dir = self.get_remote_catalogue_dir() if root==self.get_remote_root_dir() else self.get_local_catalogue_dir()
            else:
                root = self.get_dps_root_dir()
                catalogue_dir = None
            filter = "catalogue_file\s*=\s*(%s)" % ("|".join(self.__get_hierarchy()))
            for line in dcf:
               if re.search(filter,line):
                   file = root+re.sub(r"^.*?=\s*(.*?)\s*$",r"\1",line)
                   if catalogue_dir:
                       file = re.sub(r"/%s" % self.get_dps_catalogue_dir(),"/%s" % catalogue_dir,file)
                   files.append(file)
        return files 

    def get_required_remote_configuration_files_list(self):
        return self.__get_file_configuration(self.get_remote_root_dir())

    def get_required_local_configuration_files_list(self):
        return self.__get_file_configuration(self.get_local_root_dir())

    def get_required_dps_configuration_files_list(self):
        return self.__get_file_configuration()

    def ls_remote_catalogue_dirs(self):
        dirs = self.get_remote_catalogue_dirs()
        print("remote:")
        for dir in dirs:
            print(dir+":\n> "+"\n> ".join(self.listdir(dir)))

    def ls_local_catalogue_dirs(self):
        print("local:")
        for dir in self.get_local_catalogue_dirs():
            print(dir+":\n> "+"\n> ".join(os.listdir(dir)))

    def ls_dps_catalogue_dirs(self):
        print("dps:")
        for dir in self.get_dps_catalogue_dirs():
            print(dir+":\n> "+"\n> ".join(os.listdir(dir)))

    def __match_source_destination_dirs(self,source,destination):
       def get_match(key,dirs):
           for dir in dirs:
               if key in dir:
                   return dir
           return None
       return [[get_match(h_key,source),get_match(h_key,destination)] for h_key in self.__get_hierarchy()] 


    def __get_copy_tree(self,from_to_dirs,local_force_create_dir=False):
        # notes: https://realpython.com/python-thinking-recursively/
        from_files_to_dirs = list() 
        child_from_to_dirs = list()
        vos = Client() # i.e., self.listder() gives empty list on recursion.
        for from_to_dir in from_to_dirs:
            items = vos.listdir(from_to_dir[0])
            for item in items:
                node = from_to_dir[0]+item
                if self.isfile(node):
                    from_files_to_dirs.append([node,from_to_dir[1]+"."])
                elif self.isdir(node):
                    to_dir = from_to_dir[1]+item+"/"
                    if local_force_create_dir and not os.path.exists(to_dir):
                        os.makedirs(to_dir)
                    from_files_to_dirs.extend(self.__get_copy_tree([[node+"/",to_dir]]))
        return from_files_to_dirs

    def __vcp(self,source,destination):
        print("vcp {0} {1}".format(source,destination))
        try:
            self.copy(source,destination)
        except Exception as e:
            # probably a timeout error, restart vos and try again...
            print("> Whoops!")
            print("Error:",e)
            print("> We probably timed out!")
            print("> Let's try again...")
            print("> vcp {0} {1}".format(source,destination))
            self.copy(source,destination)

    def cp_catalogues_from_remote_to_local(self):
        from_to_dirs = self.__match_source_destination_dirs(self.get_remote_catalogue_dirs(),self.get_local_catalogue_dirs())
        from_files_to_dirs = self.__get_copy_tree(from_to_dirs,True)
        for from_file_to_dir in from_files_to_dirs:
            source = from_file_to_dir[0]
            destination = from_file_to_dir[1]
            self.__vcp(source,destination)


    def __re_sanitize(self,s):
        return re.sub(r"\+","\+",s)


    def cp_configuration_from_remote_to_local(self):
        source_files = self.get_required_remote_configuration_files_list()
        destination_files = self.get_required_local_configuration_files_list()
        for source_file in source_files:
            from_file = re.sub(r"^(.*?/)*","",source_file)
            for destination_file in destination_files:
                to_file = re.sub(r"^(.*?/)*","",destination_file)
                if from_file == to_file:
                    #to_dir = re.sub(r"%s$" % to_file,"",destination_file)
                    to_dir = re.sub(r"%s$" % self.__re_sanitize(to_file),"",destination_file)
                    if not os.path.isdir(to_dir):
                        print("> creating:",to_dir)
                        os.makedirs(to_dir)
                    #print("vcp {0} {1}".format(source_file,destination_file))
                    self.__vcp(source_file,destination_file)
                    break
         
        
    def cp_configuration_from_local_to_dps(self):
        source_files = self.get_required_local_configuration_files_list()
        destination_files = self.get_required_dps_configuration_files_list()
        for source_file in source_files:
            from_file = re.sub(r"^(.*?/)*","",source_file)
            for destination_file in destination_files:
                to_file = re.sub(r"^(.*?/)*","",destination_file)
                if from_file == to_file:
                    to_dir = re.sub(r"%s$" % self.__re_sanitize(to_file),"",destination_file)
                    if not os.path.isdir(to_dir):
                        print("> creating:",to_dir)
                        os.makedirs(to_dir)
                    print("cp {0} {1}".format(source_file,destination_file))
                    shutil.copyfile(source_file,destination_file)
                    break


    def flush_local_files(self):
        # NB: Paranoid method -- using realtive paths.
        cwd = os.getcwd()
        local_dirs = [re.sub(r"^%s/" % cwd,"",d) for d in self.get_local_catalogue_dirs()]
        for local_dir in local_dirs:
            cmd = "rm -rf {0}".format(local_dir)
            print(cmd)
            shutil.rmtree(local_dir)


    def flush_dps_files(self):
        # NB: Paranoid method -- using realtive paths.
        cwd = os.getcwd()
        local_dirs = [re.sub(r"^%s/" % cwd,"",d) for d in self.get_dps_catalogue_dirs()]
        for local_dir in local_dirs:
            cmd = "rm -rf {0}".format(local_dir)
            print(cmd)
            shutil.rmtree(local_dir)


if __name__ == '__main__':
    # chdir to source directory
    cwd = os.getcwd()
    if os.path.dirname(__file__):
        os.chdir(os.path.dirname(__file__))

    # get the vos i/f instance
    vos = downloader()

    # command line arguments
    parser = argparse.ArgumentParser(description="VOSpace Downloader Tool.")
    parser.add_argument('-r',action='store_true',help="set context remote")
    parser.add_argument('-s',action='store_true',help="set context local-scratch")
    parser.add_argument('-l',action='store_true',help="ls command")
    parser.add_argument('-d',action='store_true',help="root directory")
    parser.add_argument('-t',action='store_true',help="directory tree")
    parser.add_argument('-p',action='store_true',help="pull configuration files from context")
    parser.add_argument('--pull-all',action='store_true',help="pull all remote files to local-scratch")
    parser.add_argument('-f',action='store_true',help="flush dirs (-r: forbidden)")
    args = parser.parse_args()

    # exclusion filter
    xors = [
        ['r','s','pull_all'],
        ['l','d','t','p','f','pull_all'],
    ]
    for xor in xors:
        if sum([vars(args)[f] for f in xor]) > 1:
            os.chidir(cwd)
            parser.error("flags -{0} are xor'ed".format("|".join(xor)))

    if args.l:
        if args.r:
            vos.ls_remote_catalogue_dirs()
        elif args.s:
            vos.ls_local_catalogue_dirs()
        else:
            vos.ls_dps_catalogue_dirs()
    elif args.d:
        if args.r:
            print(vos.get_remote_root_dir())
        elif args.s:
            print(vos.get_local_root_dir())
        else:
            print(vos.get_dps_root_dir())
    elif args.t:
        if args.r:
            print("remote:")
            dirs = vos.get_remote_root_dirs()
        elif args.s:
            print("local:")
            dirs = vos.get_local_root_dirs()
        else:
            print("dps:")
            dirs = vos.get_dps_root_dirs()
        for dir in dirs:
            print(dir)
    elif args.p:
        if args.r:
            print("pulling configuration from remote:")
            vos.cp_configuration_from_remote_to_local()
        elif args.s:
            print("pulling configuration from local:")
            vos.cp_configuration_from_local_to_dps()
        else:
            print("pulling configuration from remote:")
            vos.cp_configuration_from_remote_to_local()
            print("pulling configuration from local:")
            vos.cp_configuration_from_local_to_dps()
    elif args.pull_all:
        vos.cp_catalogues_from_remote_to_local()
    elif args.f:
        if args.r:
            print("Connot flush remote dirs -- forbidden.")
        elif args.s:
            vos.flush_local_files()
        else:
            vos.flush_dps_files()
    else:
        parser.print_help()

    # return to calling dir
    os.chdir(cwd)
