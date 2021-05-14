import glob

class emptyError(Exception):
    pass

def pathExistCheck(check_list):
    for path in check_list:
        actual_path=glob.glob(path)
        try:
            actual_path[0]
        except IndexError:
            raise emptyError(path+" is not found!")

def setDefaultValueInConfig(cmd,config):
    if cmd=="import":
        required_parameters=["src_raw_components"]
        for i in required_parameters:
            try:
                config[i]
                if config[i]=="":
                    raise emptyError
            except (KeyError,emptyError):
                err_msg="The parameter "+i+" is not found in the setting file!"
                raise KeyError(err_msg)


        flash_options=["flash_min_overlap","flash_max_overlap"]
        for i in flash_options:
            try:
                config[i]
                if config[i]=="":
                    raise emptyError
            except (KeyError,emptyError):
                if i=="flash_min_overlap":
                    config[i]=10
                elif i=="flash_max_overlap":
                    config[i]=65

        return config