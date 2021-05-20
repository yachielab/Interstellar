import glob

class EmptyError(Exception):
    pass

def pathExistCheck(check_list):
    for path in check_list:
        actual_path=glob.glob(path)
        try:
            actual_path[0]
        except IndexError:
            raise EmptyError(path+" is not found!")

def setDefaultValueInConfig(cmd,config,distribute=False):
    if cmd=="import":
        required_parameters=["src_raw_components"]
        for i in required_parameters:
            try:
                config[i]
                if config[i]=="":
                    raise EmptyError
            except (KeyError,EmptyError):
                err_msg="The parameter "+i+" is not found in the config file."
                raise KeyError(err_msg)
        flash_options=["flash_min_overlap","flash_max_overlap"]
        for i in flash_options:
            try:
                config[i]
                if config[i]=="":
                    raise EmptyError
            except (KeyError,EmptyError):
                if i=="flash_min_overlap":
                    config[i]=10
                elif i=="flash_max_overlap":
                    config[i]=65
        return config
    
    elif cmd=="streamline":
        required_parameters=["pipeline","qoption","template_shellscript"]
        for i in required_parameters:
            if i=="qoption" and not distribute:
                continue

            try:
                config[i]
                if config[i]=="":
                    raise EmptyError
            except (KeyError,EmptyError):
                err_msg="The parameter "+i+" is not found in the config file."
                raise KeyError(err_msg)
        

        for i in ["mem_max","chunksize"]:
            try:
                config[i]
                if config[i]=="":
                    raise EmptyError
            except (KeyError,EmptyError):
                # if i== "distribute":
                #     config["distribute"]="off"
                if i=="mem_max":
                    config["mem_max"]='128'
                if i=="chunksize":
                    config["chunksize"]='2000000'
        
        for i in config["pipeline"].split(","):
            try:
                    config["mem_"+i]
                    if config["mem_"+i]=="":
                        raise EmptyError
            except (KeyError,EmptyError):
                if i in ["import","qc","to_bt","bc_sort"]:
                    config["mem_"+i]='8'
                elif i in ["mk_sval","buildTree","convert","export","tag","demultiplex"]:
                    config["mem_"+i]='16'
                elif i in ["mergeTree","correct"]:
                    config["mem_"+i]=config["mem_max"]

        return config