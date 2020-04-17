import os

def validFile(file:str):
    """[summary]
    
    :param file: [description]
    :type file: str
    :raises FileNotFoundError: [description]
    """
    if not os.path.isfile(file):
        raise FileNotFoundError(f"{file} doesn't exist")