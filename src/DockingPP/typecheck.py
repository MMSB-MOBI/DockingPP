import os

def validFile(file:str):
    if not os.path.isfile(file):
        raise FileNotFoundError(f"{file} doesn't exist")
