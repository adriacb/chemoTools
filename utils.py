import pathlib

def path_Extension(path):
    # function to return the file extension
    file_extension = pathlib.Path(path).suffix
    print("File Extension: ", file_extension)
    return file_extension
