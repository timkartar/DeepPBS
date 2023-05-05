# builtin modules
import os
import shutil

def moveFile(fileName, dest):
    path = os.path.join(dest, fileName)
    if os.path.exists(path):
        if os.path.abspath(dest) != os.getcwd():
            os.remove(path)
            shutil.move(fileName, dest)
    else:
        shutil.move(fileName, dest)
    
    return os.path.join(dest, fileName)
