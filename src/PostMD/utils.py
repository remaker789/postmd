import os
import shutil

def create_path(path, backup=False):
    """创建路径为path的文件夹。如果文件夹已存在，则将原文件夹备份，并添加.bkXXX后缀

    Args:
        path (_type_): _description_
    """    
    path += "/"
    if os.path.isdir(path) and backup:
        dirname = os.path.dirname(path)
        counter = 0
        while True:
            bk_dirname = dirname + ".bk%03d" % counter
            if not os.path.isdir(bk_dirname):
                shutil.move(dirname, bk_dirname)
                break
            counter += 1
    os.makedirs(path)
    