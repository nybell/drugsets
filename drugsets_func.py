#drugsea_func.py>

# import packages 
import os
import subprocess
from tqdm import tqdm
from subprocess import Popen, PIPE, CalledProcessError

# define function to run commands in terminal
def run_task(cmd):
    with Popen(cmd, shell=True, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for b in tqdm(p.stdout):
            # print(b, end='') # b is the byte from stdout
            next

    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)

def run_task_silent(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

# define function to create new gene set file with custom size
def setsize(path, file, size):
    
    # create name for new gene set file 
    new = file.replace('.txt', '_min'+str(size)+'.txt')
    new = "/tmp"+new

    # add path
    new = path+new

    # create file 
    with open(os.path.normpath(path+file)) as oldfile, open(os.path.normpath(new), 'w') as newfile:
        for line in oldfile:
                if len(line.split('\t')) -3 >= int(size):
                    newfile.write(line)
