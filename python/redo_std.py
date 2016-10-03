import subprocess
import os
import sys

big_bin=subprocess.Popen("export PATH=$PATH:/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/exec",shell=True)
big_bin.wait()

resubmit=True
slurm_dir="/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/work/"
slurm_min=0
do_batch=True
argv=sys.argv
if len(argv)>1: slurm_min=0 if not is_number(argv[1]) else int(argv[1]) 
if len(argv)>2: slurm_dir=argv[2]
if len(argv)>3: do_batchs=(argv[3]=="1" or argv[3].lower()=="true")

if slurm_dir[len(slurm_dir)-1:] != "/":slurm_dir+="/"
def is_number(string):
    try:
        float(string)
        return True
    except ValueError as yup:
        return False
def script_from_log(log):
    if not os.path.exists(log):
        print "You're trying to look through non-existent log file %s"%log
        return
    grep=subprocess.Popen(["grep","(script",log], stdout=subprocess.PIPE)
    line=grep.stdout.read()
    script=line[line.find("(script")+8:line.find(".sh)")+3]
    print script
    return script
        
def full_time(time):
    dummy="00:00:00"
    return dummy[0:len(dummy)-len(time)]+time

proc = subprocess.Popen(["squeue"], stdout=subprocess.PIPE)
buffer = proc.stdout.read().split("\n")
scripts=[]
kill_me=[]
old_logs=[]
files=os.listdir(slurm_dir)

for file in files:
    if "slurm-" in file and ".out" in file:
        if int(file[file.find("slurm-")+6:file.find(".out")+1]) > slurm_min:
            grep=subprocess.Popen(["grep","slurms",file], stdout=subprocess.PIPE)
            if len(grep.stdout.read())>2:
                scripts.append(script_from_log(file))

mem="8192"
precmd= ["sbatch","--mem-per-cpu",mem,"--time","4-0:0:0","-p","pleiades"] if do_batch else ["sh"]
if resubmit:
    for scr in scripts:
        bilbo=subprocess.Popen(precmd+[scr])

for log in old_logs: os.remove(slurm_dir+log)
