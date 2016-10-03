import subprocess
import os
import sys

big_bin=subprocess.Popen("export PATH=$PATH:/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/exec",shell=True)
big_bin.wait()

resubmit=True
slurm_dir="/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/work/"
do_batch=True
argv=sys.argv
if len(argv)>1: resubmit=(argv[1]=="1" or argv[1].lower()=="true")
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
for line in buffer:
    words=line.split()
    if len(words)>6:
        if is_number(words[0]) and full_time(words[5])>full_time("15:00"):
            kill_me.append(words[0])
            log="%sslurm-%s.out"%(slurm_dir,words[0])
            print "kill for time",
            scripts.append(script_from_log(log))
            old_logs.append(log)

grep=subprocess.Popen(["grep","-r","mv:",slurm_dir], stdout=subprocess.PIPE)
lines=grep.stdout.read().split("\n")
for line in lines:
    if line.find("slurm-")==-1:
        continue
    log=line[line.find("slurm-"):line.find(".out")+4]
    print "resubmit for fail ",
    scripts.append(script_from_log(log))
    old_logs.append(log)

#sys.exit()
for kill in kill_me:
    print "Killing "+kill
    yo=subprocess.Popen(["scancel",kill])
    yo.wait()

mem="2048"
precmd= ["sbatch","--mem-per-cpu",mem,"--time","4-0:0:0","-p","pleiades"] if do_batch else ["sh"]
if resubmit:
    for scr in scripts:
        bilbo=subprocess.Popen(precmd+[scr])

for log in old_logs: os.remove(slurm_dir+log)
