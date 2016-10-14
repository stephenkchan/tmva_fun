import subprocess
import sys
import os

tag="output"
slurm=False
machine=subprocess.Popen(["hostname"],stdout=subprocess.PIPE)
hostname=machine.stdout.read()
if "rc.fas.harvard.edu" in hostname:
  print "Working on a SLURM machine %s"%hostname
  slurm=True
else:
  print "Working on an LSF machine %s"%hostname
dat_base="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/" if slurm else "/afs/cern.ch/work/s/stchan/vhbb/data/tmva-data/"
code_base="/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/" if slurm else "/afs/cern.ch/work/s/stchan/vhbb/tmva_fun/"
wrk=code_base+"work/"
tmva_in_dir=dat_base+"20160419_srli/"
do_batch=True
rewrite=False
if len(sys.argv)>1: tmva_in_dir = sys.argv[1]
if len(sys.argv)>2:         tag = sys.argv[2]
if len(sys.argv)>3:    do_batch = (sys.argv[3].lower()=="true" or sys.argv[3]=="1")
if len(sys.argv)>4:    rewrite  = (sys.argv[4].lower()=="true" or sys.argv[4]=="1")

def li_tag(which,btv=True):
  nom="cut"
  if which==0:
    nom="li-only"
  elif which==1: 
    nom="li-met"
  elif which==2: 
    nom="li-ll"
  elif which==3: 
    nom="li-bb"
  elif which==4: 
    nom="li-lb0"
  elif which==5: 
    nom="li-all"
  elif which>=-1:
    nom="std"
  nom+=("-pct" if btv else "-nopct")
  return nom

def n_vars(which):
    yo=subprocess.Popen("exec_tmva_li_train RETURN_N_VARS %i"%which,shell=True,stdout=subprocess.PIPE)
    ya=yo.stdout.read()
    return int(ya[ya.rfind("\n",0,(len(ya)-2)):])

def vars(which):
    yo=subprocess.Popen("exec_tmva_li_train RETURN_VARS %i"%which,shell=True,stdout=subprocess.PIPE)
    ya=yo.stdout.read()
    return ya[ya.rfind("\n",0,(len(ya)-2)):].split()

big_bin=subprocess.Popen("export PATH=$PATH:%s/exec"%code_base,shell=True)
big_bin.wait()

if not os.path.exists(tmva_in_dir):
    print "%s for TMVA data does not exist...exiting"%tmva_in_dir
    sys.exit()
    
if tmva_in_dir[len(tmva_in_dir)-1:] != "/": tmva_in_dir+="/"
tmva_out_dir=tmva_in_dir+tag+"/"
scr_dir=code_base+"scripts/"
if not os.path.exists(scr_dir):
    os.makedirs(scr_dir)
if not os.path.exists(tmva_out_dir):
    print "Creating %s..."%tmva_out_dir
    os.makedirs(tmva_out_dir)
    os.makedirs(tmva_out_dir+"plots/")

fire_away=[]

files=os.listdir(scr_dir)
outputs=os.listdir(tmva_out_dir)
files.append(os.listdir("./"))
"""
for output in outputs:
  if ".root" in output or ".txt" in output:
    os.remove(tmva_out_dir+output)
"""
for fil in files:
  if str(fil).find("finish_")==0 and ".sh" in fil:
    os.remove(scr_dir+fil)
  elif "slurm-" in fil and ".out" in fil:
    os.remove(wrk+fil)
  elif "LSF_" in fil:
    shutil.rmtree(wrk+fil)
for tf in [1]:
  for nj in [2]:
    for lpt in [1]:
      for btv in [0]:
        for i in [0]:# range(-2,6):
        #print "\n%i for %s vars are: %s"%(len(i_vars),li_tag(i)," ".join(i_vars))
          snom="%sfinish_%s-%ijet-%spTV%s.sh"%(scr_dir,li_tag(i,btv),nj,"lo" if lpt else "hi",("_tF" if tf else "_std"))
          command="finish_li %s %s %i %i %i %i %i %i"%(tmva_in_dir,tmva_out_dir,i,nj,lpt,btv,rewrite,tf)
          fire_away.append(snom)
          script=file(snom,"w")
          script.write("#!/bin/bash\n")
          script.write("cd %swork\n"%code_base)
          script.write("echo \"About to execute (script %s) %s\"\n"%(snom,"%s \n"%(command)))
          script.write("%s\n"%(command))
          script.close()
          ok=subprocess.Popen(["chmod","777",snom])
          ok.wait()

#sys.exit()
mem="8192"
lsfq="8nh"
slrq="pleiades"
batch_cmd=["sbatch","--mem-per-cpu",mem,"--time","4-0:0:0","-p",slrq] if slurm else ["bsub","-q",lsfq]
precmd= batch_cmd if do_batch else ["sh"]

n_at=0
for scr in fire_away:
    bilbo=subprocess.Popen(precmd+[scr])
    if not do_batch: 
      n_at+=1
      print "Doing script %i of %i"%(n_at,len(fire_away))
      bilbo.wait()

