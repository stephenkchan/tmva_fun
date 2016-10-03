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
tmva_in_dir=dat_base+"20160510/"
do_batch=True
rewrite=True
if len(sys.argv)>1: tmva_in_dir = sys.argv[1]
if len(sys.argv)>2:         tag = sys.argv[2]
if len(sys.argv)>3:    do_batch = (sys.argv[3].lower()=="true" or sys.argv[3]=="1")
if len(sys.argv)>4:    rewrite  = (sys.argv[4].lower()=="true" or sys.argv[4]=="1")


big_bin=subprocess.Popen("export PATH=$PATH:%s/exec"%code_base,shell=True)
big_bin.wait()

if not os.path.exists(tmva_in_dir):
    print "%s for TMVA data does not exist...exiting"%tmva_in_dir
    sys.exit()
    
if tmva_in_dir[len(tmva_in_dir)-1:] != "/": tmva_in_dir+="/"
scr_dir=code_base+"scripts/"
if not os.path.exists(scr_dir):
    os.makedirs(scr_dir)
cuts=[60,70,77,85]
ratios=[0,10,20]
ratios=[20]
cts=[0] #truth-tagged, so it doesn't really matter
jets=[2,3]
ptv=[0,1]
fire_away=[]
def tdf_str(tdf):
  if tdf==1:
    return "_tF"
  elif tdf==2:
    return "_tD"
  return "_std"


files=os.listdir(scr_dir)
files.append(os.listdir("./"))
"""
for output in outputs:
  if ".root" in output or ".txt" in output:
    os.remove(tmva_out_dir+output)
"""
for fil in files:
  if "btag_rank_" in fil and ".sh" in fil:
    os.remove(scr_dir+fil)
  elif "slurm-" in fil and ".out" in fil:
    os.remove(wrk+fil)
  elif "LSF_" in fil:
    shutil.rmtree(wrk+fil)
for tf in [1,2]:
  for pt in ptv:
    for jet in jets:
      for cut in cuts:
        for ratio in ratios:
          for ct in cts:
            snom="%sbtag_rank_beff%i_MV2c%i_%s_%ijet_%sptv_%s.sh"%(scr_dir,cut,ratio,"cts"if ct else "dis",jet,"lo" if pt else "hi",tdf_str(tf))
            command="rank_btag %s %i %i %i %i %i %i %i"%(tmva_in_dir,cut,ratio,ct,jet,pt,rewrite,tf)
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
