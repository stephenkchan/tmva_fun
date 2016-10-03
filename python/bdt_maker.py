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
data_dir=dat_base+"20160307_srli/"
do_batch=True
rewrite=True
if len(sys.argv)>1: data_dir = sys.argv[1]
if len(sys.argv)>2:      tag = sys.argv[2]
if len(sys.argv)>3: do_batch = (sys.argv[3].lower()=="true" or sys.argv[3]=="1")
if len(sys.argv)>4: rewrite  = (sys.argv[4].lower()=="true" or sys.argv[4]=="1")
def bg_tag(bg):
  if bg==0:
      return "_zjets";
  elif bg==1:
      return "_ttbar";
  elif bg==2:
      return "_zz";
  return "_all"


def li_tag(which):
    if which==0:
        return "_li_met"
    elif which==1:
        return "_li_ll"
    elif which==2:
        return "_li_bb"
    elif which==3:
        return "_li_bb_ll"
    elif which==4:
        return "_li_only"
    elif which==5:
        return "_li_lb0"
    elif which==-1:
        return "_std";
    return "_cut";

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

if not os.path.exists(data_dir):
    print "%s for TMVA data does not exist...exiting"%data_dir
    sys.exit()
    
if data_dir[len(data_dir)-1:] != "/": data_dir+="/"
out_dir=data_dir+tag+"/"
scr_dir=code_base+"scripts/"
if not os.path.exists(scr_dir):
    os.makedirs(scr_dir)
if not os.path.exists(out_dir):
    print "Creating %s..."%out_dir
    os.makedirs(out_dir)
    os.makedirs(out_dir+"plots/")

fire_away=[]

if rewrite:
    files=os.listdir(scr_dir)
    outputs=os.listdir(out_dir)
    files.append(os.listdir("./"))
    for output in outputs:
      if ".root" in output or ".txt" in output:
        os.remove(out_dir+output)
    for fil in files:
        if "bdt_" in fil and ".sh" in fil:
            os.remove(scr_dir+fil)
        elif "slurm-" in fil and ".out" in fil:
            os.remove(scr_dir+fil)
        elif "LSF_" in fil:
            shutil.rmtree(scr_dir+fil)
    for i in range(-2,6):
        i_vars=vars(i)
        #print "\n%i for %s vars are: %s"%(len(i_vars),li_tag(i)," ".join(i_vars))
        for j in range(-1,2):
            for nv in range(len(i_vars)):
                back="_%i"%nv if nv>0 else ""
                snom="%sbdt%s%s%s.sh"%(scr_dir,li_tag(i),bg_tag(j),back)
                #print "%s "%snom,
                command="exec_tmva_li_train %s %i %i %s %i"%(data_dir,i,j,tag,nv)
                log="%stmva_out%s%s%s.txt"%(out_dir,li_tag(i),bg_tag(j),back)
                fire_away.append(snom)
                script=file(snom,"w")
                script.write("#!/bin/bash\n")
                script.write("cd %swork\n"%code_base)
                script.write("echo \"About to execute (script %s) %s\"\n"%(snom,"%s > %s\n"%(command,log)))
#                script.write("root -b -q /n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/tmva_li_train.C\(\\\"%s\\\",%i,%i,%i\\)>%stmva_out%s%s%s.txt\n"%(data_dir,i,j,nv,out_dir,li_tag(i),bg_tag(j),back))
                script.write("%s > %s\n"%(command,log))
                script.write("mv %swork/weights/TMVAClassification%s%s%s_BDT.weights.xml %s\n"%(code_base,li_tag(i),bg_tag(j),back,out_dir))
                script.close()
                ok=subprocess.Popen(["chmod","777",snom])
                ok.wait()
else:
    for i in range(-3,5):
        i_vars=vars(i)
        #print "\n%i for %s vars are: %s"%(len(i_vars),li_tag(i)," ".join(i_vars))
        for j in range(-1,2):
            for nv in range(len(i_vars)):
              back="_%i"%nv if nv>0 else ""
              file=out_dir+"TMVA"+li_tag(i)+bg_tag(j)+back+".root"
              if not os.path.exists(file):
                fire_away.append("%sbdt%s%s%s.sh"%(scr_dir,li_tag(i),bg_tag(j),back))
              elif os.stat(file).st_size < 2000:
                os.remove(file)
                fire_away.append("%sbdt%s%s%s.sh"%(scr_dir,li_tag(i),bg_tag(j),back))

#sys.exit()
mem="2048"
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

