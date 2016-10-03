import subprocess
import sys
import os

data_dir="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/20160226_srli/"
do_batch=True
rewrite=True
if len(sys.argv)>1: data_dir = sys.argv[1]
if len(sys.argv)>2: do_batch = (sys.argv[2].lower()=="true" or sys.argv[2]=="1")
if len(sys.argv)>3: rewrite  = (sys.argv[3].lower()=="true" or sys.argv[3]=="1")
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
    elif which==-1:
        return "_std";
    elif which==-2:
        return "_std_met";
    return "_cut";

def n_vars(which,stop=-1):
    yo=subprocess.Popen("exec_tmva_li_train RETURN_N_VARS %i -1 %i"%(which,stop),shell=True,stdout=subprocess.PIPE)
    ya=yo.stdout.read()
    return int(ya[ya.rfind("\n",0,(len(ya)-2)):])

def vars(which,stop=-1):
    yo=subprocess.Popen("exec_tmva_li_train RETURN_VARS %i -1 %i"%(which,stop),shell=True,stdout=subprocess.PIPE)
    ya=yo.stdout.read()
    return ya[ya.rfind("\n",0,(len(ya)-2)):].split()

big_bin=subprocess.Popen("export PATH=$PATH:/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/exec",shell=True)
big_bin.wait()

for i in range(-3,5):
  i_vars=vars(i)
  nv=len(vars(i))
  for j in range(nv):
    print "The %i vars for li%i, index %i are %s"%(n_vars(i,j),i,j," ".join(vars(i,j)))

