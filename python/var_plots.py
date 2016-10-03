import subprocess
import sys
import os

data_dir="/n/atlasfs/atlasdata/atlasdata1/stchan/vhbb/tmva-data/20160226_srli/output"
def conv_plot(plot):
  if plot[len(plot)-4:]==".eps":
    ok=subprocess.Popen(["epstopdf",plot])
    ok.wait()
  

code_dir="./"
#data_dir="./current_files/"
if len(sys.argv)>1: data_dir = sys.argv[1]
if data_dir[len(data_dir)-1:]!="/": data_dir+="/"
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
def n_vars(which):
    yo=subprocess.Popen("root -b -q /n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/tmva_li_train.C\(\\\"RETURN_N_VARS\\\",%i\)"%which,shell=True,stdout=subprocess.PIPE)
    ya=yo.stdout.read()
    return int(ya[ya.rfind("\n",0,(len(ya)-2)):])

def vars(which):
    yo=subprocess.Popen("root -b -q /n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/tmva_li_train.C\(\\\"RETURN_VARS\\\",%i\)"%which,shell=True,stdout=subprocess.PIPE)
    ya=yo.stdout.read()
    return ya[ya.rfind("\n",0,(len(ya)-2)):].split()


pdir=data_dir+"plots/"

if not os.path.exists(data_dir):
  print "%s for TMVA data does not exist...exiting"
  sys.exit()
  
if not os.path.exists(pdir):
  os.mkdir(pdir)

sig_sum=open("%ssig_summary.txt"%pdir,"w")
sig_sum.write("%15s %15s\n"%("TAG      ","SIG      "))
for i in range(-2,6):
    for j in range(4):
      tag=li_tag(i)+bg_tag(j)
      file="%sTMVA%s.root"%(data_dir,tag)
      log="%stmva_out%s.txt"%(data_dir,tag)
      if not os.path.exists(file): continue
      if os.path.exits(log):
        grep=subprocess.Popen(["grep","Significance",log],stdout=subprocess.PIPE)
        cactuar=proc.stdout.read()
        sigmund=cactuar[cactuar.find(":Significance ")+len(":Significance "):]
        sig_sum.write("%12s    %12s   \n"%(tag[1:],sigmund))
      proc=subprocess.Popen(["root","-q","-l","%smkplots_tmva.C(\"%s\")"%(code_dir,file)],stdout=subprocess.PIPE)
      proc.wait()
      for ext in [".png",".eps"]:
        for sb in ["S","B"]:
          cand_plot="./plots/CorrelationMatrix%s%s"%(sb,ext)
          if os.path.exists(cand_plot):
            new_plot="%sCorrelationMatrix%s%s%s"%(pdir,sb,tag,ext)
            mv_mtrx=subprocess.Popen(["mv",cand_plot,new_plot])
            mv_mtrx.wait()
            conv_plot(new_plot)

    bdt_cand="./plots/overtrain_BDT%s"%(ext)
    if os.path.exists(bdt_cand):
      new_plot="%sovertrain_BDT%s%s"%(pdir,tag,ext)
      mv_over=subprocess.Popen(["mv",bdt_cand,new_plot])
      conv_plot(new_plot)
    effs_cand="./plots/mvaeffs_BDT%s"%(ext)
    if os.path.exists(effs_cand):
      new_plot="%smvaeffs_BDT%s%s"%(pdir,tag,ext)
      mv_over=subprocess.Popen(["mv",effs_cand,new_plot])
      conv_plot(new_plot)

sig_sum.close()
