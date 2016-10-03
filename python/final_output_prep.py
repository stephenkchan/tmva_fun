import shutil
import sys
import os
import subprocess

data_dir="eos/atlas/user/s/stchan/tmva_fun/20160510/"
temp_dir="./ship_out"
default_tagger="MV2c20bin"
default_eff   ="beff77"

def xml_root_txt_from_log(log_dir,log_entry):
    # the ranking log name looks like "ranking_log_-MV2c20bin-beff60-2tag3jet-0_120ptv_tF.txt"
    # more generically we have "ranking_log_[variable tag]-[tag/jet]-[ptv][transformation].txt"
    # the weights file looks like "TMVAClassification-[ranking log after "log_" before ".txt"]-MV2c20bin-beff60-2tag3jet-0_120ptv-[dash separated var list]_BDT.weights.xml"
    if log_entry.find("ranking_log_")!=0 or log_entry.find(".txt")!=len(log_entry)-4 or "_tF" not in log_entry:
        return ["nope","nope","nope"]
    if log_dir[len(log_dir)-1]!="/": log_dir+="/"
    log_tag=log_entry[12:len(log_entry)-4:]
    log_name=log_dir+log_entry
    proc = subprocess.Popen(['grep',"DAISY",log_name], stdout=subprocess.PIPE)
    buffer = proc.stdout.read()
    varlist = buffer[buffer.find(": ")+2:len(buffer)-2].replace(" ","-") #last two characters are a space and a newline
    phase_space=""
    if "MV2c" in log_tag:
        phase_space = log_tag[log_tag.find("MV2c"):log_tag.find("ptv")+3]
    else:
        phase_space = default_tagger+"-"+default_eff+"-"+log_tag[log_tag.find("tag")-1:log_tag.find("ptv")+3]
    xml_file ="%sTMVAClassification-%s-%s-%s_BDT.weights.xml" % (log_dir,log_tag,phase_space,varlist)
    root_file="%sTMVA-all-%s-%s-%s.root" % (log_dir,log_tag,phase_space,varlist)
    return [log_name,xml_file,root_file]

def ranking_results_files(directory,delimiter=""):
    results=[]
    if directory[len(directory)-1]!="/": directory+="/"
    for entry in os.listdir(directory):
        xml_cand=xml_root_txt_from_log(directory,entry)
        if ".xml" in xml_cand[1] and "pct" in xml_cand[1]:
            results+=xml_cand
    return results

argv=sys.argv
flag=""
if len(argv)>1: data_dir=argv[1]
if len(argv)>2: temp_dir=argv[2]
if len(argv)>3:     flag=argv[3]

if not os.path.isdir(temp_dir): os.makedirs(temp_dir)

copy_me=ranking_results_files(data_dir,flag)
#print copy_me
#sys.exit()
for entry in copy_me:
    if not os.path.exists(entry):
        print "You're looking for some ranking result file called %s, that DOESN'T exist"%entry
        continue
    shutil.copy2(entry,temp_dir)
