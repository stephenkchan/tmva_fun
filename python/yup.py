import subprocess
import os
import sys

do_batch=True
big_bin=subprocess.Popen("export PATH=$PATH:/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/exec",shell=True)
big_bin.wait()

scripts=["/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_bb_all_6.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_bb_all_7.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_std_zjets_8.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_ll_zjets.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_met_zjets.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_met_all_10.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_only_zjets_9.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_met_zjets_5.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_std_met_zjets_6.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_ll_zjets_6.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_std_met_all_10.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_bb_all_2.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_ll_all_9.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_ll_all_7.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_bb_ll_all_2.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_ll_zjets_11.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_ll_all_8.sh",
"/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/scripts/bdt_li_ll_all_2.sh"]

mem="2048"
precmd= ["sbatch","--mem-per-cpu",mem,"--time","4-0:0:0","-p","pleiades"] if do_batch else ["sh"]

for scr in scripts:
    bilbo=subprocess.Popen(precmd+[scr])
