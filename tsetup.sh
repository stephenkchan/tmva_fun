#! /bin/sh                                                                                                                                                                                                                                   
alias e="emacs -nw"
alias p="python"
alias ll="ls -lth"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
setupATLAS
#lsetup "root 6.04.16-x86_64-slc6-gcc49-opt" #from CxAODFramework which root
export PATH=$PATH:/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/exec
rcSetup -u
rcSetup Base,2.4.28