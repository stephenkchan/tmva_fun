#! /bin/sh                                                                                                                                                                                                                                   
alias e="emacs -nw"
alias ll="ls -lth"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
setupATLAS
lsetup "root 6.02.12-x86_64-slc6-gcc48-opt"
export PATH=$PATH:/n/atlasfs/atlascode/backedup/stchan/vhbb/tmva_fun/exec
