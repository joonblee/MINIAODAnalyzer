Samples=("SingleMuon_C" "WJetsToLNu")
#Samples=("SingleMuon_C" "DYJetsToLL")
#Samples=("SingleMuon_C")

rm -rf tmp
mkdir tmp
mkdir Outputs

for Sample in ${Samples[@]}
do
  echo ""
  echo "#################################"
  echo "### Run analyzer ${Sample} ###"
  echo "#################################"

  rm -rf Outputs/${Sample}
  mkdir  Outputs/${Sample}

  ls -1 Inputs/${Sample}/*.root &> tmp/temp.txt
  sed -i 's/Inputs\/'${Sample}'\///g' tmp/temp.txt
  
  i=0
  input="tmp/temp.txt"
  
  while IFS= read -r line
  do
    i=$(($i+1))
    sed -e "s/___INPUTFILES___/${Sample}\/${line}/g" python/ConfFile_cfg.py > tmp/ConfFile_${Sample}_${i}_cfg.py
    sed -i "s/___OUTPUTFILES___/Outputs\/${Sample}\/hist_${i}.root/g" tmp/ConfFile_${Sample}_${i}_cfg.py
  done < "$input"
  
  scramv1 b
  
  
  i=0
  while IFS= read -r line
  do
    i=$(($i+1))
    echo "#!/bin/sh"                                 >  tmp/job_${Sample}_${i}.sh
    echo "cmsRun tmp/ConfFile_${Sample}_${i}_cfg.py" >> tmp/job_${Sample}_${i}.sh
  
    echo "Executable = tmp/job_${Sample}_${i}.sh"    >  tmp/submit_${Sample}_${i}.jds
    echo "Universe   = vanilla"                      >> tmp/submit_${Sample}_${i}.jds
    echo "getenv     = True"                         >> tmp/submit_${Sample}_${i}.jds
    echo "Log        = tmp/job_${Sample}_${i}.log"   >> tmp/submit_${Sample}_${i}.jds
    echo "Output     = tmp/job_${Sample}_${i}.out"   >> tmp/submit_${Sample}_${i}.jds
    echo "Error      = tmp/job_${Sample}_${i}.err"   >> tmp/submit_${Sample}_${i}.jds
    echo "Queue"                                     >> tmp/submit_${Sample}_${i}.jds
  
    chmod +x tmp/job_${Sample}_${i}.sh
    condor_submit tmp/submit_${Sample}_${i}.jds
  done < "$input"
  
done
