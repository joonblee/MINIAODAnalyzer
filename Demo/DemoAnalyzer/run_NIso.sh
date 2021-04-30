#!/bin/bash

zmass=(10 20 50 80)
nfile=100

#rm -rf tmp
mkdir tmp

rm tmp/ConfFile_NIso_*_cfg.py
rm tmp/QCD_zmass*.root
rm tmp/job_NIso_*
rm QCD_zmass*.root


cd ../..
eval `scram runtime -sh`
cd Demo/DemoAnalyzer

scramv1 b

for z in ${zmass[@]}
do
  #for i in {1..${nfile}}
  for ((i=1;i<=${nfile};i++))
  do
    echo ""
    echo "#####################################"
    echo "### M(Z) = $z ### run number = $i ###"
    echo "#####################################"
  
    sed -e "s/!RUN!/${i}/g" python/ConfFile_NIso_cfg.py > tmp/ConfFile_NIso_zmass${z}_${i}_cfg.py
    sed -i "s/!ZMASS!/${z}/g" tmp/ConfFile_NIso_zmass${z}_${i}_cfg.py

    echo "#!/bin/sh"                                      >  tmp/job_NIso_zmass${z}_${i}.sh
    echo "cmsRun tmp/ConfFile_NIso_zmass${z}_${i}_cfg.py" >> tmp/job_NIso_zmass${z}_${i}.sh
  
    echo "Executable = tmp/job_NIso_zmass${z}_${i}.sh"    >  tmp/submit_NIso_zmass${z}_${i}.jds
    echo "Universe   = vanilla"                           >> tmp/submit_NIso_zmass${z}_${i}.jds
    echo "getenv     = True"                              >> tmp/submit_NIso_zmass${z}_${i}.jds
    echo "Log        = tmp/job_NIso_zmass${z}_${i}.log"   >> tmp/submit_NIso_zmass${z}_${i}.jds
    echo "Output     = tmp/job_NIso_zmass${z}_${i}.out"   >> tmp/submit_NIso_zmass${z}_${i}.jds
    echo "Error      = tmp/job_NIso_zmass${z}_${i}.err"   >> tmp/submit_NIso_zmass${z}_${i}.jds
    echo "Queue"                                          >> tmp/submit_NIso_zmass${z}_${i}.jds
  
    chmod +x tmp/job_NIso_zmass${z}_${i}.sh
    condor_submit tmp/submit_NIso_zmass${z}_${i}.jds
  done

done

nRootFiles=`ls tmp/QCD_zmass*.root | grep -c "root"`
while [ $nRootFiles -ne $((${#zmass[@]} * ${nfile})) ]
do
  sleep 10
  nRootFiles=`ls tmp/QCD_zmass*.root | grep -c "root"`
  echo $nRootFiles
done

for z in ${zmass[@]}
do
  hadd QCD_zmass${z}.root tmp/QCD_zmass${z}_*.root 
done




#sed -e "s/___INPUTFILES___/NIso_zmass${z}\/${line}/g" python/ConfFile_cfg.py > tmp/ConfFile_NIso_zmass${z}_${i}_cfg.py
#sed -i "s/___OUTPUTFILES___/Outputs\/NIso_zmass${z}\/hist_${i}.root/g" tmp/ConfFile_NIso_zmass${z}_${i}_cfg.py

