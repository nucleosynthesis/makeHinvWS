#!/bin/bash 

DATE=$(date +"%m-%d-%y")
mkdir -p  ./vbfinv/2017/$DATE

text2workspace.py all_percategory.txt --channel-masks
combine all_percategory.root --run blind | tee limit.log

combine all_percategory.root -M FitDiagnostics  --saveShapes --saveWithUncertainties --setParameters mask_SR=1 --plots  --saveNormalizations


combineTool.py -M Impacts -d all_percategory.root -m 125 --doInitialFit --robustFit 1 -t -1 --rMin -0.5 --rMax 0.5   
combineTool.py -M Impacts -d all_percategory.root -m 125  --robustFit 1 -t -1 --rMin -0.5 --rMax 0.5  --doFits --parallel 4 
combineTool.py -M Impacts -d all_percategory.root -m 125  --robustFit 1 -t -1 --rMin -0.5 --rMax 0.5  -o impacts.json
plotImpacts.py -i impacts.json -o impacts

combine all_percategory.root -M MultiDimFit --algo grid --rMax 0.8 -t -1 --expectSignal 0 --points 20 


plot1DScan.py --main-label "Expected" higgsCombineTest.AsymptoticLimits.mH120.root
python ~/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py fitDiagnostics.root -g plots.root --skipFitB

#mv fitDiagnostics.root *.png impacts.pdf ./vbfinv/2017/$DATE
#cp fitDiagnostics.root *.png impacts.pdf ~/www/private/higgs/hinv/vbfinv/2017/$DATE/
