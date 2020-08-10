cmsenv
LD_LIBRARY_PATH=~/work/private/Hinv/workspace/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/lib:$LD_LIBRARY_PATH
cd test_df_MTR_2017_2020v1
root -b '../makeallWS.C("2017","MTR")'
cd ../test_df_MTR_2018_2020v1
root -b '../makeallWS.C("2018","MTR")'
cd ../test_df_VTR_2017_2020v1
root -b '../makeallWS.C("2017","VTR")'
cd ../test_df_VTR_2018_2020v1
root -b '../makeallWS.C("2018","VTR")'
cd ../
