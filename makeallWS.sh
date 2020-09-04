cmsenv
cd test_df_MTR_2017_2020v1
root -b '../makeallWS.C("2017","MTR")'
root -b '../makeSWS.C("2017","MTR")'
cd ../test_df_MTR_2018_2020v1
root -b '../makeallWS.C("2018","MTR")'
root -b '../makeSWS.C("2018","MTR")'
cd ../test_df_VTR_2017_2020v1
root -b '../makeallWS.C("2017","VTR")'
root -b '../makeSWS.C("2017","VTR")'
cd ../test_df_VTR_2018_2020v1
root -b '../makeallWS.C("2018","VTR")'
root -b '../makeSWS.C("2018","VTR")'
cd ../
