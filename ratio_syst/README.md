The main script is this one `plotTF.py` (copy this and also the other one here, wherever you need it)

You need to change a path to wherever your `all_percategory.root` workspace file is. This is the one that gets created from the `makeWS_percategory.C` script 
Modify this line to the correct path `tf = transferFactorSys.TFSystematics("./all_percategory.root")`  
Even better, if you separate into folders, you can pass in the `sys.argv[1]` argument eg  `"./%s/all_percategory.root"%sys.argv[1]`

then, to run just do :
`python plotTF.py MTR_2017 5 0 2  "41.5 fb^{-1} (13 TeV, 2017)" "MTR"  "WJETS_WMUNU" "WJETS_SR"`

The last two are the processes to make the ratio for (in the example it looks for the Wjets in the single muon divided by the Wjets in the signal region)
