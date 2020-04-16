import pandas as pd

regions = ['Zmumu', 'Zee', 'Wmunu', 'Wenu', 'SR']


for region in regions:
   df = pd.read_csv(region+'/QCD.csv', sep = '\t')
   df['content'] = 0
   df['error'] =  0
   df.to_csv(region+'/QCD.csv', sep = '\t')
