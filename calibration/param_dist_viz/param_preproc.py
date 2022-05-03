# -*- coding: utf-8 -*-
from pathlib import Path
import pandas as pd
import seaborn as sns

def read_param_file(paths, currents, idx_dict):
    plist = [list()]*len(currents)    
    for p in paths:
        param_df = pd.read_excel(p)
        for n,c in enumerate(currents):
            plist[n].append(gen_param_dict(param_df[c], idx_dict[c]))

    plist = [pd.DataFrame(x) for x in plist]
    return plist

def gen_param_dict(p, idx):
    param_dict = {}
    for i in idx:
        param_dict[str(i)] = p.iloc[i-1]
    return param_dict

def compare_param_dist(pdf_wt, pdf_ko, idx):
    pdf_wt['Group'] = 'WT'
    pdf_ko['Group'] = 'Mgat1KO'
    pdf = pd.concat([pdf_wt, pdf_ko], ignore_index=True)

    p = sns.displot(pdf, x=str(idx), hue='Group', kind='kde', fill=True)
    p.set_xlabels('Parameter '+str(idx), fontsize=12)

def melt_param(pdf_wt, pdf_ko):
    pdf_wt['Group'] = 'WT'
    pdf_ko['Group'] = 'Mgat1KO'
    pdf = pd.concat([pdf_wt, pdf_ko], ignore_index=True)
    return pdf.melt(id_vars='Group', var_name='param')
    
exp_num =  '45'
currents = ['iktof', 'ikslow1', 'ikslow2', 'ikss']
idx_dict = {'iktof': [1, 2, 4, 5, 7, 11, 13],
    'iktos': [2, 3],
    'ikslow1': [1, 2, 4, 5, 9, 10, 11],
    'ikslow2': [2, 3],
    'ikss': [1, 2, 3, 4]}

base_dir = Path('/Users/haedongkim/Documents/Kv_data_assim/calibration')
wt_dir = base_dir / ('calib_exp'+exp_num+'_wt')
ko_dir = base_dir / ('calib_exp'+exp_num+'_ko')
fpaths_wt = wt_dir.glob('*.xlsx')
fpaths_ko = ko_dir.glob('*.xlsx')

plist_wt = read_param_file(fpaths_wt, currents, idx_dict)
plist_ko = read_param_file(fpaths_ko, currents, idx_dict)

pktof = melt_param(plist_wt[0], plist_ko[0])
pkslow1 = melt_param(plist_wt[1], plist_ko[1])
pkslow2 = melt_param(plist_wt[2], plist_ko[2])
pkss = melt_param(plist_wt[3], plist_ko[3])

pktof.to_csv('exp'+exp_num+'_pktof.csv', index=False)
pkslow1.to_csv('exp'+exp_num+'_pkslow1.csv', index=False)
pkslow2.to_csv('exp'+exp_num+'_pkslow2.csv')
pkss.to_csv('exp'+exp_num+'_pkss.csv', index=False)

# for v in idx_dict['kto']:
#     compare_param_dist(pkto_wt, pkto_ko, v)
