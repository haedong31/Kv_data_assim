# -*- coding: utf-8 -*-
from pathlib import Path
import pandas as pd
import seaborn as sns

def read_param_file(paths, currents, idx_dict):
    plist = []
    for _ in range(len(currents)):
        plist.append([])
        
    for p in paths:
        param_df = pd.read_excel(p)
        for n,c in enumerate(currents):
            param_dict = gen_param_dict(p, param_df[c], idx_dict[c])
            plist[n].append(param_dict)

    plist = [pd.DataFrame(x) for x in plist]
    return plist

def gen_param_dict(p, pdf, idx):
    param_dict = {}
    for i in idx:
        param_dict[str(i)] = pdf.iloc[i-1]
    
    p_wo_ext = p.with_suffix('')
    param_dict['file'] = p_wo_ext.name
    return param_dict

def melt_param(pdf_wt, pdf_ko):
    pdf_wt['Group'] = 'WT'
    pdf_ko['Group'] = 'Mgat1KO'
    pdf = pd.concat([pdf_wt, pdf_ko], ignore_index=True)
    return pdf.melt(id_vars=['Group','file'], var_name='param')

def compare_param_dist(pdf_wt, pdf_ko, idx):
    pdf_wt['Group'] = 'WT'
    pdf_ko['Group'] = 'Mgat1KO'
    pdf = pd.concat([pdf_wt, pdf_ko], ignore_index=True)

    p = sns.displot(pdf, x=str(idx), hue='Group', kind='kde', fill=True)
    p.set_xlabels('Parameter '+str(idx), fontsize=12)
    
exp_num =  '45'
currents = ['iktof','ikslow1','ikslow2','ikss']
idx_dict = {'iktof': [1, 2, 4, 5, 7, 11, 13],
    'iktos': [2, 3],
    'ikslow1': [1, 2, 4, 5, 9, 10, 11],
    'ikslow2': [2, 3],
    'ikss': [1, 2, 3, 4]}

file_names = pd.read_excel('file_names_4half4.xlsx', index_col=None)

wt_dir = Path.cwd() / ('calib_exp'+exp_num) / 'wt'
ko_dir = Path.cwd() / ('calib_exp'+exp_num) / 'mgat1ko'

fpaths_wt = [Path(wt_dir/(p+'.xlsx')) for p in file_names['WT'].tolist()]
fpaths_ko = [Path(ko_dir/(p+'.xlsx')) for p in file_names['MGAT1KO'].tolist()]

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
