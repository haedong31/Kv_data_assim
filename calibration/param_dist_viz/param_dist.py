# -*- coding: utf-8 -*-
from pathlib import Path
import pandas as pd
import seaborn as sns


def read_param_file(paths, idx_dict):
    kto_row_list = []
    kslow1_row_list = []
    kslow2_row_list = []
    kss_row_list = []
    for p in paths:
        param_df = pd.read_excel(p)
        kto_row_list.append(gen_param_dict(param_df['ikto'], idx_dict['kto']))
        kslow1_row_list.append(gen_param_dict(param_df['ikslow1'], idx_dict['kslow1']))
        kslow2_row_list.append(gen_param_dict(param_df['ikslow2'], idx_dict['kslow2']))
        kss_row_list.append(gen_param_dict(param_df['ikss'], idx_dict['kss']))

    pkto = pd.DataFrame(kto_row_list)
    pkslow1 = pd.DataFrame(kslow1_row_list)
    pkslow2 = pd.DataFrame(kslow2_row_list)
    pkss = pd.DataFrame(kss_row_list)

    return (pkto, pkslow1, pkslow2, pkss)

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
    

exp_num =  '16'
idx_dict = {'kto': [1, 2, 6, 10, 13, 14, 15, 16, 17],
    'kslow1': [1, 2, 3, 4, 5, 8, 9, 11, 12, 13],
    'kslow2': [1, 3],
    'kss': [3, 4],
    'kur': [1, 3],
    'k1': [1, 3, 5, 7]}

file_dir_wt = Path('./calib_exp'+exp_num+'_wt')
file_dir_ko = Path('./calib_exp'+exp_num+'_ko')
file_paths_wt = file_dir_wt.glob('*.xlsx')
file_paths_ko = file_dir_ko.glob('*.xlsx')

(pkto_wt,pkslow1_wt,pkslow2_wt,pkss_wt) = read_param_file(file_paths_wt, idx_dict)
(pkto_ko,pkslow1_ko,pkslow2_ko,pkss_ko) = read_param_file(file_paths_ko, idx_dict)

# for v in idx_dict['kto']:
#     compare_param_dist(pkto_wt, pkto_ko, v)

pkto = melt_param(pkto_wt, pkto_ko)
pkslow1 = melt_param(pkslow1_wt, pkslow1_ko)
pkslow2 = melt_param(pkslow2_wt, pkslow2_ko)
pkss = melt_param(pkss_wt, pkss_ko)

# pkto.to_csv('pkto.csv', index=False)
pkslow1.to_csv('pkslow1.csv', index=False)
pkslow2.to_csv('pkslow2.csv', index=False)
pkss.to_csv('pkss.csv', index=False)
