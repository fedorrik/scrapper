import numpy as np
import pandas as pd
from FlowCytometryTools import FCMeasurement
from os import listdir
from sys import exit


def process_datadir_fc(drug, datadir, tubes, threshold):
    # take all fcs files in datadir except F8 and E2 (PNC1)
    fcs_files = [fcs for fcs in listdir(datadir) if fcs[-4:] == '.fcs' and fcs.split('.')[0].split('-')[-1] in tubes]
    # take all tubes protein_id in datadir
    protein_ids = sorted(list(set([fcs_file.split('-')[2][:-4] for fcs_file in fcs_files])))
    # init df with all fc
    drug_fc_df = pd.DataFrame(index=protein_ids, columns=['{} MIC/2'.format(drug), '{} MIC'.format(drug)])
    # read autoflu fsc files
    autoflu_fcs_df_zero = FCMeasurement(ID='1', datafile=datadir+'01-Well-F7.fcs').data.rename({'PE-A': 'PE'}, axis='columns').query('PE < {}'.format(threshold))
    autoflu_fcs_df_small = FCMeasurement(ID='2', datafile=datadir+'02-Well-F7.fcs').data.rename({'PE-A': 'PE'}, axis='columns').query('PE < {}'.format(threshold))
    autoflu_fcs_df_big = FCMeasurement(ID='3', datafile=datadir+'03-Well-F7.fcs').data.rename({'PE-A': 'PE'}, axis='columns').query('PE < {}'.format(threshold))
    # calculate autoflu ratios
    autoflu_ratio_zero = np.median(autoflu_fcs_df_zero['FITC-A'] / autoflu_fcs_df_zero['FSC-A'])
    autoflu_ratio_small = np.median(autoflu_fcs_df_small['FITC-A'] / autoflu_fcs_df_small['FSC-A'])
    autoflu_ratio_big = np.median(autoflu_fcs_df_big['FITC-A'] / autoflu_fcs_df_big['FSC-A'])
    # memorising problemic tubes
    empty_control = []
    empty_small = []
    empty_big = []
    # iterating
    for protein_id in protein_ids:
        # creating names of 3 files with protein_id (3 concentrations)
        file_zero = '01-Well-{}.fcs'.format(protein_id)
        file_small = '02-Well-{}.fcs'.format(protein_id)
        file_big = '03-Well-{}.fcs'.format(protein_id)
        # reading these files using FlowCytometryTools library
        fcs_df_zero = FCMeasurement(ID='1', datafile=datadir+file_zero).data.rename({'PE-A': 'PE'}, axis='columns').query('PE < {}'.format(threshold))
        fcs_df_small = FCMeasurement(ID='2', datafile=datadir+file_small).data.rename({'PE-A': 'PE'}, axis='columns').query('PE < {}'.format(threshold))
        fcs_df_big = FCMeasurement(ID='3', datafile=datadir+file_big).data.rename({'PE-A': 'PE'}, axis='columns').query('PE < {}'.format(threshold))
        # calculating fitca/fsca ratios for each cell (smth like concentration of protein)
        ratio_zero = np.median(fcs_df_zero['FITC-A'] / fcs_df_zero['FSC-A'])
        ratio_small = np.median(fcs_df_small['FITC-A'] / fcs_df_small['FSC-A'])
        ratio_big = np.median(fcs_df_big['FITC-A'] / fcs_df_big['FSC-A'])
        # comparing ratios with autoflu level
        # if autoflu is close to signal or autoflu > signal => remove value (grey in heatmap)
        protein_name = df_proteins.loc[protein_id, 'protein']
        if autoflu_ratio_zero > ratio_zero * 0.9:
            empty_control.append(protein_name)
        if autoflu_ratio_small > ratio_small * 0.9:
            empty_small.append(protein_name)
        if autoflu_ratio_big > ratio_big * 0.9:
            empty_big.append(protein_name)
        # calculating fc using
        fc_small = ratio_small / ratio_zero
        fc_big = ratio_big / ratio_zero
        # append fc to drug_fc_df
        drug_fc_df.loc[protein_id] = fc_small, fc_big
    return drug_fc_df, empty_control, empty_small, empty_big

# Z-score instead of fold change
def ztest_medians(a, b, a_autoflu, b_autoflu):
    x1 = np.median(a) - np.median(a_autoflu)
    x2 = np.median(b) - np.median(b_autoflu)
    q1 = (np.std(a)+np.std(a_autoflu))/np.sqrt(len(a))
    q2 = (np.std(b)+np.std(b_autoflu))/np.sqrt(len(b))
    z = (x1 - x2) / np.sqrt(q1**2 + q2**2)
    return z

def process_datadir_z(drug, datadir, tubes, threshold):
    # take all fcs files in datadir except F8 and E2 (PNC1)
    fcs_files = [fcs for fcs in listdir(datadir) if fcs[-4:] == '.fcs' and fcs.split('.')[0].split('-')[-1] in tubes]
    # take all tubes protein_id in datadir
    protein_ids = sorted(list(set([fcs_file.split('-')[2][:-4] for fcs_file in fcs_files])))
    # init df with all zscores
    drug_zscores_df = pd.DataFrame(index=protein_ids, columns=['{} MIC/2'.format(drug), '{} MIC'.format(drug)])
    # read autoflu fsc files
    autoflu_fcs_df_zero = FCMeasurement(ID='1', datafile=datadir+'01-Well-F7.fcs').data.rename({'PE-A': 'PE'}, axis='columns').query('PE < {}'.format(threshold))
    autoflu_fcs_df_small = FCMeasurement(ID='2', datafile=datadir+'02-Well-F7.fcs').data.rename({'PE-A': 'PE'}, axis='columns').query('PE < {}'.format(threshold))
    autoflu_fcs_df_big = FCMeasurement(ID='3', datafile=datadir+'03-Well-F7.fcs').data.rename({'PE-A': 'PE'}, axis='columns').query('PE < {}'.format(threshold))
    # calculate autoflu ratios
    autoflu_ratio_zero = autoflu_fcs_df_zero['FITC-A'] / autoflu_fcs_df_zero['FSC-A']
    autoflu_ratio_small = autoflu_fcs_df_small['FITC-A'] / autoflu_fcs_df_small['FSC-A']
    autoflu_ratio_big = autoflu_fcs_df_big['FITC-A'] / autoflu_fcs_df_big['FSC-A']
    # iterating
    for protein_id in protein_ids:
        # creating names of 3 files with protein_id (3 concentrations)
        file_zero = '01-Well-{}.fcs'.format(protein_id)
        file_small = '02-Well-{}.fcs'.format(protein_id)
        file_big = '03-Well-{}.fcs'.format(protein_id)
        # reading these files using FlowCytometryTools library
        fcs_df_zero = FCMeasurement(ID='1', datafile=datadir+file_zero).data.rename({'PE-A': 'PE'}, axis='columns').query('PE < {}'.format(threshold))
        fcs_df_small = FCMeasurement(ID='2', datafile=datadir+file_small).data.rename({'PE-A': 'PE'}, axis='columns').query('PE < {}'.format(threshold))
        fcs_df_big = FCMeasurement(ID='3', datafile=datadir+file_big).data.rename({'PE-A': 'PE'}, axis='columns').query('PE < {}'.format(threshold))
        # calculating fitca/fsca ratios for each cell (smth like concentration of protein)
        ratio_zero = fcs_df_zero['FITC-A'] / fcs_df_zero['FSC-A']
        ratio_small = fcs_df_small['FITC-A'] / fcs_df_small['FSC-A']
        ratio_big = fcs_df_big['FITC-A'] / fcs_df_big['FSC-A']
        # calculating z-score using ztest_medians function
        zscore_small = ztest_medians(ratio_small, ratio_zero, autoflu_ratio_small, autoflu_ratio_zero)
        zscore_big = ztest_medians(ratio_big, ratio_zero, autoflu_ratio_big, autoflu_ratio_zero)
        # append zscore to drug_zscores_df
        drug_zscores_df.loc[protein_id] = zscore_small, zscore_big
    return drug_zscores_df

def process_datadir_dead(drug, datadir, tubes, threshold):
    # take all fcs files in datadir except F8 and E2 (PNC1)
    fcs_files = [fcs for fcs in listdir(datadir) if fcs[-4:] == '.fcs' and fcs.split('.')[0].split('-')[-1] in tubes]
    # take all tubes protein_id in datadir
    protein_ids = sorted(list(set([fcs_file.split('-')[2][:-4] for fcs_file in fcs_files])))
    # df
    live_percent_df = pd.DataFrame()
    # iterating
    for protein_id in protein_ids:
        # creating names of 3 files with protein_id (3 concentrations)
        file_zero = '01-Well-{}.fcs'.format(protein_id)
        file_small = '02-Well-{}.fcs'.format(protein_id)
        file_big = '03-Well-{}.fcs'.format(protein_id)
        # reading these files using FlowCytometryTools library
        fcs_df_zero = FCMeasurement(ID='1', datafile=datadir+file_zero).data.rename({'PE-A': 'PE'}, axis='columns')
        fcs_df_small = FCMeasurement(ID='2', datafile=datadir+file_small).data.rename({'PE-A': 'PE'}, axis='columns')
        fcs_df_big = FCMeasurement(ID='3', datafile=datadir+file_big).data.rename({'PE-A': 'PE'}, axis='columns')
        # count percent of live cells and put in df
        if len(fcs_df_zero) > 0:
            live_percent_df.loc[protein_id, '{} C0'.format(drug)] = len(fcs_df_zero.query('PE >= {}'.format(threshold))) / len(fcs_df_zero)
        if len(fcs_df_small) > 0:
            live_percent_df.loc[protein_id, '{} MIC/2'.format(drug)] = len(fcs_df_small.query('PE >= {}'.format(threshold))) / len(fcs_df_small)
        if len(fcs_df_big) > 0:
            live_percent_df.loc[protein_id, '{} MIC'.format(drug)] = len(fcs_df_big.query('PE >= {}'.format(threshold))) / len(fcs_df_big)
    return live_percent_df

def delete_from_data(data, drugs_to_delete):
    if len(drugs_to_delete) > 0:
        data = data.drop([drug + ' MIC' for drug in drugs_to_delete])
        data = data.drop([drug + ' MIC/2' for drug in drugs_to_delete])
        # drop C0 from dead df
        if  'C0' in ' '.join(data.index):
            data = data.drop([drug + ' C0' for drug in drugs_to_delete])
    return data
    

# read PE-A threshold
with open('parameters.txt') as f:
    threshold = [line.strip() for line in f if line[0] != '#'][-1]

# associate proteins with plate tubes
df_proteins = pd.read_csv('plate.txt', sep=' ', index_col=0, names=['tube', 'protein'])
tubes = df_proteins.index

# read all_drugs.txt
with open('all_drugs.txt') as f:
	all_drugs = [' '.join(line.strip().split()[:-1]) for line in f]

# read drugs_to_add.txt
with open('drugs_to_add.txt') as f:
    drugs_to_add = {}
    for line in f:
        if line[0] == '#' or line == '':
            continue
        drug_name, drug_path = line.strip().split('\t')
        drugs_to_add[drug_name] = drug_path

# check if any new drag already exist in a data
drugs_to_delete = []
drugs_to_skip = []
for drug_name in drugs_to_add:
    if drug_name in all_drugs:
        print('{} already exist in the data files. Would you like to a.overwrite or b.skip? a/b'.format(drug_name))
        while True:
            answer = input('Type a or b and press Enter: ')
            if answer == 'a':
                # add drug to list to remove it from data
                drugs_to_delete.append(drug_name)
                break
            elif answer == 'b':
                # add drug to list to remove it from drugs_to_add dict
                drugs_to_skip.append(drug_name)
                break
            else:
                pass

# remove drug which want to skip from drugs_to_add dict
if len(drugs_to_skip) > 0:
    for drug in drugs_to_skip:
        drugs_to_add.pop(drug)

# exit script if nothing to add
if len(drugs_to_add) == 0:
    print('Nothing to add')
    exit()

### FOLD CHANGE ###
# read databases
fc_data = pd.read_csv('data/data.fc.tsv', sep='\t').set_index('drug')
problems_data = pd.read_csv('data/data.fc-prblms.tsv', sep='\t').set_index('drug')
# delete drugs
fc_data = delete_from_data(fc_data, drugs_to_delete)
problems_data = delete_from_data(problems_data, drugs_to_delete)
# problemic dicts
empty_control_dict = {}
empty_small_dict = {}
empty_big_dict = {}
# counting
for drug_name in drugs_to_add:
    drug_path = drugs_to_add[drug_name]
    # add slash at the end of path
    if drug_path[-1] != '/':
        drug_path += '/'
    # processing files
    fc_data_to_add = df_proteins.copy()
    fc_df, empty_control, empty_small, empty_big = process_datadir_fc(drug_name, drug_path, tubes, threshold)
    fc_data_to_add = pd.concat([fc_data_to_add, fc_df], axis=1)
    fc_data_to_add = fc_data_to_add.set_index('protein').T
    # log fc data
    fc_data_to_add = fc_data_to_add.fillna(0).apply(np.log2)
    # add processed data to database
    fc_data = pd.concat([fc_data, fc_data_to_add])
    print('{} added to data.fc'.format(drug_name))
    # fill problemic dicts
    empty_control_dict[drug_name] = empty_control
    empty_small_dict[drug_name] = empty_small
    empty_big_dict[drug_name] = empty_big

# problems df
problems_df = fc_data_to_add.copy()
problems_df.loc[:] = np.nan
for drug in empty_control_dict:
    for protein_name in empty_control_dict[drug]:
        if protein_name in empty_small_dict[drug]: # both control and experement lack protein
            problems_df.loc[drug+' MIC/2', protein_name] = 0
        else: # only control lacks protein
            problems_df.loc[drug+' MIC/2', protein_name] = 1
        if protein_name in empty_small_dict[drug]: # both control and experement lack protein
            problems_df.loc[drug+' MIC', protein_name] = 0
        else: # only control lacks protein
            problems_df.loc[drug+' MIC', protein_name] = 1
            
    for protein_name in empty_small_dict[drug]:
        if protein_name not in empty_control_dict[drug]: #  only experement lacks protein
            problems_df.loc[drug+' MIC/2', protein_name] = -1
    for protein_name in empty_small_dict[drug]:
        if protein_name not in empty_control_dict[drug]: #  only experement lacks protein
            problems_df.loc[drug+' MIC', protein_name] = -1

# add problem data to database
problems_data = pd.concat([problems_data, problems_df])
# write to databases
fc_data.to_csv('data/data.fc.tsv', sep='\t', index_label='drug')
problems_data.to_csv('data/data.fc-prblms.tsv', sep='\t', index_label='drug')


### Z-SCORE ###
# read data.tsv
data = pd.read_csv('data/data.z.tsv', sep='\t').set_index('drug')
# delete drugs
data = delete_from_data(data, drugs_to_delete)
# counting
for drug_name in drugs_to_add:
    drug_path = drugs_to_add[drug_name]
    # add slash at the end of path
    if drug_path[-1] != '/':
        drug_path += '/'
    # processing files
    data_to_add = df_proteins.copy()
    drug_zscores_df = process_datadir_z(drug_name, drug_path, tubes, threshold)
    data_to_add = pd.concat([data_to_add, drug_zscores_df], axis=1)
    data_to_add = data_to_add.set_index('protein').T
    # merge to data
    data = pd.concat([data, data_to_add])
    print('{} added to data.z'.format(drug_name))
# write to data.tsv
data.to_csv('data/data.z.tsv', sep='\t', index_label='drug')


### DEAD CELLS ###
# read data.tsv
data = pd.read_csv('data/data.dead.tsv', sep='\t').set_index('drug')
# delete drugs
data = delete_from_data(data, drugs_to_delete)
# counting
for drug_name in drugs_to_add:
    drug_path = drugs_to_add[drug_name]
    # add slash at the end of path
    if drug_path[-1] != '/':
        drug_path += '/'
    # processing files
    data_to_add = df_proteins.copy()
    drug_zscores_df = process_datadir_dead(drug_name, drug_path, tubes, threshold)
    data_to_add = pd.concat([data_to_add, drug_zscores_df], axis=1)
    data_to_add = data_to_add.set_index('protein').T
    # merge to data
    data = pd.concat([data, data_to_add])
    print('{} added to data.dead'.format(drug_name))
# write to data.tsv
data.to_csv('data/data.dead.tsv', sep='\t', index_label='drug')
# write drug names
with open('all_drugs.txt', 'w') as f:
    f.write('\n'.join(data.index))

