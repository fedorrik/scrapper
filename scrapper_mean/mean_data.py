import pandas as pd
from collections import Counter


def get_mean_df(df_path):
    data = pd.read_csv(df_path, sep='\t').set_index('drug')
    drugs = [' '.join(i.split()[:-2]) for i in data.index]
    drugs = list(Counter(drugs).keys())
    mean_data = pd.DataFrame()
    for drug in drugs:
        drug_data_big = data.loc[['{} rep{} MIC'.format(drug, str(i)) for i in range(1, 4)]]
        drug_data_small = data.loc[['{} rep{} MIC/2'.format(drug, str(i)) for i in range(1, 4)]]
        mean_big = pd.DataFrame(drug_data_big.mean(), columns=[drug + ' MIC']).T
        mean_small = pd.DataFrame(drug_data_small.mean(), columns=[drug + ' MIC/2']).T
        mean_data = pd.concat([mean_data, mean_big, mean_small])
    mean_data.index.name = 'drug'
    mean_data.to_csv('data/{}'.format(df_path.split('/')[-1]), sep='\t')

def get_mean_df_dead(df_path):
    data = pd.read_csv(df_path, sep='\t').set_index('drug')
    drugs = [' '.join(i.split()[:-2]) for i in data.index]
    drugs = list(Counter(drugs).keys())
    mean_data = pd.DataFrame()
    for drug in drugs:
        drug_data_zero = data.loc[['{} rep{} C0'.format(drug, str(i)) for i in range(1, 4)]]
        drug_data_big = data.loc[['{} rep{} MIC'.format(drug, str(i)) for i in range(1, 4)]]
        drug_data_small = data.loc[['{} rep{} MIC/2'.format(drug, str(i)) for i in range(1, 4)]]
        mean_zero = pd.DataFrame(drug_data_zero.mean(), columns=[drug + ' C0']).T
        mean_big = pd.DataFrame(drug_data_big.mean(), columns=[drug + ' MIC']).T
        mean_small = pd.DataFrame(drug_data_small.mean(), columns=[drug + ' MIC/2']).T
        mean_data = pd.concat([mean_data, mean_zero, mean_big, mean_small])
    mean_data.index.name = 'drug'
    mean_data.to_csv('data/{}'.format(df_path.split('/')[-1]), sep='\t')

get_mean_df('../data/data.z.tsv')
get_mean_df('../data/data.fc.tsv')
#???('data.fc-prblms.tsv')
get_mean_df_dead('../data/data.dead.tsv')
