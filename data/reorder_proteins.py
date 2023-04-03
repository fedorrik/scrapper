import pandas as pd


def reorder(df_path):
    df = pd.read_csv(df_path, sep='\t').set_index('drug')
    new_order = ['Gcv2', 'His7', 'Trp3', 'Aro2', 'Lys1', 'Hom2', 'Aro3', 'Car1', 'Sam4', 'Sui2', 'Gcv3', 'Rnr4', 'Pdc1', 'Pyc2', 'Tps3', 'Tsl1', 'Pgm2', 'Ald4', 'Hxk1', 'Pda1', 'Idp1', 'Tal1', 'Hxt3', 'Pgk1', 'Pmi40', 'Hog1', 'Cwp1', 'Mdh2', 'Erg3', 'Erg10', 'Lsc1', 'Pdr5', 'Dur1', 'Uga1', 'Gpd2', 'Gpd1', 'Ahp1', 'Grx1', 'Uth1', 'Yhb1', 'Tsa1', 'Sod1', 'Trx2', 'Ssa2', 'Ssa1', 'Ssa4', 'Ypk1', 'Sti1', 'Bat2', 'Tdh1', 'Hsp12', 'Hsp104', 'Hsp78', 'Nth1', 'Sse1', 'Hsc82', 'Tps1', 'Nsr1', 'Hsp42', 'Aha1', 'Hsp26', 'Ydj1', 'Pre4', 'Htb2', 'Bmh1', 'BY4741']
    df = df.loc[:, new_order]
    df.to_csv(df_path, sep='\t')

reorder('data.z.tsv')
reorder('data.fc.tsv')
reorder('data.fc-prblms.tsv')
reorder('data.dead.tsv')
