import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date


def draw_fc_heatmap(drugs_to_draw, scale):
    # remove C0
    drugs_to_draw = [i for i in drugs_to_draw if 'C0' not in i and 'control' not in i]
    n_drugs = len(drugs_to_draw)
    # read data.tsv, take selected drugs
    fc_data = pd.read_csv('data/data.fc.tsv', sep='\t').set_index('drug').loc[drugs_to_draw]
    problems_data = pd.read_csv('data/data.fc-prblms.tsv', sep='\t').set_index('drug').loc[drugs_to_draw]
    # plot heatmap
    pic_height = n_drugs*0.29
    datestamp = str(date.today()).replace('-', '')
    cbar_dict = {'location': 'top', 'pad': 0, 'aspect': 30, 'shrink': 0.5, 'label': 'Fold Change'}
    plt.figure(figsize=(16, pic_height))
    custom_palette = sns.color_palette([(0, 1, 1), (0.6, 0.6, 0.6), (1, 0, 1)])
    ax = sns.heatmap(fc_data, vmin=-scale, vmax=scale, center=0, cmap='coolwarm', cbar_kws=cbar_dict).set_facecolor('xkcd:grey')
    ax = sns.heatmap(problems_data.fillna(np.nan), vmin=-1, vmax=1, center=0, cmap=custom_palette, linewidths=0.5, linecolor='grey', ax=ax, cbar=False)
    plt.ylabel('')
    plt.xlabel('')
    plt.savefig('heatmaps/fc-heatmap.height{}.scale{}.{}.png'.format(n_drugs, scale, datestamp), facecolor='white', bbox_inches='tight', pad_inches=0.1)

def draw_z_heatmap(drugs_to_draw, scale):
    # remove C0
    drugs_to_draw = [i for i in drugs_to_draw if 'C0' not in i and 'control' not in i]
    n_drugs = len(drugs_to_draw)
    # read data.tsv, take selected drugs
    data = pd.read_csv('data/data.z.tsv', sep='\t').set_index('drug').loc[drugs_to_draw]
    # plot heatmap
    pic_height = n_drugs*0.29
    datestamp = str(date.today()).replace('-', '')
    cbar_dict = {'location': 'top', 'pad': 0, 'aspect': 30, 'shrink': 0.5, 'label': 'z-score'}
    plt.figure(figsize=(16, pic_height))
    heatmap = sns.heatmap(data.fillna(0), mask=data.isnull(), vmin=-scale, vmax=scale, center=0, cmap='coolwarm', linewidths=0.5, linecolor='grey', cbar_kws=cbar_dict).set_facecolor('xkcd:grey')
    plt.ylabel('')
    plt.xlabel('')
    plt.savefig('heatmaps/z-heatmap.height{}.scale{}.{}.png'.format(n_drugs, scale, datestamp), facecolor='white', bbox_inches='tight', pad_inches=0.1)

def draw_dead_heatmap(drugs_to_draw, scale):
    n_drugs = len(drugs_to_draw)
    # read data.tsv, take selected drugs
    data = pd.read_csv('data/data.dead.tsv', sep='\t').set_index('drug').loc[drugs_to_draw]
    # plot heatmap
    pic_height = n_drugs*0.29
    datestamp = str(date.today()).replace('-', '')
    cbar_dict = {'location': 'top', 'pad': 0, 'aspect': 30, 'shrink': 0.5, 'label': 'dead cells, %'}
    plt.figure(figsize=(16, pic_height))
    heatmap = sns.heatmap(data*100, vmin=0, vmax=scale, center=scale/2, cmap='rocket_r', linewidths=0.5, linecolor='grey', cbar_kws=cbar_dict)
    plt.ylabel('')
    plt.xlabel('')
    plt.savefig('heatmaps/dead-heatmap.height{}.scale{}.{}.png'.format(n_drugs, scale, datestamp), facecolor='white', bbox_inches='tight', pad_inches=0.1)


with open('drugs_to_draw.txt') as f:
    # list with drugs to draw
    drugs_to_draw = [line.strip() for line in f if line[0] != '#']
with open('draw_parameters.txt') as f:
    draw_parameters = [line.strip() for line in f if line[0] != '#']
    # list with heatmaps to draw
    heatmaps_to_draw = draw_parameters[:-4]
    # dict with heatmap's scales
    heatmaps_scales = {i.split()[0]: int(i.split()[1]) for i in draw_parameters[-4:-1]}

# draw heatmaps
if 'fold change' in heatmaps_to_draw:
    print('drawing fold change heatmap...')
    draw_fc_heatmap(drugs_to_draw, heatmaps_scales['fc'])

if 'z-score' in heatmaps_to_draw:
    print('drawing z-score heatmap...')
    draw_z_heatmap(drugs_to_draw, heatmaps_scales['z'])

if 'dead percent' in heatmaps_to_draw:
    print('drawing dead percent heatmap...')
    draw_dead_heatmap(drugs_to_draw, heatmaps_scales['dead'])

print('Heatmaps are drawn. Check heatmaps folder')
