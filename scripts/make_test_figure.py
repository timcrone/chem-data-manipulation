import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

cols = ["mol1", "mol2", "sim", "qed1", "qed2"]
covid_df = pd.read_csv("data/covid.csv")
covid_model_df = pd.read_csv("model/test/covid_model_10_target_results.txt", sep=" ")
covid_model_df.columns = cols
covid_model_df = covid_model_df.merge(covid_df.copy(), left_on="mol1", right_on="smiles")
orgnl_model_df = pd.read_csv("model/test/newmodel_1_target_results.txt", sep=" ")
orgnl_model_df.columns = cols

covid_model_df['model'] = 'Team5'

orgnl_model_df = orgnl_model_df.merge(covid_df.copy(), left_on="mol1", right_on="smiles")
orgnl_model_df['model'] = 'Native'
covid_models = pd.concat([covid_model_df, orgnl_model_df])
covid_models['sim'] = covid_models.apply(lambda x: np.nan if x['mol1'] == x['mol2'] else x['sim'], axis=1)
covid_models['qed_delta'] = covid_models['qed2'] - covid_models['qed1']
covid_models = covid_models.query("qed_delta > 0.0")
data=covid_models.sort_values(by=['sim'], ascending=False)
data.drop(columns=['smiles']).to_csv("figures/test_results.csv")
import seaborn as sns
from matplotlib import gridspec
fig = plt.figure(figsize=(8, 6)) 
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3]) 
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])

sns.set(style="whitegrid")
ax1 = sns.boxplot(x="drug_name", y="sim", hue="model",data=covid_models.fillna(0.0), palette="Set3", ax=ax1)
ax2 = sns.boxplot(x="model", y="sim", hue="model",data=covid_models.fillna(0.0), palette="Set3", ax=ax0)
ax2.legend_.remove()
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90)
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=0)
plt.subplots_adjust(bottom=0.25)
ax1.set_title("Per drug target similarity")
ax2.set_title("Per model similarity")
fig.suptitle("Similarity comparison between original and optimized model")
fig.savefig('figures/target_by_drug_model.png')
