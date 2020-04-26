import pandas as pd
import matplotlib.pyplot as plt
import os

# parse train log

dfs_log = []
i = 0
with open("model/train_logs/hg2g_v3.log", 'r') as fh:
    lines = fh.readlines()
    lines = [line for line in lines if line.startswith('[') or line.startswith('learning')]
    chunked_lines = []
    chunked = []
    for line in lines:
        line = line.strip()
        if line.startswith('['):
            batch, rec = line.split("] ")
            batch = batch.replace("[","")
            rec = rec.split(', ')
            rec = [x for x in rec if ':' in x]
            rec_d = {x:y for x,y in [l.split(':') for l in rec]}
            rec_d['batch'] = batch
            chunked.append(rec_d)
            
        else:
            df_log = pd.DataFrame.from_dict(chunked)
            df_log['epoch'] = i
            df_log['learning_rate'] = line.split(':')[-1]
            dfs_log.append(df_log)
            chunked=[]
            i+=1

df_log = pd.concat(dfs_log)
df_log['loss'] = df_log['loss'].astype(float)
df_log['KL'] = df_log['KL'].astype(float)
df_log['epoch'] = df_log['epoch'].astype(int)
log_grp_stats = df_log.groupby("epoch")['loss','KL'].agg(['mean','std'])
print(log_grp_stats)   
log_grp_stats = log_grp_stats.reset_index(drop=True)
dfs = []
filepath = "model/validation/covid_results_v3/"
for file in os.listdir(filepath):
    df = pd.read_csv(os.path.join(filepath, file), sep=' ', header=None)
    df.columns = ['mol1','mol2','sim','qed1','qed2']
    df['epoch'] = file.split('.')[-1]
    dfs.append(df)

df_val = pd.concat(dfs)
df_val['epoch'] = df_val['epoch'].astype(int)
df_val['qed_delta'] = df_val['qed2'] - df_val['qed1']
val_grp_stats = df_val.groupby("epoch")['sim','qed_delta'].agg(['median','std'])
print(val_grp_stats)
val_grp_stats = val_grp_stats.reset_index(drop=True)
#print(val_grp_stats.sort_values(by=['epoch'], ascending=True)) 

plt.style.use('seaborn-whitegrid')
fig, ax = plt.subplots( nrows=1, ncols=2 ) 

#x =  val_grp_stats.epoch.values # np.linspace(0, 10, 50)
loss_mean = log_grp_stats['loss']['mean'].values
loss_std = log_grp_stats['loss']['std'].values
kl_mean = log_grp_stats['KL']['mean'].values
kl_std = log_grp_stats['KL']['std'].values
qed_std = val_grp_stats['qed_delta']['std'].values
qed_mean = val_grp_stats['qed_delta']['median'].values
epoch = val_grp_stats.index.values
sim_std = val_grp_stats['sim']['std'].values
sim_mean = val_grp_stats['sim']['median'].values
#y = np.sin(x) + dy * np.random.randn(50)
ax[0].errorbar(epoch, loss_mean, yerr=loss_std, fmt='.--r')
ax[0].errorbar(epoch, kl_mean, yerr=kl_std, fmt='.--k')
ax[1].errorbar(epoch, qed_mean, yerr=qed_std, fmt='.--r')
ax[1].errorbar(epoch, sim_mean, yerr=sim_std, fmt='.--k')
ax[0].set_title("mean loss (red), KL (black)\n per epoch")
ax[1].set_title("median qed_delta (red) \n similarity (black) per epoch")
#ax.plot([0,1,2], [10,20,3])

fig.savefig('figures/validation.png')   # save the figure to file
plt.close(fig)    # close the figure window