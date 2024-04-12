#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import seaborn as sns; import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from IPython.display import display_html
from itertools import chain, cycle
import scipy.stats as stats
import statsmodels.api as sm
import matplotlib.patches as patches


class Plotting:
    
    def __init__(self, study: str):
        self.parameter_labels = ['reward sensitivity','loss sensitivity','reward learning rate','loss learning rate', 'appetitive Pavlovian bias', 'aversive Pavlovian bias', 'noise', 'bias']
        self.parameter_labels_split = ['reward\nsensitivity','loss\nsensitivity','reward\nlearning rate','loss\nlearning rate', 'appetitive\nPavlovian bias', 'aversive\nPavlovian bias', 'noise', 'bias']
        if study == 'panda': self.time_label = ['baseline', 'week 2', 'week 6']
    
    
    def correlation_matrix_plot(self, data: pd.DataFrame, ax, fig, xticklabels=None, yticklabels=None):
        rho = data.corr()
        pval = data.corr(method=lambda x, y: stats.pearsonr(x, y)[1]) - np.eye(*rho.shape)
        p = pval.applymap(lambda x: ''.join(['*' for t in np.array([0.01,0.05,0.1]) /((len(data.columns)**2-1)/2) if x<=t]))

        matrix = np.triu(rho,k=1)
        annot = (rho.round(2).astype(str) + p).values
        for k in range(len(data.columns)): annot[k,k] = '1'
        if xticklabels == None: xticklabels = self.parameter_labels_split
        if yticklabels == None: yticklabels = self.parameter_labels_split
        for k in range(len(data.columns)): 
            sns.heatmap(rho,ax=ax, annot = annot, fmt='', cbar=False, mask=matrix,  vmin=-1, vmax=1, \
                        cmap="YlGnBu", xticklabels=xticklabels, yticklabels=yticklabels)
            

    def confusion_matrix_plot(self, matrix: np.array, ax, fig, xticklabels=None, yticklabels=None):
        if len(matrix.shape)==3:
            annot = np.diag([str(np.round(np.mean(matrix[i,i,:]),2)) + '\n±' +  str(np.round(np.std(matrix[i,i,:]),2)) for i in range(matrix.shape[0])])
            matrix = np.mean(matrix,axis=2)
        else:
            annot = np.round(np.diag(matrix),2).astype('str')
            annot[annot=='0.0']=''
        if xticklabels == None: xticklabels = self.parameter_labels_split
        if yticklabels == None: yticklabels = self.parameter_labels_split
        sns.heatmap(matrix, annot = annot, fmt='', ax=ax,vmin=-1, vmax=1,cmap="YlGnBu",center=0, xticklabels=xticklabels, yticklabels=yticklabels)
            
            
    def error_shaded_line_plot(self, data_to_plot, opt_plot):
        if "excluding" in opt_plot: ex = opt_plot['excluding']
        else: ex = np.full([np.shape(data_to_plot)[2]], False)
        fig, ax = plt.subplots(1,np.shape(data_to_plot)[1],figsize=(20, 10))
        clrs = sns.color_palette("husl", 5)
        with sns.axes_style("darkgrid"):
            for i in range(4):
                epochs = np.arange(np.shape(data_to_plot)[0])
                meanst = np.nanmean(data_to_plot[:,i,~ex],axis=1)
                sdt = np.nanstd(data_to_plot[:,i,~ex],axis=1)/np.sqrt(np.shape(data_to_plot)[2])
                ax[i].plot(epochs, meanst, c=clrs[i],linewidth=3)
                ax[i].fill_between(epochs, meanst - sdt, meanst + sdt, alpha=0.3, facecolor=clrs[i])
                if "surr_to_plot" in opt_plot:
                    ax[i].plot(epochs, np.nanmean(opt_plot["surr_to_plot"][:,i,~ex],axis=1),linewidth=3,color=clrs[i],linestyle='--')
                ax[i].set_title(opt_plot["title"][i],size=18,weight='bold')
                ax[i].set_xlabel(opt_plot["xlabel"],size=18)
                ax[i].set(ylim=(0, 1))
                if opt_plot["format"]:
                    ax[i].axvline(24,0,1,linewidth=1,color='k',linestyle='--',alpha=0.3)
                    ax[i].axvline(48,0,1,linewidth=1,color='k',linestyle='--',alpha=0.3)
                if i == 0:
                    ax[i].set_ylabel(opt_plot["ylabel"],size=18)
        if "legend" in opt_plot: plt.legend(opt_plot["legend"], fontsize=12);
            

    def group_line_plot(self, data_to_plot, opt_plot):
        if "excluding" in opt_plot: ex = opt_plot['excluding']
        else: ex = np.full([np.shape(data_to_plot)[2]], False)
        fig, ax = plt.subplots(1,np.shape(data_to_plot)[1],figsize=(20, 10))
        clrs = sns.color_palette("husl", 5)
        with sns.axes_style("darkgrid"):
            for i in range(4):
                epochs = np.arange(np.shape(data_to_plot)[0])
                for g in range(2):
                    ax[i].plot(epochs, np.nanmean(data_to_plot[:,i,~ex & (opt_plot["group"] == g + 1)],axis=1), c=clrs[g*3])
                    if "surr_to_plot" in opt_plot:
                        ax[i].plot(epochs, np.nanmean(opt_plot["surr_to_plot"][:,i,~ex & (opt_plot["group"] == g + 1)],axis=1),linewidth=3,c=clrs[g*3],linestyle='--')
                    ax[i].set_title(opt_plot["title"][i],size=18,weight='bold')
                    ax[i].set_xlabel(opt_plot["xlabel"],size=18)
                    ax[i].set(ylim=(0, 1))
                if opt_plot["format"]:
                    ax[i].axvline(24,0,1,linewidth=1,color='k',linestyle='--',alpha=0.3)
                    ax[i].axvline(48,0,1,linewidth=1,color='k',linestyle='--',alpha=0.3)
                if i == 0: ax[i].set_ylabel(opt_plot["ylabel"],size=18)
        if "legend" in opt_plot: plt.legend(opt_plot["legend"], fontsize=12);


    def goprobability_scatter_plot(self, data, opt_plot):
        if "excluding" in opt_plot: ex = opt_plot['excluding']
        else: ex = np.full([np.shape(data['a_go'])[2]], False)
        fig, ax = plt.subplots(1, 4, figsize=(20, 5))
        for i in range(4):
            ax[i].scatter(np.nanmean(data['a_go'][:, i, ~ex], axis=0),
                          np.nanmean(opt_plot["surr_to_plot"][:, i, ~ex], axis=0))
            ax[i].set_title(opt_plot["title"][i], size=18, weight='bold')
            ax[i].set_xlabel('go probability empirical', size=18)
            ax[i].set_ylabel('go probability model', size=18)
            ax[i].set(ylim=(0, 1), xlim=(0, 1))
            

    def over_trials_plot(self, fig, ax, data, surrdata, opt):
        cmapmine = opt['colormap']

        Nsj = opt['Nsj']
        for i in range(len(ax)):
            if opt['plot_individual_choices']:
                im = ax[i].imshow((data[:, i, :] == 1).astype(int).T, cmap=cmapmine, aspect='auto', origin='lower')
            epochs = np.arange(np.shape(data)[0])
            meanst = np.nanmean(data[:, i, :], axis=1)
            sdt = np.nanstd(data[:, i, :], axis=1) / np.sqrt(Nsj)
            p, = ax[i].plot(epochs, meanst * Nsj, linewidth=1, label=opt['legend'][0])
            ax[i].fill_between(epochs, (meanst - sdt) * Nsj, (meanst + sdt) * Nsj, alpha=0.3)
            ax[i].plot(epochs, np.nanmean(surrdata[:, i, :], axis=1) * Nsj, color=p.get_color(), linewidth=2, \
                       linestyle='--', label=opt['legend'][1])
            if i == 0:
                ax[i].set_yticks([1, np.round(Nsj / 4), np.round(Nsj / 2), np.round(Nsj / 4 * 3), Nsj])
                ax[i].set_yticklabels(opt['yticklabels'])
                ax[i].set_ylabel(opt['ylabel'])
            else:
                ax[i].set_yticks([1, np.round(Nsj / 4), np.round(Nsj / 2), np.round(Nsj / 4 * 3), Nsj])
                ax[i].set_yticks([])
            ax[i].set_title(opt['title'][i], weight='bold')
            ax[i].set_xlabel(opt['xlabel'])
            if opt['plot_chance_line']:
                ax[i].axhline(np.round(Nsj / 2), color='black', linewidth=2, linestyle='--', alpha=0.3)
        fig.tight_layout()

        if opt['plot_individual_choices']:
            values = {'go': 1, 'nogo': 0}
            colors = [im.cmap(im.norm(value)) for value in values.values()]
            pt = [patches.Patch(facecolor=colors[i], edgecolor='lightgrey', label=list(values.keys())[i]) for i in
                  range(len(values))]
            l0 = plt.legend(handles=pt, bbox_to_anchor=(opt['legend_location'][0], 1), loc=2, borderaxespad=0.2)
            plt.gca().add_artist(l0);

        return p
    

    def boxplot_over_time(self, fig, ax, data, opt):

        vals, names, xs = [], [], []
        for i, col in enumerate(data.columns):
            vals.append(data[col].dropna().values)
            names.append(col)
            xs.append(np.random.normal(i + 1, 0.04, data[col].dropna().values.shape[0]))

        p = ax.boxplot(vals)
        ax.set_xticks(opt['xticks'])
        ax.set_xticklabels(opt['xlabels'])
        palette = opt['colors']
        for x, val, c in zip(xs, vals, palette):
            ax.scatter(x, val, alpha=0.4, color=c)
        ax.tick_params(labelsize=18)
        ax.set_title(opt['title'], fontsize=20, weight='bold');
        ax.set_ylabel(opt['ylabel'], fontsize=18);
        for t in opt['vlines']:
            ax.axvline(t, 0, 1, linewidth=1, color='black', linestyle='--', alpha=0.1);
        fig.tight_layout()

        return p


    def display_side_by_side(self, *args, titles=cycle([''])):
        html_str = ''
        for df, title in zip(args, chain(titles, cycle(['</br>']))):
            html_str += '<th style="text-align:center"><td style="vertical-align:top">'
            html_str += f'<h2>{title}</h2>'
            html_str += df.to_html().replace('table', 'table style="display:inline"')
            html_str += '</td></th>'
        display_html(html_str, raw=True)
        
            
    def check_mle_residuals(self, model):
        fig,axes = plt.subplots(2,2, figsize = (16, 9))

        sns.distplot(model.resid, kde_kws = {"shade" : True, "lw": 1}, fit = stats.norm, ax=axes[0,0])
        axes[0,0].set_title("KDE Plot of Model Residuals (Blue) and Normal Distribution (Black)")
        axes[0,0].set_xlabel("Residuals")

        sm.qqplot(model.resid, dist = stats.norm, line = 's', ax=axes[0,1])
        axes[0,1].set_title("Q-Q Plot")

        sns.scatterplot(y = model.resid, x = model.fittedvalues, ax=axes[1,0])
        axes[1,0].set_title("RVF Plot")
        axes[1,0].set_xlabel("Fitted Values")
        axes[1,0].set_ylabel("Residuals")

        sns.boxplot(x = model.model.groups, y = model.resid, ax=axes[1,1])
        axes[1,1].set_title("Distribution of Residuals for parameter values")
        axes[1,1].set_ylabel("Residuals")
        axes[1,1].set_xlabel("parameter values")
        fig.tight_layout()
        

    def convert_pvalue_to_asterisks(self, pvalue: float) -> str:
        if pvalue <= 0.0001:
            return "****"
        elif pvalue <= 0.001:
            return "***"
        elif pvalue <= 0.01:
            return "**"
        elif pvalue <= 0.05:
            return "*"
        return "ns"
    

    def plot_significance(self, ax, x: np.array, text: str):
        ylim = ax.get_ylim()
        # ax.set_ylim(ylim[0], ylim[1] + max(ylim)/10)
        for i in range(len(text)):
            if '*' in text[i]:
                text_len = len(text[i])
            else:
                text_len = 2
            ax.text(x=x[i] - text_len * 0.05, y=max(ylim), s=text[i])
            ax.hlines(max(ylim), x[i] - 0.5, x[i] + 0.5, 'black')
            