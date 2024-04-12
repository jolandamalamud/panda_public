import numpy as np
import pandas as pd
import scipy.io
import statsmodels.formula.api as smf
from scipy.stats import stats
from scipy.stats import chi2_contingency
from statsmodels.stats.diagnostic import het_white
from typing import Optional
import pingouin as pg


class GNGstats:

    def __init__(self):
        pass
    

    def ttest2(self, group1, group2, labels, grouplabels=None):
        mm = np.round(np.array([np.nanmean(group1, axis=1), np.nanmean(group2, axis=1)]), 2).T
        ss = np.round(np.array([np.nanstd(group1, axis=1), np.nanstd(group2, axis=1)]), 2).T
        if grouplabels==None: grouplabels = ['group1', 'group2']
        for i in range(len(labels)):
            stats = scipy.stats.ttest_ind(group1[i, :], group2[i, :], equal_var=False)
            print(grouplabels[0] + ': ' + str(mm[i, 0]) + ' (' + str(ss[i, 0]) + \
                  '), ' + grouplabels[1] + ': ' + str(mm[i, 1]) + ' (' + str(
                ss[i, 1]) + '), pval: ' + str(round(stats.pvalue, 3)) + ' -> ' + labels[i])
        pass
    

    def correlation(self, X, Y, xlabel, ylabel):
        for i in range(len(X[:, 0])):
            x = X[i,:]
            mask = ~np.isnan(x) & ~np.isnan(Y)
            stats = scipy.stats.linregress(x[mask], Y[mask])
            print(ylabel + ' and ' + xlabel[i] + ': ' + str(round(stats.slope, 2)) + ', pval: ' + str(
                round(stats.pvalue, 3)))

        pass
    
    
    # function comparing a variable between sertraline and placebo group (2-sample ttest, printing group mean (SD))
    def group_ttest(self, variable, group, labels):
        mask = ~variable.isna()
        group1 = variable[(group==1) & mask]
        group2 = variable[(group==0) & mask]
        stat = scipy.stats.ttest_ind(group1, group2, equal_var=False)
        tab = pd.concat((group1, group2), axis=1, ignore_index=True).mean().T
        tab = pd.concat((tab, pd.DataFrame([stat.statistic, stat.pvalue]))).round(3)
        tab.index = labels + ['t','p']      
        return tab.T
    

    def chi2_test(self, variable, group, labels):
        tab = pd.concat((variable[group==1].value_counts(), \
               variable[group==0].value_counts()),axis=1)
        tab.columns = labels
        stat, p, dof, expected = chi2_contingency(tab)
        tab.loc[:,'X2'] = stat
        tab.loc[:,'p'] = p
        return tab.round(3)
    

    # function running GLM and printing desired output (beta and pvalue) or whole summary
    def glm(self, formula, data, var_of_interest):
        model = smf.glm(formula, data, missing='drop').fit()
        if var_of_interest:
            for i in var_of_interest:
                print(i + ': \tbeta: ' + str(np.round(model.params[i], 2)) \
                      + ',\tCI: [' + str(np.round(model.conf_int().loc[i][0], 2)) \
                      + ',' + str(np.round(model.conf_int().loc[i][1], 2)) + ']' \
                      + ',\tpvalue: ' + str(np.round(model.pvalues[i], 4)))

        return model
    

    # function running MLE and printing desired output (beta and pvalue) or whole summary
    def mle(self, formula, data, var_of_interest, random_effects):
        if random_effects:
            model = smf.mixedlm(formula, data, re_formula='1', vc_formula={random_effects: '0 + ' + random_effects}, \
                                groups='subject', missing='drop').fit()
            if model.converged == False: 
                print('model not converged -> random effects removed:')
                model = smf.mixedlm(formula, data, groups='subject', missing='drop').fit()
        else:
            model = smf.mixedlm(formula, data, groups='subject', missing='drop').fit()
        if var_of_interest:
            for i in var_of_interest:
                print(i + ': \tbeta: ' + str(np.round(model.params[i], 2)) \
                      + ',\tCI: [' + str(np.round(model.conf_int().loc[i][0], 2)) \
                      + ',' + str(np.round(model.conf_int().loc[i][1], 2)) + ']' \
                      + ',\tpvalue: ' + str(np.round(model.pvalues[i], 4)))

        return model
    

    def model_comparison(self, datapath: str, models: list) -> list:
        ibic = []
        for m in models:
            modelling = self.load_modelling_results(datapath, [m])
            ibic.append(modelling['ibic'])
        return ibic
    
    
    def check_mle_residuals(self, model):
        print('normality:')
        labels = ["Statistic", "p-value"]
        norm_res = stats.shapiro(model.resid)
        for key, val in dict(zip(labels, norm_res)).items():
            print(key, val)
            
        print('homoskedasticity:')
        het_white_res = het_white(model.resid, model.model.exog)
        labels = ["LM Statistic", "LM-Test p-value", "F-Statistic", "F-Test p-value"]
        for key, val in dict(zip(labels, het_white_res)).items():
            print(key, val)
            
            
    def calculate_trt_reliability(self, df: pd.DataFrame, columns: list, time: list, condition: Optional[list]=None) -> pd.DataFrame:
        trt_reliability = pd.DataFrame()
        for j in columns:
            x = [j + str(i) for i in time]
            if condition is not None:
                y = df.loc[condition, x].corr().to_numpy()
            else:
                y = df[x].corr().to_numpy()
            trt_reliability[j] = [y[0,1],y[0,2], y[1,2]]
        return trt_reliability.rename(index={0: 'baseline to T = ' + str(time[1]) + 'w', 1: 'baseline to T = ' + str(time[2]) + 'w', 2: 'T = ' + str(time[1]) + 'w to T = ' + str(time[2]) + 'w'})
    
    
    def normalize(self, df):
        return (df-df.mean())/df.std(ddof=0)
    
    
    def calculate_icc(self, df: pd.DataFrame, columns: list, time: list) -> pd.DataFrame:
        tmp = []
        for i in columns:
            Br = pd.melt(df[[i + str(t) for t in time]].reset_index(), id_vars='index')
            icc_df = pg.intraclass_corr(data=Br, targets='index', raters='variable', ratings='value', nan_policy='omit')
            tmp.append(icc_df.iloc[2])
        return pd.DataFrame(tmp, index=columns).drop(['Type', 'Description'], axis=1)


    def calculate_preregistered_hypothesis(self, mle_df: pd.DataFrame(), df: pd.DataFrame(), subsample: pd.DataFrame(), covariates: list):
    
        included_subjects = mle_df['exclusion']==0
        timing = [mle_df['time']<2, (mle_df['time']==0)|(mle_df['time']==2), mle_df['time']<3]
        timing_label = ['2','6','over time', 'time x group']
        
        print('\nH1: aversive Pav related to sertraline?')
        p = 'av_Pav'
        for i in range(3):
            tmp = (df['exclusion'+str(i)]==0)&subsample[:len(df)]
            print('T = '+ (['0'] + timing_label)[i]+', placebo: ' + str(np.round(df[p + str(i)][tmp&(df['group']==0)].mean(),2)) + '(' \
            + str(np.round(scipy.stats.sem(df[p + str(i)][tmp&(df['group']==0)]),2)) + ') sertraline: ' \
            + str(np.round(df[p + str(i)][tmp&(df['group']==1)].mean(),2)) + '(' \
            + str(np.round(scipy.stats.sem(df[p + str(i)][tmp&(df['group']==1)]),2)) +')', end='\t')
            if i > 0:      
                self.mle(p + ' ~ group + time' + covariates, mle_df[timing[i-1]&included_subjects&subsample], \
                            ['group'],[])

        print('over time', end='\t\t\t\t\t')
        self.mle(p + ' ~ group + time' + covariates, \
                mle_df[included_subjects&subsample], ['group'],[])

        print('time group interaction', end='\t\t\t\t')
        self.mle(p + ' ~ group * time' + covariates, \
                mle_df[included_subjects&subsample], ['group:time'],[]);
        
        print('\nH2: aversive Pav related to anxiety?')
        for i in range(3):
            tmp = mle_df[timing[i]&included_subjects&subsample]
            print(timing_label[i] + ':', end='\t')
            self.mle('gadlog ~ ' + p + ' + group + time' + covariates, \
                    tmp, [p], p) 
            
        print('\nH4: can baseline aversive Pav predict anxiety at week 12?')
        p = 'av_Pav0'
        covariates_reduced = covariates.replace("+ morisky + tablet","")
        self.glm('gad3log ~ ' + p + ' + gad0log + group ' + covariates_reduced, \
                df[(df['exclusion0']==0)&subsample[:len(df)]].astype(float), [p]);
        
        print('\nH5: appetitive Pav related to anxiety?')
        p = 'app_Pav'
        for i in range(3):
            tmp = mle_df[timing[i]&included_subjects&subsample]
            print(timing_label[i] + ':', end='\t')
            self.mle('phqlog ~ ' + p + ' + group + time' + covariates, \
                    tmp, [p], p)

        print('\nH6: reward sensitivity related to anhedonia?')
        p = 'rew_se'
        for i in range(3):
            tmp = mle_df[timing[i]&included_subjects&subsample]
            print(timing_label[i] + ':', end='\t')
            self.mle('anhlog ~ ' + p + ' + group + time' + covariates, \
                tmp, [p], p)


    def calculate_exploratory_stats(self, mle_df: pd.DataFrame(), df_panda, subsample: pd.DataFrame(), covariates: list):
    
        included_subjects = mle_df['exclusion']==0
        
        print('\n1: drug effect on loss learning rate at week 2')
        self.mle('loss_LR ~ group + time' + covariates, \
                mle_df[(mle_df['time']<2)&included_subjects&np.tile(subsample,3)], ['group'],[])
        
        print('\n2: drug effect on change in loss learning rate between baseline and week 2')
        covariates_adjust = covariates.replace("+ morisky + tablet","+ adaptmorisky1 + tablet1")
        self.glm('loss_LR_slope01 ~ group' + covariates_adjust, \
                                df_panda[(df_panda[['exclusion0',\
                                                    'exclusion1']] == 0).all(axis=1)&subsample], ['group'])
        
        print('\n3: loss learing rate related to anxiety over time')
        self.mle('gadlog ~ loss_LR + group + time' + covariates, \
                mle_df[included_subjects&np.tile(subsample,3)], ['loss_LR'],[])
        
        print('\n4: change in aversive Pavlovian bias predicting treatment outcome (at week 12)')
        p = 'av_Pav'
        j = 'bdi'
        self.glm(j + '3log ~' + p + '_slope01 + ' + j + '0log + group' + covariates_adjust, \
                    df_panda[(df_panda[['exclusion0', 'exclusion1']] == 0).all(axis=1)&subsample], \
                [p + '_slope01'])
        self.glm(j + '3log ~' + p + '_slope01 * group +' +  j + '0log' + covariates_adjust, \
                    df_panda[(df_panda[['exclusion0', 'exclusion1']] == 0).all(axis=1)&subsample], \
                [p + '_slope01:group'])
        print('-' * 100)
        for i in range(2):
            print(['placebo','sertraline'][i], end=': ')
            self.glm(j + '3log ~' + p + '_slope01 +' +  j + '0log' + covariates_adjust, \
                        df_panda[(df_panda['group'] == i) & \
                                (df_panda[['exclusion0', 'exclusion1']] == 0).all(axis=1)&subsample], \
                    [p + '_slope01'])


                



