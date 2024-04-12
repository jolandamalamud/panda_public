import numpy as np
import pandas as pd
import scipy.io
import numpy as np
from typing import Optional

class PrepData:

    def __init__(self, study: str, filepath: str, psychiatric_questionnaire: list):
        self.filepath = filepath
        self.psychiatric_questionnaire = psychiatric_questionnaire
        self.study = study
        self.gng_conditions = ['g2w', 'g2a', 'ng2w', 'ng2a']
        self.gng_variables = ['exclusion', 'iL', 'goprotot', 'acctot'] + \
                ['gopro_' + i for i in self.gng_conditions] + ['acc_' + i for i in self.gng_conditions] + \
                ['switch_' + i for i in self.gng_conditions] + ['stay_' + i for i in self.gng_conditions]
        self.parameter_labels = ['rew_se', 'loss_se', 'rew_LR', 'loss_LR', 'app_Pav', 'av_Pav', 'noise', 'bias', 'rbias']
        if study == 'panda': self.fu_weeks = [0,2,6,12]


    def load_data(self, gngfilepath: Optional[str] = None) -> pd.DataFrame:
        if gngfilepath == None:
            gngfilepath = self.filepath + 'data/gng_data/' + self.study + '_datafile_date.mat'
        mat = scipy.io.loadmat(gngfilepath)
        D = mat["D"][0, :]
        return D
    

    def extract_data(self, D: np.array) -> dict:
        Nruns = len(D)
        correct_action = np.array([1, 1, 2, 2])
        outcome = np.array([1,0,1,0])
        data = {
            "subids": np.full(Nruns, np.nan),
            "sess": np.full(Nruns, np.nan),
            "stimuli": np.full([96, Nruns], np.nan),
            "a": np.full([96, Nruns], np.nan),
            "r": np.full([96, Nruns], np.nan),
            "a_correct": np.full([24, 4, Nruns], np.nan),
            "a_go": np.full([24, 4, Nruns], np.nan),
            "to_late": np.full([24,4, Nruns], np.nan),
            "rt": np.full([96, Nruns], np.nan),
            "rt_ss": np.full([24,4, Nruns], np.nan),
            "switch": np.full([4, Nruns], np.nan),
            "stay": np.full([4, Nruns], np.nan)
        }
        for sj in range(Nruns):
            data["subids"][sj] = D[sj][1]
            data["sess"][sj] = D[sj][0]
            s = D[sj][4]
            a = D[sj][5]
            r = D[sj][6]
            rt = D[sj][7]
            data["stimuli"][:, sj] = np.squeeze(s)
            data["a"][:, sj] = np.squeeze(a)
            data["r"][:, sj] = np.squeeze(r)
            data["rt"][:, sj] = np.squeeze(rt)
            for ss in range(4):
                ass = a[s == ss + 1]
                rss = r[s == ss + 1]
                i = ass != 3
                data["a_correct"][i, ss, sj] = (ass[i] == correct_action[ss]).astype(int)
                data["a_go"][i, ss, sj] = (ass[i] == 1).astype(int)
                data["to_late"][:, ss, sj] = ~i
                data["rt_ss"][i, ss, sj] = rt[(s == ss + 1)][i]
                ass = ass[:8]
                rss = rss[:8]
                if sum(rss[0:-1]==(outcome[ss]-1)) > 0:
                    data["switch"][ss, sj] = (ass[1:] != ass[0:-1])[rss[0:-1]==(outcome[ss]-1)][:8].mean()
                if sum(rss[0:-1]==(outcome[ss])) > 0:
                    data["stay"][ss, sj] = (ass[1:] == ass[0:-1])[rss[0:-1]==outcome[ss]][:8].mean()
        return data
    

    def calculate_outcome_probabilities(self, D: np.array, ex: Optional[np.array] = None):
        if ex == None:
            ex = np.ones((np.shape(D),1))
        rr = D['r'][:, ex]
        ss = D['stimuli'][:, ex]
        aa = D['a'][:, ex]
        Nsj = np.shape(rr)[1]
        nn = np.zeros([Nsj, 4])
        dd = np.zeros([Nsj, 4])
        cr = [1, 1, 2, 2]
        rew = [1, 0, 1, 0]
        for i in range(Nsj):
            for s in range(4):
                nn[i, s] = np.nansum((rr[ss[:, i] == s + 1, i] == rew[s]) & \
                                     (aa[ss[:, i] == s + 1, i] == cr[s]))
                dd[i, s] = np.sum(aa[ss[:, i] == s + 1, i] == cr[s])
        return nn, dd
    

    def merge_over_sessions(self, data: dict) -> dict:
        ids = data["subids"]
        sess = data["sess"]
        Nsj = len(np.unique(ids))
        data_merged = {
            "subids": np.unique(ids),
            "stimuli": np.full([288, Nsj], np.nan),
            "a": np.full([288, Nsj], np.nan),
            "a_correct": np.full([72, 4, Nsj], np.nan),
            "a_go": np.full([72, 4, Nsj], np.nan)
        }
        if 'surr' in data: data_merged['surr'] = np.full([72, 4, Nsj], np.nan)
        for sj, subid in enumerate(np.unique(ids)):
            subsess = sess[ids == subid]
            if subsess.size > 0:
                for c, i in enumerate(subsess.astype(int)):
                    idx = (ids == subid) & (sess == subsess[c])
                    data_merged["stimuli"][(i - 1) * 96:i * 96, sj] = np.squeeze(data["stimuli"][:, idx])
                    data_merged["a"][(i - 1) * 96:i * 96, sj] = np.squeeze(data["a"][:, idx])
                    data_merged["a_correct"][(i - 1) * 24:i * 24, :, sj] = np.squeeze(data["a_correct"][:, :, idx])
                    data_merged["a_go"][(i - 1) * 24:i * 24, :, sj] = np.squeeze(data["a_go"][:, :, idx])
                    if 'surr' in data_merged: data_merged['surr'][(i - 1) * 24:i * 24, :, sj] = np.squeeze(data["surr"][:, :, idx])
        return data_merged
    

    def create_questionnaire_item_list(self, questionnaires: list, df: pd.DataFrame):
        y = []
        new_item_list = []
        for q in questionnaires:
            for i in range(1, questionnaires[q] + 1):
                x = [x for x in df.columns if q + str(i) in x]
                y.append(x)
            item_list = [x for xs in y for x in xs]
            item_list = [x for x in item_list if 'bin' not in x]

        for i in item_list:
            if '12wk' in i: new_item_list.append(i[:len(q) + 1] + '_4')
            elif '6wk' in i: new_item_list.append(i[:len(q) + 1] + '_3')
            elif '2wk' in i: new_item_list.append(i[:len(q) + 1] + '_2')
            else: new_item_list.append(i[:4] + '_1')

        return item_list, new_item_list
    

    def create_df(self, data: pd.DataFrame, modelling: np.array, rct_data: pd.DataFrame, log_transfer: Optional[list] = None) -> pd.DataFrame:
        # create dataframe of all gng task variables
        number_of_sessions = max(data['sess']).astype(int)
        
        gng_columns = self.parameter_labels + self.gng_variables
        gng_data = np.vstack((modelling['most_parsimonious']['emmap'], modelling['random']['emmap'], modelling['most_parsimonious']['excluding'], modelling['most_parsimonious']['iL'], np.nanmean(data['a_go'],axis=(0,1)), \
                      np.nanmean(data['a_correct'],axis=(0,1)), np.nanmean(data['a_go'],axis=0), \
                      np.nanmean(data['a_correct'],axis=0), data['switch'], data['stay']))
        
        df_gng = pd.DataFrame()
        df_gng['subids'] = np.unique(rct_data['ID'])
        
        for t in range(1,number_of_sessions+1):
            session_columns = [c + str(t-1) for c in gng_columns]
            df_gng[session_columns] = ""
            for sj, subid in enumerate(df_gng['subids']):
                id_match = data['subids'] == subid
                idx = np.where(id_match & (data['sess']==t))
                if len(idx[0]) > 0:
                    df_gng.loc[sj, session_columns] = np.squeeze(gng_data[:, idx])
        
        # create dataframe of all rct variables
        df_rct = pd.DataFrame(columns = rct_data.columns)
        rct_id_match = np.full_like(df_gng['subids'], np.nan)
        
        for sj, subid in enumerate(df_gng['subids']):
            idx = np.where(rct_data['ID'] == subid)
            if len(idx[0]) > 0:
                df_rct.loc[sj] = \
                    rct_data.iloc[idx].values[0]
        
        if log_transfer != None:
            for i in log_transfer:
                df_rct[i + 'log'] = np.log10(df_rct[i]+1)
        
        # combine gng task and rct data
        df = pd.concat([df_gng, df_rct], axis=1)
        df = df.replace(r'^\s*$', np.nan, regex=True)
        no_group_data = df_rct['group'].isna()
        df = df[~no_group_data].reset_index()
        if self.study == 'panda': df['group'] = (df['group']==2).astype(int)
        return df
    

    def load_modelling_results(self, datapath: str, model_name: list) -> dict:
        complex_model = model_name[0]
        mat = scipy.io.loadmat(self.filepath + datapath + complex_model + '.mat')
        modelling = {
            "emmap": mat["E"],
            "emmap_variance": mat["V"],
            "ml": mat["stats"]["EML"]
        }
        if "bf" in mat:
            modelling["iL"] = mat["bf"][0, 0][1][0]
            modelling["ibic"] = mat["bf"][0, 0][-2][0]
        if len(model_name) > 1:
            simple_model = model_name[1]
            mat = scipy.io.loadmat(self.filepath + datapath + simple_model + '.mat')
            simple_iL = mat["bf"][0, 0][1][0]
            modelling["excluding"] = (modelling["iL"] - simple_iL) < 3
        return modelling
    
    
    def load_modelfits(self, models: Optional[list] = None) -> dict:
        if models == None: models = ['ll2b2a2epxb', 'llb']
        M = {'most_parsimonious': dict(), 'random': dict()}
        M['most_parsimonious'] = self.load_modelling_results('results/model_fits/', models)
        M['random'] = self.load_modelling_results('results/model_fits/', [models[1]])
        return M
    
    
    def transform_parameters(self, params: np.array) -> np.array:
        transform_params = np.zeros_like(params[:7,:])
        transform_params[:2,:] = np.exp(params[:2,:])
        transform_params[2:4,:] = 1/(1+np.exp(-params[2:4,:]))
        transform_params[4:6,:] = np.exp(params[4:6,:])
        transform_params[6,:] = 1/(1+np.exp(-params[6,:]))
        return transform_params
    

    def load_surrogate_data(self, datapath: str, model_name: str) -> np.array:
        mat = scipy.io.loadmat(self.filepath + datapath +  'SurrogateData_' + model_name + '.mat')
        surrogate_data = mat["surr"]
        return surrogate_data
    

    def extract_surrogate_data(self, data: dict, surr_data: np.array) -> dict:
        Nsj = np.shape(surr_data)[0]
        Nsamples = np.shape(surr_data)[1]
        Ntrials = np.shape(data["stimuli"])[0]
        surr_a = np.full([Ntrials, Nsamples, Nsj], np.nan)
        mean_surr_a = np.full([Ntrials, Nsj], np.nan)
        data["surr"] = np.full([int(Ntrials / 4), 4, Nsj], np.nan)
        for sj in range(Nsj):
            s = data["stimuli"][:, sj]
            a = data["a"][:, sj]
            if Ntrials == 96:
                i = (a != 3)
            else:
                i = ~np.isnan(s)
            for k in range(100):
                surr_a[i, k, sj] = np.squeeze(surr_data[sj, k][0])
            mean_surr_a[i, sj] = np.nanmean(surr_a[i, :, sj] == 1, axis=1)
            for ss in range(4):
                ii = s == (ss + 1)
                if Ntrials == 96:
                    data["surr"][:, ss, sj] = mean_surr_a[ii, sj]
                else:
                    data["surr"][0:sum(ii), ss, sj] = mean_surr_a[ii, sj]

        return data
    

    def load_rctdata(self, rctfilepath: Optional[str] = None) -> pd.DataFrame:
        if rctfilepath == None:
            rctfilepath = self.filepath + 'data/trial_data/gng_' + self.study + '_data.csv'
        dfRCT = pd.read_csv(rctfilepath)
        return dfRCT.iloc[:, 1:]
    

    def match_ids(self, data: dict, ids1: np.array, ids2: np.array) -> np.array:
        data_matched = np.full(len(ids1), np.nan)
        for sj, subid in enumerate(ids1):
            if sum(ids2 == subid) > 0:
                data_matched[sj] = float(data[ids2 == subid].iloc[0])
        return data_matched
    

    def create_mle_df(self, df: pd.DataFrame, time_varying_columns: list, time_constant_columns: list, weeks: Optional[list] = None) -> pd.DataFrame:
        mle_df = pd.DataFrame(columns=time_varying_columns + time_constant_columns)

        for c in time_varying_columns:
            mle_df[c] = pd.to_numeric(
                pd.concat([df[c + '0'], df[c + '1'], df[c + '2']], axis=0, ignore_index=True))
        for c in time_constant_columns:
            mle_df[c] = pd.concat([df[c], df[c], df[c]], axis=0, ignore_index=True)

        mle_df['subject'] = np.tile(np.arange(len(df)),3)
        mle_df['group'] = np.hstack((np.tile(0, len(df)), df['group'], df['group']))
        mle_df['time'] = np.hstack((np.tile(0, len(df)), np.tile(1, len(df)), np.tile(2, len(df))))
        
        if weeks != None: 
            mle_df['weeks'] = np.hstack((np.tile(weeks[0], len(df)), \
                                         np.tile(weeks[1], len(df)), \
                                         np.tile(weeks[2], len(df))))
        
        for i in self.psychiatric_questionnaire:
            mle_df[i + 'log'] = np.hstack((df[i + '0log'],df[i + '1log'],df[i + '2log']))
            
        return mle_df

    def parameter_change(self, df: pd.DataFrame()) -> pd.DataFrame():
        df_change_list = []
        x = np.arange(3)
        for i in self.parameter_labels:
            # between baseline and week 2
            df_change_list.append(pd.DataFrame({i + '_slope01': df[i + '1'] - df[i + '0']}))
            # between week 2 and week 6
            df_change_list.append(pd.DataFrame({i + '_slope12': df[i + '2'] - df[i + '1']}))
            # Compute slope over time points for each subject
            slope_list=[]
            for sj in range(len(df)):
                y = df[[i + str(0), i + str(1), i + str(2)]].iloc[sj]
                idx = ~y.isna()
                if sum(idx) > 1:
                    slope_list.append(np.polyfit(x[idx],y[idx],1)[0])
                else:
                    slope_list.append(np.nan)
            slope_column = pd.DataFrame({i + '_slope': slope_list})
            df_change_list.append(slope_column)

        # # Concatenate all DataFrames in the list
        df = pd.concat([df] + df_change_list, axis=1)
        return df.astype(float)
        
            

                
                



