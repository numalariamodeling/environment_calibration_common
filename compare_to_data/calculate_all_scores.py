##### Import required packages #####
# standard packages
import sys
from pathlib import Path
from math import sqrt, exp
import os
from statistics import mean
import numpy as np
import pandas as pd
import warnings
from datetime import datetime
import math
from scipy.special import gammaln
from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=SettingWithCopyWarning)
pd.options.mode.chained_assignment = None  # default='warn'
# from within environment_calibration_common submodule
# sys.path.append("../")
# from load_inputs import load_sites
# from helpers import load_coordinator_df
# from source 'simulations' directory
sys.path.append("../../simulations")

import manifest

def load_coordinator_df(characteristic=False, set_index=True):
    csv_file = manifest.sweep_sim_coordinator_path if characteristic else manifest.simulation_coordinator_path
    coord_df = pd.read_csv(csv_file)
    if set_index:
        coord_df = coord_df.set_index('option')
    return coord_df

coord_df = load_coordinator_df(characteristic=False, set_index=True)

#site = coord_df.at['site','value']

start_year = int(coord_df.at['simulation_start_year','value'])

def load_case_data(site):
    case_df= pd.read_csv(os.path.join(manifest.base_reference_filepath,
                                      coord_df.at['incidence_comparison_reference','value']))
    case_df=case_df[case_df['site']==site]                  
    return(case_df)

def load_prevalence_data(site):
    ### Load reference PCR prevalence data
    refprev = pd.read_csv(os.path.join(manifest.base_reference_filepath,
                                      coord_df.at['prevalence_comparison_reference','value']))
    refprev = refprev[refprev['site']==site]
    return(refprev)

def load_prevalence_data_U2(site):
    ### Load reference PCR prevalence data
    refprev = pd.read_csv(os.path.join(manifest.base_reference_filepath,
                                      coord_df.at['prevalence_comparison_reference_U2','value']))
    refprev = refprev[refprev['site']==site]
    return(refprev)

def prepare_inset_chart_data_PCR(site):
    ### Load InsetChart.csv
    ic = pd.read_csv(Path(os.path.join(manifest.simulation_output_filepath,site, "InsetChart_PCR.csv")))
    # convert Time to year and month
    ic['year'] = [np.trunc(t/365) for t in ic['time']]
    ic['month'] = ic.apply(lambda row: np.trunc((row['time'] - (row['year']*365))/30.001)+1, axis=1)
    ic = ic[ic['month']<=12]
    return(ic)
  
def prepare_inset_chart_data_EIR(site):
    ### Load InsetChart.csv
    ic = pd.read_csv(Path(os.path.join(manifest.simulation_output_filepath,site, "InsetChart_EIR.csv")))
    # convert Time to year and month
    ic['year'] = [np.trunc(t/365) for t in ic['time']]
    ic['month'] = ic.apply(lambda row: np.trunc((row['time'] - (row['year']*365))/30.001)+1, axis=1)
    ic = ic[ic['month']<=12]
    return(ic)  
  
def compare_all_age_PCR_prevalence(site):
    #### SCORE - Monthly PCR Parasite Prevalence ####
    #################################################
    ic = prepare_inset_chart_data_PCR(site)
    refpcr = load_prevalence_data(site)
    refpcr_U2 = load_prevalence_data_U2(site)
    # get average prevalence from insetchart for each year/month/sample_id 
    score3 = ic.groupby(['Sample_ID','year','month']).agg(prevalence=('PCR Parasite Prevalence', 'mean')).reset_index()
    # convert reference_pcr 'year' to start at 0, like simulations
    refpcr['year'] = [(y - start_year) for y in refpcr['year']]
    refpcr_U2['year'] = [(y - start_year) for y in refpcr_U2['year']]
    # merge with reference data on ['year','month']
    score3 = score3[['year', 'month', 'Sample_ID', 'prevalence']]
    score3 = score3.merge(refpcr,on=["year","month"], how='left')
    score4 = score3.merge(refpcr_U2,on=["year","month"], how='left')
    # score as sqrt(|diff|)
    score3['prevalence_score'] = score3.apply(lambda row: sqrt(abs(row['prevalence']-row['ref_prevalence'])**2), axis=1)
    score4['prevalence_U2_score'] = score4.apply(lambda row: sqrt(abs(row['prevalence']-row['ref_prevalence'])**2), axis=1)
    # get average score by sample_id
    score3 = score3.groupby('Sample_ID').agg(prevalence_score=('prevalence_score','mean')).reset_index()
    score4 = score4.groupby('Sample_ID').agg(prevalence_score=('prevalence_U2_score','mean')).reset_index()
    return score3, score4
  
def compare_PfPR_prevalence(site,agebin):
    #### SCORE - Monthly PCR Parasite Prevalence ####
    #################################################
    # load reference data
    if agebin <= 2:
        refpfpr = load_prevalence_data_U2(site)
        fname = "PfPR_monthly_U2.csv"
        score_name = "prevalence_U2_score"
    else:
        refpfpr = load_prevalence_data(site)
        fname = "PfPR_monthly_U5.csv"
        score_name = "prevalence_score"
    refpfpr = refpfpr[refpfpr['age'] == agebin]
    # convert reference_pcr 'year' to start at 0, like simulations
    refpfpr['year'] = [(y) for y in refpfpr['year']]
    sim_pfpr = pd.read_csv(Path(os.path.join(manifest.simulation_output_filepath,site,fname)))
    # filter to age of interest
    sim_pfpr = sim_pfpr[sim_pfpr['agebin']==agebin]
    # get mean PfPR by month, year, and Sample_ID across runs
    sim_pfpr = sim_pfpr.groupby(['Sample_ID', 'month','Year'])['PfPR'].agg(np.nanmean).reset_index()
    sim_pfpr['year']=sim_pfpr['Year']
    sim_pfpr = sim_pfpr[['year','month','Sample_ID','PfPR']]

    # merge simulated normalized monthly incidence with reference data on ['month']
    score = sim_pfpr.merge(refpfpr, on =['month','year'], how="left")
    # score as sqrt(|diff|)
    score['prevalence_score'] = score.apply(lambda row: sqrt(abs(row['PfPR']-row['prevalence'])**2), axis=1)
    # get average score by sample_id
    score = score.groupby('Sample_ID').agg(prevalence_score=('prevalence_score','mean')).reset_index()
    score = score.rename(columns={'prevalence_score': score_name})
    return score


def check_EIR_threshold(site):
    ic = prepare_inset_chart_data_EIR(site)
    # get maximum year of simulation
    last_year= max(ic['year'])
    
    #### SCORE - Maximum monthly EIR threshold ####
    ###############################################
    # filter insetchart to last 10 years of simulation
    score4 = ic[ic['year'] >= last_year-10]
    # sum daily to monthly EIR
    score4 = score4.groupby(['Sample_ID','Run_Number','month','year']).agg(EIR=('Daily EIR','sum')).reset_index()
    # Average across Run Numbers
    score4 = score4.groupby(['Sample_ID','year','month']).agg(EIR=('EIR','mean')).reset_index()
    # Find min/max for parameter set
    score4 = score4.groupby('Sample_ID').agg({'EIR': ['min', 'max']}).droplevel(axis=1, level=0).reset_index()
    score4['eir_score'] = 1*((score4['max']>=100) | (score4['min']==0))
    #print(score4.to_string())
    #print(score4)
    return score4

def compute_incidence_likelihood(combined_df):
    df = combined_df
    df['ref.Trials'] = 1000
    df['sim.Trials'] = 1000
    df['ref.Observations'] = df['ref.Trials'] * df['norm_repincd']
    df['sim.Observations'] = df['sim.Trials'] * df['norm_simincd']
    df['ll'] = gammaln(df['ref.Trials'] + 1) + gammaln(df['sim.Trials'] + 2) - gammaln(df['ref.Trials'] + df['sim.Trials'] + 2) + gammaln(
        df['ref.Observations'] + df['sim.Observations'] + 1) + gammaln(
        df['ref.Trials'] - df['ref.Observations'] + df['sim.Trials'] - df['sim.Observations'] + 1) - gammaln(df['ref.Observations'] + 1) - gammaln(
        df['ref.Trials'] - df['ref.Observations'] + 1) - gammaln(df['sim.Observations'] + 1) - gammaln(df['sim.Trials'] - df['sim.Observations'] + 1)
    #print(ll)
    return df #mean

def compare_incidence_shape(site,agebin):
    #### Load incidence data
    case_df = load_case_data(site)
    # filter to DS_Name
    case_df = case_df[case_df['site']==site]
    # filter to age of interest
    case_df = case_df[case_df['age']==agebin]
    # convert case_df 'year' to start at 0, like simulations
    case_df['year'] = [(y - start_year) for y in case_df['year']]
    #print(case_df)
    # find max incidence by year    
    case_df=case_df.merge(case_df.groupby('year')['case'].agg(np.nanmax).reset_index(name='max_incd'), on='year',how='left')
    # normalize monthly incidence so each year maxes at 1
    case_df['norm_repincd'] = case_df.apply(lambda row: (row['case'] / row['max_incd']), axis=1)
    # get average normalized incidence per month
    rcases = case_df.groupby(['month','site'])['norm_repincd'].agg(np.nanmean).reset_index()
    
    sim_cases = pd.read_csv(Path(os.path.join(manifest.simulation_output_filepath,site,"ClinicalIncidence_monthly.csv")))
    # filter to age of interest
    sim_cases = sim_cases[sim_cases['agebin']==agebin]
    # get mean population and clinical cases by month, year, and Sample_ID
    sim_cases['Inc'] = sim_cases['Cases'] / sim_cases['Pop']
    sim_cases=sim_cases.merge(sim_cases.groupby(['Sample_ID','Year'])['Inc'].agg(np.nanmax).reset_index(name='max_simincd'), on=['Sample_ID',"Year"],how='left').reset_index()
    sim_cases['norm_simincd'] = sim_cases.apply(lambda row: (row['Inc'] / row['max_simincd']), axis=1)
    #print(sim_cases)
    # get mean normalized incidence by month/sample_id (across years)
    sim_cases = sim_cases.groupby(['Sample_ID', 'month'])['norm_simincd'].agg(np.nanmean).reset_index()
    # merge simulated normalized monthly incidence with reference data on ['month']
    score1 = sim_cases.merge(rcases, on ='month')
    score1 = score1.dropna(subset=['norm_simincd']).reset_index()

    ll = compute_incidence_likelihood(score1)
    ll = ll.groupby(['Sample_ID'])['ll'].agg(np.nansum).reset_index()
    score1 = score1.drop(columns=['ll']).merge(ll,on='Sample_ID')
    #print(score1)
    score1=score1.groupby(['Sample_ID'])['ll'].agg(np.nanmean).reset_index()
    score1.rename(columns={'ll': 'shape_score'}, inplace=True)
    score1['shape_score'] = score1['shape_score'].abs()
    #print(score1)
    #print(score1)
    return score1
  
def compare_annual_incidence(site,agebin):
    ##### Score - Mean annual cases per person #####
    ################################################
    site_df = pd.read_csv(manifest.site_coordinator_path)
    site_df = site_df[site_df['site'] == site]
    # Reference data
    rcases = load_case_data(site)
    rcases = rcases[rcases['site']==site]
    rcases['incidence']=rcases['case'] / site_df['population']      # assuming population of 10,0000
    rcases = rcases.groupby(['age','year'])[['incidence']].agg(np.nansum).reset_index()
    target=rcases[rcases['age']==agebin]
    target=target['incidence'].mean()
    
    ### Load analyzed monthly MalariaSummaryReport from simulation
    sim_cases = pd.read_csv(Path(os.path.join(manifest.simulation_output_filepath,site,"ClinicalIncidence_monthly.csv")))
    # filter to age
    sim_cases = sim_cases[sim_cases['agebin']==agebin]
    sim_cases['Inc'] = sim_cases['Cases'] #/ sim_cases['Pop']
    # Average annual incidence in each year present
    sim_cases = sim_cases.groupby(['Sample_ID', 'Year','agebin'])['Inc'].agg(np.nanmean).reset_index()
    # Average annual incidence across all years
    score2 = sim_cases.groupby(['Sample_ID','agebin'])['Inc'].agg(np.nanmean).reset_index()
    # Score e^(absolute difference) vs. target
    score2['intensity_score'] = score2.apply(lambda row: exp(abs(row['Inc']-target)/target), axis=1)
    return score2

def compute_all_scores(site,incidence_agebin=100,prevalence_agebin=5, prevalence_agebin_U2=2):
    # merge unweighted scores into one dataframe, and return
    
    scores = check_EIR_threshold(site)
    scores = scores[['Sample_ID','eir_score']]
    if(coord_df.at['incidence_comparison','value']):
        score1 = compare_incidence_shape(site,agebin=incidence_agebin)
        scores = scores.merge(score1[['Sample_ID','shape_score']], how='outer', on='Sample_ID')
        # score2 = compare_annual_incidence(site,agebin=incidence_agebin)
        # scores = scores.merge(score2[['Sample_ID','intensity_score']], how='outer', on='Sample_ID')
    if(coord_df.at['prevalence_comparison','value']):
        if(coord_df.at['prevalence_comparison_diagnostic','value']=="PCR"):
            score3 = compare_all_age_PCR_prevalence(site)
            scores = scores.merge(score3[['Sample_ID','prevalence_score']], how='outer', on='Sample_ID')
        if(coord_df.at['prevalence_comparison_diagnostic','value']=="Microscopy"):
            score3 = compare_PfPR_prevalence(site, agebin=prevalence_agebin)
            score4 = compare_PfPR_prevalence(site,agebin=prevalence_agebin_U2)
            scores = scores.merge(score3[['Sample_ID','prevalence_score']], how='outer', on='Sample_ID')
            scores = scores.merge(score4[['Sample_ID','prevalence_U2_score']], how='outer', on='Sample_ID')
    return scores


if __name__ == '__main__':

  coord_df = pd.read_csv(manifest.simulation_coordinator_path)
  site = coord_df.at['site', 'value']
  print(site)
  
  print(compare_incidence_shape(site,agebin=100))
#   print(compare_annual_incidence(site,agebin=100))
  print(compare_all_age_PCR_prevalence(site))
  print(check_EIR_threshold(site))
  
  
  