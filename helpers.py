####################################
##### Import required packages #####
####################################

##### Standard packages #####
import os
import json
import itertools
import pandas as pd
import warnings
import sys
import math
from pathlib import Path
import numpy as np
from datetime import datetime
from functools import partial
import emod_api.config.default_from_schema_no_validation as dfs
import emod_api.demographics.Demographics as Demog
import emodpy_malaria.demographics.MalariaDemographics as Demographics
from emod_api.demographics.DemographicsTemplates import CrudeRate   
import emodpy_malaria.malaria_config as conf
import emodpy_malaria.malaria_config as malaria_config
from emodpy_malaria.malaria_config import set_drug_param
from emodpy_malaria.interventions.treatment_seeking import add_treatment_seeking
from emodpy_malaria.interventions.drug_campaign import add_drug_campaign
from emodpy_malaria.interventions.vaccine import add_scheduled_vaccine, add_triggered_vaccine
from emodpy_malaria.interventions.common import add_triggered_campaign_delay_event
from emodpy_malaria.interventions.usage_dependent_bednet import add_scheduled_usage_dependent_bednet 
from emodpy_malaria.reporters.builtin import (
    add_event_recorder,
    add_malaria_summary_report,
    add_report_node_demographics,
    add_report_vector_stats
)
##### from within environment_calibration_common submodule ##### 

##### from malaria_vaccdrug_campaigns import add_vaccdrug_campaign ##### 
# from source 'simulations' directory
sys.path.append("../simulations")
import manifest


##################################
##### EMODpy building blocks #####
##################################

##### EMODpy building blocks #####

def set_param_fn(config):
    coord_df=load_coordinator_df()
    
    config = conf.set_team_defaults(config, manifest)
    # Add vector species to simulation
    set_species_from_file(config, vector_file = coord_df.at['vector_filepath','value'])
    # Turn on/off InsetChart
    config.parameters.Enable_Default_Reporting = 1
    config.parameters.Simulation_Duration = 70 * 365 + 1
    config.parameters.x_Temporary_Larval_Habitat = 1
    config.parameters.Age_Initialization_Distribution_Type = "DISTRIBUTION_COMPLEX"
    config.parameters.Climate_Model = "CLIMATE_BY_DATA"
    config.parameters.Air_Temperature_Filename = os.path.join(
        "climate", "air_temperature_daily.bin"
    )
    config.parameters.Land_Temperature_Filename = os.path.join(
        "climate", "air_temperature_daily.bin"
    )
    config.parameters.Rainfall_Filename = os.path.join(
        "climate", "rainfall_daily.bin"
    )
    config.parameters.Relative_Humidity_Filename = os.path.join(
        "climate", "relative_humidity_daily.bin"
    )
    config.parameters.Enable_Demographics_Birth = 1
    config.parameters.Enable_Natural_Mortality = 1
    config.parameters.Custom_Individual_Events = ["Received_Treatment","Received_Vehicle",
                                                  "Received_SMC_Vaccine",
                                                  "Bednet_Using","Bednet_Discarded","Bednet_Got_New_One"]
    
    # SMC parameters
    drug_param_dict = {'drug_box_day': 2.0, 'drug_irbc_killing': 10.8, 'drug_hep_killing': 3.64}
    drug_box_day = drug_param_dict['drug_box_day']
    drug_irbc_killing = drug_param_dict['drug_irbc_killing']
    drug_hep_killing = drug_param_dict['drug_hep_killing']
    make_vehicle_drug(config,drug_box_day=drug_box_day, drug_irbc_killing=drug_irbc_killing,
                      drug_hep_killing=drug_hep_killing)
    return config


def set_param(simulation, param, value):
    """
    Set specific parameter value
    Args:
        simulation: idmtools Simulation
        param: parameter
        value: new value
    Returns:
        dict
    """
    return simulation.task.set_parameter(param, value)


def update_sim_random_seed(simulation, value):
    simulation.task.config.parameters.Run_Number = value
    return {"Run_Number": value}


def add_outputs(task, site):
    """
    Requesting reports/outputs to the task.
    """
    coord_df = load_coordinator_df(characteristic=False, set_index=True)
    simulation_years = int(coord_df.at['simulation_years','value'])
    sim_start_year = int(coord_df.at['simulation_start_year','value'])
    
    # By default, inset_chart is included 
    # change the value 'config.parameters.Enable_Default_Reporting = 1' to 0 to turn it OFF #
    # This covers the case where PCR prevalence is desired and is currently used to check EIR
    
    ### Incidence-related reports ###     
    if(coord_df.at['incidence_comparison','value']):
        incidence_df = pd.read_csv(os.path.join(manifest.base_reference_filepath,coord_df.at['incidence_comparison_reference','value']))
        incidence_agebins = sorted([float(a) for a in incidence_df['age'].unique()])
        first_year = int(incidence_df['year'].min()) - sim_start_year
        last_year = int(incidence_df['year'].max()) - sim_start_year + 1
        
        # Logical bounds on years to request
        if first_year < 0:
            first_year=0
        if last_year > simulation_years:
            last_year=simulation_years
        
        if(coord_df.at['incidence_comparison_frequency','value']=='monthly'):
            for year in range(first_year, last_year):
                start_day = 0 + 365 * year
                sim_year = sim_start_year + year
                add_malaria_summary_report(task,manifest,
                                           start_day=start_day,end_day=365 + year * 365,
                                           reporting_interval=30,age_bins=incidence_agebins,
                                           max_number_reports=13,pretty_format=True,
                                           filename_suffix=f"Monthly_incidence_{sim_year}")
                                           
        if(coord_df.at['incidence_comparison_frequency','value']=='annual'):
            add_malaria_summary_report(task,manifest,
                                       start_day=first_year*365,end_day=last_year*365+365,
                                       reporting_interval=365,age_bins=incidence_agebins,
                                       max_number_reports=last_year-first_year+1,
                                       pretty_format=True,
                                       filename_suffix=f"Yearly_incidence_{first_year+sim_start_year}_to_{last_year+sim_start_year}")
    ### Prevalence-related reports ### 
    if(coord_df.at['prevalence_comparison','value']):
        if(coord_df.at['prevalence_comparison_diagnostic','value']!="PCR"):
            prevalence_df = pd.read_csv(os.path.join(manifest.base_reference_filepath,
                                                    coord_df.at['prevalence_comparison_reference','value']))
            prevalence_agebins =  sorted([float(a) for a in prevalence_df['age'].unique()])
            first_year = int(prevalence_df['year'].min()) - sim_start_year
            last_year = int(prevalence_df['year'].max()) - sim_start_year + 1
            
            # Logical bounds on years to request
            if first_year < 0:
                first_year=0
            if last_year > simulation_years:
                last_year=simulation_years
                
            if(coord_df.at['prevalence_comparison_frequency','value']=='monthly'):
                for year in range(first_year, last_year):
                    start_day = 0 + 365 * year
                    sim_year = sim_start_year + year
                    add_malaria_summary_report(task,manifest,
                                               start_day=start_day,end_day=365 + year * 365,
                                               reporting_interval=30,age_bins=prevalence_agebins,
                                               max_number_reports=13,pretty_format=True,
                                               filename_suffix=f"Monthly_prevalence_{sim_year}")
                                               
            if(coord_df.at['prevalence_comparison_frequency','value']=='annual'):
                add_malaria_summary_report(task,manifest,
                                           start_day=first_year*365,end_day=last_year*365,
                                           reporting_interval=365,age_bins=prevalence_agebins,
                                           max_number_reports=last_year-first_year+1,
                                           pretty_format=True,
                                           filename_suffix=f"Yearly_prevalence_{first_year+sim_start_year}_to_{last_year+sim_start_year}")
        
    return
  

def build_camp(site, coord_df=None):
    """
    Build a campaign input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    """
    # create campaign object
    if coord_df is None:
        coord_df = pd.read_csv(manifest.simulation_coordinator_path)
        coord_df = coord_df.set_index('option')

    camp = build_standard_campaign_object(manifest)
    # === INTERVENTIONS === #

    # health-seeking
    if (not pd.isna(coord_df.at['CM_filepath','value'])) and (coord_df.at['CM_filepath','value'] != ''):
        hs_df = pd.read_csv(manifest.input_files_path / Path(coord_df.at['CM_filepath','value']))
    else:
        hs_df = pd.DataFrame()
  
    if not hs_df.empty:
        # case management for malaria
        add_health_seeking(camp,hs_df)
    
    # NMFs
    if (not pd.isna(coord_df.at['NMF_filepath','value'])) and (coord_df.at['NMF_filepath','value'] != ''):
        nmf_df = pd.read_csv(manifest.input_files_path / Path(coord_df.at['NMF_filepath','value']))
    else:
        nmf_df = pd.DataFrame()
    
    if (not pd.isna(coord_df.at['NMF_filepath','value'])) and (coord_df.at['NMF_filepath','value'] != ''):
        if not hs_df.empty:
            add_nmf_hs(camp, hs_df, nmf_df)
    
    # SMC
    if (not pd.isna(coord_df.at['SMC_filepath','value'])) and (coord_df.at['SMC_filepath','value'] != ''):
        smc_df = pd.read_csv(manifest.input_files_path / Path(coord_df.at['SMC_filepath','value']))
    else:
        smc_df = pd.DataFrame()
    if not smc_df.empty:
        add_smc(camp,smc_df)

    # ITNS
    itn_df = pd.DataFrame()
    if (not pd.isna(coord_df.at['ITN_filepath','value'])) and (coord_df.at['ITN_filepath','value'] != ''):
        if (not pd.isna(coord_df.at['ITN_age_filepath','value'])) and (coord_df.at['ITN_age_filepath','value'] != ''):
            if (not pd.isna(coord_df.at['ITN_season_filepath','value'])) and (coord_df.at['ITN_season_filepath','value'] != ''):
                itn_df = pd.read_csv(manifest.input_files_path / Path(coord_df.at['ITN_filepath','value']))
                itn_age = pd.read_csv(manifest.input_files_path / Path(coord_df.at['ITN_age_filepath','value']))
                itn_season = pd.read_csv(manifest.input_files_path / Path(coord_df.at['ITN_season_filepath','value']))
    
    if not itn_df.empty:
        # Distribute ITNs with age- and season-based usage patterns
        add_itns(camp,itn_df,itn_age,itn_season)

    return camp


def build_standard_campaign_object(manifest):
    import emod_api.campaign as campaign
    campaign.set_schema(manifest.schema_file)
    return campaign


def set_simulation_scenario(simulation, site, csv_path):
    # get information on this simulation setup from coordinator csv
    coord_df = pd.read_csv(csv_path)
    coord_df = coord_df.set_index('option')

    # === set up config === #
    # simulation duration
    simulation_duration = int(coord_df.at['simulation_years','value'])*365
    simulation.task.config.parameters.Simulation_Duration = simulation_duration
    simulation.task.config.parameters.Enable_Vital_Dynamics = 1
    demographics_filename = str(coord_df.at['demographics_filepath','value'])
    #print(demographics_filename)
    if demographics_filename and demographics_filename != 'nan':
        simulation.task.transient_assets.add_asset(manifest.input_files_path / Path(demographics_filename))
        simulation.task.config.parameters.Demographics_Filenames = [demographics_filename.rsplit('/',1)[-1]]
    simulation.task.config.parameters.Age_Initialization_Distribution_Type = 'DISTRIBUTION_COMPLEX'

    # === set up campaigns === #
    build_camp_partial = partial(build_camp, site=site, coord_df=coord_df)
    simulation.task.create_campaign_from_callback(build_camp_partial)
 
    return {"Site": site, 'csv_path': str(csv_path)}


set_simulation_scenario_for_matched_site = partial(set_simulation_scenario, csv_path=manifest.simulation_coordinator_path)
set_simulation_scenario_for_characteristic_site = partial(set_simulation_scenario, csv_path=manifest.sweep_sim_coordinator_path)


def build_demog():
    """
    This function builds a demographics input file for the DTK using emod_api.
    """
    coord_df = load_coordinator_df()
    
    demog= Demographics.from_template_node(lat=float(coord_df.at['lat','value']), 
                                           lon=float(coord_df.at['lon','value']), 
                                           pop=int(coord_df.at['pop','value']), 
                                           forced_id=1, 
                                           init_prev=float(coord_df.at['prev0','value']), 
                                           include_biting_heterogeneity=True)
    
    return demog


################################
##### Intervention Helpers #####
################################

##### Treatment Seeking - Symptomatic Malaria #####
def add_health_seeking(camp,hs_df):
    coord_df=load_coordinator_df()
    for r, row in hs_df.iterrows():
        sim_year=int(row['year'])-int(coord_df.at['simulation_start_year','value'])
        if sim_year >=0:
            sim_day=sim_year*365
            sim_day=sim_day+day_of_year(int(row['month']),int(row['day']),int(row['year']))
            
            cm_coverage_by_age =  {'trigger': str(row['trigger']),                                         
                                   'coverage': float(row['coverage']),                                     
                                   'agemin': float(row['age_min']),                                        
                                   'agemax': float(row['age_max']),
                                   'rate': float(row['rate'])}
            add_treatment_seeking(camp,                                                  
                                  start_day = sim_day,                                        
                                  duration = int(row['duration']),                                          
                                  drug=["Artemether","Lumefantrine"],                                             
                                  targets=[cm_coverage_by_age],                                                                
                                  broadcast_event_name="Received_Treatment")


##### Treatment Seeking - Nonmalaria Causes #####
def add_nmf_hs(camp, hs_df, nmf_df):
    # if no NMF rate is specified, assume all age groups have 0.0038 probability each day
    if nmf_df.empty:
        nmf_df = pd.DataFrame({'U5_nmf': [0.0038], 'adult_nmf': [0.0038]})
    elif nmf_df.shape[0] != 1:
        warnings.warn('The NMF dataframe has more than one row. Only values in the first row will be used.')
    nmf_row = nmf_df.iloc[0]
    # apply the health-seeking rate for clinical malaria to NMFs
    for r, row in hs_df.iterrows():
        add_nmf_hs_from_file(camp, row, nmf_row)


def add_nmf_hs_from_file(camp, row, nmf_row):
    if row['trigger'] == "NewClinicalCase":
        coord_df=load_coordinator_df()
        sim_year=int(row['year'])-int(coord_df.at['simulation_start_year','value'])
        if sim_year >=0:
            sim_day=sim_year*365
            sim_day=sim_day+day_of_year(int(row['month']),int(row['day']),int(row['year']))
            hs_coverage = float(row['coverage'])
            duration = int(row['duration'])
            hs_agemin = int(row['age_min'])
            hs_agemax = int(row['age_max'])
            if 'drug' in row.index:
                drug_code = "AL"
            else:
                drug_code = "AL"
            if sim_day == 0:  # due to dtk diagnosis/treatment configuration, a start day of 0 is not supported
                sim_day = 1  # start looking for NMFs on day 1 (not day 0) of simulation
                if duration > 1:
                    duration = duration - 1
            nmf_rate=0
            if hs_agemax > 5:
                nmf_rate = nmf_row['adult_nmf']
            if hs_agemax <= 5:
                nmf_rate = nmf_row['U5_nmf']
            # workaround for maximum duration of 1000 days is to loop, creating a new campaign every 1000 days
            separate_durations = [1000] * int(np.floor(duration/1000))  # create a separate campaign for each 1000 day period
            if (duration - np.floor(duration/1000) > 0):  # add final remaining non-1000-day duration
                separate_durations = separate_durations + [int(duration - np.floor(duration/1000) * 1000)]
            separate_start_days = sim_day + np.array([0] + list(np.cumsum(separate_durations)))
            for dd in range(len(separate_durations)):
                if nmf_rate > 0:
                    add_drug_campaign(camp, 'MSAT', drug_code=drug_code, start_days=[separate_start_days[dd]],
                                      target_group={'agemin': hs_agemin, 'agemax': hs_agemax},
                                      coverage=nmf_rate*hs_coverage,
                                      repetitions=separate_durations[dd], tsteps_btwn_repetitions=1,
                                      diagnostic_type='PF_HRP2', diagnostic_threshold=5,
                                      receiving_drugs_event_name='Received_NMF_Treatment')


##### SMC using vaccine-drug formulation #####
def make_vehicle_drug(config,drug_box_day: float = 0, drug_irbc_killing: float = 0, drug_hep_killing: float = 0):
    if drug_box_day:
        set_drug_param(config,drug_name="Vehicle",parameter="Drug_Decay_T1",value=drug_box_day)
        set_drug_param(config,drug_name="Vehicle",parameter="Drug_Decay_T2",value=drug_box_day)
    if drug_irbc_killing:
        set_drug_param(config,drug_name="Vehicle",parameter="Max_Drug_IRBC_Kill",value=drug_irbc_killing)
    if drug_hep_killing:
        set_drug_param(config,drug_name="Vehicle",parameter="Drug_Hepatocyte_Killrate",value=drug_hep_killing)
    return {'drug_box_day': drug_box_day,
            'drug_irbc_killing': drug_irbc_killing,
            'drug_hep_killing': drug_hep_killing}


def add_vaccdrug_smc(campaign,start_days: list, coverages: list,
                     target_group: dict = None,
                     node_ids: list = None,
                     receiving_vaccine_event: str = None, receiving_drugs_event: str = None,
                     listening_duration: int = -1, trigger_condition_list: list = None,
                     ind_property_restrictions: dict = None,
                     target_residents_only: int = 1,
                     check_eligibility_at_trigger: bool = False):
    """
        Add a vaccine + vehicle drug intervention to approximate efficacy of SMC. This intervention uses default
        parameters corresponding to SMC with SPAQ, if not otherwise specified via vaccine_param_dict and drug_param_dict.
        The vehicle drug instantly clear parasites (blood stage + liver stage) and the prophylactic effect is added
        by the vaccine event.

        Campaign type specifications:
        For SMC, the drug event (MDA drug campaign) is initiated at the specified simdays and generates a broadcast event
        that is used to trigger the vaccine event without delay.

    Args:
        campaign: campaign object to which the intervention will be added, and schema_path container
        start_days: list of days on which to run the interventions
        coverages: list of coverages for each day in start_days
        target_group: A dictionary of to specify age range for SMC
            Default is Everyone.

             Example::

                 {'agemin': x, 'agemax': y} for campaign_type = SMC

        node_ids: The list of nodes to apply this intervention to (**Node_List**
            parameter). If not provided, set value of NodeSetAll.
        receiving_vaccine_event:  Event to send out when person received vaccine.
            Default: 'Received_<campaign_type>_VaccDrug'
        receiving_drugs_event: Event to send out when person received drugs.
            Event name needs to include 'Received_Vehicle' in it, as otherwise overwritten in drug_campaigns function
            (see drug_campaigns.py L247)
            Default: SMC: 'Received_Vehicle'
        listening_duration: Length of time, in days, for which the triggered event will be listening for the triggers
        trigger_condition_list: List of events that will begin a triggerable
            campaign if campaign_type is SMC. If campaign_type is PMC, campaign is triggered by birth per default.
        ind_property_restrictions:  List of IndividualProperty key:value pairs that
            individuals must have to receive the diagnostic intervention.
            For example, ``[{"IndividualProperty1":"PropertyValue1"},
            {"IndividualProperty2":"PropertyValue2"}]``. Default is no restrictions.
        target_residents_only: When set to True the intervention is only distributed to individuals that began the
            simulation in that node.
        check_eligibility_at_trigger: If triggered event is delayed, you have an
            option to check individual/node's eligibility at the initial trigger
            or when the event is actually distributed after delay. (for example, a person might've aged out of the
            intervention before the initial trigger and the intervention distribution)
        receiving_drugs_event_name: Event to send out when person received drugs.
            Event name needs to include 'Received_Vehicle' in it, as otherwise overwritten in drug_campaigns function
            (see drug_campaigns.py L247)
            Default: SMC: 'Received_Vehicle'; PMC: 'Received_Vehicle_X' with X being number of PMC dose
        num_iiv_groups: Number of individual drug response groups.
            If >1, ind_property_restrictions is set to {'DrugResponseGroup': val} if campaign_type = PMC, not used for
            SMC. Default: SMC: Not used; PMC: 1, IIV only acts on the vaccine event, not the drug event
        receiving_drugs_event: Specify whether to deploy the parasite clearing drug event or the vaccine event only.
            Default: True
            Exception: for PMC, set to False

    Returns:
        dictionary of tags
    """
    if len(start_days) != len(coverages):
        raise ValueError(f"Length of start_days - {len(start_days)}, should be equal to length of coverages - "
                         f"{len(coverages)}, but it's not.\n")
    
    vaccine_param_dict = {'vacc_initial_effect': 0.598, 'vacc_box_duration': 21.7, 'vacc_decay_duration': 1.18}
    vaccine_initial_effect = vaccine_param_dict['vacc_initial_effect']
    vaccine_box_duration = vaccine_param_dict['vacc_box_duration']
    vaccine_decay_duration = vaccine_param_dict['vacc_decay_duration']
    target_age_min = 0
    target_age_max = 125
    if target_group:
        target_age_min = target_group['agemin']
        target_age_max = target_group['agemax']
    for (d, cov) in zip(start_days, coverages):
        add_drug_campaign(campaign, campaign_type='MDA',
                          drug_code='Vehicle',
                          start_days=[d + 1],  # when triggered and triggering interventions are deployed on the
                          # same day, that-day events are not "heard" by the triggered intervention, so I'm
                          # off-setting them by one day
                          coverage=cov,
                          repetitions=-1,
                          tsteps_btwn_repetitions=-1,
                          listening_duration=listening_duration,
                          target_group=target_group,
                          ind_property_restrictions=ind_property_restrictions,
                          receiving_drugs_event_name=receiving_drugs_event,
                          trigger_condition_list=trigger_condition_list,
                          target_residents_only=target_residents_only,
                          node_ids=node_ids,
                          check_eligibility_at_trigger=check_eligibility_at_trigger)
    add_triggered_vaccine(campaign,
                          start_day=start_days[0],  # otherwise it won't "hear" the first round of drugs
                          demographic_coverage=1,
                          trigger_condition_list=[receiving_drugs_event],
                          listening_duration=listening_duration,
                          target_age_min=target_age_min,
                          target_age_max=target_age_max,
                          node_ids=node_ids,
                          ind_property_restrictions=ind_property_restrictions,
                          intervention_name="RTSS",
                          vaccine_type="AcquisitionBlocking",
                          vaccine_initial_effect=vaccine_initial_effect,
                          vaccine_box_duration=vaccine_box_duration,
                          vaccine_decay_time_constant=vaccine_decay_duration / math.log(2),
                          efficacy_is_multiplicative=True,
                          broadcast_event=receiving_vaccine_event)
    return {'smc_cov': sum(coverages) / len(coverages),
            'total_smc_rounds': len(coverages)}
            
            
def add_smc(camp,smc_df):
    coord_df=load_coordinator_df(characteristic=False, set_index=True)
    sim_start_yr = int(coord_df.at['simulation_start_year','value'])
    for r, row in smc_df.iterrows():
         smc_year=int(row['year']) - sim_start_yr
         smc_month=int(row['month'])
         smc_day=int(row['day'])
         smc_start = smc_year*365 + day_of_year(smc_month,smc_day,smc_year)
         add_vaccdrug_smc(camp,start_days=[smc_start],
                          coverages=[row['coverage']],
                          target_group={'agemin': row['agemin'],
                                        'agemax': row['agemax']},
                          receiving_vaccine_event="Received_SMC_Vaccine", receiving_drugs_event="Received_Vehicle")                 


##### ITNs with seasonal and age-dependent usage #####
def add_itns(camp,itn_df,itn_age,itn_season):
    coord_df=load_coordinator_df(characteristic=False, set_index=True)
    sim_start_yr = int(coord_df.at['simulation_start_year','value'])
    itn_seasonal_usage = {"Times": list(itn_season['season_time']),
                          "Values":list(itn_season['season_usage'])}
    for r, row in itn_df.iterrows():
        itn_year = int(row['year'])
        itn_month = int(row['month'])
        itn_day = int(row['day'])
        itn_start = itn_year-sim_start_yr
        itn_start = itn_start*365 + day_of_year(itn_month,itn_day,itn_year)
        itn_discard_config = {"Expiration_Period_Distribution": "WEIBULL_DISTRIBUTION",                        
                              "Expiration_Period_Kappa": float(row['discard_k']),                        
                              "Expiration_Period_Lambda": float(row['discard_l'])}                       
        itn_age_year = itn_age[itn_age['year']==itn_year]                                                          
        itn_age_bins = itn_age_year['age']                                                                
        itn_age_usage =itn_age_year['age_usage']                                                            
        add_scheduled_usage_dependent_bednet(camp, intervention_name = "UsageDependentBednet",                   
                                             start_day = itn_start,                      
                                             demographic_coverage = float(row['coverage']),          
                                                   killing_initial_effect = float(row['kill_effect']),    
                                                   killing_decay_time_constant = int(row['kill_decay']),   
                                                   blocking_initial_effect = float(row['block_effect']),   
                                                   blocking_decay_time_constant=int(row['block_decay']),   
                                             age_dependence = {"Times": list(itn_age_bins),                
                                                               "Values": list(itn_age_usage)},             
                                             seasonal_dependence = itn_seasonal_usage,                     
                                             discard_config = itn_discard_config)
                                             
                                             
####################################        
##### General helper functions #####
####################################

def day_of_year(month, day, year):
    """
    Converts a given month and day to the day of the year.

    Args:
    month: The month (1-12).
    day: The day of the month (1-31).
    year: The year.

    Returns:
    The day of the year (1-366).
    """
    return datetime(year, month, day).timetuple().tm_yday


def get_comps_id_filename(site: str, level: int = 0):
    folder_name = manifest.comps_id_folder
    if level == 0:
        file_name = folder_name / (site + '_exp_submit')
    elif level == 1:
        file_name = folder_name / (site + '_exp_done')
    elif level == 2:
        file_name = folder_name / (site + '_analyzers')
    else:
        file_name = folder_name / (site + '_download')
    return file_name.relative_to(manifest.CURRENT_DIR).as_posix()


def load_coordinator_df(characteristic=False, set_index=True):
    csv_file = manifest.sweep_sim_coordinator_path if characteristic else manifest.simulation_coordinator_path
    coord_df = pd.read_csv(csv_file)
    if set_index:
        coord_df = coord_df.set_index('option')
    return coord_df


def get_suite_id():
    if os.path.exists(manifest.suite_id_file):
        with open(manifest.suite_id_file, 'r') as id_file:
            suite_id = id_file.readline()
        return suite_id
    else:
        return 0
      
      
def generate_demographics():
    coord_df=load_coordinator_df(characteristic=False, set_index=True)
    site = coord_df.at['site','value']
    latitude=coord_df.at['lat','value']
    longitude=coord_df.at['lon','value']
    population=coord_df.at['pop','value']
    prev0 = coord_df.at['prev0','value']
    BR = coord_df.at['birth_rate','value']
    new_nodes = [Demog.Node(lat=float(latitude),
                            lon=float(longitude),
                            pop=int(population),
                            name=str(site),
                            forced_id=1)]
    demog = Demographics.MalariaDemographics(nodes=new_nodes,
                                             idref=str(site),
                                             init_prev=float(prev0),
                                             include_biting_heterogeneity=True)
    print("Setting Equilibrium Vital Dynamics")
    demog.SetEquilibriumVitalDynamics(CrudeRate(float(BR)))
    print("Getting Equilibrium Age Distribution")
    demog.SetEquilibriumAgeDistFromBirthAndMortRates(CrudeRate(float(BR)),
                                                     CrudeRate(float(BR)))
    print("Amending Birth Rate")
    demog.SetBirthRate(CrudeRate(float(BR) * int(population)))
    with open(f"{manifest.input_files_path}/demographics_files/{site}_demographics.json", "w") as outfile:
        json.dump(demog.to_dict(), outfile, indent=3, sort_keys=True)
        print(f"Saved to {manifest.input_files_path}/demographics_files/{site}_demographics.json")
    return demog
  

def extract_climate(flatten_temp=True):
    import time
    from emodpy_malaria.weather import (generate_weather, weather_to_csv, WeatherVariable, 
                                        csv_to_weather)
    coord_df=load_coordinator_df(characteristic=False, set_index=True)
    coord_df=coord_df[['value']]
    #print(coord_df)
    # ---| Request weather files |---
    site = coord_df.at['site','value']
    start_yr = int(coord_df.at['climate_start_year','value'])
    length = int(coord_df.at['climate_year_dur','value'])
    extractdir = f'{manifest.input_files_path}/tmp/'
    outdir = os.path.join(f'{manifest.input_files_path}/site_climate', site)
    if not os.path.exists(extractdir):
        os.makedirs(extractdir)
    site_climate=coord_df.transpose().reset_index()
    site_climate = site_climate[['lon','lat','nodes']]
    #print(site_climate)
    site_climate.to_csv(f"{manifest.simulation_input_filepath}/{site}_climate.csv")
    weather_dir = extractdir
    startdate = start_yr * 1000 + 1
    enddate = (start_yr + length - 1) * 1000 + 365
    wr = generate_weather(platform="Calculon",
                          site_file=f"{manifest.simulation_input_filepath}/{site}_climate.csv",
                          start_date=startdate,
                          end_date=enddate,
                          node_column="nodes",
                          id_reference=site,
                          local_dir=weather_dir,
                          force=True)
    time.sleep(10)
    df, wa = weather_to_csv(weather_dir)
    weather_columns = {
        WeatherVariable.AIR_TEMPERATURE: "airtemp",
        WeatherVariable.RELATIVE_HUMIDITY: "humidity",
        WeatherVariable.RAINFALL: "rainfall"
    }
    weather_filenames = {
        WeatherVariable.AIR_TEMPERATURE: "air_temperature_daily.bin",
        WeatherVariable.RELATIVE_HUMIDITY: "relative_humidity_daily.bin",
        WeatherVariable.RAINFALL: "rainfall_daily.bin"
    }
    # Remove extra day in 2016
    df = df[df.steps != 1096]
    df.steps = [x if x < 1096 else x - 1 for x in df.steps]
    df1 = df.copy()
    print(df1['airtemp'])
    # Set constant air temperature to the mean
    print("flattening temp")
    if flatten_temp:
        airtempMean = df1["airtemp"].mean()
        df1.loc[:, "airtemp"] = airtempMean
    print(df1['airtemp'])
    csv_to_weather(df1, attributes=wa, weather_columns=weather_columns,
                   weather_dir=outdir,
                   weather_file_names=weather_filenames)


###########################################        
##### Calibration parameter "plug-in" #####
###########################################
    
def add_calib_param_func(simulation, calib_params, sets, hab_base = 1e8, const_base = 1e6):
    X = calib_params[calib_params['param_set'] == sets]
    X = X.reset_index(drop=True)
    # Temperature Shift: Ensure climate model is enabled before setting temperature offsets
    if simulation.task.config.parameters.Climate_Model == "CLIMATE_BY_DATA":
        simulation.task.config.parameters.Air_Temperature_Offset = float(X['emod_value'][0])
        simulation.task.config.parameters.Land_Temperature_Offset = float(X['emod_value'][0])
    else:
        warnings.warn("Climate model is not enabled; skipping temperature offsets.")
    # Vectors
    coord_df = load_coordinator_df(characteristic=False, set_index=True)
    # Load vector file
    vdf = pd.read_csv(os.path.join(manifest.input_files_path,coord_df.at['vector_filepath','value']))
    # Get list of species
    s = [species for species in vdf['species']]
    for r in range(len(s)):
        # Update habitat available to species
        habitat1 = dfs.schema_to_config_subnode(
            manifest.schema_file, ["idmTypes", "idmType:VectorHabitat"]
        )
        # constantly available habitat
        habitat1.parameters.Habitat_Type = "CONSTANT"
        malaria_config.set_species_param(
            simulation.task.config, vdf['species'][r], "Habitats", habitat1.parameters, overwrite=True
        )
        #habitat1.parameters.Max_Larval_Capacity = int(const_base * (vdf['fraction'][r] * vdf['constant'][r]) * float(X['emod_value'][1]))
        # temporarily available habitat
        habitat2 = dfs.schema_to_config_subnode(
            manifest.schema_file, ["idmTypes", "idmType:VectorHabitat"]
        )
        habitat2.parameters.Habitat_Type = "TEMPORARY_RAINFALL"
        malaria_config.set_species_param(simulation.task.config, vdf['species'][r], "Habitats", habitat2.parameters, overwrite=False)
        # semipermanent habitat
        habitat3 = dfs.schema_to_config_subnode(
            manifest.schema_file, ["idmTypes", "idmType:VectorHabitat"]
        )
        habitat3.parameters.Habitat_Type = "WATER_VEGETATION"
        malaria_config.set_species_param(simulation.task.config, vdf['species'][r], "Habitats", habitat3.parameters, overwrite=False)
        malaria_config.set_max_larval_capacity(
            simulation.task.config, vdf['species'][r], "CONSTANT", const_base * (vdf['fraction'][r] * vdf['constant'][r]) * float(X['emod_value'][1])
        )
        malaria_config.set_max_larval_capacity(
            simulation.task.config, vdf['species'][r], "TEMPORARY_RAINFALL", hab_base * (vdf['fraction'][r] * vdf['temp_rain'][r]) * float(X['emod_value'][2])
        )
        malaria_config.set_max_larval_capacity(
            simulation.task.config, vdf['species'][r], "WATER_VEGETATION", hab_base * (vdf['fraction'][r] * vdf['water_veg'][r]) * float(X['emod_value'][3])
        )
    return {'Sample_ID':sets}


def set_climate(config, shift):
    """
    Set climate to specific files, currently hardcoded to apply 
    a change of <shift> degrees to each daily land & temperature value
    in the climate input files
    """
    config.parameters.Air_Temperature_Offset = shift
    config.parameters.Land_Temperature_Offset = shift
    return {"Temperature_Shift": shift}


def set_species_from_file(config, vector_file):
    """
    Set up mosquito species and behavioral parameters
    """
    # Load input vector file
    vdf = pd.read_csv(os.path.join(manifest.input_files_path,vector_file))   
    # Get list of species
    s = [species for species in vdf['species']]
    # Add species to simulation
    conf.add_species(config, manifest, s)
    # For each species...
    for r in range(len(s)):    
        # Set vector parameter - anthropophily
        conf.set_species_param(config, 
                               species = vdf['species'][r],
                               parameter='Anthropophily',
                               value=vdf['anthropophily'][r],
                               overwrite=True)
        # Set vector parameter - indoor feeding fraction
        conf.set_species_param(config, 
                               species = vdf['species'][r],
                               parameter='Indoor_Feeding_Fraction',
                               value=vdf['indoor_feeding'][r],
                               overwrite=True)                                                                                     
    return
