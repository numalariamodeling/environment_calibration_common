##### Import required packages #####
# standard packages
import argparse
import os
import sys
from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.core import ItemType
import pandas as pd
# from within analyzers/
sys.path.append(os.path.dirname(__file__))
from .analyzer_collection import (
    EventReporterAnalyzer,
    MonthlyPfPRAnalyzer,
    MonthlyIncidenceAnalyzer,
    AnnualPfPRAnalyzer,
    AnnualIncidenceAnalyzer,
    InsetChartAnalyzer,
    EventReporterSummaryAnalyzer,
    NodeDemographicsAnalyzer,
    VectorStatsAnalyzer,
    EIRAnalyzer,
    PCRAnalyzer
)
# from within environment_calibration_common submodule
sys.path.append("../")
from helpers import load_coordinator_df

# from source 'simulations' directory
sys.path.append("../../simulations")
import manifest


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--site", dest="site", type=str, required=True)
    parser.add_argument("--expid", dest="expid", type=str, required=True)

    return parser.parse_args()


def analyze_experiment(platform, expid, wdir):
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    coord_df = load_coordinator_df()
    sim_start_year = int(coord_df.at['simulation_start_year','value'])
    sim_years = int(coord_df.at['simulation_years','value'])
    analyzers = []
    # custom analyzers
    sweep_variables = ['Run_Number', 'Sample_ID']
    # EIR analyzer
    analyzers.append(EIRAnalyzer(sweep_variables=sweep_variables,
                                                working_dir=wdir,
                                                start_day=min(0,(sim_years-10)*365),
                                                end_day = sim_years*365,
                                                channels=["Daily EIR"]))
    if coord_df.at['prevalence_comparison','value']:
        prevalence_df_U5 = pd.read_csv(os.path.join(manifest.base_reference_filepath,
                                                    coord_df.at['prevalence_comparison_reference','value']))

        prevalence_df_U2 = pd.read_csv(os.path.join(manifest.base_reference_filepath,
                                                     coord_df.at['prevalence_comparison_reference_U2','value']))

        if(coord_df.at['prevalence_comparison_diagnostic','value']=="PCR"):
            prevalence_start_u5 = int(prevalence_df_U5['year'].min()) - sim_start_year
            prevalence_end_u5 = int(prevalence_df_U5['year'].max()) - sim_start_year
            analyzers.append(PCRAnalyzer(sweep_variables=sweep_variables,
                                                working_dir=wdir,
                                                start_day=prevalence_start_u5*365,
                                                end_day = 365+prevalence_end_u5*365,
                                                channels=["PCR Parasite Prevalence"]))
        
        if coord_df.at['prevalence_comparison_diagnostic','value']=="Microscopy":
            prevalence_start_u5 = int(prevalence_df_U5['year'].min())
            prevalence_end_u5 = int(prevalence_df_U5['year'].max())
            if coord_df.at['prevalence_comparison_frequency','value']=="monthly":
                print(f"Prevalence: {prevalence_start_u5} to {prevalence_end_u5}")
                analyzers.append(MonthlyPfPRAnalyzer(sweep_variables=sweep_variables,
                                                     working_dir=wdir,
                                                     start_year=prevalence_start_u5,
                                                     end_year= prevalence_end_u5,
                                                     grp = 'U5'))
                                                     
            prevalence_start_u2 = int(prevalence_df_U2['year'].min())
            prevalence_end_u2 = int(prevalence_df_U2['year'].max()+1)
            if coord_df.at['prevalence_comparison_frequency','value']=="monthly":
                print(f"Prevalence: {prevalence_start_u2} to {prevalence_end_u2}")
                analyzers.append(MonthlyPfPRAnalyzer(sweep_variables=sweep_variables,
                                                     working_dir=wdir,
                                                     start_year=prevalence_start_u2,
                                                     end_year=prevalence_end_u2,
                                                     grp = 'U2'))

        if(coord_df.at['prevalence_comparison_diagnostic','value']=="PCR"):
            prevalence_start_u2 = int(prevalence_df_U2['year'].min()) - sim_start_year
            prevalence_end_u2 = int(prevalence_df_U2['year'].max()) - sim_start_year
            analyzers.append(PCRAnalyzer(sweep_variables=sweep_variables,
                                                working_dir=wdir,
                                                start_day=prevalence_start_u2*365,
                                                end_day = 365+prevalence_end_u2*365,
                                                channels=["PCR Parasite Prevalence"]))
        
        if coord_df.at['prevalence_comparison_diagnostic','value']=="Microscopy":

            if coord_df.at['prevalence_comparison_frequency','value']=="annual":
                analyzers.append(AnnualPfPRAnalyzer(sweep_variables=sweep_variables,
                                                    working_dir=wdir,
                                                    start_year=prevalence_start_u5,
                                                    end_year=prevalence_end_u5))

            if coord_df.at['prevalence_comparison_frequency','value']=="annual":
                analyzers.append(AnnualPfPRAnalyzer(sweep_variables=sweep_variables,
                                                    working_dir=wdir,
                                                    start_year=prevalence_start_u2,
                                                    end_year=prevalence_end_u2))
    # Don't change these - used for fitting #
    if coord_df.at['incidence_comparison','value']:
        incidence_df = pd.read_csv(os.path.join(manifest.base_reference_filepath,
                                                    coord_df.at['incidence_comparison_reference','value']))
        incidence_start = int(incidence_df['year'].min())
        incidence_end = int(incidence_df['year'].max())
        if(coord_df.at['incidence_comparison_frequency','value']=="monthly"):
                print("pass")
                analyzers.append(MonthlyIncidenceAnalyzer(sweep_variables=sweep_variables,
                                                         working_dir=wdir,
                                                         start_year=incidence_start,
                                                         end_year=incidence_end+1))
        if(coord_df.at['incidence_comparison_frequency','value']=="annual"):        
                analyzers.append(AnnualIncidenceAnalyzer(sweep_variables=sweep_variables,
                                                         working_dir=wdir,
                                                         start_year=incidence_start,
                                                         end_year=incidence_end+1))                                    
    # analyzers.append(EventReporterAnalyzer(sweep_variables=sweep_variables,
    #                                        working_dir=wdir,
    #                                        time_cutoff=int(report_start_day),
    #                                        event_list=["Received_Treatment"],
    #                                        output_filename="events"))
    # analyzers.append(VectorStatsAnalyzer(sweep_variables=sweep_variables,
    #                                      working_dir=wdir,
    #                                      start_time=int(report_start_day/365),
    #                                      end_time=int(simulation_duration/365),
    #                                      ))
    
    manager = AnalyzeManager(platform=platform,
                             configuration={},
                             ids=[(expid, ItemType.EXPERIMENT)],
                             analyzers=analyzers,
                             partial_analyze_ok=True,
                             max_workers=16)
    
    manager.analyze()


if __name__ == "__main__":
    
    from idmtools.core.platform_factory import Platform
    import manifest
    args = parse_args()
    platform = Platform('SLURM_LOCAL', job_directory=manifest.job_directory)
    outdir = args.site
    analyze_experiment(platform, 
                       args.expid,
                       os.path.join(manifest.output_dir, outdir))
