##### Import required packages #####
# standard packages
import argparse
import os
from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform
import pandas as pd

# from within analyzers/
from .analyzer_collection import (EventReporterAnalyzer, 
                                MonthlyPfPRAnalyzer,
                                InsetChartAnalyzer)
# from source 'simulations' directory
sys.path.append("../../")
from environment_calibration.helpers import load_coordinator_df
sys.path.append("/simulations")
import manifest

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', dest='expt_name', type=str, required=False)
    parser.add_argument('-id', dest='expt_id', type=str, required=True)

    return parser.parse_args()

if __name__ == "__main__":
    simyears = 30
    args = parse_args()
    wdir = os.path.join(manifest.output_dir, args.expt_name)

    if not os.path.exists(wdir):
        os.makedirs(wdir)

    sweep_variables = ['Run_Number', 'Habitat_Type', 'Habitat_Multiplier', 
                       'Constant_Multiplier', 'CM_Coverage']
    channels = ['Adult Vectors', 'Daily EIR', 'Daily Bites per Human', 'Rainfall']

    coord_df = load_coordinator_df(characteristic=False, set_index=True)
    simulation_years = int(coord_df.at['simulation_years','value'])
    sim_start_year = int(coord_df.at['simulation_start_year','value'])

    prevalence_df_U5 = pd.read_csv(os.path.join(manifest.base_reference_filepath,
                                                     coord_df.at['prevalence_comparison_reference','value']))

    prevalence_df_U2 = pd.read_csv(os.path.join(manifest.base_reference_filepath,
                                                     coord_df.at['prevalence_comparison_reference_U2','value']))

    prevalence_agebins = sorted([float(a) for a in prevalence_df_U5['age'].unique()])
    prevalence_agebins = sorted([float(a) for a in prevalence_df_U2['age'].unique()])

    first_year_U5 = int(prevalence_df_U5['year'].min()) - sim_start_year
    last_year_U5 = int(prevalence_df_U5['year'].max()) - sim_start_year

    first_year_U2 = int(prevalence_df_U2['year'].min()) - sim_start_year
    last_year_U2 = int(prevalence_df_U2['year'].max()) - sim_start_year

    # Logical bounds on years to request
    first_year = max(first_year, 0)
    last_year = min(last_year, simulation_years)

    first_year_U5 = max(first_year_U5, 0)
    last_year_U5 = min(last_year_U5, simulation_years)

    first_year_U2 = max(first_year_U2, 0)
    last_year_U2 = min(last_year_U2, simulation_years)



    with Platform('SLURM_LOCAL', job_directory=manifest.job_directory) as platform:
        analyzers = []
        analyzers.append(MonthlyPfPRAnalyzer(sweep_variables=sweep_variables,
                                            start_year=2015,
                                            end_year=2019,
                                            working_dir=wdir))
                                            
        analyzers.append(MonthlyPfPRAnalyzer(sweep_variables=sweep_variables,
                                            start_year=2015,
                                            end_year=2019,
                                            grp="U5",  # Explicitly setting U5
                                            working_dir=wdir))

        analyzers.append(MonthlyPfPRAnalyzer(sweep_variables=sweep_variables,
                                            start_year=2015,
                                            end_year=2019,
                                            grp="U2",  # Explicitly setting U2
                                            working_dir=wdir))
                                            
        analyzers.append(EventReporterAnalyzer(sweep_variables=sweep_variables,
                                               working_dir=wdir))
        analyzers.append(InsetChartAnalyzer(channels=channels,
                                            sweep_variables=sweep_variables,
                                            working_dir=wdir,
                                            start_day=(simyears-4)*365,
                                            end_day=99999))

        manager = AnalyzeManager(configuration={},
                                 ids=[(args.expt_id, ItemType.EXPERIMENT)],
                                 analyzers=analyzers, partial_analyze_ok=True)
        
        manager.analyze()
        print(f"\nSaving outputs to: {wdir}")
