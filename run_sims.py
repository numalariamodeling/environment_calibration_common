##### Import required packages #####
# standard packages
import argparse
import numpy as np
import os
from pathlib import Path
from functools import \
    partial 
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from emodpy.emod_task import EMODTask
import sys
# from within environment_calibration_common submodule
from helpers import set_param_fn, update_sim_random_seed, build_demog, set_simulation_scenario_for_characteristic_site, \
    set_simulation_scenario_for_matched_site, get_comps_id_filename, add_outputs, add_calib_param_func,extract_climate,\
    generate_demographics
from utils_slurm import submit_scheduled_analyzer
# from source 'simulations' directory
sys.path.append("../simulations")
import manifest as manifest

def submit_sim(site=None, nSims=1, characteristic=False, priority=manifest.priority, my_manifest=manifest,
               not_use_singularity=False, X=None):
    """
    This function is designed to be a parameterized version of the sequence of things we do 
    every time we run an emod experiment. 
    """
    # Create a platform
    # Show how to dynamically set priority and node_group
    platform = Platform("SLURM_LOCAL", job_directory=manifest.job_directory, partition='b1139', time='6:00:00',
                            account='b1139', modules=['singularity'], max_running_jobs=250, mem=2500)
    
    experiment = create_exp(characteristic, nSims, site, my_manifest, not_use_singularity,platform, X)

    # The last step is to call run() on the ExperimentManager to run the simulations.
    experiment.run(wait_until_done=False, platform=platform)

    # Additional step to schedule analyzer to run after simulation finished running
    submit_scheduled_analyzer(experiment, platform, site, analyzer_script='../environment_calibration_common/run_analyzers.py', mem=60000) #, mem=10000)

    # Save experiment id to file
    comps_id_file = get_comps_id_filename(site=site)
    with open(comps_id_file, "w") as fd:
        fd.write(str(experiment.uid))
    print()
    print(str(experiment.uid))
    return str(experiment.uid)

def add_calib_params(task, param, value, ptype):
    if ptype == "integer": #and (param != "Max_Individual_Infections" or (param == "Max_Individual_Infections" and value > 1)):
        task.set_parameter(param, int(value))
    #elif ptype == "integer" and param == "Max_Individual_Infections" and value <= 1:
     #   task.set_parameter('Enable_Superinfection', int(0))
    elif ptype in ["float", "double"]:
        task.set_parameter(param, float(value))
    elif ptype in ['string']:
        task.set_parameter(param, str(value))
        
    return {param: value}

def create_exp(characteristic, nSims, site, my_manifest, not_use_singularity, platform, X):
    task = _create_task(my_manifest, site)

    if not not_use_singularity:
        task.set_sif(manifest.SIF_PATH, platform)
        #task.set_sif(my_manifest.sif_id.as_posix())
    builder, exp_name = _create_builder(task,characteristic, nSims, site, X)
    # create experiment from builder
    print("created builder")
    experiment = Experiment.from_builder(builder, task, name=exp_name)
    print("created experiment")
    return experiment


def _create_builder(task,characteristic, nSims, site, X):
    # Create simulation sweep with builder
    builder = SimulationBuilder()
    exp_name = "validation_" + site
    # Sweep run number
    builder.add_sweep_definition(update_sim_random_seed, range(nSims))
    print("sweep run_number")
    # Sweep sites and seeds - based on values in simulation_coordinator csv
    if characteristic:
        builder.add_sweep_definition(set_simulation_scenario_for_characteristic_site, [site])
    else:
        builder.add_sweep_definition(set_simulation_scenario_for_matched_site, [site])
    print("setting scenario")
    builder.add_sweep_definition(partial(add_calib_param_func, calib_params=X), np.unique(X['param_set']))
    print("sweep builder params")
    print("made builder")
    return builder, exp_name


def _create_task(my_manifest,site):
    # create EMODTask
    print("Creating EMODTask (from files)...")
    task = EMODTask.from_default2(config_path="my_config.json",
                                  eradication_path=str(my_manifest.eradication_path),
                                  ep4_custom_cb=None,
                                  campaign_builder=None,
                                  schema_path=str(my_manifest.schema_file),
                                  param_custom_cb=set_param_fn,
                                  demog_builder=build_demog,
                                  )

    task.config.parameters.Birth_Rate_Dependence = "FIXED_BIRTH_RATE"
    task.common_assets.add_directory(
       Path(os.path.join(manifest.input_files_path, 'site_climate', site)), relative_path="climate"
    )
    add_outputs(task,site)
    return task


if __name__ == "__main__":
    # TBD: user should be allowed to specify (override default) erad_path and input_path from command line 
    # plan = EradicationBambooBuilds.MALARIA_LINUX
    # print("Retrieving Eradication and schema.json from Bamboo...")
    # get_model_files( plan, manifest )
    # print("...done.")

    parser = argparse.ArgumentParser(description='Process site name')
    parser.add_argument('--site', '-s', type=str, help='site name',
                        default="test_site")
    parser.add_argument('--nSims', '-n', type=int, help='number of simulations', default=1)
    parser.add_argument('--characteristic', '-c', action='store_true', help='site-characteristic sweeps')
    parser.add_argument('--not_use_singularity', '-i', action='store_true',
                        help='not using singularity image to run in Comps')
    parser.add_argument('--priority', '-p', type=str,
                        choices=['Lowest', 'BelowNormal', 'Normal', 'AboveNormal', 'Highest'],
                        help='Comps priority', default=manifest.priority)
    parser.add_argument('--calib_params', '-X', type=str, help='calib parameter set')
    args = parser.parse_args()
    generate_demographics(site=args.site)
    submit_sim(site=args.site, nSims=args.nSims, characteristic=args.characteristic, priority=args.priority,
               not_use_singularity=args.not_use_singularity, X=args.calib_params)
