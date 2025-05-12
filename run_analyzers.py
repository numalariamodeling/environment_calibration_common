##### Import required packages #####
# standard packages
import sys
import os
import argparse
from pathlib import Path
from idmtools.core.platform_factory import Platform
# from within environment_calibration_common submodule
sys.path.append("/home/upf3610/b1139/ipti_pmc/environment_calibration/environment_calibration_common")
from helpers import get_comps_id_filename, load_coordinator_df

# from .helpers import get_comps_id_filename, load_coordinator_df
# from helpers import get_comps_id_filename, load_coordinator_df

sys.path.append("/home/upf3610/b1139/ipti_pmc/environment_calibration/environment_calibration_common")
from analyzers.analyze import analyze_experiment

# from analyzers.analyze import analyze_experiment
# from source 'simulations' directory

# sys.path.append("../simulations")
sys.path.append("/home/upf3610/b1139/ipti_pmc/environment_calibration/simulations")
import manifest as manifest

def run_analyzers(site: str, expid: str = None, characteristic: bool = False) -> (bool, str):
    """
    Wait for experiment to be done and run relevant analyzers for site on Comps with SSMT
    Args:
        site ():
        characteristic ():

    Returns: If experiment is succeeded, returns analyzer work item status and id,
             if not, return experiment status and id.

    """
    #platform = Platform(manifest.platform_name)
    platform = Platform('SLURM_LOCAL',job_directory=manifest.job_directory, mem=80000)
    comps_id_file = get_comps_id_filename(site=site)
    if expid:
        exp_id = expid
    else:
        with open(comps_id_file, 'r') as id_file:
            exp_id = id_file.readline()
    # Wait for experiment to be done
    a_ok=True
    if a_ok:#check_experiment(site, platform): #checks if succeeded, removed for now to allow for partial analysis
        coord_df = load_coordinator_df(characteristic = characteristic, set_index = True)
        # for expt_name, id in exp_name_id.items():
        # site = expt_name.replace('validation_', '')
        # determine the analyzers to run for each site
        wdir=manifest.simulation_output_filepath
        wdir=Path(os.path.join(wdir,site))

        if not os.path.exists(wdir):
            os.makedirs(wdir)

        analyzers_id_file = get_comps_id_filename(site=site, level=2)
        print(analyzers_id_file)
        

        analyze_experiment(platform, exp_id, wdir)

        with open(analyzers_id_file, 'w') as id_file:
            id_file.write(site)
        with open(os.path.join(manifest.simulation_output_filepath,site,'finished.txt')  , 'w') as f: 
            f.write("I'm done running :]") 
        #return wi.succeeded, wi.uid
        return
    else:
        return False, exp_id


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process site name')
    parser.add_argument('--site', '-s', type=str, help='site name',
                        default="")  # not sure if we want to make this required argument
    parser.add_argument(
        "-i",
        "--expid",
        type=str,
        default=None
    )
    args = parser.parse_args()
    sites=[]
    coord_df=load_coordinator_df()
    my_site=coord_df.at['site','value']
    sites=[my_site]
    for site in sites:
        run_analyzers(site,args.expid)

