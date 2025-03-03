##### Import required packages #####
# standard packages
import copy
import os
import sys
from typing import List, Dict
import subprocess
import re
import pandas as pd
from COMPS import AuthManager
from COMPS.Data import QueryCriteria, Simulation
from idmtools.core import ItemType
from idmtools.entities.iplatform import IPlatform
from idmtools_platform_slurm.slurm_operations.operations_interface import SlurmOperations
# from source 'simulations' directory
sys.path.append("../simulations")
import manifest


def _get_serialized_filenames(num_cores, timesteps):
    if num_cores == 1:
        serialized_list = [f"state-{str(timesteps).zfill(5)}.dtk"]
    else:
        serialized_list = [f"state-{str(timesteps).zfill(5)}-{str(core_count).zfill(3)}.dtk"
                           for core_count in range(num_cores)]
    return serialized_list


def _get_core_counts(sim_id, platform):
    # TODO, get num_cores from simulation in slurm
    # sim = platform.get_item(sim_id, ItemType.SIMULATION, raw=True)
    # sim.refresh(QueryCriteria().select_children('hpc_jobs'))
    # num_cores = int(sim.hpc_jobs[-1].configuration.max_cores)
    num_cores = 1
    return num_cores

def get_workdir_from_simulations(platform: 'IPlatform', comps_simulations: List[Simulation]) -> Dict[str, str]:
    """
    Get COMPS filepath
    Args:
        platform: idmtools Platform
        comps_simulations: COMPS Simulations
    Returns: dictionary with simid as key and filepath as value
    """

    if platform.environment.upper() == "SLURMSTAGE" or platform.environment.upper() == "CALCULON":
        mounts = AuthManager.get_environment_macros(platform.environment)['DOCKER_MOUNTS']
        mounts = {v[0]: v[1:4] for v in [m.split(';') for m in mounts.split('|')]}
        # pretend I'm on Linux and set the Linux mapping environment variables
        for k, v in mounts.items():
            os.environ[k] = ';'.join([v[0], v[2]])
    sim_work_dir = {str(sim.id): sim.hpc_jobs[-1].working_directory for sim in comps_simulations if sim.hpc_jobs}

    return sim_work_dir


def create_sim_directory_map(exp_id: str, platform: 'IPlatform'):
    """
        Return a dataframe which contains simulation's working_path and tags.
    Args:
        exp_id: experiment id
        platform: idmtools platform
    Returns:
        dataframe
    """
    # Get idmtools Experiment
    exp = platform.get_item(exp_id, ItemType.EXPERIMENT, raw=False)

    tags_list = []
    for sim in exp.simulations:
        tags = copy.deepcopy(sim.tags)
        tags.update(dict(simid=sim.id))
        tags_list.append(tags)
    df = pd.DataFrame(tags_list)

    if len(df) != len(exp.simulations):
        print(f'Warning: not all jobs in {exp_id} succeeded', ':', )
        print(len(exp.simulations) - len(df), 'unsucceeded')

    simulations = exp.simulations
    dir_list = []
    for sim in simulations:
        dir_dict = {"simid": str(sim.id), "serialized_file_path": platform._op_client.get_directory(sim)}
        dir_list.append(dir_dict)

    df_dir = pd.DataFrame(dir_list)
    df = pd.merge(left=df, right=df_dir, on='simid')

    return df


def build_burnin_df(exp_id: str, platform,serialize_days):
    """
    return dataframe which contains serialized_file_path, serialized_population_filenames
    Args:
        exp_id:
        platform:
    Returns:
        dataframe:
        Run_Number | Serialization_Time_Steps | task_type | sweep_tag | simid | serialized_file_path|Num_Cores|Serialized_Population_Filenames
    Note, Serialized_Population_Filenames depends on n_cores. if n_cores = 2, Serialized_Population_Filenames look
    like these: state-00050-000.dtk, state-00050-001.dtk
    """

    df = create_sim_directory_map(exp_id, platform)
    # add Num_Cores to df
    df["Num_Cores"] = df["simid"].apply(_get_core_counts, platform=platform)
    #print(list(df.columns))
    #print(df.head())

    #try:
    burnin_length_in_days = serialize_days#int(df["Serialization_Time_Steps"].iloc[0].strip('[]'))
    #except AttributeError:
        # different versions of pandas save this as either a string or a list
        #burnin_length_in_days = df["Serialization_Time_Steps"].iloc[0][-1]

    df["Serialized_Population_Filenames"] = df["Num_Cores"].apply(_get_serialized_filenames, timesteps=burnin_length_in_days)
    df = df.reset_index(drop=True)
    return df
    
### TODO idmtools probably has a predefined shell template   
def shell_header_quest(A='b1139', p='b1139', t='02:00:00', N=1, 
                       ntasks_per_node=1, mem=8000, job_name='myjob',
                       arrayJob=None, c=1, node_spec=None):
    """Requires a 'logs' subfolder to write in .err and .out files, alternatively logs/ needs to be removed"""
    if not os.path.exists('../simulations/log'):
        os.makedirs(os.path.join('../simulations/log'))
   
    header = f'#!/bin/bash\n' \
             f'#SBATCH -A {A}\n' \
             f'#SBATCH -p {p}\n' \
             f'#SBATCH -t {t}\n' \
             f'#SBATCH -N {N}\n' \
             f'#SBATCH --ntasks-per-node={ntasks_per_node}\n' \
             f'#SBATCH --mem={mem}\n' \
             f'#SBATCH --job-name="{job_name}"\n'\
             f'#SBATCH -c {c}\n'
    if node_spec is not None:
        constraint = f'#SBATCH --constraint={node_spec}\n'
    if arrayJob is not None:
        array = arrayJob
        err = '#SBATCH --error=log/slurm_%A_%a.err\n'
        out = '#SBATCH --output=log/slurm_%A_%a.out\n'
        header = header + array + err + out #+ constraint
    else:
        err = f'#SBATCH --error=log/{job_name}.%j.err\n'
        out = f'#SBATCH --output=log/{job_name}.%j.out\n'
        header = header + err + out #+ constraint
    return header
        
def submit_scheduled_analyzer(experiment, platform, site, analyzer_script, mem=20000):
    #wdir = os.path.abspath(os.path.dirname(__file__))
    wdir =  os.path.abspath(os.path.join(__file__ ,"../../simulations"))
    ## Write bash file to submit
    header_post = shell_header_quest(job_name=f'analyze_exp', t='06:00:00', mem=mem, c='8')
    pymodule = '\n\nmodule purge all' \
    ### pycommand(s) - additional python or R scripts to directly run after analyzer can be added below
    pycommand = f'\n{manifest.VENV_PATH}/bin/python {analyzer_script} --site {site} --expid {experiment.id}' 
    
    # if not os.path.exists('analyzers/batch'):
    #     os.makedirs(os.path.join('analyzers/batch'))    
    
    file = open(os.path.join(wdir,f'run_analyzer_{experiment.uid}.sh'), 'w') #'analyzers','batch',
    file.write(header_post + pymodule + pycommand)
    file.close()
    
    header_post_wait = shell_header_quest(job_name=f'wait_{site}', t='00:05:00', mem=50)
    batchscript = f'run_analyzer_{experiment.uid}.sh'
    batchcommand = f'\ncd {wdir}' \
                   f'\nsbatch {batchscript}'
                   ##os.path.join(wdir,"analyzers","batch"
    file = open(os.path.join(wdir,f'wait_analyzer_{experiment.uid}.sh'), 'w') #,'analyzers','batch'
    file.write(header_post_wait + pymodule + batchcommand)
    file.close()
    
    job_id = platform._op_client.get_job_id(experiment.id, experiment.item_type)
    job_id = job_id[0]
    script_path = os.path.join(wdir,f'wait_analyzer_{experiment.uid}.sh') #,'analyzers','batch' # save under different names, will require cleanup
    print(script_path)
    result = subprocess.run([f'sbatch --dependency=afterok:{job_id} {script_path}'], shell=True, stdout=subprocess.PIPE)
    print(result,flush=True)
    result = result.stdout.decode('utf-8').strip()
    print(result,flush=True)
    SBATCH_REGEX = re.compile('^[a-zA-Z ]+(?P<id>\d+)$')
    # print(SBATCH_REGEX.match(result),flush=True)
    job_id_analyzer = SBATCH_REGEX.match(result).group('id')
    #Regular expression to match the job ID
    regex = r'Submitted batch job (\d+)'
    # Using re.search to find the match
    match = re.search(regex, result)
    if match:
        job_id = match.group(1)
        print(job_id)  
    else:
        print("Problem with scheduled analyzer job id")
        exit(1)

    print(f"Experiment job id: {job_id}")
    print(f"Analyzer job id: {job_id_analyzer}")
    
    #return (job_id_analyzer)  # if needed to add another dependency slurm submission
    return ()