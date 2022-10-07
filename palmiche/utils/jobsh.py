#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Fri Jun 25 21:00:10 2021
Author        : Alejandro Martínez León
Mail          : [alejandro.martinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
import socket, glob
sbatch_translation = {
    #https://slurm.schedmd.com/sbatch.html
    'A':'account',
    'account':'account',
    'acctg-freq':'acctg-freq',
    'a':'array','array':'array',
    'batch':'batch',
    'bb':'bb',
    'bbf':'bbf',
    'b':'begin', 'begin': 'begin',
    'D':'chdir', 'chdir':'chdir',
    'cluster-constraint':'cluster-constraint',
    'M':'clusters', 'clusters':'clusters',
    'comment':'comment',
    'C':'constraint', 'constraint':'constraint',
    'container':'container',
    'contiguous':'contiguous',
    'S':'core-spec', 'core-spec':'core-spec',
    'cores-per-socket':'cores-per-socket',
    'cpu-freq':'cpu-freq',
    'cpus-per-gpu':'cpus-per-gpu',
    'c':'cpus-per-task', 'cpus-per-task':'cpus-per-task',
    'deadline':'deadline',
    'delay-boot':'delay-boot',
    'd':'dependency', 'dependency':'dependency',
    'm':'distribution', 'distribution':'distribution',
    'e':'error', 'error':'error',
    'x':'exclude', 'exclude':'exclude',
    'exclusive':'exclusive',
    'export':'export',
    'export-file':'export-file',
    'B':'extra-node-info', 'extra-node-info':'extra-node-info',
    'get-user-env':'get-user-env',
    'gid':'gid',
    'gpu-bind':'gpu-bind',
    'gpu-freq':'gpu-freq',
    'G':'gpus', 'gpus':'gpus',
    'gpus-per-node':'gpus-per-node',
    'gpus-per-socket':'gpus-per-socket',
    'gpus-per-task':'gpus-per-task',
    'gres':'gres',
    'gres-flags':'gres-flags',
    'h':'help', 'help':'help',
    'hint':'hint',
    'H':'hold', 'hold':'hold',
    'ignore-pbs':'ignore-pbs',
    'i':'input', 'input':'input',
    'J':'job-name', 'job-name':'job-name',
    'kill-on-invalid-dep':'kill-on-invalid-dep',
    'L':'licenses', 'licenses':'licenses',
    'mail-type':'mail-type',
    'mail-user':'mail-user',
    'mcs-label':'mcs-label',
    'mem':'mem',
    'mem-bind':'mem-bind',
    'mem-per-cpu':'mem-per-cpu',
    'mem-per-gpu':'mem-per-gpu',
    'mincpus':'mincpus',
    'network':'network',
    'nice':'nice',
    'k':'no-kill', 'no-kill':'no-kill',
    'no-requeue':'no-requeue',
    'F':'nodefile', 'nodefile':'nodefile',
    'w':'nodelist','nodelist':'nodelist',
    'N':'nodes', 'nodes':'nodes',
    'n':'ntasks', 'ntasks':'ntasks',
    'ntasks-per-core':'ntasks-per-core',
    'ntasks-per-gpu':'ntasks-per-gpu',
    'ntasks-per-node':'ntasks-per-node',  
    'ntasks-per-socket':'ntasks-per-socket',    
    'open-mode':'open-mode',
    'o':'output', 'output':'output',
    'O':'overcommit', 'overcommit':'overcommit',
    's':'oversubscribe', 'oversubscribe':'oversubscribe',
    'parsable':'parsable',
    'p':'partition', 'partition':'partition',
    'power':'power',    
    'priority':'priority',    
    'profile':'profile',
    'propagate':'propagate',
    'q':'qos', 'qos':'qos',
    'Q':'quiet','quiet':'quiet',
    'reboot':'reboot',
    'requeue':'requeue',
    'reservation':'reservation',
    'signal':'signal',  
    'sockets-per-node':'sockets-per-node',
    'spread-job':'spread-job',
    'switches':'switches',
    'test-only':'test-only',
    'thread-spec':'thread-spec',
    'threads-per-core':'threads-per-core',
    't':'time', 'time':'time', #Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds".
    'tmp':'tmp',
    'uid':'uid',
    'usage':'usage',
    'use-min-nodes':'use-min-nodes',
    'v':'verbose', 'verbose':'verbose',
    'V':'version', 'version':'version',
    'W':'wait', 'wait':'wait',
    'wait-all-nodes':'wait-all-nodes',
    'wckey':'wckey',
    'wrap':'wrap',
}
class JOB:
    def __init__(self, sbatch_keywords = None, GROMACS_version = 2021.4, mdrun_keywords = None, build_GROMACS_section = None) -> None:
        self.hostname = socket.gethostname()
        self.GROMACS_version = GROMACS_version
        self.build_GROMACS_section = build_GROMACS_section
        self.default_times = {
            'deflt': '2-00:00:00',
            'long': '5-00:00:00',
            'debug': '30:00',
            'medium': '2-00:00:00',
            'fat+': '2-00:00:00',
            'fat': '2-00:00:00',
            'int': '2-00:00:00',
            'gpu': '2-00:00:00',
            'gpu:test': '30:00',
            'gpu-hub':'2-00:00:00',
        }
        self.sbatch_keywords = {
            'job-name': 'run',
            'output': 'myjob.out',
            'error': 'myjob.err',
            'nice': 0,
            'nodes': 1,
            'gpus':1,
        }
        # This part need to be improved
        if self.hostname == 'smaug' or self.hostname.startswith('fang'):
            self.sbatch_keywords['cpus-per-task'] = 12
            self.sbatch_keywords['partition'] = 'deflt'
        elif self.hostname.startswith('gwd'): # This is quiet lazy, but is because I don't know the names of the gwdg nodes
            self.sbatch_keywords['cpus-per-task'] = 10
            self.sbatch_keywords['partition'] = 'gpu-hub'#'medium'#'gpu-hub'
            self.sbatch_keywords['exclude'] = 'gwdo[161-180]'
            # self.sbatch_keywords['account'] = 'all'
            # self.sbatch_keywords['constraint'] = 'scratch'
        else:
            self.sbatch_keywords['cpus-per-task'] = 12
            self.sbatch_keywords['partition'] = 'deflt'
         
        # Translating the user sbatch_keywords in case that something was provided
        if sbatch_keywords:       
            sbatch_keywords_translated = {}
            for key in sbatch_keywords:
                try:
                    sbatch_keywords_translated[sbatch_translation[key]] = sbatch_keywords[key]
                except:
                    raise Exception(f'The keyword {key} of sbatch_keywords is not a valid sbatch keyword argument')       
            
            # Update with the user sbatch keword provided
            self.sbatch_keywords.update(sbatch_keywords_translated)
            
            # Setting the correct time
            if 'time' not in sbatch_keywords_translated:
                self.sbatch_keywords['time'] = self.default_times[self.sbatch_keywords['partition']]

        # Deffining default options. In case that flags was provided, the value must be a bool type
        self.mdrun_keywords = {
            'nt': self.sbatch_keywords['cpus-per-task'], # This is not so straightforwart. Because ntasks and cpus-per-task, but is a default option, you could just provided a new value
            'cpi': True,
            'stepout': 5000,
            'v':True,
        }
        # Updating the values of self.mdrun_keywords with the provided by the user
        if mdrun_keywords:        
            self.mdrun_keywords.update(mdrun_keywords)
        
    # Now we will start to construct the file
    def string(self):
        job_str = '#!/bin/bash\n'
        if self.hostname == 'smaug' or self.hostname.startswith('fang') or self.hostname.startswith('gwd'):
            # SBATCH section
            for key in self.sbatch_keywords:
                job_str += f"#SBATCH --{key}={self.sbatch_keywords[key]}\n"
            # Echo some important variables
            job_str += (
                "\n\n# This block is echoing some SLURM variables\n"\
                "echo \"Job execution start: $(date)\"\n"\
                "echo \"JobID = $SLURM_JOBID\"\n"\
                "echo \"Host = $SLURM_JOB_NODELIST\"\n"\
                "echo \"Jobname = $SLURM_JOB_NAME\"\n"\
                "echo \"Subcwd = $SLURM_SUBMIT_DIR\"\n"\
                "echo \"SLURM_TASKS_PER_NODE = $SLURM_TASKS_PER_NODE\"\n"\
                "echo \"SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK\"\n"\
                "echo \"SLURM_CPUS_ON_NODE = $SLURM_CPUS_ON_NODE\"\n\n"\
                )
            # Loading modules
            if self.hostname == 'smaug' or self.hostname.startswith('fang'):
                job_str += (
                        f"source /data/shared/opt/gromacs/{self.GROMACS_version}/bin/GMXRC\n"\
                        f"PATH=$PATH:\"/data/shared/opt/gromacs/{self.GROMACS_version}/bin\"\n"\
                        f"#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:\n\n"
                        )
            # This could give some problems for testing GROMACS versions
            elif self.hostname.startswith('gwd'):
                job_str += 'source /usr/users/cmb/shared/spack-0.17/shared.bash\n'\
                'if [[ $HOSTNAME == gwdo* ]]; then\n'\
                f'  module load  ivybridge/gromacs@{self.GROMACS_version}\n'\
                'else\n'\
                f'  module load gromacs@{self.GROMACS_version}\n'\
                'fi\n\n'\
        
        job_str += 'cd $(pwd)\n\n'
        # Constructing GROMACS command
        if self.build_GROMACS_section:
            job_str += self.build_GROMACS_section
        else:
            if 's' not in self.mdrun_keywords and 'deffnm' not in self.mdrun_keywords:
                tpr_name = glob.glob('*.tpr')[0].split('.')[0]
                CMD_GROMACS = f'gmx mdrun -v -deffnm {tpr_name}'
            else:
                CMD_GROMACS = 'gmx mdrun'
                try:
                    tpr_name = self.mdrun_keywords['s'][0].split('.')[0]
                except:
                    try:
                        tpr_name = self.mdrun_keywords['deffnm'][0].split('.')[0]
                    except:
                        tpr_name = 'log'

            END_CMD_GROMACS = f' >> {tpr_name}.lis 2>&1'

            for key in self.mdrun_keywords:
                if type(self.mdrun_keywords[key]) == bool:
                    CMD_GROMACS += f' -{key}'
                else:
                    CMD_GROMACS += f' -{key} {self.mdrun_keywords[key]}'
            job_str += CMD_GROMACS + END_CMD_GROMACS

        return job_str

    def write(self, name = 'job.sh'):
        with open(name, "w") as out:
            out.write(self.string())

if __name__ == '__main__':
    
    pass
    G = JOB(mdrun_keywords={'deffnm': 'mine'}, build_GROMACS_section='gmx mdrun -nt 12 -cpi -stepout 5000 -v -deffnm taka >& tak.lis')
    print(G.string())