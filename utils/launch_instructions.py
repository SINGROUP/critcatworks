# This script instructs the remote worker to launch jobs
import fireworks
print("Hello world")
import logging

from fireworks import Firework, FWorker, LaunchPad, Workflow
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fireworks.queue.queue_launcher import launch_rocket_to_queue
from fireworks.user_objects.queue_adapters.common_adapter import *
from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

import os, time, errno
import pathlib

def _recover_offline(lp, fworker_name = None, ignore_errors = False, print_errors = False):
    """
    Slightly modified command line function (recover_offline) from fireworks. v 1.8.4
    """
    lp
    failed_fws = []
    recovered_fws = []

    for l in lp.offline_runs.find({"completed": False, "deprecated": False},
                                  {"launch_id": 1, "fw_id":1}):
        if fworker_name and lp.launches.count({"launch_id": l["launch_id"],
                                               "fworker.name": fworker_name}) == 0:
            continue
        fw = lp.recover_offline(l['launch_id'], ignore_errors, print_errors)
        if fw:
            failed_fws.append(l['fw_id'])
        else:
            recovered_fws.append(l['fw_id'])

    lp.m_logger.info("FINISHED recovering offline runs. {} job(s) recovered: {}".format(
        len(recovered_fws), recovered_fws))
    if failed_fws:
        lp.m_logger.info("FAILED to recover offline fw_ids: {}".format(failed_fws))
    return

def _alt_recover_offline(lp):
    # recover offline fireworks
    ids =lp.get_fw_ids()
    for idx in ids:
        lp.recover_offline(launch_id = idx)
    return



if __name__ == "__main__":
    IS_QUEUE = True
    IS_OFFLINE = True
    # path used to log info

    logpath = ".fireworks/fw_logs"
    try:
        os.makedirs(logpath)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    abspath = str(pathlib.Path(logpath).resolve())

    if IS_QUEUE:
        logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.INFO)
    else:
        logging.basicConfig(filename = abspath + "/logfile_ranked_adsites.log", level=logging.INFO)

    # set up the LaunchPad 
    #launchpad = LaunchPad(logdir=".", strm_lvl='INFO')
    launchpad = LaunchPad(host = "austerity-shard-00-00-hgeov.mongodb.net:27017",
        port = 27017,
        name = "fireworks",
        username = "mjcritcat",
        password = "heterogeniuscatalysis",
        logdir =  ".",
        strm_lvl = "INFO",
        ssl =  True,
        authsource = "admin")

    # run the jobs found in the launchpad
    # (added previously by other routine)
    if IS_OFFLINE:
        rlaunch_command = "rlaunch  singleshot --offline"
    else:
        rlaunch_command = "rlaunch  singleshot"

    if IS_QUEUE:

        dft = CommonAdapter(
            q_type="SLURM",
            queue="test",
            nodes= None,
            ntasks= 48,
            walltime= '00:02:00',
            constraint='hsw',
            account= None,
            job_name= 'dfttestrun',
            pre_rocket= "module load cp2k-env/4.1-hsw",
            post_rocket= "echo 'current fashion: post-modern rocket after running dft'",
            logdir= abspath,
            rocket_launch= rlaunch_command)
        lightweight = CommonAdapter(
            q_type="SLURM",
            queue="test",
            nodes= 1,
            ntasks= 8,
            walltime= '00:00:30',
            constraint='hsw',
            account= None,
            job_name= 'light',
            pre_rocket= None,
            post_rocket= None,
            logdir= abspath,
            rocket_launch= rlaunch_command)
        medium = CommonAdapter(
            q_type="SLURM",
            queue="test",
            nodes= 1,
            ntasks= 8,
            walltime= '00:01:00',
            constraint='hsw',
            account= None,
            job_name= 'medium',
            pre_rocket= None,
            post_rocket= None,
            logdir= abspath,
            rocket_launch= rlaunch_command)

        for i in range(0, 400):
            for category, adapter in zip(['dft', 'medium', 'lightweight'], [dft, medium, lightweight]):
                launch_rocket_to_queue(launchpad, FWorker(category=category), adapter,
                    launcher_dir=abspath , create_launcher_dir=True, reserve=True)
            launch_rocket_to_queue(launchpad, FWorker(), lightweight,
                launcher_dir=abspath , create_launcher_dir=True, reserve=True)
            time.sleep(3)

            _recover_offline(lp = launchpad, fworker_name = None)


    else:
        rapidfire(launchpad, FWorker(category=['dft', 'medium', 'lightweight']))


