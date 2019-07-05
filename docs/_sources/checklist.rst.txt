Checklist
=========

Obligatory
----------

These points must be met:

coming soon



- Did the test run with a LJ-potential terminate correctly?
- Did you change the database to "ncdb"? (Argument in worfklow: extdb_connect = {"db_name": "ncdb"},)
- Is the information for the launchpad correct? (username, password) in both workflow submission and launching script?
- Are the walltimes in the launching script sensible?
- Does worker_target_path (absolute path) exist and is it void of previous runs? (older folders and files might interfere with the CP2K calculation)

Recommended
-----------

These points are not critical, 
yet it is recommended to check them
for the best user experience:

coming soon


Troubleshooting
---------------

If your workflow stops ahead of time, it might have fizzled or been defused. 

If it fizzled, it means there is an unexpected failure. Investigate the error in the corresponding launch directory. Once you have fixed the issue (e.g. missing file, typo, walltime, bug in the code, ...) you can rerun a workflow from a specific firework through

```sh
lpad -l <LAUNCHPAD> rerun_fws -i <FIREWORKID>
```

If the workflow has been defused, there is usually an expected reason for it. A set convergence criterion might have been met or a loop has become static (no more changes detected) ...

In most cases you can choose to keep on running it until the next stop point. Use the command

```sh
lpad -l reignite_wflows -i <WORKFLOWID>
```


If you have found a bug, please raise an issue. You can also make a pull request if you have found a solution to the problem.


