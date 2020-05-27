Checklist
=========

Before a production run of a workflow, it is recommended to go once through this checklist.

Obligatory
----------

These points must be met:


- Is the information for the launchpad correct? (username, password) in both workflow submission and launching script?
- Are the arguments *username* and *password* correct in the workflow function?
- Does *worker_target_path* on the computing ressource exist?
- Did the test run with a LJ-potential terminate correctly?
- Did you switch to the correct DFT template (*template_path*)?
- Did you change the database to the production run database? (Argument in worfklow: e.g. :code:`extdb_connect = {"db_name": "ncdb"}`)
- Are the walltimes and number of cores in the QueueAdapters of the launching script sensible?

Recommended
-----------

These points are not critical, 
yet it is recommended to check them
for the best user experience:

- Did you set a walltime in the CP2K template? If no walltime is set or the walltime is not a little bit lower than the walltime specified during submission, Fireworks might not mark jobs FIZZLED but categorize them as RUNNING which were cancelled due to runtime. This halts the workflow and requires manual input (see Troubleshooting below).
- Did you provide a consistent (same DFT-level) *reference_energy* or *atomic_energies*? This may reduce your post-processing work.
- Is the location of your logpath ok? Some larger data files and output files are dumped there.


Troubleshooting
---------------

If your workflow stops ahead of time, it might have fizzled or been defused. 

If it fizzled, it means there is an unexpected failure. Investigate the error in the corresponding launch directory. Once you have fixed the issue (e.g. missing file, typo, walltime, bug in the code, ...) you can rerun a workflow from a specific firework through


:code:`lpad -l <LAUNCHPAD> rerun_fws -i <FIREWORKID>`

If the workflow has been defused, there is usually an expected reason for it. A set convergence criterion might have been met or a loop has become static (no more changes detected) ...

In most cases you can choose to keep on running it until the next stop point. Use the command

:code:`lpad -l <LAUNCHPAD> reignite_wflows -i <WORKFLOWID>`


If jobs are marked as RUNNING but have clearly finished, they might be lost. Use :code:`lpad -l <LAUNCHPAD> detect_lostruns`
to look for them. Either add the argument :code:`--fizzle` or
:code:`--rerun`. For more options use the :code:`--help` argument.


If you have found a bug, please raise an issue on github. You can also make a pull request if you have found a solution to the problem.


