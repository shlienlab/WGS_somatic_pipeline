# WGS_somatic_pipeline
A WDL workflow and wrapper for executing tumour + normal alignment, somatic small variant discovery, somatic structural variant and somatic copy number alteration discovery.

### Installing
1. Clone this repository to you running enviroment.
2. Alter the config.ini file to provide paths for:
   1. Cromwell jar file.
   2. Java path to use with cromwell.
   3. Path to use as the cromwell execution directory.
   4. Cloud Cromwell client API URL.
   5. Number of resource clones to use.
3. Alter the config/cancer_analysis_cromwell_replacement.conf file.
   * This file is the cromwell execution config and is currently set to run for the developers CentOS7 Toque/Moab cluster.
   * Alter this file to work with your system.


### Resources
The resources directory within this repo has empty place holder files to show where to place files for this workflow to run.  The directory labelled zero can be copied and renumbered to provide more resource files.  Cromwell can make extensive use of hard links when running hundreds of samples through the workflows.  On the developers system, CentOS7, there is a hard link limit of 1000 links, and thus it was found necessary to have duplicated resource files.  The number of duplicates can be given within the config.ini file to spread hard links across the duplicated resources if required.
[runtime_parameters]: resource_element_to_use = 0, resource_to_use = 0


### Running the Pipeline
1. Compose your TSV file (use the -h option to find an example TSV file)
    * There is a different format if you are supplying FASTQs or BAMs
2. Start the pipeline:
    * If you have the shlienlab module loaded: Execute `$ start_genome samples.tsv ...`
    * Otherwise you can call it from its full path: Execute `$ /hpf/largeprojects/adam/local/wdl_pipelines/wdl_somatic_cancer_genome_pipeline/start_genome samples.tsv ...`
3. Provide output directories for the parts you want to run:
    * `-a-out <alignment output directory>`: If provided, TSV must be FASTQ type; if not, TSV must be BAM type
    * `-ssm-out <SSM output directory>`
    * `-sv-out <SV output directory>`
    * `-cnv-out <CNV output directory>`
4. If you want to run with non-default inputs, provide the `--custom-inputs` option and follow the instructions
    * Look in the input section of `workflows/somatic_cancer_genome_analysis_workflow.wdl` to see what inputs you can override


### How `start_genome` works
* `start_genome` reads your TSV file and creates the required files and directory structure:
    1. A patient directory in each analysis output directory
    2. Sample directories in each patient directory, where output will be saved:
        * Alignment: `N_-_<id>` or `T_-_<id>` for each tumour or normal
        * SSM/SV/CNV: `N_-_<normal id>+T_-_<tumour id>` for each normal/tumour pair
    3. `wdl_logs+run_info` directory in each sample directory
    4. `progress` directory in each log directory, which will contain info for the resume mechanism
    5. `wdl_run_info` directory in each log directory, which will contain:
        * `ran_from.txt` specifying where run info can be found, or
        * Once per patient: `<patient>_analysis_inputs.json` which contains the inputs for the pipeline, `<patient>_analysis.sh` which is the script to run Cromwell, and `cromwell_logs` which will hold the cromwell log files
        * Each `wdl_run_info` directory will also contain a compressed copy of the WDL file used
* Unless the `--custom-inputs` option is provided, `start_genome` also submits `<patient>_analysis.sh` to the cluster
* Note: A separate Cromwell instance is used for each patient


##### How the Workflows Work
* Task names are defined at the beginning so that they can be passed to the `CheckProgress` task, which returns a `Map[task name, true/false]` specifying which tasks have been completed
* Each task call is then placed in an if block, and will only run if `Map[task name] == "false"` (ie. that task hasn't been completed yet)
* Because the outputs from task calls become optional values (undefined if the task doesn't run), each task's outputs are defined outside of the task
* When specifying the input to the next task, we select either the optional output, or the pre-defined value; ie:
    * `task2_input = select_first([task2.output, task2_output])`
    * This ensures that the execution flow of the workflow remains intact (since this is inferred by output/input chains), while still being able to run `task2` if `task1` has previously been completed (therefore having a null [undefined] output)
    * A special case is after a scatter, where we must specify `task2`'s input as: `if (length(select_all(task1.output)) == length(whatever_was_scattered_on)) then select_all(task1.output) else task1_output`, since `task1.output` is an array of optional values
* When the command of a task completes successfully, a file is written so that if the workflow is restarted, `Map[task name] == "true"` and that task will not re-run

##### Structure of the Repositories
* **bin**: Contains code used by the pipeline
* **config**: Contains a configuration file for Cromwell)
* **dockerfiles**: Contains the files and scripts required to make the necessary containers for the pipeline.
* **example_files**: Contains example TSV files
* **resources**: Should contain all the resources required by tasks within the workflow, e.g. reference genome, GATK processes reference, annotations.  
* **templates**: Contains template files for the inputs json file and the execution bash script
* **workflows**: Contains WDL files:
    * align_workflow.wdl: the alignment pipeline
    * ssm_workflow.wdl: the SSM pipeline
    * sv_workflow.wdl: the SV pipeline
    * cnv_workflow.wdl: the CNV pipeline
    * somatic_cancer\_genome\_analysis\_workflow.wdl: runs the above workflows
    * common_tasks.wdl: contains tasks used by multiple workflows
    * default_values.wdl: contains default input values
* **`start_genome`**: Python script which parses the supplied TSV file and writes the inputs json and execution script to disk, as well as submits the execution script to HPF.
