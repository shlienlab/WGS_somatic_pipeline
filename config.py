from configparser import ConfigParser
import os


def verify_resources(resource_dir):
    output = []
    if not os.path.isdir(RESOURCE_DIR):
        raise NotADirectoryError(
            "The resources directory is not present.  Ensure the resources directory "
            "is created and formated properly. "
        )
    first_dir_list = os.listdir(os.path.join(resource_dir, RESOURCE_DIRS_TO_CHECK[0]))
    for element in first_dir_list:
        present_in_all = True
        for dir in RESOURCE_DIRS_TO_CHECK:
            if not os.path.isdir(os.path.join(RESOURCE_DIR, dir, element)):
                present_in_all = False
        if present_in_all:
            output.append(element)
    return output


def save_resource_to_use(new_value):
    config.set("runtime_parameters", "RESOURCE_TO_USE", str(new_value))
    with open("config.ini", "w") as configfile:
        config.write(configfile)


config = ConfigParser()
config.read("config.ini")

CROMWELL_JAR_PATH = config.get("cromwell_paths", "CROMWELL_JAR_PATH")
CROMWELL_JAVA_VERSION = config.get("cromwell_paths", "CROMWELL_JAVA_VERSION")
CROMWELL_EXECUTION_DIRECTORY = config.get(
    "cromwell_paths", "CROMWELL_EXECUTION_DIRECTORY"
)
CROMWELL_CLOUD_URL = config.get("cromwell_paths", "CROMWELL_CLOUD_URL")
RESOURCE_ELEMENT_TO_USE = config.getint("runtime_parameters", "resource_element_to_use")

REPO_PATH = os.path.dirname(os.path.realpath(__file__))
CONFIG_FILE = REPO_PATH + "/config/cancer_analysis_cromwell_replacement.conf"
ENTRY_WORKFLOW_FILE = (
    REPO_PATH + "/workflows/somatic_cancer_genome_analysis_workflow.wdl"
)
ALIGN_WORFKLOW_FILE = REPO_PATH + "/workflows/align_workflow.wdl"
SSM_WORKFLOW_FILE = REPO_PATH + "/workflows/ssm_workflow.wdl"
SV_WORKFLOW_FILE = REPO_PATH + "/workflows/sv_workflow.wdl"
CNV_WORKFLOW_FILE = REPO_PATH + "/workflows/cnv_workflow.wdl"
GRIDSS_WORKFLOW_FILE = REPO_PATH + "/workflows/gridss_workflow.wdl"
INPUTS_TEMPLATE_FILE_PATH = (
    REPO_PATH + "/templates/cancer_analysis_workflow_inputs.json.template"
)
EXECUTE_TEMPLATE_FILE_PATH = REPO_PATH + "/templates/execute_analysis.sh.template"
EXAMPLE_FASTQ_TAB_FILE = REPO_PATH + "/example_files/example_fastq_input.tsv"
EXAMPLE_BAM_TAB_FILE = REPO_PATH + "/example_files/example_bam_input.tsv"
RESOURCE_DIR = REPO_PATH + "/resources"

LOG_DIR_RELATIVE = "/wdl_logs+run_info"
PROGRESS_DIR_RELATIVE = LOG_DIR_RELATIVE + "/progress"
RUN_DIR_RELATIVE = LOG_DIR_RELATIVE + "/wdl_run_info"
CROMWELL_LOG_DIR_NAME = "cromwell_logs"
RESOURCE_DIRS_TO_CHECK = [
    "bwa_reference",
    "cnv",
    "dummy_key",
    "ssm",
    "ssm_and_sv_reference",
    "sv",
]
RESOURCE_ELEMENTS = verify_resources(RESOURCE_DIR)
if len(RESOURCE_ELEMENTS) == 0:
    raise ResourceWarning(
        "Resource directory is not implemented corrected.  See README for instructions."
    )
