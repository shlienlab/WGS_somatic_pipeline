version 1.0

task CheckProgress {
    input {
        File sa_key
        String sample
        Array[String] task_names
        String progress_dir
        String log_dir
        String environment
    }

    String sa_key_base = basename(sa_key)

    command {
        set -eo pipefail

        my_sa_key=~{sa_key}

        if [ "~{environment}" = "LOCAL" ]; then
            cp ~{sa_key} /scratch/
            my_sa_key=/scratch/~{sa_key_base}
            rm ~{sa_key}

            mkdir -p ~{progress_dir}
            cd ~{progress_dir}

            for i in ~{sep=" " task_names}; do
                if [ ! -f $i ]; then
                    printf "$i\tfalse\n" > $i
                    chmod g+w $i
                fi
            done
            cat ~{sep=" " task_names}
        elif [ "~{environment}" = "CLOUD" ]; then
            set +e
            gcloud auth activate-service-account --key-file ~{sa_key}
            for task in ~{sep=" " task_names}; do
                gsutil -q stat ~{progress_dir}/$task
                if [ $? -eq 1 ]; then
                    printf "$task\tfalse\n" > $task
                    gsutil mv $task ~{progress_dir}/
                fi
            done
            set -e
            # this command line may become too long
            gsutil cat ~{sep=" " prefix(progress_dir + "/", task_names)}
        else
            echo "ERORR: environment must be one of ['LOCAL', 'CLOUD']. ~{environment} is not a valid value." >&2
            exit 1
        fi
    }

    output {
        Map[String, String] task_completion = read_map(stdout())
    }

    runtime {
        docker: "gcr.io/dnastack-dropbox-157402/cloud-sdk-base:1"
        memory: "4 GB"
        walltime: "02:30:00"
        task_name: "CheckProgress.~{sample}"
        log_dir: log_dir
        output_dir: progress_dir
    }
}

task WriteManifest {
    input {
        File sa_key
        String patient
        String sample
        String experiment = "wgs"
        String tumor_normal
        String file_path
        String file_type
        String environment
        String task_name
        String log_dir
    }

    String manifest = "~{file_path}.irodsmanifest"

    command {
        set -eo pipefail

        if [ "~{environment}" = "LOCAL" ]; then
            rm ~{sa_key}

            if [ ! -f "~{manifest}" ]; then
                echo -e patient_id'\t'sample_id'\t'experiment_type'\t'tumor_normal'\t'hpf_file_path'\t'file_type >> ~{manifest}
                echo -e ~{patient}'\t'~{sample}'\t'~{experiment}'\t'~{tumor_normal}'\t'~{file_path}'\t'~{file_type} >> ~{manifest}
            fi
        elif [ "~{environment}" = "CLOUD" ]; then
            set +e
            gcloud auth activate-service-account --key-file ~{sa_key}

        else
            echo "ERROR: environment must be one of ['LOCAL', 'CLOUD']. ~{environment} is not a valid value." >&2
            exit 1
        fi
    }

    output {
        String ManifestWritten = "~{manifest}"
    }

    runtime {
        memory: "4 GB"
        walltime: "02:30:00"
        task_name: "WriteManifest.~{sample}"
        log_dir: log_dir
    }
}

task FindWorkDir{
    input {
        File sa_key
        String patient
        String sample
        String my_workflow
        String environment
        String log_dir
        String output_dir
    }

    command <<<
        set -eo pipefail
        if [ "~{environment}" = "LOCAL" ]; then
            rm ~{sa_key}
            MY_PATH=$PWD
            DELIMITER="/~{my_workflow}/"
            TO_DEL=${MY_PATH%"$DELIMITER"*}
            echo "$TO_DEL" >> "~{output_dir}/execution_dirs.txt"
        fi
    >>>

    output {
        String exe_dir = "~{output_dir}/execution_dirs.txt"
    }

    runtime {
        memory: "4 GB"
        walltime: "02:30:00"
        task_name: "FindWorkDir.~{patient}.~{sample}.~{my_workflow}"
        log_dir: log_dir
        output_dir: output_dir
    }
}

