version 1.0

import "common_tasks.wdl" as CommonTasks

workflow GRIDSS {
	input {
		String output_dir
		String patient
		String sample
		File tumor_bam
		File tumor_bai
		File normal_bam
		File normal_bai
		File snv_vcf
		File gridss_purple_linx_data
		File reference
		File reference_index
		File reference_dict
		File reference_amb
		File reference_ann
		File reference_bwt
		File reference_sa
		File reference_pac

		String repo_dir
		File sa_key
		String environment
	}

	# Directories
	String log_dir = "~{output_dir}/wdl_logs+run_info"
	String progress_dir = "~{log_dir}/progress"

	# Task names
	String filterBAM_Normal_task_name = "~{patient}.~{sample}.filterBAM_normal"
	String filterBAM_Tumor_task_name = "~{patient}.~{sample}.filterBAM_tumor"
	String GRIDSS_task_name = "~{patient}.~{sample}.gridss"
	# String Remove_Files_task_name = "~{patient}.~{sample}.Remove_Files.GRIDSS"

	Array[String] task_names = [filterBAM_Normal_task_name, filterBAM_Tumor_task_name, GRIDSS_task_name]

	if (environment=="LOCAL") {
		call CommonTasks.FindWorkDir as GRIDSSFindWorkDir {
			input:
				sa_key = sa_key,
				patient=patient,
				sample=sample,
				my_workflow='GRIDSS',
				environment=environment,
				log_dir=log_dir,output_dir = output_dir
		}
	}

	# Task Calls and Output Filenames
	call CommonTasks.CheckProgress {
		input:
			sa_key = sa_key,
			sample = sample,
			task_names = task_names,
			progress_dir = progress_dir,
			log_dir = log_dir,
			environment = environment
	}

	String filterBAM_Tumor_output_bam = "~{output_dir}/~{sample}_T.filtered.bam"
	String filterBAM_Tumor_output_bai = "~{output_dir}/~{sample}_T.filtered.bam.bai"
	if (CheckProgress.task_completion[filterBAM_Tumor_task_name] == "false") {
		call filterBAM as filterBAM_Tumor {
			input:
				bam = tumor_bam,
				bam_index = tumor_bai,
				out_bam = filterBAM_Tumor_output_bam,
				out_bai = filterBAM_Tumor_output_bai,
				log_dir = log_dir,
				task_name = filterBAM_Tumor_task_name,
				environment = environment,
				sa_key = sa_key,
				progress_dir = progress_dir,
				output_dir = output_dir
		}
	}

	String filterBAM_Normal_output_bam = "~{output_dir}/~{sample}_N.filtered.bam"
	String filterBAM_Normal_output_bai = "~{output_dir}/~{sample}_N.filtered.bam.bai"
	if (CheckProgress.task_completion[filterBAM_Normal_task_name] == "false") {
		call filterBAM as filterBAM_Normal {
			input:
				bam = normal_bam,
				bam_index = normal_bai,
				out_bam = filterBAM_Normal_output_bam,
				out_bai = filterBAM_Normal_output_bai,
				log_dir = log_dir,
				task_name = filterBAM_Normal_task_name,
				environment = environment,
				sa_key = sa_key,
				progress_dir = progress_dir,
				output_dir = output_dir
		}
	}

	String run_GRIDSS_output_archive = "~{output_dir}/~{sample}.gridss.tar.gz"
	Array [String] run_GRIDSS_files_to_delete = [filterBAM_Tumor_output_bam, filterBAM_Tumor_output_bai, filterBAM_Normal_output_bam, filterBAM_Normal_output_bai, tumor_bam, tumor_bai]
	if (CheckProgress.task_completion[GRIDSS_task_name] == "false") {
		call run_GRIDSS {
			input:
				sample = sample,
				tumor_bam = select_first([filterBAM_Tumor.output_bam, filterBAM_Tumor_output_bam]),
				tumor_bai = select_first([filterBAM_Tumor.output_bai, filterBAM_Tumor_output_bai]),
				normal_bam = select_first([filterBAM_Normal.output_bam, filterBAM_Normal_output_bam]),
				normal_bai = select_first([filterBAM_Normal.output_bai, filterBAM_Normal_output_bai]),
				snv_vcf = snv_vcf,
				gridss_purple_linx_data = gridss_purple_linx_data,
				reference = reference,
				reference_index = reference_index,
				reference_dict = reference_dict,
				reference_amb = reference_amb,
				reference_ann = reference_ann,
				reference_bwt = reference_bwt,
				reference_sa = reference_sa,
				reference_pac = reference_pac,
				out_archive = run_GRIDSS_output_archive,
				log_dir = log_dir,
				task_name = GRIDSS_task_name,
				environment = environment,
				sa_key = sa_key,
				progress_dir = progress_dir,
				output_dir = output_dir,
				files_to_delete = run_GRIDSS_files_to_delete
		}
	}

	# if (environment=="LOCAL") {
	# 	if (CheckProgress.task_completion[Remove_Files_task_name] == "false") {
	# 		call Remove_Files {
	# 			input:
	# 				patient = patient,
	# 				sample = sample,
	# 				GRIDSS_output_tar = select_first([run_GRIDSS.output_archive, run_GRIDSS_output_archive]),
	# 				files_to_delete = run_GRIDSS_files_to_delete,
	# 				log_dir = log_dir,
	# 				task_name = Remove_Files_task_name,
	# 				environment = environment,
	# 				output_dir = output_dir,
	# 				progress_dir = progress_dir
	# 		}
	# 	}
	# }

	output {
		String output_archive = select_first([run_GRIDSS.output_archive, run_GRIDSS_output_archive])
	}
}

task filterBAM {
	input {
		File bam
		File bam_index
		String out_bam
		String out_bai
		String log_dir
		String task_name
		String environment
		String output_dir
		String progress_dir
		File sa_key
	}

	String bam_base = basename(bam)
	String bam_index_base = basename(bam_index)
	String sa_key_base = basename(sa_key)

	String out_bam_base = basename(out_bam)
	String out_bai_base = basename(out_bai)

	Int threads = 4
	Int disk_size = ceil(size(bam, "GB") * 4 + 100)

	command <<<
		set -eo pipefail

		echo "In filterBAM"
		env

		my_bam=~{bam}
		my_bai=~{bam_index}
		my_sa_key=~{sa_key}

		my_out_bam=~{out_bam_base}
		my_out_bai=~{out_bai_base}

		if [ "~{environment}" = "LOCAL" ]; then

			cp ~{bam} /scratch/
			my_bam=/scratch/~{bam_base}
			rm ~{bam}

			cp ~{bam_index} /scratch/
			my_bai=/scratch/~{bam_index_base}
			rm ~{bam_index}

			cp ~{sa_key} /scratch/
			my_sa_key=/scratch/~{sa_key_base}
			rm ~{sa_key}

			my_out_bam=/scratch/~{out_bam_base}
			my_out_bai=/scratch/~{out_bai_base}
		fi

		samtools view \
		$my_bam \
		-b \
		-o $my_out_bam \
		-@ ~{threads} \
		1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y

		samtools index -@ ~{threads} $my_out_bam

		# Cleanup
		EXIT_STATUS=$?
		echo EXIT STATUS $EXIT_STATUS

		if [ $EXIT_STATUS -eq 0 ]; then
			if [ "~{environment}" = "LOCAL" ]; then
				cp $my_out_bam ~{out_bam}
				cp $my_out_bai ~{out_bai}
				printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
			elif [ "~{environment}" = "CLOUD" ]; then
				gcloud auth activate-service-account --key-file ~{sa_key}
				printf "~{task_name}\ttrue\n" > ./~{task_name}
				gsutil cp ~{out_bam_base} ~{output_dir}
				gsutil cp ~{out_bai_base} ~{output_dir}
				gsutil cp ./~{task_name} ~{progress_dir}
			fi
		else
			exit 1
		fi
	>>>

	output {
		String output_bam = "~{out_bam}"
		String output_bai = "~{out_bai}"
	}

	runtime {
		docker: "gcr.io/dnastack-dropbox-157402/private/profyle/samtools:1.9"
		cpu: threads
		memory: "4 GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 2
		walltime: "48:00:00"
		task_name: task_name
		log_dir: log_dir
		output_dir: output_dir
	}
}

task run_GRIDSS {
	input {
		String sample
		File tumor_bam
		File tumor_bai
		File normal_bam
		File normal_bai
		File snv_vcf
		File gridss_purple_linx_data
		File reference
		File reference_index
		File reference_dict
		File reference_amb
		File reference_ann
		File reference_bwt
		File reference_sa
		File reference_pac
		String out_archive
		String log_dir
		String task_name

		String environment
		String output_dir
		String progress_dir
		File sa_key

		Array [String] files_to_delete
	}

	String linx_ref_genome_version = "HG37"

	String tumor_bam_base = basename(tumor_bam)
	String tumor_bai_base = basename(tumor_bai)
	String normal_bam_base = basename(normal_bam)
	String normal_bai_base = basename(normal_bai)
	String snv_vcf_base = basename(snv_vcf)
	String gridss_purple_linx_data_base = basename(gridss_purple_linx_data)
	String reference_index_base = basename(reference_index)
	String reference_dict_base = basename(reference_dict)
	String reference_amb_base = basename(reference_amb)
	String reference_ann_base = basename(reference_ann)
	String reference_bwt_base = basename(reference_bwt)
	String reference_sa_base = basename(reference_sa)
	String reference_pac_base = basename(reference_pac)
	String sa_key_base = basename(sa_key)

	String out_archive_base = basename(out_archive)
	String out_base = sub(basename(out_archive), ".tar.gz", "")

	String reference_base = basename(reference)
	String reference_name = sub(reference_base, "(\\.fasta|\\.fa)", "")

	Int disk_size = ceil((size(tumor_bam, "GB") + size(normal_bam, "GB") + size(gridss_purple_linx_data, "GB") + size(reference, "GB")) * 1.3 + 50)
	Int threads = 8

	command <<<
		set -eo pipefail

		echo "In GRIDSS"
		env

		my_tumor_bam=~{tumor_bam}
		my_tumor_bai=~{tumor_bai}
		my_normal_bam=~{normal_bam}
		my_normal_bai=~{normal_bai}
		my_snv_vcf=~{snv_vcf}
		my_gridss_purple_linx_data=~{gridss_purple_linx_data}
		my_reference=~{reference}
		my_reference_index=~{reference_index}
		my_reference_dict=~{reference_dict}
		my_reference_amb=~{reference_amb}
		my_reference_ann=~{reference_ann}
		my_reference_bwt=~{reference_bwt}
		my_reference_sa=~{reference_sa}
		my_reference_pac=~{reference_pac}
		my_sa_key=~{sa_key}

		my_out_archive=~{out_archive_base}
		my_out=~{out_base}

		refdata_dir=gridss-dbs
		reference_to_ref_data_dir="$(pwd)/$refdata_dir"

		if [ "~{environment}" = "LOCAL" ]; then
			cp ~{tumor_bam} /scratch/
			my_tumor_bam=/scratch/~{tumor_bam_base}
			rm ~{tumor_bam}

			cp ~{tumor_bai} /scratch/
			my_tumor_bai=/scratch/~{tumor_bai_base}
			rm ~{tumor_bai}

			cp ~{normal_bam} /scratch/
			my_normal_bam=/scratch/~{normal_bam_base}
			rm ~{normal_bam}

			cp ~{normal_bai} /scratch/
			my_normal_bai=/scratch/~{normal_bai_base}
			rm ~{normal_bai}

			cp ~{snv_vcf} /scratch/
			my_snv_vcf=/scratch/~{snv_vcf_base}
			rm ~{snv_vcf}

			cp ~{gridss_purple_linx_data} /scratch/
			my_gridss_purple_linx_data=/scratch/~{gridss_purple_linx_data_base}
			rm ~{gridss_purple_linx_data}

			cp ~{reference} /scratch/
			my_reference=/scratch/~{reference_base}
			rm ~{reference}

			cp ~{reference_index} /scratch/
			my_reference_index=/scratch/~{reference_index_base}
			rm ~{reference_index}

			cp ~{reference_dict} /scratch/
			my_reference_dict=/scratch/~{reference_dict_base}
			rm ~{reference_dict}

			cp ~{reference_amb} /scratch/
			my_reference_amb=/scratch/~{reference_amb_base}
			rm ~{reference_amb}

			cp ~{reference_ann} /scratch/
			my_reference_ann=/scratch/~{reference_ann_base}
			rm ~{reference_ann}

			cp ~{reference_bwt} /scratch/
			my_reference_bwt=/scratch/~{reference_bwt_base}
			rm ~{reference_bwt}

			cp ~{reference_sa} /scratch/
			my_reference_sa=/scratch/~{reference_sa_base}
			rm ~{reference_sa}

			cp ~{reference_pac} /scratch/
			my_reference_pac=/scratch/~{reference_pac_base}
			rm ~{reference_pac}

			cp ~{sa_key} /scratch/
			my_sa_key=/scratch/~{sa_key_base}
			rm ~{sa_key}

			my_out_archive=/scratch/~{out_archive_base}
			my_out=/scratch/~{out_base}

			refdata_dir=/scratch/gridss-dbs
			reference_to_ref_data_dir=$refdata_dir
		fi

		mkdir $my_out

		mkdir $refdata_dir

		tar -zxvf $my_gridss_purple_linx_data --directory $refdata_dir

		mkdir -p $refdata_dir/refgenomes/~{reference_name}
		cp $my_reference $my_reference_index $my_reference_dict $my_reference_amb $my_reference_ann $my_reference_bwt $my_reference_sa $my_reference_pac $refdata_dir/refgenomes/~{reference_name}/

		time \
		/opt/gridss-purple-linx/gridss-purple-linx.sh \
			--threads ~{threads} \
			--jvmheap 60g \
			--sample ~{sample} \
			--tumour_bam $my_tumor_bam \
			--normal_bam $my_normal_bam \
			--ref_dir $reference_to_ref_data_dir \
			--ref_genome_version ~{reference_name}/~{reference_base} \
			--linx_ref_genome_version ~{linx_ref_genome_version} \
			--output_dir $my_out \
			--nosnvvcf

		tar -zcvf $my_out_archive --exclude=*.bam --exclude=*.bai $my_out
		# Adding excludes for bam files and bam indexes, I don't think we need these and they are quite large

		# Cleanup
		EXIT_STATUS=$?
		echo EXIT STATUS $EXIT_STATUS

		if [ $EXIT_STATUS -eq 0 ]; then
			if [ "~{environment}" = "LOCAL" ]; then
				# Don't copy the tar.gz to the output directory, extract it there in place
				#cp $my_out_archive ~{output_dir}
				# untar the archive, but strip the first directory since we know it will be called 'scratch'
				tar -xzvf $my_out_archive --strip 1 -C ~{output_dir}
				printf "~{task_name}\ttrue\n" > ~{progress_dir}/~{task_name}
				#rm ~{sep=' ' files_to_delete}
				for f in ~{sep=' ' files_to_delete}; do [ ! -e $f ] || rm $f; done
			elif [ "~{environment}" = "CLOUD" ]; then
				gcloud auth activate-service-account --key-file ~{sa_key}
				printf "~{task_name}\ttrue\n" > ./~{task_name}
				gsutil cp ~{out_archive_base} ~{output_dir}
				gsutil cp ./~{task_name} ~{progress_dir}
				gsutil rm ~{sep=' ' files_to_delete}
			fi
		else
			exit 1
		fi
	>>>

	output {
		String output_archive = "~{out_archive}"
	}

	runtime {
		docker: "gcr.io/dnastack-dropbox-157402/private/profyle/gridss-purple-linx:1.3.2_1.17_beta"
		cpu: threads
		memory: "80 GB"
		disks: "local-disk " + disk_size + " HDD"
		walltime: "48:00:00"
		log_dir: log_dir
		task_name: task_name
		output_dir: output_dir
	}
}

task Remove_Files{
	input {
		String patient
		String sample
		File GRIDSS_output_tar
		Array [String] files_to_delete
		String log_dir
		String task_name

		String environment
		String output_dir
		String progress_dir
	}

	command {
		set -eo pipefail
		if [ "~{environment}" = "LOCAL" ]; then
			#rm ~{sep=' ' files_to_delete}
			for f in ~{sep=' ' files_to_delete}; do [ ! -e $f ] || rm $f; done
		fi
		# Cleanup
		EXIT_STATUS=$?
		echo EXIT STATUS $EXIT_STATUS
	}

	output {
		String exe_dir = "~{output_dir}/execution_dirs.txt"
	}

	runtime {
		memory: "4 GB"
		walltime: "00:30:00"
		task_name: task_name
		log_dir: log_dir
		output_dir: output_dir
	}
}
