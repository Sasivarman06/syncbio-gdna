
# <code>syncbio-gdna <b>unlock</b></code>

## 1. About 
The `syncbio-gdna` executable is composed of several inter-related sub commands. Please see `syncbio-gdna -h` for all available options.

This part of the documentation describes options and concepts for <code>syncbio-gdna <b>unlock</b></code> sub command in more detail. With minimal configuration, the **`unlock`** sub command enables you to unlock a pipeline output directory. 

If the pipeline fails ungracefully, it may be required to unlock the working directory before proceeding again. Snakemake will inform a user when it may be necessary to unlock a working directory with an error message stating: `Error: Directory cannot be locked`. 

Please verify that the pipeline is not running before running this command. If the pipeline is currently running, the workflow manager will report the working directory is locked. This is the default behavior of snakemake, and it is normal. Do NOT run this command if the pipeline is still running! Please kill the master job and its child jobs prior to running this command.

Unlocking syncbio-gdna pipeline output directory is fast and easy! In its most basic form, <code>syncbio-gdna <b>unlock</b></code> only has *one required input*.

## 2. Synopsis
```text
$ syncbio-gdna unlock [--help] [--output OUTPUT]
```

The synopsis for this command shows its parameters and their usage. Optional parameters are shown in square brackets.

A user **must** provide an output directory to unlock via `--output` argument. After running the unlock sub command, you can resume the pipeline from where it left off by re-running it. 

You can always use the `-h` option for information on a specific command. 

### 2.1 Required Arguments  

  `--output OUTPUT` 
> **Output directory to unlock.**  
> *type: path*
> 
> Path to a previous run's output directory. This will remove a lock on the working directory. Please verify that the pipeline is not running before running this command.  
> ***Example:*** `--output /data/$USER/syncbio-gdna_results`

### 2.2 Options

Each of the following arguments are optional and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

## 3. Example
```bash 
# Step 1.) Unlock a pipeline output directory
syncbio-gdna unlock --output /data/$USER/syncbio-gdna_results

# Step 2.) Resume the pipeline after unlocking
syncbio-gdna run \
    --samplesheet my_samples.csv \
    --output /data/$USER/syncbio-gdna_results \
    --rerun-incomplete
```

