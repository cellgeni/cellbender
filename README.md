# Cellbender
Repo containing scripts for running cellbender

### Directory structure should be as followed:

top level directories - actions, data, work

actions - contains all the scripts which are used to run souporcell

data - contains subdirectories for each sample with sampleid as name, output is written to these subdirectories

work - contains logs subdirectory which stores all the log files for running cellbender

### Getting Data

Ensure the file `irods.txt` contains 2 tab separated fields. The first contains a sample IDs and the second contains a path to the location on irods. 
Each row contains a single sample ID and its corresponding irods path as shown below:

```bash
HCA_A_LNG11599017_and_HCA_A_LNG11986504	/seq/illumina/cellranger-arc/cellranger-arc201_count_308b45a587ebfa6f443e3bb754bb625f
HCA_A_LNG12177503_and_HCA_A_LNG11986505	/seq/illumina/cellranger-arc/cellranger-arc201_count_aa07d3102ccae90a449dca434e47b104
```

Run the script `get-cellbender.sh` inside the `data` directory and it will download all the relevant data. If you have different file names to the default edit `get-cellbender.sh` as the comments instruct. I recommend using a screen session for this process.
To run the script simply:

```bash
cd data
../actions/get-cellbender.sh
```

### Running cellbender

Ensure your sample ID list is contained within the `samples.txt` file. It should be a single field file containing one sample ID per row. I.e.

```bash
HCA_A_LNG11599017_and_HCA_A_LNG11986504
HCA_A_LNG12177503_and_HCA_A_LNG11986505
```

Run the script `bsub-cellbender.sh` inside the `work` directory. 

```bash
cd work
../actions/bsub-cellbender.sh
```

If you need to change the resources used then edit the `bsub-cellbender.sh` script as appropriate. If you want to change the cellbender parameters
then edit the `cellbender-0-2-0.sh` script.

Jobs will be submitted to the FARM, they can be monitored with the command `bjobs` or you can look at the log and error files within the `logs` subdirectory.

### Quality Control

Running the script `cellbender.qc.R` will assess the quality of the output. It runs on a single cellbender output directory at a time. I.e.

```bash
cd data/sampleID
Rscript cellbender.qc.R cellbender
```

Please note this script installs some packages so there needs to be a writable path in r libs path.
