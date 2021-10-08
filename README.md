# A deep learning based pipeline for TSS coverage denoising and feature extraction from shallow cell free DNA sequencing
AECT (Autoencoder on Cell-free DNA TSS coverage profile) is an autoencoder based method to denoise the TSS coverage profiles generated by shallow cfDNA sequencing. 
A set of pre-processing steps on cfDNA sequencing data, including GC bias adjustment, copy number variation normalization were also integrated.
AECT improves robustness of TSS coverage quantification, and improved sensitivity and specificity of TSS profile based classifier.

## Installation  
	
#### install from GitHub
	Using setuptools:
	git clone git://github.com/hanbw0120/AECT
	cd AECT
	python setup.py install
	
	Alternatively, you can install all dependent packages and directly run the script as "python AECT.py <parameters>":
	conda install numpy pandas keras
	conda install tensorflow=1.14.0

## Quick Start

#### Usage
AECT.py -i <input_data> -f [file_type] -o [output_data]

#### Input
* a matrix file in either csv (default) or tab-seperated format
* row is TSS, and column is sample
* an example data is provided in [改]

#### Other options  
* size of the hidden layers, separated with commas, default is 128,64,32,64,128: [--hidden_size]
* batch size, default is 32: [-b] or [--batch_size]
* epochs fr training, default is 500: [-e] or [epochs]
change iterations by watching the convergence of loss, default is 30000: [-i] or [--max_iter]  
* reduces learning rate if validation loss does not improve in a given number of epochs, default is 10: [--reduce_lr]
* stops training if validation loss does not improve in a given number of epochs, default is 15: [--early_stop]

#### Run use Docker image
* a Docker file was upload to: [改]
* load the Docker file: 
	gunzip -c AECT_v0.9.tar.gz | docker load
* run AECT:
	docker run --ipc=host -v "/PATH/TO/DATA/":"/PATH/TO/DATA/" aect:v0.9 AECT.py -i /PATH/TO/DATA/example/example_data.csv -o out.csv

## Other Scripts

#### GC adjustment
Usage: python run_gc_adj.py -b <input bam> -d <reference bed> -g <genome2bit>

* Deeptools (https://github.com/deeptools/deepTools) and pybedtools (https://pypi.org/project/pybedtools) are required 
* input bam should have been removed PCR duplicates
* Reference bed file contains information of TSSs, in our study, we used -1~+1kbp of TSS, such as "chr17	7589868	7591868	NM_000546	0	-	TP53"
* Genome 2bit file can be downloaded from http://hgdownload.cse.ucsc.edu/gbdb/
* Raw counts (*.raw.bed) and gc-adjust counts (*.gc.bed) are generated
* You can modify other parameters refer to instruction of deepTools.

#### CNA normalization
Usage: python cnv_norm.py <input bed> <input cnv>  
* input cnv is generated by DNAcopy package 
