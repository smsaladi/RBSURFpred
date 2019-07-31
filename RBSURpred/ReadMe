RBSURpred
==========
RBSURpred is a extension of the REGAd3p predictor for predicting absolute accessible surface area of residues 
in protein. The predictor is  based on regularized exact regression 
method using 3rd order polynomial as kernel and  with different metaheuristics for the
optimization of weights. It takes protein sequence (standard FASTA format) 
as input and outputs per residue predicted ASA values(both absolute and binary values of ASA). This software 
also predicts secondary structure internally to use the probabilities
of three types of secondary structure (helix, coil and beta).

By Sumit Tarafder, Toukir Ahmed Real, Sumaiya Iqbal, Md Tamjidul Hoque and Dr. M. Sohel Rahman.
Availability
============
RBSURpred is available at:

https://github.com/Sumit46/RBSURFpred


***QUICK RUN***
===============
See 'QuickRunGuidelines' for quick and complete run of RBSURpred.


***The rest of this file contains detail description of the software.

Table of Contents
=================
- How To Run RBSURpred
	- Prerequisite 
	- Input feature format
	- Brief of source codes
- Software material guideline
- Examples
- Quick Run Guidelines	

How To Run RBSURpred
=====================

Prerequisite
============

1) PSI-BLAST with NR database (from NCBI toolkit) 
Availability > ftp://ftp.ncbi.nih.gov/blast/
Output > PSSM
		 
3) IUPred (Installed in 'AdditionalFiles' directory of Software)
Availability > http://iupred.enzim.hu/
Output > prediction of intrinsically unstructured protein 
Note: 
	> The output (short and long disorder) 	
	
4) libSVM
Availability > http://www.csie.ntu.edu.tw/~cjlin/libsvm
Note:
	> used to generate 3 class classification model for secondary structure
	> need to predict secondary structure probabilities per residue which are used as features to predict ASA
	
5) GCC
Availability > http://gcc.gnu.org/
The source codes are written in C/C++. To compile and execute, GCC is needed.

6) DAVAR
Included in AdditionalFiles folder.

7) SPINE-X
Also included in AdditionalFiles folder.


Input feature format
====================

For a protein sequence, following per residue features must be collected in space separated columns:

1) Column 1: actual annotation
				- actual value, collected from DSSp
				- for a test protein sequence, use dummy: 1 for SS and 0 for ASA		
2) Column 2: amino acid type
				- one numerical value out of 20
				-  (1, 2, 3, ..., 20) for (A R N D C Q E G H I L K M F P S T W Y V)
				
3) Column 3 - 9: 7 physical parameters
					- available at RBSURpred/Software/AdditionalFiles/physiochemical_properties.txt	

4) Column 10 - 29: 20 PSSMs
					- generate with PSI-BLAST and normalize using 9.0	
										
5) Column 30 - 50: 1 monogram and 20 bigrams
				- generate using source code provided with RBSURpred
				- normalize the values with exponent (6.0)
				
6) Column 51 - 52: long and short disorder probabilities
				- collected by running IUpred

7) Column 53 - 55: 3 secondary structure (helix, strand and coil) probabilities
					- generate using SVM model integrated within RBSURpred

8) Column 56     : 1 Position Specific Estimated Energy indicator
					- generate using RBSURpred model
9) Column 57     : 1 angle fluctuation feature namely delta phi
                   - generated using DAVAR software

10) Column 58    :  1 angle fluctuation feature namely delta psi
                   - generated using DAVAR software

11) Column 59    : 1 terminal indicator
				- 5 residues from N-terminal: -1.0, -0.8, -0.6, -0.4, -0.2
				- 5 residues from C-terminal: +1.0, +0.8, +0.6, +0.4, +0.2
				- rest: 0.0
				
Brief of source codes
======================

1) processInput.cpp
	- creates a subdirectory for each 'id' in Features/
	- copy .fasta file into that directory
	- If fasta file not found in Input, exit.
	- log files are generated in the Output/log directory

2) Generate PSSM 
	- Integrated within the script run_REGAd3p+.sh
		- SET path of PSI-BLAST (BLAST/bin) and NR database within the script
	- creates the output with name 'id.out', 'id.pssm' in RBSURpred/Software/Features/id

3) Execute IUPred 
	- Integrated within the script run_REGAd3p.sh
		- SET path of installed IUPred within the script
	- creates the output with name 'id.iupredS', 'id.iupredL' in RBSURpred/Software/Features/id
	
4) computeMG_BG.cpp
	- create two file in Features/id directory 
		- name format: 'id.monogram' and 'id.bigram'
	- log files are generated in the Output/log directory

5) collectFeatures_SS.cpp
	- collects 52 features pre residue for secondary structure prediction
	- Keep the output with name 'id.initialSS.features' in RBSURpred/Software/Features/id
	
6) prepareInput_SS.cpp
	- apply windowing on 52 features for SS prediction and generate final formatted input for libsvm
	- Keep the output with name 'id.initialSS.input' in RBSURpred/Software/Features/id
	
	
7) predictSS_processOutput.cpp	
	- Predict SS by SVM model and process (format) the output
		- SET path of installed libsvm within the script
	- creates the output with name 'id.ss.svm.scale', 'id.ss.svm.predict' in RBSURpred/Software/Features/id
	- log files are generated in the Output/log directory
	- final formatted output is kept in /Output/prediction/$id/SS/'id.SSp'
	
8) collectFeatures_ASA.cpp
	- collects 55 features pre residue for ASA prediction
	- Keep the output with name 'id.ASA.features' in RBSURpred/Software/Features/id
	
9) prepareInput_ASA.cpp
	- apply windowing on 55 features for ASA prediction with 
	an extra bias column of value 1 and generate final formatted input for exact regression
	- Keep the output with name 'id.ASA.input' in RBSURpred/Software/Features/id
		
10) predictASA_processOutput.cpp	
	- Predict ASA by exact method and process (format) the output
	- log files are generated in the Output/log directory
	- final formatted output is kept in /Output/prediction/$id/ASA/'id.ASAp'

11) PSEE_calculation.cpp

    - Calculate sequence composition based Position Specific Estimated Energy (sPSEE)
    - Neighborhood size can be varried by specifying the start and end the region around both side of the target residue(not recommended).

12) collect_davar_features.cpp

    - Collect davar features from SPINEX and other features for DAVAR software in order
    - to do angle fluctuations prediction

13) collect_phi_psi_fluctuations.cpp

    - Collect delta phi and  delta psi fluctuation output from DAVAR Software

Software material guideline
================================

Data
====
- SSD_TR1001.list
	- list of protein ids for training 
- SSD_TS295.list
	- list of protein ids for testing
	
Additional Files
=================
- physiochemical_properties.txt 
	- contains physical parameters for 20 amino acids
- IUPred
	- intrinsically unstructured protein prediction used as features
- EASA_new.txt
	- Contains ASA normalizing values in Gly - X -Gly reference state to calculate RSA 
- DT_predicted_energy.txt
	- Energy matrix file for the calculation of PSEE
		
Models
=======
- SVM model for secondary structure prediction
	- RBSURpred/Software/Models/SS_SVM/SSD_TR1001.libsvm.scale.model
- Weights for ASA prediction (55 features)
	- RBSURpred/Software/Models/ASA_weight/weight1.txt
- Weights for ASA prediction (58 features)
	- RBSURpred/Software/Models/ASA_WEIGHT_REGAd3p_plus/weight1.txt
			
Features
========
- A protein ('id') specific feature files and intermediate outputs throughout the prediction:
- id.fasta
	- FASTA formatted file
- id.pssm
	- PSSM
- id.iupredS
	- IUPred (short) output
- id.iupredL
	- IUPred (short) output
- id.monogram
	- Monograms
- id.bigram
	- Bigrams
- id.PSEE
	- Newly added PSEE feature
- id.dphi
	- Newly added angle fluctuation feature named delta phi
- id.dpsi
	- Newly added angle fluctuation feature named delta psi
- id.initialSS.features
	- per residue 52 features for SS prediction
- id.intialSS.input
	- final input with windowing for SS prediction
- id.ss.svm.scale
	- scaled feature file for SS prediction
- id.ss.svm.predict
	- output of SS prediction from libsvm (not formatted), used as feature	
- id.initialASA.features
	- per residue 55 features for ASA prediction
- id.ASA.input
	- final input with windowing for ASA prediction
	
Codes
======
- Source code (in C/C++) necessary to use RBSURpred+(guidelines mentioned above)

Scripts
=======
- predict.sh
	- script for complete run of ASA prediction
	- run_REGAd3p+.sh is used to extract all the features necessary for ASA prediction.
	
		
Input
=====
- FASTA directory contains all proteins to be predicted in standard FASTA format with .fasta extension
- id_list.txt in "list" foldercontains id list of proteins to be predicted without any extension
	
Output
======
- log
	- contains log of intermediate steps throughout prediction

- prediction
	- SS
		- id.SSp ===> secondary structure prediction output by SVM
		- Format: residue serial number, amino acid, predicted secondary structure and three probabilities
	- ASA
		- id.ASApnew ==> ASA prediction output by RBSURpred
		- Format: residue serial number, amino acid, predicted ASA
	- RSA 
	    - Prediction outputs of ASA(RSA or binary value) will be at "Output/prediction/$id/ASA"
		- 'id.RSA'
Example
=======
- Sample files are given in Input folder.
	- use it for sample run and find the outputs in Output folder.
	
Quick Run Guidelines
=====================
Available in Software/QuickRunGuidelines
	- Use only those short number of steps for complete and quick run!!


Thanks!!
ENJOY!!
