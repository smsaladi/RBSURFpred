RBSURpred
==========
RBSURpred is a extension of the REGAd3p predictor for predicting absolute accessible surface area of residues 
in protein. The predictor is  based on regularized exact regression 
method using 3rd order polynomial as kernel and  with different metaheuristics for the
optimization of weights. It takes protein sequence (standard FASTA format) 
as input and outputs per residue predicted ASA values (both absolute and binary values of ASA). This software 
also predicts secondary structure internally to use the probabilities
of three types of secondary structure (helix, coil and beta).

By Sumit Tarafder, Toukir Ahmed Real, Sumaiya Iqbal, Md Tamjidul Hoque and Dr. M. Sohel Rahman.

Availability
============
RBSURpred is available at:
https://github.com/Sumit46/RBSURFpred


Quick Run Guide
===============

1)     SET input
	- Redirect into RBSURpred/Software/Input
	- Follow the instructions in "ReadMe_Input.txt" to complete "id_list.txt"
	
2)  - SET path variables within script 'run_REGAd3p+.sh'
	- SET path of PSI-BLAST (BLAST/bin) and NR database
	- SET path of IUPred source codes (given within the software in AdditionalFiles directory)
	- SET path of libSVM installation directory 

3) - Make sure to run "make" command in libsvm package installation directory.
    - Make sure to install mysql-server using linux ommand - "sudo apt-get install mysql-server"	
	
4) Run prediction
	- Redirect into RBSURpred/Software/Scripts	
	- Execute 
		- 'predict.sh'
		- run predict.sh (SET the permission if required)
	- Prediction outputs of secondary structure will be at "Output/prediction/$id/SS"
		- 'id.SSp'
	- Prediction outputs of ASA(real value) will be at "Output/prediction/$id/ASA"
		- 'id.ASApnew'
	- Prediction outputs of ASA(RSA or binary value) will be at "Output/prediction/$id/ASA"
		- 'id.RSA'		
For detail description of Software, users are requested to read 'ReadMe' in 'RBSURpred' package.
		
Thanks!!
ENJOY!!		












