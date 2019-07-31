#!/bin/bash

#purpose:Prepare and collect 58 features and make final feature file.
#author: Sumit Tarafder

export IUPred_PATH=../AdditionalFiles/iupred
iupred_path="/home/sumit/RBSURpred/RBSURpred/Software/AdditionalFiles/iupred"
input_path="../Input";
id_path="../Input/list";
source_code_path="../Codes";
libsvm_path="../AdditionalFiles/libsvm";

sequence="../Features";
feature="../../../../Features";

spinex_path="../AdditionalFiles/spineXpublic";
davar_path="../AdditionalFiles/davar";

ext=".fasta"

#The ncbi-blast folder is in(should be)the same directory as libsvm...

psiblast_path="../../../ncbi-blast-2.4.0+/bin";
nr_database_path="../../../ncbi-blast-2.4.0+/bin/nr/nr";

fastapath="../Input/FASTA";

cd ../Codes
g++ -o inputasa prepareInput_ASA.cpp
g++ -o asaf collectFeatures_ASA.cpp

cd ../Scripts

printf "\nStarting------------------\n";
file="$1"
printf "ID: "$file"\n";
	#===========================================================================================================
	# check and process input
	if [ -f $input_path/FASTA/$file.fasta ];
	then
		if [ -f ../Features/$file/$file.fasta ];
		then
			printf "Directory and fasta file already exists!!\n";
		else	
			cd $source_code_path;
			g++ -o start processInput.cpp
			./start $file;
			printf "processing Input...";
			printf "...DONE!!!\n";
			rm start;
			cd ../Scripts
		fi	
	else
		printf "Error --- FASTA file not found in Input/FASTA directory!!!\n";
		exit 1;
	fi
	
	#===========================================================================================================
	#Run IUPred (short & long)
	if [ -f ../Features/$file/$file.iupredS ] && [ -f ../Features/$file/$file.iupredL ];
	then
		printf "IUPred short & long already exists!!\n";
	else
		printf "running IUPred (short & long)...";
		cc $iupred_path/iupred.c -o $iupred_path/iupred;
		$iupred_path/iupred $sequence/$file/$file.fasta short > $sequence/$file/$file.iupredS;
		$iupred_path/iupred $sequence/$file/$file.fasta long > $sequence/$file/$file.iupredL;
		printf "...DONE!!!\n";
	fi	
	
	#===========================================================================================================
	#===========================================================================================================
	#Collect PSSM
	if [ -f ../Features/$file/$file.pssm ];
	then
		printf "PSSM already exists!!\n";
	else
		printf "Collecting PSSM...";
		filename="$file$ext"

        #printf "running PSI-BLAST For ID : $line ...";
		
		#$psiblast_path/psiblast -query $fastapath/$filename -db $nr_database_path -out $sequence/$file/$file.out -num_iterations 3 -num_threads 16 -out_ascii_pssm $sequence/$file/$file.pssm;
		cp ../$file.pssm ../Features/$file
		printf "...DONE!!!\n";
	fi	
	
	#===========================================================================================================
	#Run Monogram, Bigram computation
	if [ -f ../Features/$file/$file.monogram ] && [ -f ../Features/$file/$file.bigram ];
	then
		printf "Monorgam & Bigram exists!!\n";
	else
		cd $source_code_path;
		printf "generating monogram and bigram...";
		g++ -o mg_bg computeMG_BG.cpp
		./mg_bg $file
		printf "...DONE!!!\n";
		rm mg_bg;
	fi	
	#===========================================================================================================
	#===========================================================================================================
	#Gather 52 features for SS prediction
	
	if [ -f ../Features/$file/$file.initialSS.features ]
	then
			printf "Initial SS feature exists!!\n";
	else 
	    cd $source_code_path;
	   # printf "collecting features for SS prediction...";
	    g++ -o ssf collectFeatures_SS.cpp
	    ./ssf $file
	    #printf "...DONE!!!\n";
	    rm ssf;
	fi
	
	#===========================================================================================================
	#==========================================================================================================
	#apply windowing and prepare final input for initial SS prediction by libsvm
	if [ -f ../Features/$file/$file.initialSS.input ]
	then
		printf "SS input file exists!!\n";
	else 
	    cd $source_code_path;
	    printf "generating input for SS prediction...";
	    g++ -o inputss prepareInput_SS.cpp
	    ./inputss $file
	    printf "...DONE!!!\n";
	    rm inputss;
	fi
	#===========================================================================================================
	
	#===========================================================================================================
	#predict secondary structure and process output
	cd $source_code_path;
	if [ -f ../Output/prediction/$file/SS/$file.SSp ];
	then
		printf "Predicted Secondary structure output already exists!!\n";
	else
		mkdir ../Output/prediction/$file;
		mkdir ../Output/prediction/$file/SS;
		printf "predicting secondary structure and processing output...";
		g++ -o predictss predictSS_processOutput.cpp
		./predictss $file $libsvm_path
		printf "...DONE!!!\n";
		rm predictss;
	fi	
	#===========================================================================================================
	#===========================================================================================================
	#Gather 55 features for ASA prediction
	cd $source_code_path;
	printf "collecting features for ASA prediction...";
	g++ -o asaf collectFeatures_ASA.cpp
	./asaf $file 0 52
	rm asaf
	printf "...DONE!!!\n";
	#===========================================================================================================
	#===========================================================================================================
	#Prepare exact input for ASA prediction
	cd $source_code_path;
	printf "generating input for ASA prediction...";
	
	./inputasa $file 55 9
	printf "...DONE!!!\n";
	#===========================================================================================================

	#===========================================================================================================
	# Predict ASA
	cd $source_code_path;
	#mkdir ../Output/prediction/$file;
	mkdir ../Output/prediction/$file/ASA;
	
	printf "predicting ASA...";
	g++ -o predictasa predictASA_processOutput.cpp
	./predictasa $file 55 2 1 1
    mv ../Output/prediction/$file/ASA/$file.ASApnew ../Features/$file/$file.ASAp;
	printf "...DONE!!!\n";
	rm predictasa;
	#===========================================================================================================
	
	#===========================================================================================================
	#Run PSEE computation
	
	cd $source_code_path;
	printf "generating PSEE...";
	g++ -o psee PSEE_calculation.cpp
	./psee $file 1 9
	printf "...DONE!!!\n";
	rm psee;
	#===========================================================================================================
	
	#delete previous id.ASA.features and id.ASA.input
	
	rm ../Features/$file/$file.ASA.features
	rm ../Features/$file/$file.ASA.input
	
	#===========================================================================================================
	#gather angle fluctuatuons
	cd ../Scripts
	#===========================================================================================================
	#Run SPINE-X
	if [ -f ../Features/$file/$file.spXout ];
	then
		printf "SPINE X output exists!!\n";
	else
		printf "running SPINE-X...";
		truncate $spinex_path/test/list1 --size 0; 
		echo $file > $spinex_path/test/list1;
		mv ../Features/$file/$file.fasta ../Features/$file/$file
		cp ../Features/$file/$file $spinex_path/test/sequence;
		cp ../Features/$file/$file.pssm ../Features/$file/$file.mat
		mv ../Features/$file/$file.mat $spinex_path/test/profile;
		chmod +x $spinex_path/spX.pl;
		$spinex_path/spX.pl $spinex_path/test/list1 $spinex_path/test/profile/ > ../Output/log/log_spinex.txt;
		cp $spinex_path/test/spXout/$file.spXout ../Features/$file;
		rm -r spXout;
		printf "...DONE!!!\n";
		mv ../Features/$file/$file ../Features/$file/$file.fasta
	fi
	#===========================================================================================================
	
	#===========================================================================================================
	# run DAVAR for angle fluctuation
	if [ -f ../Features/$file/$file.dphi ] && [ -f ../Features/$file/$file.dpsi ];
	then
		printf "phi fluctuation & psi fluctuation exists!!\n";
	else
		# collecting features for DAVAR
		printf "collecting features for angle fluctuation prediction...";
		cd $source_code_path;
		g++ -o davar_features collect_davar_features.cpp
		./davar_features $file;
		printf "...DONE!!!\n";
		rm davar_features;
		
		#Run DAVAR
		cd $davar_path/
		printf "running angle (phi,psi) fluctuation prediction...";
		
		#echo "Round 1...";
		cd dphi/round1/;
		truncate genn_test.db --size 0;
		echo "$feature/$file/$file.davar_features" >> genn_test.db;
		chmod +x genn_test.pl;
		./genn_test.pl > mytmp.txt;
	
		#echo "Round 2...";
		cd ../round2/
		truncate genn_test.db --size 0;
		echo "$feature/$file/$file.davar_features" >> genn_test.db;
		chmod +x genn_test.pl;
		./genn_test.pl > mytmp.txt;
		
		#echo "Round 3...";
		cd ../round3/
		truncate genn_test.db --size 0;
		echo "$feature/$file/$file.davar_features" >> genn_test.db;
		chmod +x genn_test.pl;
		./genn_test.pl > mytmp.txt;
		
		#echo "Round 4...";
		cd ../round4/
		truncate genn_test.db --size 0;
		echo "$feature/$file/$file.davar_features" >> genn_test.db;
		chmod +x genn_test.pl;
		./genn_test.pl > mytmp.txt;
		
		#echo "Round 5...";
		cd ../round5/
		truncate genn_test.db --size 0;
		echo "$feature/$file/$file.davar_features" >> genn_test.db;
		chmod +x genn_test.pl;
		./genn_test.pl > mytmp.txt;
		
		cd ../../
	
		#echo "Round 1...";
		cd dpsi/round1/;
		truncate genn_test.db --size 0;
		echo "$feature/$file/$file.davar_features" >> genn_test.db;
		chmod +x genn_test.pl;
		./genn_test.pl > mytmp.txt;
		
		#echo "Round 2...";
		cd ../round2/
		truncate genn_test.db --size 0;
		echo "$feature/$file/$file.davar_features" >> genn_test.db;
		chmod +x genn_test.pl;
		./genn_test.pl > mytmp.txt;
		
		#echo "Round 3...";
		cd ../round3/
		truncate genn_test.db --size 0;
		echo "$feature/$file/$file.davar_features" >> genn_test.db;
		chmod +x genn_test.pl;
		./genn_test.pl > mytmp.txt;
		
		#echo "Round 4...";
		cd ../round4/
		truncate genn_test.db --size 0;
		echo "$feature/$file/$file.davar_features" >> genn_test.db;
		chmod +x genn_test.pl;
		./genn_test.pl > mytmp.txt;
	
		#echo "Round 5...";
		cd ../round5/
		truncate genn_test.db --size 0;
		echo "$feature/$file/$file.davar_features" >> genn_test.db;
		chmod +x genn_test.pl;
		./genn_test.pl > mytmp.txt;
		
		cd ../../../../Codes;
		g++ -o phi_psi collect_phi_psi_fluctuations.cpp
		./phi_psi
		rm phi_psi
		
		cp test.dphi ../Features/$file/$file.dphi;
		cp test.dpsi ../Features/$file/$file.dpsi;
		
		rm test.dphi;
		rm test.dpsi;
		printf "...DONE!!!\n";
	fi


	#===========================================================================================================
	#Gather 58 features for ASA prediction(includes PSEE)
	cd $source_code_path;
	#printf "collecting features for ASA prediction...";
	g++ -o asaf collectFeatures_ASA.cpp
	#52/31
	./asaf $file 1 52
	rm asaf
	#printf "...DONE!!!\n";
	
	#===========================================================================================================
	#===========================================================================================================
	#Prepare exact input for ASA prediction
	cd $source_code_path;
	#printf "generating input for ASA prediction...";
	
	./inputasa $file 58 9
	
	#printf "$file...DONE!!!\n";
	cd ../Scripts
	

cd ../Codes
rm inputasa

    #===========================================================================================================
printf "\nFeature Extraction Complete.\n";
