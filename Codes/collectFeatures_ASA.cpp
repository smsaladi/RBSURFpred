/*
Author: Sumaiya Iqbal
Edited By : Sumit Tarafder
Part of RBSURpred software
Collect 55 features for initial ASA prediction by RBSURpred
*/

#define _CRT_SECURE_NO_DEPRECATE
#define LINE_SIZE 30000 //Maximum line size that can be read from the input file
#define BIGRAM_MONOGRAM_LOG_MEDIAN 6.0

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <malloc.h>
#include<ctype.h>
#include<math.h>

FILE *nfp;
FILE *fasta;
FILE *svm_predict;
FILE *feature2;
FILE *feature1;
FILE *dphi;
FILE *dpsi;
FILE *spsee;
FILE *count_p;

char make_d[100];
char terminal[10];
char normalOuputFile[] = "../Output/log/log_feature_collection_for_ASA.txt";
char ppFile[] = "../AdditionalFiles/physiochemical_properties.txt";
char id[LINE_SIZE];
char svm_predict_file_name[LINE_SIZE];
char svm_prediction[LINE_SIZE];
char feature_set2_file_name[LINE_SIZE];
char feature_set1_file_name[LINE_SIZE];
char sPSEE_file_name[LINE_SIZE];
char energy_line[LINE_SIZE];
char dphiFile[300];
char dpsiFile[300];
char philine[100];
char psiline[100];
char cline[100];


char rline[LINE_SIZE];
char fseq[LINE_SIZE];
char all_feature[LINE_SIZE];
char feature_part_1[LINE_SIZE];
char all_feature_filename[300];
char fastaFileName[300];
char cflag[100];
char PSEE_val[1000];
char count_line[100];

int aa_counter = 0;
int pp_counter = 0;
int lineLength = 0;
int seqNumber = 0;
int seqLength = 0;
int seqNum = 0;
int local_counter = 0;
int flag = 0;
int i =0;
int l=0;
int count1=0;
int count2=0;

/*Start of SUB-ROUTINE: substring - It returns a pointer to the substring */
char *substring(char *string, int position, int length)
{
	char *pointer;
	int c;

	pointer = (char *)malloc(length + 1);

	if (pointer == NULL)
	{
		printf("Unable to allocate memory.\n");
		exit(EXIT_FAILURE);
	}

	for (c = 0; c < position - 1; c++)
		string++;

	for (c = 0; c < length; c++)
	{
		*(pointer + c) = *string;
		string++;
	}

	*(pointer + c) = '\0';

	return pointer;
}
/*End of SUB-ROUTINE: substring - It returns a pointer to the substring */


/* SUB-ROUTINE: Trim */
char* ltrim(char *s)
{
	while (isspace(*s)) s++;
	return s;
}

char *rtrim(char *s)
{
	char* back = s + strlen(s);
	while (isspace(*--back));
	*(back + 1) = '\0';
	return s;
}

char *trim(char *s)
{
	return rtrim(ltrim(s));
}

int main(int argc, char *argv[])
{
	// collect ID
	strcpy(id, argv[1]);

	strcpy(cflag, argv[2]);
	flag = atoi(cflag);
	count1 = atoi(argv[3]);
	count2 = count1+1;

	nfp = fopen(normalOuputFile, "wb+");
	if (nfp == NULL) {
		fprintf(stderr, "Can't log output file!\n");
		exit(1);
	}

	//====================================================================================================================================
	// open svm predict file for ss probabilities
	strcpy(svm_predict_file_name, "../Features/");
	strcat(svm_predict_file_name, id);
	strcat(svm_predict_file_name, "/");
	strcat(svm_predict_file_name, id);
	strcat(svm_predict_file_name, ".svm.ss.predict");

	svm_predict = fopen(svm_predict_file_name, "r");					   // Open svm predict from step1 file for reading
	if (svm_predict == NULL) {
		fprintf(stderr, "Can't open svm predict file!\n");
		exit(1);
	}

	// skip header
	fgets(svm_prediction, sizeof svm_prediction, svm_predict);	// skipp header

	//====================================================================================================================================

	//====================================================================================================================================
	// prepare fasta file name, open and collect sequence length
	strcpy(fastaFileName, "../Features/");
	strcat(fastaFileName, id);
	strcat(fastaFileName, "/");
	strcat(fastaFileName, id);
	strcat(fastaFileName, ".fasta");


	//printf("Fasta File name: %s\n", fastaFileName);

	fasta = fopen(fastaFileName, "r");
	if (fasta == NULL) {
		fprintf(stderr, "Can't open fasta File!\n");
		exit(1);
	}

	fgets(fseq, sizeof fseq, fasta);	// skip header
	fgets(fseq, sizeof fseq, fasta);	// collect sequence
	fclose(fasta);

	int tr = 0;
	seqLength = 0;
	while ((fseq[tr] >= 'A') && (fseq[tr] <= 'Z'))
	{
		seqLength++;
		tr++;
	}
	//====================================================================================================================================

	//====================================================================================================================================


	strcpy(feature_set2_file_name,"../Features/");
	strcat(feature_set2_file_name, id);
	strcat(feature_set2_file_name, "/");
	strcat(feature_set2_file_name, id);
	strcat(feature_set2_file_name, ".ASA.features");
	fprintf(nfp, "feature_set2_file_name: %s\n", feature_set2_file_name);

	feature2 = fopen(feature_set2_file_name, "wb+");					// Open all feature2 file for writing
	if (feature2 == NULL) {
		fprintf(stderr, "Can't open rsa feature file!\n");
		exit(1);
	}
	//====================================================================================================================================

	//====================================================================================================================================
	// prepare seature set 1 file name and read on that file
	strcpy(feature_set1_file_name, "../Features/");
	strcat(feature_set1_file_name, id);
	strcat(feature_set1_file_name, "/");
	strcat(feature_set1_file_name, id);
	strcat(feature_set1_file_name, ".initialSS.features");


	feature1 = fopen(feature_set1_file_name, "r");			// Open feature step 1 file for reading
	if (feature1 == NULL) {
		fprintf(stderr, "Can't open feature1 file!\n");
		exit(1);
	}


	if(flag ==1)///Extra 3 features
    {

        strcpy(sPSEE_file_name, "../Features/");
        strcat(sPSEE_file_name, id);
        strcat(sPSEE_file_name, "/");
        strcat(sPSEE_file_name, id);
        strcat(sPSEE_file_name, ".PSEE");

        spsee = fopen(sPSEE_file_name, "r");							// open annotated fasta file and read annotation
        if (spsee == NULL) {
            fprintf(nfp,"%s\n", id);
            fprintf(stderr, "Can't open energy file!\n");
            exit(1);
            //continue;
        }
        fgets(energy_line, sizeof energy_line, spsee);
        fgets(energy_line, sizeof energy_line, spsee);


        // prepare dphi output file, open to read
        strcpy(dphiFile, "../Features/");
        strcat(dphiFile, id);
        strcat(dphiFile, "/");
        strcat(dphiFile, id);
        strcat(dphiFile, ".dphi");

        dphi = fopen(dphiFile, "r");							// open fasta file and read sequence
        if (dphi == NULL) {
            fprintf(stderr, "Can't open dphi file!\n");
            exit(1);
        }

        //====================================================================================================================================
        // prepare dpsi output file, open to read
        strcpy(dpsiFile, "../Features/");
        strcat(dpsiFile, id);
        strcat(dpsiFile, "/");
        strcat(dpsiFile, id);
        strcat(dpsiFile, ".dpsi");

        dpsi = fopen(dpsiFile, "r");							// open fasta file and read sequence
        if (dpsi == NULL) {
            fprintf(stderr, "Can't open dpsi file!\n");
            exit(1);
        }
		// Take psi fluctuation value

    }
	//====================================================================================================================================

	//====================================================================================================================================
	// Merger all features
	// Feature serial is: class, residue indication, pp, pssm, ss, sa, torsion angle, monogram, bigram, terminal
	int res_pos = 0,counter=0;
	while (res_pos < seqLength)
	{
		//====================================================================================================================================
		// Take ASA annottaion --- dummy 0
		strcpy(all_feature, "0 ");
		//====================================================================================================================================

		//====================================================================================================================================
		// take feature set by step 1
		fgets(rline, sizeof rline, feature1);

		int tk_count = 0;
		char *tk;

		tk = strtok(rline, " ");
		while (tk != NULL)
		{
			tk_count++;
			if ((tk_count >= 2) && (tk_count <= count1))
			{
				strcat(all_feature, tk);						// take aa + pp + pssm + bi/mono + iu == 51 features
				strcat(all_feature, " ");
			}
			if (tk_count == count2)
			{
				strcpy(terminal, tk);					// save terminal
			}
			tk = strtok(NULL, " ");
		}
		//====================================================================================================================================

		//====================================================================================================================================

		// collection svm prediction probabilities -- SS0 information
		fgets(svm_prediction, sizeof svm_prediction, svm_predict);

		// collect SS0 probabilities

		char *ss_prob;
		int ss_count = 0;
		ss_prob = strtok(svm_prediction, " ");
		while (ss_prob != NULL)
		{
			ss_count++;
			if ((ss_count >= 2) && (ss_count <= 4))
			{
				strcat(all_feature, trim(ss_prob));						// take ss0 probabilities == 3 features
				strcat(all_feature, " ");
			}
			ss_prob = strtok(NULL, " ");
		}

		///+++++++++++++++++++++++++Edit Starts Collect PSEE features+++++++++++++++++++++

		if(flag == 1)//then collect PSEE else don't
        {
            // prepare sPSEE file name and open to read
            fgets(energy_line, sizeof energy_line, spsee);
            i = 0;
            char *te;
            te = strtok(energy_line, " ");
            while (te != NULL)
            {
                i++;
                if (i == 6)
                {
                    strcpy(PSEE_val, te);
                    l = strlen(PSEE_val);
                    l--;
                    PSEE_val[l] = '\0';
                }

                te = strtok(NULL, " ");
            }
            //printf("\n %d. PSEE val = %s\n",counter+1,PSEE_val);
            strcat(all_feature, PSEE_val);
            strcat(all_feature, " ");

            fgets(philine, sizeof philine, dphi);
            fgets(psiline, sizeof psiline, dpsi);

            // Take phi fluctuation value
            int phi_l = strlen(philine);
            phi_l = phi_l - 2;
            char* phi_sub = substring(philine, 1, phi_l);
            strcat(all_feature, phi_sub);
            strcat(all_feature, " ");
		//====================================================================================================================================

		//====================================================================================================================================
		// Take psi fluctuation value
            int psi_l = strlen(psiline);
            psi_l = psi_l - 2;
            char* psi_sub = substring(psiline, 1, psi_l);
            strcat(all_feature, psi_sub);
           strcat(all_feature, " ");


        }


		///+++++++++++++++++++++++++++++++++Edit ends Collect PSEE features++++++++++++++++++++++++++
		//====================================================================================================================================

		//====================================================================================================================================
		// take terminal information
		strcat(all_feature, trim(terminal));
		//====================================================================================================================================

		//====================================================================================================================================
		fprintf(feature2, "%s\n", all_feature);
		//====================================================================================================================================

		res_pos++;
	}

	fclose(nfp);
	fclose(feature2);
	fclose(feature1);
	fclose(svm_predict);
	//fclose(spsee);
	//fclose(dphi);
	//fclose(dpsi);
}
