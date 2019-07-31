/*
Author: Sumaiya Iqbal
Edited By : Sumit Tarafder
Part of RBSURpred sofware
--- Predict ASA by exact method with loading weights
--- Process formatted output
*/

#define _CRT_SECURE_NO_DEPRECATE
#define LINE_SIZE 40000					//Maximum line size that can be read from the input file

#define KERNEL 3

#include<stdio.h>
#include<stdlib.h>
#include<fstream>
#include<string.h>
#include <malloc.h>
#include<ctype.h>
#include<math.h>

using namespace std;
FILE *nfp;
FILE *X;
FILE *wgt;
FILE *fasta;
FILE *asap;
FILE *idlist;
FILE* read1;
FILE* read2;
FILE* write1;

char X_filename[200];
char weight_filename[100];
char weight_filename2[200] ;
char fastaFileName[200];
char ASApFileName[200];
char aa[] = "ARNDCQEGHILKMFPSTWYV";

char id[10];
char fsequence[LINE_SIZE];
char featureLine[LINE_SIZE];
char wghtLine[100];
char ASAp_str[500];
char wline[500];
char feat[100];
char no[10];
char fold[10];
char wind_arr[100];
char filename[100];
char readline1[100];
char readline2[100];

double featureVal = 0.0;
double wt_par = 1.15;
double wghtVal = 0.0;
int seqLength = 0;
int resPos = 0;
double ASAp_val = 0.0;
int spc = 0;
int constl = 0;
int FEATURES;
int flag=0,weight_flag=0;

double reference_state[20];
double ASAr[1000];
int res_count=0;


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
	//=============================================================================================================
	// collect id
	int window = 9;
	sprintf(wind_arr,"%d",window);
	strcpy(id, argv[1]);///file
	strcpy(feat, argv[2]);///55
	strcpy(no,argv[4]);///1
	FEATURES = atoi(feat);
	FEATURES = FEATURES*window*3 +1;
	strcpy(feat,argv[3]);///1
	flag = atoi(feat);///1
	strcpy(feat,argv[5]);///1
	weight_flag=atoi(feat);

	//==========================================================================================================================

	//==========================================================================================================================
        //==========================================================================================================================
        // Prepare fasta file name, open and collect fasta sequence for formatted output
        strcpy(fastaFileName, "../Features/");
        strcat(fastaFileName, id);
        strcat(fastaFileName, "/");
        strcat(fastaFileName, id);
        strcat(fastaFileName, ".fasta");

        fasta = fopen(fastaFileName, "r");
        if (fasta == NULL) {
            fprintf(stderr, "Can't open fasta File!\n");
            exit(1);
        }

        fgets(fsequence, sizeof fsequence, fasta);	// skip header
        fgets(fsequence, sizeof fsequence, fasta);	// collect sequence
        int tr = 0;
        seqLength = 0;
        while ((fsequence[tr] >= 'A') && (fsequence[tr] <= 'Z'))
        {
            seqLength++;
            tr++;
        }
        //==========================================================================================================================

        //=========================================================================================================================
        strcpy(X_filename, "../Features/");
        strcat(X_filename, id);
        strcat(X_filename, "/");
        strcat(X_filename, id);
        strcat(X_filename, ".ASA.input");

        X = fopen(X_filename, "r");
        if (X == NULL) {
            fprintf(stderr, "Can't open ASA input File!\n");
            exit(1);
        }

        //==========================================================================================================================

        //==========================================================================================================================
        // prepare predicted ASA file name and write
        strcpy(ASApFileName, "../Output/prediction/");
        strcat(ASApFileName, id);
        strcat(ASApFileName, "/ASA/");
        strcat(ASApFileName, id);
        strcat(ASApFileName, ".ASApnew");

        asap = fopen(ASApFileName, "wb+");
        if (asap == NULL) {
            fprintf(stderr, "Can't open feature File!\n");
            exit(1);
        }
        fprintf(asap, "Absolute Accessible Surface Area (ASA) Prediction Output by REGAd3p\n");
        fprintf(asap, "TARGET: %s, ", id);
        fprintf(asap, "Length: %d\n", seqLength);
        fprintf(asap, "\n");
        fprintf(asap, "SR#  AA  ASAp\n");
        //==========================================================================================================================

        //==========================================================================================================================
        // extend kernel to degree 3 polynomial with feature set internally and multiply with weight to get predicted ASA
        resPos = 0;
        while (resPos < seqLength)
        {
            ASAp_val = 0;
            fgets(featureLine, sizeof featureLine, X);						// collect feature

            if(flag == 1)///58
            {
                strcpy(weight_filename,"../Models/ASA_WEIGHT_RBSURpred_plus");
                strcat(weight_filename,"/weight1.txt");
                wgt = fopen(weight_filename, "r");
                								// open weights for reading
                if (wgt == NULL) {
                    fprintf(stderr, "Can't feature extended File!\n");
                    exit(1);
                }
            }

            else if(flag == 2)///55
            {
                if(weight_flag == 1)
                {
                    strcpy(weight_filename2,"../Models/ASA_WEIGHT/weight");
                    strcat(weight_filename2,no);
                }
                else strcpy(weight_filename2,"../Models/pca_weights1/weight1");
                strcat(weight_filename2,".txt");
                wgt = fopen(weight_filename2, "r");								// open weights for reading
                if (wgt == NULL) {
                    fprintf(stderr, "Can't open weight File!\n");
                    exit(1);
                }
            }
            char *tf;
            tf = strtok(featureLine, " ");
            int token_track = 1;
            //int feat_val =0;
            //feat_val = atoi(no);
            while (tf != NULL)
            {

                if(token_track == 1)
                {
                    featureVal = atof(tf);
                    fgets(wghtLine, sizeof wghtLine, wgt);						// collect feature
                    wghtVal = atof(wghtLine);
                    //printf("\nweight for token track %d is = %lf",token_track,wghtVal);
                    ASAp_val = ASAp_val + (featureVal * wghtVal);
                }
                if((token_track >= 2) && (token_track <= FEATURES))
                {
                    featureVal = atof(tf);

                    fgets(wghtLine, sizeof wghtLine, wgt);									// collect wight
                    wghtVal = atof(wghtLine);
                    ASAp_val = ASAp_val + (featureVal * wghtVal * wt_par);							// 1st order
                    //printf("\nweight for token track %d is = %lf",token_track,wghtVal);
                    fgets(wghtLine, sizeof wghtLine, wgt);									// collect wight
                    wghtVal = atof(wghtLine);
                    ASAp_val = ASAp_val + (featureVal * featureVal * wghtVal* wt_par);				// 2nd order
                    //printf("\nweight for token track %d is = %lf",token_track,wghtVal);
                    fgets(wghtLine, sizeof wghtLine, wgt);									// collect wight
                    wghtVal = atof(wghtLine);
                    //printf("\nweight for token track %d is = %lf",token_track,wghtVal);
                    ASAp_val = ASAp_val + (featureVal * featureVal * featureVal *wghtVal* wt_par);	// 3rd order
                }


                tf = strtok(NULL, " ");
                token_track++;

            }
            //printf("\ntrack = %d",token_track);

            fclose(wgt);

            //==========================================================================================================================
            // write serial number
            sprintf(wline, "%d", resPos + 1);
            spc = 0;
            constl = strlen(wline);
            while (spc < 5 - constl)
            {
                strcat(wline, " ");
                spc++;
            }
            //==========================================================================================================================

            //==========================================================================================================================
            // write amino acid
            char* fresidue1 = substring(fsequence, resPos + 1, 1);		// residues
            strcat(wline, fresidue1);
            strcat(wline, "   ");
            resPos++;
            //==========================================================================================================================

            // write predicted ASA
            if (ASAp_val < 0.0)
            {
                ASAp_val = 0.0;
            }
            else ASAp_val = floor(ASAp_val);
            sprintf(ASAp_str, "%.0lf", ASAp_val);
            strcat(wline, ASAp_str);
            strcat(wline, "   ");

            fprintf(asap, "%s\n", trim(wline));				// write predicted asa
        }
        fclose(asap);
	fclose(X);

	fclose(fasta);


	///=======================Calculating RSA====================================


	strcpy(filename,"../AdditionalFiles/EASA_new.txt");
	read1 = fopen(filename, "r");
    if (read1 == NULL) {
        fprintf(stderr, "Can't open easa file!\n");
        exit(1);
    }
    int index1=0;
    while(fgets(readline1,sizeof readline1,read1))
    {
        reference_state[index1++] = atof(substring(readline1,3,5));
    }
    fclose(read1);

    //for(int i=0;i<index1;i++)cout<<reference_state[i]<<endl;

    //strcpy(id,substring(readline3,0,5));

    strcpy(filename, "../Output/prediction/");
    strcat(filename,id);
    strcat(filename,"/ASA/");
    strcat(filename,id);
    strcat(filename,".ASApnew");

    read2 = fopen(filename, "r");
    if (read2 == NULL) {
        fprintf(stderr, "Can't open reg file for id %s!\n",id);
        exit(1);
    }

    fgets(readline2,sizeof readline2,read2);
    fgets(readline2,sizeof readline2,read2);
    fgets(readline2,sizeof readline2,read2);
    fgets(readline2,sizeof readline2,read2);

    int index=0;
    res_count=0;

    while(fgets(readline2,sizeof readline2,read2))
    {
        char *val1 = substring(readline2,10,7);
        char * residue = substring(readline2,6,1);

        ///Substring method's 2nd Parameter defines [(how much forward it will go)-1] before extracting
        while(1)
        {
            if(!strcmp(residue,substring(aa,index+1,1)))
            {
                break;
            }
            else index++;
        }
        ///Now I know the index
       // cout<<index<<endl;
       // cout<<atof(val1)<<" "<<atof(val2)<<endl;
        ASAr[res_count] = atof(val1)/reference_state[index];

        res_count++;
        index=0;
    }


    strcpy(filename, "../Output/prediction/");
    strcat(filename,id);
    strcat(filename,"/ASA/");
    strcat(filename,id);
    strcat(filename,".RSA");

    write1 = fopen(filename, "w");

    if (write1 == NULL) {
        fprintf(stderr, "Can't create file!\n");
        exit(1);
    }

    for(int i=0;i<res_count;i++)
    {
        fprintf(write1,"%lf\n",ASAr[i]);
    }

    fclose(write1);
    fclose(read2);
	//=============================================================
}
