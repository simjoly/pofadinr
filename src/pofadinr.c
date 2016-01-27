/*

  pofadinr.c

  This file is part of the pofad R package. 

  It is distributed under the GPL 2.

  (c) Simon Joly, Emmanuel Paradis, 2016

*/


#include <R.h>

/* returns 8 if the base is known surely, 0 otherwise */
#define KnownBase(a) (a & 8)

/* Check if the base contains at least one nucleotide */
#define KnownBase_pofad(a) (a & 240) >= 16

/* returns 4 if the base has a gap, 0 otherwise */
#define IsGap(a) (a & 4)

/* returns 1 if both bases are different surely, 0 otherwise */
#define DifferentBase(a, b) (a & b) < 16

#define CHECK_PAIRWISE_DELETION_POFAD\
    if (x[s1] >= 16 && x[s2] >= 16) L++;\
    else continue;

/* 
	Function: count_nucleotides

	Function that counts the number of nucleotides
	present in one sample
*/

int count_nucleotides(int nucleotide) {
    int count = 0;
    int cur = nucleotide >> 4; // Ignore the last four bits (extra info)
    while (cur != 0) {
        cur &= cur - 1;
        count++;
    }
    return count;
}

/*  genpofad distance */

double genpofad_dist(int a, int b){
	int a_union_b = count_nucleotides(a & b);
	int only_a = count_nucleotides(a);
	int only_b = count_nucleotides(b);
	int max_ab = (only_a >= only_b ? only_a : only_b);
	return (1-((double) a_union_b / max_ab));
}

/*  matchstates distance */

double matchstate_dist(int a, int b){
	int a_exclusiveor_b = count_nucleotides(a ^ b);
	int total_nuc = count_nucleotides(a) + count_nucleotides(b);
	return ((double) a_exclusiveor_b / total_nuc );
}


void distDNA_genpofad(unsigned char *x, int *n, int *s, double *d) {
	int i1, i2, s1, s2, L, target;
	double dist;
	L = *s;
	target = 0;
	for (i1=1; i1<*n; i1++) { //loop over the sequences
		for (i2 = i1 + 1; i2 <= *n; i2++) {
			dist = 0;  //total distance between sequences
			for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1 += *n, s2 += *n) {
				if (x[s1] == x[s2]) continue;
				if (DifferentBase(x[s1],x[s2])) {dist+=1; continue;}
				else dist+=genpofad_dist(x[s1],x[s2]);
			}
			d[target] = ((double) dist/L);
			target++;
		}
	}
}

void distDNA_genpofad_pairdel(unsigned char *x, int *n, int *s, double *d) {
	int i1, i2, s1, s2, L, target;
	double dist;
	target = 0;
	for (i1=1; i1<*n; i1++) { //loop over the sequences
		for (i2 = i1 + 1; i2 <= *n; i2++) {
			dist = 0;  //total distance between sequences
			L = 0; //sequence length
			for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1 += *n, s2 += *n) {
				CHECK_PAIRWISE_DELETION_POFAD
				if (x[s1] == x[s2]) continue;
				if (DifferentBase(x[s1],x[s2])) {dist+=1; continue;}
				else dist+=genpofad_dist(x[s1],x[s2]);
			}
			d[target] = ((double) dist/L);
			target++;
		}
	}
}


void distDNA_matchstates(unsigned char *x, int *n, int *s, double *d) {
	int i1, i2, s1, s2, L, target;
	double dist;
	L = *s;
	target = 0;
	for (i1=1; i1<*n; i1++) { //loop over the sequences
		for (i2 = i1 + 1; i2 <= *n; i2++) {
			dist = 0;  //total distance between sequences
			for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1 += *n, s2 += *n) {
				if (x[s1] == x[s2]) continue;
				if (DifferentBase(x[s1],x[s2])) {dist+=1; continue;}
				else dist+=matchstate_dist(x[s1],x[s2]);
			}
			d[target] = ((double) dist/L);
			target++;
		}
	}
}

void distDNA_matchstates_pairdel(unsigned char *x, int *n, int *s, double *d) {
	int i1, i2, s1, s2, L, target;
	double dist;
	target = 0;
	for (i1=1; i1<*n; i1++) { //loop over the sequences
		for (i2 = i1 + 1; i2 <= *n; i2++) {
			dist = 0;  //total distance between sequences
			L = 0; //sequence length
			for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1 += *n, s2 += *n) {
				CHECK_PAIRWISE_DELETION_POFAD
				if (x[s1] == x[s2]) continue;
				if (DifferentBase(x[s1],x[s2])) {dist+=1; continue;}
				else dist+=matchstate_dist(x[s1],x[s2]);
			}
			d[target] = ((double) dist/L);
			target++;
		}
	}
}


void GlobalDeletionDNA(unsigned char *x, int *n, int *s, int *keep)
{
    int i, j;

    for (j = 0; j < *s; j++) {
        i = *n * j;
	while (i < *n * (j + 1)) {
	    if (KnownBase_pofad(x[i])) i++;
	    else {
	        keep[j] = 0;
			break;
	    }
	}
    }
}


void pofad(unsigned char *x, int *n, int *s, int *model, double *d,
	      int *pairdel)
{

	/*
		*x:		the data matrix (DNAbin)
		*n:		number of sequences
		*s:		number of nucleotides
		*model:	distance method
		*d:		the distance matrix (most often asymetrix)
		*pairdel 	pairwise deletion (boolean)
	*/

    switch (*model) {
    case 1 : if (*pairdel) distDNA_genpofad_pairdel(x, n, s, d);
             else distDNA_genpofad(x, n, s, d); break;
    case 2 : if (*pairdel) distDNA_matchstates_pairdel(x, n, s, d);
             else distDNA_matchstates(x, n, s, d); break;
    }

}


void consensusDNA(unsigned char *x, unsigned char *y, int *n, int *s) {

	/*
		*x:		the data matrix
		*y:		the consensus sequence
		*n:		number of sequences
		*s:		number of nucleotides
	*/

	// Note: presently treats gaps as missing data

	int s1, n1;
	unsigned char temp;
	for (s1 = 1; s1 <= *s; s1++) {				//loop over the nucleotides
		temp = 0;
		for (n1 = *n*(s1-1); n1 < *n*(s1-1) + *n; n1++) { //loop over the sequences
			temp = (temp | x[n1]);
			if (count_nucleotides(temp) >= 1) 	//IF at least one nucleotide,
				temp = temp & 248; 					//make sure to remove gap and unknown base
			if (count_nucleotides(temp) > 1) 	//IF more than two nucleotides,
				temp = temp & 246; 					//remove 'known base'
			if (IsGap(temp)) 					//IF nucleotide is a gap,
				temp = temp & 252; 					//make sure to remove unkown base
		}
		y[s1-1] = (unsigned char) temp;
	}

}
