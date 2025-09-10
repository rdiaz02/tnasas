//   Copyright (C) 2003, 2004,2005, 2006  Ramón Díaz-Uriarte

//  This program is free software; you can redistribute it and/or
//  modify it under the terms of the GNU General Public License
//  as published by the Free Software Foundation; either version 2
//  of the License, or (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


///////////////////////////////////////////////////////////////////////////////



// Code for Diagonal Linear Discriminant Analysis.
// Two functions:
//     - dlda: returns the predictions, given a
//             training set and a test set.
//     - dldaPredError: returns relative prediction
//             error given training and test sets (and
//             true classes of test set cases).

// DLDA is like typical discriminant analysis, except we assume
// diagonal var-cov matrices (i.e., covariance of variables = 0)
// and, for each variable, the same variance in all classes.

// The code is originally based on the stat.diag.da code from Dudoit and
// Fridlyand, in the sma package for R. A very direct translation of the
// the R code to C++ is the dlda_original. I then tried to make it faster
// and the result is dlda. I also searched for C/C++ code for DLDA and LDA
// in statlib and netlib and I didn find anything (there is code in Fortran,
// but for general discriminant analysis).

// For now, this code is supposed to be called as dynamically loaded
// compiled code from R, using ".C".


// 2025-09-18: the code is failing with segfault.



#include<cmath>
//#include<iostream>

//xx: use R_alloc instead of new and delete.


extern "C" {

// dlda: Diagonal Linear Discriminat Analysis. Arguments are:

//  ls: learning set matrix of covariates. In ls subjects change faster
//      than genes; in R ls is a matrix with subjects as rows
//      (and genes as columns).

//  cll: the class to which each subject belongs. cll goes from 0
//      to (num_classes - 1) in steps of 1.

//  nk: number of subjects in each class.

//  ts: testing set, covariates. In ts genes change faster than subjects;
//      in R ts is a matrix with genes as rows. Makes life simpler.
//      For testing, obviously, just transpose ls and pass as ts.

//  num_classes: number of classes (or length of cll + 1).

//  nsubjects_ls: the number of subjects in the learning set.

//  num_genes: yeap, number of genes (or columns in ls).

//  nsubjects_ts: number of subjects in the training set.

//  predictions: the vector of predicted classes for the training set.
//       This is what we want.


// The basic steps of dlda are:
//     a) Get means and variances of every variable (gene)
//        for every class.
//     b) Get the pooled variance for every gene (pooled over
//        classes).
//     c) For every test case, find (squared) diference between observed
//        value of variable (gene) and mean of variable of each
//        class, divided by variance of variable.
//     d) The predicted class is the one with smaller value in c).


void dlda(double *ls, int *cll, int *nk, double *ts, int *num_classes,
	  int *nsubjects_ls,
	  int *num_genes, int *nsubjects_ts, int *predictions) {

  int n_total = (*num_genes) * (*nsubjects_ls);

  // mean and var arrays: genes vary faster within classes
  double *sum_x = new double[(*num_classes) * (*num_genes)];
  double *sum_x2  = new double[(*num_classes) * (*num_genes)];
  double *vars  = new double[*num_genes];


  //xx: add error checking of new allocation

  // Zero the allocated memory
  for (int i = 0; i < ((*num_classes) * (*num_genes)); ++i) {
    sum_x[i] = 0.0;
    sum_x2[i] = 0.0;
  }
  for (int i = 0; i < (*num_genes); ++i) vars[i] = 0.0;


  // Get means and variances of each class for each gene, and pooled variance.
  int k = 0;
  int m = 0;
  int current_class;
  double current_val;
  int pos_g_by_c_array;
  for (int gene_index = 0; gene_index < *num_genes; ++gene_index) {
    for (int subj_ind = 0; subj_ind < *nsubjects_ls; ++subj_ind) {
      current_class = cll[subj_ind];
      current_val = ls[k++];
      sum_x[(current_class * (*num_genes)) + gene_index] += current_val;
      sum_x2[(current_class * (*num_genes)) + gene_index] +=
	(current_val * current_val);
    }

    for (int class_index = 0; class_index < *num_classes; ++class_index) {
      pos_g_by_c_array = class_index * (*num_genes) + gene_index;
      sum_x2[pos_g_by_c_array] -=
	(sum_x[pos_g_by_c_array] * sum_x[pos_g_by_c_array]) / nk[class_index];
      //variance (without divisor)
      vars[gene_index] += sum_x2[pos_g_by_c_array];//pooled variance;
      //no need to divide by ((*nsubjects_ls) - (*num_classes)) since
      // common to all genes.
      sum_x[pos_g_by_c_array] /= nk[class_index];
    }
  }

  // Discriminant scores and prediction avoiding looping several times.
  int predicted_class;
  double actual_value;
  double past_score;
  double *score_class = new double[(*num_classes)];
  double tmp;
  k = 0;

  for(int subject_ts = 0; subject_ts < *nsubjects_ts; ++subject_ts) {
    for (int m = 0; m < (*num_classes); ++m) score_class[m] = 0.0;

    for (int gene_index = 0; gene_index < *num_genes; ++gene_index) {
      actual_value = ts[k++];
      for (int class_index = 0; class_index < *num_classes; ++class_index) {
	tmp = actual_value - sum_x[(class_index * (*num_genes)) + gene_index];
	score_class[class_index] += tmp * tmp / vars[gene_index];
// 	score_class[class_index] +=
// 	  pow(actual_value - sum_x[(class_index *(*num_genes))+gene_index], 2)
// 	  /vars[gene_index];
      }
    }

    //Now, get smallest score
    predicted_class = 0;
    past_score = score_class[0];
    for(int class_index = 1; class_index < *num_classes; ++class_index) {
      if(score_class[class_index] < past_score) {
	past_score = score_class[class_index];
	predicted_class = class_index;
      }
    }
    predictions[subject_ts] = predicted_class;
  }

  delete [] score_class;
  delete [] sum_x;
  delete [] sum_x2;
  delete [] vars;
}





//dldaPredError returns the relative error

  void dldaPredError(double *ls, int *cll, int *nk, double *ts, int *tcll,
		     int *num_classes,
		     int *nsubjects_ls, int *num_genes, int *nsubjects_ts,
		     double *PredError) {
    int *predictions = new int[*nsubjects_ts];
    dlda(ls, cll, nk, ts, num_classes, nsubjects_ls, num_genes, nsubjects_ts,
	 predictions);

    int PredNeReal = 0;
    *PredError = 0.0;
    for(int i = 0; i < (*nsubjects_ts); ++i) {
      if(predictions[i] != tcll[i]) ++PredNeReal;
    }
    *PredError = static_cast<double>(PredNeReal)/(*nsubjects_ts);
    delete [] predictions;

  }


/* void dlda_original(double *ls, int *cll, double *ts,  */
/* 		   int *num_classes, int *nsubjects_ls,  */
/* 	  int *num_genes, int *nsubjects_ts, int *predictions) { */
/* // A (very direct) translation to C++ of the code in the R package  */
/* // sma, function stat.diag.da, */
/* // by Sandrine Dudoit & Jane Fridlyand, version 0.5.10. The sma package */
/* // has GNU GPL license. */
/* // I implement it here using common variance. */
/* // This code is kept here for historical reasons. The function */
/* // can be optimized, as done above. */

/*   int n_total = (*num_genes) * (*nsubjects_ls); */

/*   int *nk = new int[*num_classes]; */
/*   for (int i = 0; i < (*num_classes); ++i) nk[i] = 0; */

/*   // Number of subjects in each class */
/*   for (int i = 0; i < *nsubjects_ls; ++i) { */
/*     ++nk[cll[i]]; */
/*   } */


/*   // genes vary faster within classes */
/*   double *sum_x = new double[(*num_classes) * (*num_genes)]; */
/*   double *sum_x2  = new double[(*num_classes) * (*num_genes)]; */
/*   double *vars  = new double[*num_genes];  */
/*   //since using DLDA, pooled over classes. */

/*   //xx: add error checking of new allocation */

/*   // Zero the allocated memory */
/*   for (int i = 0; i < ((*num_classes) * (*num_genes)); ++i) { */
/*     sum_x[i] = 0.0; */
/*     sum_x2[i] = 0.0; */
/*   } */
/*   for (int i = 0; i < (*num_genes); ++i) vars[i] = 0.0; */



/*   // Get means and variances of each class for each gene */
/*   int current_class; */
/*   int current_gene; */
/*   for (int i = 0; i < n_total; ++i) { */
/*     current_class = cll[i % (*nsubjects_ls)]; */
/*     current_gene = i / (*nsubjects_ls); */
/*     sum_x[(current_class * (*num_genes)) + current_gene] += ls[i]; */
/*     sum_x2[(current_class * (*num_genes)) + current_gene] += (ls[i] * ls[i]); */
/*   } */


/*   int k = 0; */
/*   for (int class_index = 0; class_index < *num_classes; ++class_index) { */
/*     for (int gene_index = 0; gene_index < *num_genes; ++gene_index) { */
/*       sum_x2[k] = (sum_x2[k] - ((sum_x[k] * sum_x[k])/nk[class_index])); */
/*       //	/(nk[class_index] - 1); //variance */
/*       sum_x[k] /= nk[class_index]; //mean */
/* 	++k; */
/*     } */
/*   } */

/*   // Pooled variance */
/*   k = 0; */
/*   for (int class_index = 0; class_index < *num_classes; ++class_index) { */
/*     for (int gene_index = 0; gene_index < *num_genes; ++gene_index) { */
/*       vars[gene_index] += sum_x2[k]; */
/*       ++k; */
/*     } */
/*   } */
/*   for (int gene_index = 0; gene_index < *num_genes; ++gene_index) { */
/*     vars[gene_index] /= ((*nsubjects_ls) - (*num_classes)); */
/*   } */

/*   // Discriminants scores */
/*   double *discr = new double[(*nsubjects_ts) * (*num_classes)];  */
/*   for (int i = 0; i < ((*nsubjects_ts) * (*num_classes)); ++i) discr[i] = 0.0; */

/*   k = 0; */
/*   for (int subject_ts = 0; subject_ts < *nsubjects_ts; ++subject_ts) { */
/*     for (int gene_index = 0; gene_index < *num_genes; ++gene_index) { */
/*       for (int class_index = 0; class_index < *num_classes; ++class_index) { */
/*       discr[ (subject_ts * (*num_classes)) + class_index] +=  */
/* 	pow(ts[k] -  sum_x[(class_index * (*num_genes)) + gene_index], 2) */
/* 	/vars[gene_index]; */
/*       } */
/*       ++k; */
/*     } */
/*   } */


/*   // Predicted class */
/*   int predicted_class; */
/*   double min_score; */

/*   for(int subject_ts = 0; subject_ts < *nsubjects_ts; ++subject_ts) { */
/*     min_score = discr[(subject_ts * (*num_classes))]; */
/*     predicted_class = 0; */
/*      for (int class_index = 1; class_index < *num_classes; ++class_index) { */
/*        if(discr[(subject_ts * (*num_classes)) + class_index] < min_score) { */
/* 	 predicted_class = class_index; */
/* 	 min_score = discr[(subject_ts * (*num_classes)) + class_index]; */
/*        } */
/*      } */
/*      predictions[subject_ts] = predicted_class; */
/*   } */


/*    delete [] nk; */
/*    delete [] sum_x; */
/*    delete [] sum_x2; */
/*    delete [] vars; */
/*    delete [] discr; */
/* } */

} // extern "C"
