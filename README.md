![GitHub release (latest by date)](https://img.shields.io/github/v/release/jjsayleraxio/JTIMLmaster?label=current%20release)
![GitHub Release Date](https://img.shields.io/github/release-date/jjsayleraxio/JTIMLmaster)

# JTIMLmaster

__Author__: Sangsoon Woo, Consultant, Cytel

__Maintainer__: Joseph Sayler, Analyst, Cytel

### INTRODUCTION

This vignette is created to show how to apply machine learning (ML) approaches in building a predictive model. We will use RNA-seq data set to train and evaluate the performance of different ML models. 
The ML approaches considered in this training are:

  * Penalized Regression model (PR)
  * Random Forest (RF)
  * Support Vector Machine (SVM)
  * Neural Network (NN)
  * Gredient Boosting Model (GBM)

The workflow of the ML application to omics data will be following:
  1) Data pre-processing : filtering, normalization
  2) External feature selection based on univariate modeling
  3) Application of ML to pre-processed data 
  4) Selection of predictors 
  5) Model comparison based on accurasy measures or misclassification

### STUDY
The report was based on two Post-traumatic stress disorder (PTSD) data sets. The PTSD affects 7~8% of the general US population and is higher(up to 20%) among troops returned from the wars in Iraq and Afghanistan. This large difference in incidence rates indicates that life-threatening life experience may perturbate molecular level functions related to mentality which may induce the development of PTSD. Therefore, understanding the molecular mechanisms involved might help to reduce the morbidity and mortality associated with PTSD.

#### WTC Data 
The data is RNA-seq data available in Gene Expression Omnibus under the GEO extension GSE97356. Transcriptome-wide expression study using RNA sequencing of whole blood was conducted in 324 World Trade Center responders. This data was used to build a predictive model for PTSD status and top ranked genes obtained from Penalized regression model are included in the final model which is also used to estimate polygenic score for an independent samples (Marine data). 

#### Marine Data 
The data is RNA-seq data available in Gene Expression Omnibus under the GEO extension GSE64813. All subjects in the data were males. Whole blood samples were obtained from 124 MRS II US Marine participants who served a seven month deployment. For each patient, blood was drawn 1 month prior to deployment and again at 3 months after deployment. RNA-seq data was generated from the whole blood samples. The study also includes PTSD status of the patients which will be compared to predicted PTSD status using the predictive model trained in WTC data set. The polygenic score based on predictors was estimated in Marine data sets. 

### Inspiration
Can you train a machine learning model using WTC data set to accurately predict whether or not the people in Marine data have PTSD or not?
