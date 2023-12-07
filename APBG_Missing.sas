* Purpose: Create simulated data for a presentation;
* Author:  Michael Fitzgerald and Tiffany Evans;
* Date:    26 Nov 2023;
/*****************************************
* Change log
* Name: Tiffany Evans
* Date: 29NOV2023
* Reason: Added single imputation analyses and tipping point analysis using J2R method. 
          Specified reference for Primary analysis and diff and exp options to lsmeans, 
          updated the NegBin_PMI30 macros to use time not log of time;
*****************************************/

%bootloader;

%inc "&g_pathname.Macros\NegBin_PMI30.sas";
%inc "&g_pathname.Macros\TPA_singleImp_wDelta.sas";

%let n=300;
%let seed=727098;

ods word file="\\pha-fs01-uk\Training Area\PhaPharma\PRAXIS2\User_Collab_MLT\Output\Output\APBG_Missing_3TE.docx";

*** Simulate data ***;
data dat1;
  call streaminit(&seed.);
  do subjid=1 to &n.;
    * Tretament;
    trt=rand('BERN', .5);
    * Strata;
    region=rand('BERN', .25);
    * Intercurrent events;
    peot=rand('BERN', .1);
    ttie=rand('UNIFORM',1,51); * Time to intercurrent event;
    if peot=0 then ttie=52;
    ln_ttie=log(ttie);
    * Count data;
    resp1=round(rand('negbinomial',0.7+trt*0.08,2)*ttie/52,1); 

    output;
  end;
run;


*** Event counts by treatement ***;
proc means data=dat1 n mean std median q1 q3 min max;
  class trt;
  var resp1;
run;


*** Review the missing data pattern ***;
proc freq data=dat1;
  tables trt*peot;
run;

proc means data=dat1 n mean std median q1 q3 min max;
  where peot=1;
  class trt;
  var ttie resp1;
run;

proc freq data=dat1;
  where peot=1;
  tables trt*resp1;
run;


**********************************************************************************************;
*** Primary analysis - Negbin ***;
**********************************************************************************************;
proc genmod data=dat1 ;
  class trt region / ref=first;
  model resp1 = trt region / dist=negbin offset=ln_ttie ;
  lsmeans trt / diff ilink cl exp ;
  ods output lsmeans=NBPMI_LSMP diffs=NBPMI_LSMDP ParameterEstimates=parms;
run;


*#################################################################################################;
**** Analyses from AZALEA SAP, 4.2.5.5 first paragraph. 
*#################################################################################################;
*       "First a MAR analysis will be performed where for each participant the rate after withdrawal 
        lambda1 is assumed to be the same as their rate before withdrawal labmda2, which itself is calculated 
        based on their randomized treatment group and baseline covariates."
      We are not certain of the meaning because the first half of the sentence sounds like estimates
      for each individual but second half sounds like model based estimates. We assume model based 
      estimates are required because lambdas do not have the term i for subject, so Lambda 1 is assumed
      to be the model based rate. Also think that the numbering of lambdas should be the other way aroud
      because it is chronologic and that notation has been used in the last paragraph of 4.2.5.5.
*##################################################################################################;



    **********************************************************************************************;
    *** Sensitivity analyses 1 and 2 (MAR, single imputation, model based both pre and post withdrawal)
          Sens analysis 1 is the upper left corner, Sens analysis 2 is the TPA with non-zero deltas ***;
    **********************************************************************************************;
    *     This form of imputation will decrease p-values because all patients who drop out have the 
          same imputed value based on their covariate pattern. Artificially reducing the variance. 
          But the model based imputation will yeild non-zero events (though very small and non-
          integer) for all who dropped out...;
    **********************************************************************************************;
    
    **** Macro takes parameter estimates from the primary analysis above and calculates each subjects 
      linear predictor (corresponding to the predicted mean rate for the covariate pattern);
    %*TPA(assump=1);
      
    
    
    
    **********************************************************************************************;
    *** Sensitivity analysis 1a and 2a (MAR, single imputation, observed pre withdrawal rate plus 
               model-based post-withdrawal rate)
          Sens analysis 1a is the upper left corner, Sens analysis 2a is the TPA with non-zero deltas ***;
    **********************************************************************************************;
    *    This form of imputation will be similar to above but artificial variance reduction might not be as severe?  ;
    **********************************************************************************************;
    **** Macro takes parameter estimates from the primary analysis above and calculates each subjects 
      linear predictor to use for the post-withdrawal imputation of number of events;
    %TPA(assump=2);
    
    
    **********************************************************************************************;
    *** Sensitivity analysis 1b and 2b (MAR, single imputation, observed both pre and post-withdrawal )
          Sens analysis 1a is the upper left corner, Sens analysis 2a is the TPA with non-zero deltas ***;
    **********************************************************************************************;
    *    ;
    **********************************************************************************************;
    **** Macro takes parameter estimates from the primary analysis above and calculates each subjects 
      linear predictor to use for the post-withdrawal imputation of number of events;
    %TPA(assump=3);



        
**********************************************************************************************;
*** Sensitivity analysis 3 - MAR Multiple Imputation ***;
**********************************************************************************************;
* MAR;
%NegBinMI(ds=dat1, 
          treat=trt, 
          withdraw=peot, 
          covariates=region, 
          class=region, 
          response=resp1, 
          time=ttie, 
          Maxtime=52, 
          method=MAR, 
          DeltaV= ,
          ref=0, 
          lsmopt=, 
          label=,
          methodv=, 
          refv=, 
          dependence=COND, 
          alpha=0.05, 
          bayesseed=27438, 
          randseed=78550, 
          out=NegBinMI_MAR_Out, 
          debug=0, 
          nbi=1000, 
          nimpute=100, 
          thin=20);
proc print data=NegBinMI_MAR_Out;
run;

*** Sensitivity analysis 2 - MI Tipping point ***;
* Tipping;
%macro tip1;

  proc datasets lib=work nolist;
    delete _:;
  quit;

  %let cnt=0;
  %do dp=0 %to 6; * Delta placebo - multiply by -0.25 later to get the correct increments;
    %do da=0 %to 6; * Delta active - multiply by 0.25 later to get the correct increments;
      %let cnt=%eval(&cnt.+1);
      %put cnt=&cnt. dp=&dp. da=&da.;
      
      data _dat1;
        set dat1;
        if trt=0 then delta=exp(&dp.*-0.25);
        else if trt=1 then delta=exp(&da.*0.25);
      run;
      
      %NegBinMI(ds=_dat1, 
        treat=trt, 
        ref=0,
        withdraw=peot, 
        covariates=region, 
        class=region, 
        response=resp1, 
        time=ttie, 
        Maxtime=52, 
        method=MAR, 
        DeltaV=delta,
        dependence=COND, 
        alpha=0.05, 
        bayesseed=27438, 
        randseed=78550, 
        out=_NegBinMI_Out, 
        debug=0, 
        nbi=1000, 
        nimpute=100, 
        thin=20);
      
        data _tip&cnt.;
          set _NegBinMI_Out;
          dp=&dp.*-0.25;
          da=&da.*0.25;
        run;
            
    %end;
  %end;

  data _out1(keep=probt est dp da);
    set _tip:;
    where _trt ne .;
    length est $200;
    est=put(estimate, 5.2)||" ("||put(uppcl-lowcl,5.2)||")";
  run;
  
  proc sort data=_out1;
    by dp;
  run;
  
  * Table of p-values;
  proc transpose data=_out1 out=outp(drop=_:) prefix=trt1_;
    by dp;
    id da;
    idlabel da;
    var probt;
  run;
  
  * Format significant p-values as green and not significant as red;
  proc format;
    value back low-0.05 ='green'
               >0.05-high ='lightred';
  run;

  * Output to table;
  ods html style=journal; 
  proc report data=outp split="~";
    columns dp ("Active" trt1:);

    define dp /order order=internal descending "Placebo";
    define trt1: /display style(column)={background=back.};
  run;  
  
  * Table of estimates;
  proc transpose data=_out1 out=oute(drop=_:) prefix=trt1_;
    by dp;
    id da;
    idlabel da;
    var est;
  run;
  
  * Output to table;
  ods html style=journal; 
  proc report data=oute split="~";
    columns dp ("Active" trt1:);

    define dp /order order=internal descending "Placebo";
    define trt1: /display ;
  run;  
  

%mend tip1;
%tip1;


**********************************************************************************************;
***  Controlled MI ***;
**********************************************************************************************;
* Jump to reference - rate after withdrawal is shifted to that of placebo arm, but imputation is conditional on each patients observed rate before withdrawal;

%*NegBinMI(ds=dat1, 
          treat=trt, 
          withdraw=peot, 
          covariates=region, 
          class=region, 
          response=resp1, 
          time=ttie, 
          Maxtime=52, 
          method=J2R, 
          DeltaV= ,
          ref=0, 
          lsmopt=, 
          label=,
          methodv=, 
          refv=, 
          dependence=COND, 
          alpha=0.05, 
          bayesseed=27438, 
          randseed=78550, 
          out=NegBinMI_J2R_Out, 
          debug=0, 
          nbi=1000, 
          nimpute=100, 
          thin=20) ;


* Copy reference - rate before and after withdrawal are shifted to that of reference, but imputation is still conditional on each patients observed rate before withdrawal;
%*NegBinMI(ds=dat1, 
          treat=trt, 
          withdraw=peot, 
          covariates=region, 
          class=region, 
          response=resp1, 
          time=ttie, 
          Maxtime=52, 
          method=CR, 
          DeltaV= ,
          ref=0, 
          lsmopt=, 
          label=,
          methodv=, 
          refv=, 
          dependence=COND, 
          alpha=0.05, 
          bayesseed=27438, 
          randseed=78550, 
          out=NegBinMI_CR_Out, 
          debug=0, 
          nbi=1000, 
          nimpute=100, 
          thin=20) ;

* Tipping point analysis using controlled MI and shift;
%macro tip2;

  proc datasets lib=work nolist;
    delete _:;
  quit;

  %let cnt=0;
  %do dp=0 %to 6; * Delta placebo - multiply by -0.25 later to get the correct increments;
    %do da=0 %to 6; * Delta active - multiply by 0.25 later to get the correct increments;
      %let cnt=%eval(&cnt.+1);
      %put cnt=&cnt. dp=&dp. da=&da.;
      
      data _dat2;
        set dat1;
        if trt=0 then delta=exp(&dp.*-0.25);
        else if trt=1 then delta=exp(&da.*0.25);
      run;
      
      %NegBinMI(ds=_dat2, 
        treat=trt, 
        ref=0,
        withdraw=peot, 
        covariates=region, 
        class=region, 
        response=resp1, 
        time=ttie, 
        Maxtime=52, 
        method=J2R, 
        DeltaV=delta,
        dependence=COND, 
        alpha=0.05, 
        bayesseed=27438, 
        randseed=78550, 
        out=_NegBinMI_Out, 
        debug=0, 
        nbi=1000, 
        nimpute=100, 
        thin=20);
      
        data _tip&cnt.;
          set _NegBinMI_Out;
          dp=&dp.*-0.25;
          da=&da.*0.25;
        run;
            
    %end;
  %end;

  data _out1(keep=probt dp da);
    set _tip:;
    where _trt ne .;
  run;
  
  proc sort data=_out1;
    by dp;
  run;
  
  proc transpose data=_out1 out=Tip2out2(drop=_:) prefix=trt1_;
    by dp;
    id da;
    idlabel da;
    var probt;
  run;
  
  * Print into a table;
  proc format;
    value back low-0.05 ='green'
               >0.05-high ='lightred';
  run;

  * Output to table;
  ods html style=journal; 
  proc report data=Tip2out2 split="~";
    columns dp ("Active" trt1:);

    define dp /order order=internal descending "Placebo";
    define trt1: /display style(column)={background=back.};
  run;  
proc datasets;
  delete _tip:;
quit;
  
%mend tip2;
%*tip2;


**********************************************************************************************;
* Format Analyses for output;
**********************************************************************************************;
data Prim(keep = expestimate probz lowerexp upperexp trt _trt rename=(expestimate=estimate probz=prob lowerexp=lowcl upperexp=uppcl));
  set NBPMI_LSMP(drop=probz) NBPMI_LSMDP ;
run;

/* data MAR_sing1(keep = expestimate probz lowerexp upperexp trt _trt rename=(expestimate=estimate probz=prob lowerexp=lowcl upperexp=uppcl)); */
/*   set NBPMAR_LSMP(drop=probz) NBPMAR_LSMDP ; */
/* run; */

data comb(keep=estimate trt _trt prob probt lowcl uppcl Analysis);
  length Analysis $60;
  set prim(in=a) /*MAR_sing1(in=b)*/ negbinmi_mar_out(in=c);
  if a then Analysis = "Primary analysis";
/*   if b then Analysis = "Sensitivity analysis 1: MAR Single imputation using model based estimates"; */
  if c then Analysis = "Sensitivity analysis 3: MAR Multiple imputation";
run;
  
/*   ods html style=journal;  */
title "Primary analysis and Sensitivity analysis 3 - MAR MI";
proc print data=comb;run;

title "Sensitivity analysis 1 and 2 - MAR, single imputation, model based pre and post, Tipping point ";
  proc report data=TPA1_out2 split="~";
    columns dp ("Active" trt1:);
    define dp /order order=internal descending "Placebo";
    define trt1: /display style(column)={background=back.};
  run;


title "MF - Sensitivity analysis 3 (upper left cell)  and further Tipping point - MAR, Multiple imputation ";
  proc report data=out2 split="~";
    columns dp ("Active" trt1:);
    define dp /order order=internal descending "Placebo";
    define trt1: /display style(column)={background=back.};
  run;  

title "Tipping point, jump to reference, Multiple imputation";
  proc report data=Tip2out2 split="~";
    columns dp ("Active" trt1:);
    define dp /order order=internal descending "Placebo";
    define trt1: /display style(column)={background=back.};
  run;  
title;


ods word close;




/* title "Sensitivity analysis 1a and 2a Tipping point MAR, single imputation, observed pre and model based and post";   */
/*     proc report data=TPA2_out2 split="~"; */
/*     columns dp ("Active" trt1:); */
/*     define dp /order order=internal descending "Placebo"; */
/*     define trt1: /display style(column)={background=back.}; */
/*   run; */
/*  */
/* title "Sensitivity analysis 1b and 2b Tipping point MAR, single imputation, observed rate pre and post";   */
/*     proc report data=TPA3_out2 split="~"; */
/*     columns dp ("Active" trt1:); */
/*     define dp /order order=internal descending "Placebo"; */
/*     define trt1: /display style(column)={background=back.}; */
/*   run; */
  
