**********************************************************************************************;
*** Sensitivity analysis 1 and 2 - MAR and Tipping point analysis - single imputation ***;
**********************************************************************************************;


%macro TPA( 
            assump= 1 /* assumption for the MAR imputation: 
                            1 = model-based estimate of rate used pre and post withdrawal
                            2 = observed rate used for pre-withdrawal, model-based used for post withdrawal
                            3 = observed rate used for both pre and post withdrawal */  
          );
          
          
********************************************************************************;
***** First transpose parameter estimates to merge onto original data;
*******************************************************************************;
proc transpose data=parms(keep=parameter estimate where=(estimate ne 0)) out=parmtest;
  id parameter;
run;

***** Merge. Then calculate linear predictor (log of rate), add maxtime of observation;
proc sql;
  create table MAR_out1 as
  select a.*, b.*, 
         intercept + Btrt*trt + Breg*region  as pred ,                            /* log of the rate*/
         (case when peot=0 then ttie else 52 end) as maxtime,                     /* update time for dropouts to 52 weeks*/
         log(calculated maxtime) as ln_maxtime                                   /* log of time for offset*/
  from dat1 as a full join parmtest(rename=(trt=Btrt region=Breg)) as b on a.subjid ne b.intercept
  ;
  quit;



*******************************************************************************;
***** Cycle through to impute outcome based on assumption and apply delta to the post-witdrawal rate, then fit each model for delta;
*******************************************************************************;
  %let cnt=0;
  %do dp=0 %to 6; * Delta placebo - multiply by -0.25 later to get the correct increments;
    %do da=0 %to 6; * Delta active - multiply by 0.25 later to get the correct increments;
      %let cnt=%eval(&cnt.+1);
      %put cnt=&cnt. dp=&dp. da=&da.;
      
      data _MAR_out1;
        set MAR_out1;
        if trt=0 then delta=exp(&dp.*-0.25);
        else if trt=1 then delta=exp(&da.*0.25);
        
        *For those who withdraw, calc number of events before and after withdrawal ;
  
        if peot=1 then do;
          *calc number of events pre-withdrawal - model or observed;
            b1t1o= resp1 ;
            b1t1m= exp(pred)*ttie;
            
            *Calc number of events post-withdrawal - model or observed;
            b2t2m= exp(pred)*(maxtime-ttie) * delta;
            b2t2o= resp1/ttie *(maxtime-ttie) * delta;
            
            
           *impute number of events for each assumption ;
           resp_a1 = b1t1m + b2t2m;
           resp_a2 = b1t1o + b2t2m;
           resp_a3 = b1t1o + b2t2o;
        end;    
            
       *generate outcome variable ;
       if peot ne 1 then resp_a = resp1;
        *set the imputed outcome based on chosen assumption;
        else             
          %if &assump = 1 %then %do;
           resp_a = resp_a1;
          %end;
            %else %if &assump = 2 %then %do;
              resp_a = resp_a2;
            %end;
              %else %if &assump = 3 %then %do;
               resp_a = resp_a3;
              %end;   
      
      run;

    **********************;
    ***** Fit model with imputed outcome, offset is the maximum time - 52 weeks;

      ODS RESULTS OFF;
      ODS select none;
      option nonotes;
      proc genmod data=_MAR_out1 ;
        class trt region / ref=first;
        model resp_a = trt region / dist=negbin offset=ln_maxtime ;
        lsmeans trt / diff ilink cl exp ;
        ods output  diffs=_NBPTPA_LSMDP_&cnt ;
      run; 
      option notes;
      ods Results on;
      ODS select all;


        data _TPA&assump._&cnt;*(keep=probz dp da) ;
          set _NBPTPA_LSMDP_&cnt ;
          dp=&dp.*-0.25;
          da=&da.*0.25;
        run;
     %end;
  %end;

*******************************************************************************;
***** Format output;
*******************************************************************************;
  
    data _out1(keep=probz est dp da);
    set _TPA&assump._:;
        est=put(expestimate, 5.2)||" ("||put(upperexp-lowerexp,5.2)||")";
  run;
  
  proc sort data=_out1;
    by dp;
  run;
  
  proc transpose data=_out1 out=TPA&assump._outp(drop=_:) prefix=trt1_;
    by dp;
    id da;
    idlabel da;
    var probz;
  run;
  
  * Print into a table;
  proc format;
    value back low-0.05 ='green'
               >0.05-high ='lightred';
  run;

  * Output to table;
  ods html style=journal; 
  proc report data=TPA&assump._outp split="~";
    columns dp ("Active" trt1:);

    define dp /order order=internal descending "Placebo";
    define trt1: /display style(column)={background=back.};
  run;  
  
    * Table of estimates;
  proc transpose data=_out1 out=TPA&assump.oute(drop=_:) prefix=trt1_;
    by dp;
    id da;
    idlabel da;
    var est;
  run;
  
  * Output to table;
  ods html style=journal; 
  proc report data=TPA&assump.oute split="~";
    columns dp ("Active" trt1:);

    define dp /order order=internal descending "Placebo";
    define trt1: /display ;
  run;  

  
  proc datasets noprint;
    *delete _NBPTPA: _tpa: _out: _mar_out1;
  quit;
  
%mend TPA;


