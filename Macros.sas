/******************************************************************************/
/**This code contains the macros for: *****************************************/
/**%SIMULATE: creates simulated cohort studies*********************************/
/**%ABM: runs an ABM taking parameter inputs from any specified source*********/
/**%METHOD_COMPARISON: runs %simulate, %abm, and %gformula under given inputs**/
/**%SEEDGEN: generates matrix of random number seeds***************************/
/**%LABELS2: creates labels for ABM results************************************/
/**%PRINTRESULTS: labels and outputs results of each method********************/
/**%CREATEPARAMS: creates input parameter datasets for creation of the cohorts*/
/**%FINAL_RESULTS: combines macros into a single run of the simulation study***/
/******************************************************************************/



/******************************************************************************/
/* model: logit(L) = sgamma0 + sgamma1*U +sgamma2*L_1 + sgamma3*L_2 +sgamma4*A_1 + sgamma5*A_2*/
/* model: logit(A) = salpha0 + salpha1*L + salpha2*L_1+ salpha3*L_2 + salpha4*A_1 + salpha5*A_2 */
/******************************************************************************/
/*Macro SIMULATE: Simulate data under a SNAFTM with exponential failure times*/
%macro simulate(
   subjects= ,     /* number of subjects in the simulated data sets */
   tpoints=  ,     /* maximum number of measurements for each subject in the simulated data */           
   psi1= ,         /* true psi = AFT parameter for causal effect of treatment on survival */
	
   sgamma0 = , sgamma1 = , sgamma2 = , sgamma3 = , sgamma4 = , sgamma5 = , sgamma6 = ,	
   kappa = , tau = , 

   salpha0 = , salpha1 = , salpha2 = , salpha3 = , 


   cut= , 	/*t for int(t) cutpoint, for use with exact algortihm*/
   lam= ,	/*scale parameter for weibull*/
   shape=1,	/*shape parameter for weibull; set to 1 for T~exp(lam)*/

   dataout = Simdata_all,
   datasum = Simdata_surv, 
   contU = 0, 
	pCD0=
);

   data &dataout;
      length  A A_l1 A_l2 CD CD_l1 CD_l2 Y /*Y_l1*/ T 3;

        sT0 = %sysevalf(&sT0) ;
        sCD = %sysevalf(&sCD) ;
        sA = %sysevalf(&sA) ;
        sY = %sysevalf(&sY);

        %let tpoints_sim = %sysevalf(&tpoints + 2); /*simulate two additional time points to generate pre-baseline CD4. Times 1 and 2 are thus pre-baseline and have A = 0 and Y=0*/
   
        cut = %sysevalf(&cut);

	do id=1 to %sysevalf(&subjects); 
		call ranexp(sT0,expvar); 
               	T0=(expvar/(&lam**&shape))**(1/&shape); /*T0 ~exp(lambda)*/
		
		if &contU = 0 then do;			
			if T0 < cut then IT0=1; /*U = 1 if T0 < x*/
			else IT0=0;    	       /*U = 0 if T0>=x*/			
		end;
		else if &contU = 1 then do;
			IT0=T0; 	
		end;
			
		maxT= &tpoints_sim;
		
		initA = .;
				
	        array CDarray CDp2 CDp1 CD1-CD&tpoints;
        	array pCDarray pCDp2 pCDp1 pCD1-pCD&tpoints;
		array Aarray Ap2 Ap1 A1-A&tpoints;
		array pAarray pAp2 pAp1 pA1-pA&tpoints;
		array Yarray Yp2 Yp1 Y1-Y&tpoints;
		array Ymarray Ymp2 Ymp1 Ym1-Ym&tpoints;      
		array tarray tp2 tp1 t1-t&tpoints;

		time = 1.0;
	

 	do j=1 to &tpoints_sim;

		if j = maxT +1 then leave;
		
		if j = 1 then do;	
			A_l1 = 0; A_l2 = 0; CD_l1 = 0; CD_l2 = 0;
		
		end;
		
		else if j = 2 then do;
			A_l1 = Aarray(j-1); A_l2 = 0; CD_l1 = CDarray(j-1); CD_l2 = 0;

		end;
		
		else do;
			A_l1 = Aarray(j-1); A_l2 = Aarray(j-2); CD_l1 = CDarray(j-1); CD_l2 = CDarray(j-2);

		end;


		if j le 2 then do; 

			pCDarray(j) = 1;
			CDarray(j) = 1;
		end;
                     
		else if j = 3 then do;

			pCDarray(j) = 1-&pCD0;
          			call ranbin(sCD,1,pCDarray(j),CDarray(j));                 
		end;
		else if j > 3 then do;

			pCDarray(j) = (&sgamma4)*A_l1 + (1-&sgamma4)*(exp(&sgamma1 + &kappa*(log(IT0 + &tau))))/(1+(exp(&sgamma1 + &kappa*(log(IT0 + &tau)))));
	
			if pCDarray(j)=1 or pCDarray(j) = 0 then CDarray(j) = pCDarray(j);
 			else do;                      
          			call ranbin(sCD,1,pCDarray(j),CDarray(j));                 
	 		end; 	
		end;

		CD = CDarray(j);
		
		if A_l1 = 1 then do; 
			Aarray(j) = 1; 
		end;

		else do;
			if j in (1,2) then do;
				Aarray(j) = 0; pAarray(j) = 0;
			end;
			else do;
			logitA = &salpha0 + &salpha1*(1-CD) + &salpha2*(1-CD_l1) +&salpha3*(1-CD_l2);
			pAarray(j) = 1/(1+exp(-logitA)); 
			call ranbin(sA,1,pAarray(j),Aarray(j)); 
			end;
		end;
		
	    	A = Aarray(j);
			

		if j = maxT then do;

			Tcounter = 0;
			Tholder = 0;
			TAholder = exp(&Psi1*Aarray(1));
			
			do while (Tholder + TAholder <= IT0 and (Tcounter+1) < maxT);

				Tcounter + 1;
				Tholder = Tholder + TAholder;
				TAholder = exp(&Psi1*Aarray(Tcounter+1));
			end;
			Tfunc = Tcounter +(1/TAholder)*(IT0-Tholder);
		end;
	end; /*end j to &tpoints*/


	do tloop =3 to &tpoints_sim; /*baseline is at tpoint = 3*/
		tpoint = tloop - 2;				  
        	tpoint2 = tpoint - 1;
                         
		CD  = CDarray(tloop);
		pCD = pCDarray(tloop);

		CD_0 = CDarray(3);
		CD_l1 = CDarray(tloop - 1);
		CD_l2 = CDarray(tloop - 2);
					 
		A_0 = Aarray(3);         
		A  = Aarray(tloop);
		pA = pAarray(tloop);
		A_l1 = Aarray(tloop - 1);
		A_l2 = Aarray(tloop - 2);

		T = Tfunc - 3;
		
		if tpoint < Tfunc <= tpoint +1 and tpoint < maxT then Y = 1;
		else Y = 0;
                                       
		output;

	end; /* end of tloop loop */
   end;  /* end of id loop */
     
       	keep id  tpoint tpoint2 time T0 IT0  Y A A_l1 A_l2 A_0 CD CD_0 CD_l1 CD_l2 pA T maxT pCD;
   run;	
   	
   	data &dataout;
		set &dataout;
		if tpoint ne . and y ne .;
	run;

	proc means data = &dataout;
	var T0 IT0;
	run;
		
	proc sort data=&dataout;
		by tpoint;
	data temp1;
		set &dataout;
		where tpoint=1;
	proc means data=temp1 /*noprint*/;
  		var A CD;
		output out=temp  mean=pA_0 pCD_0 ;
	run;

	proc means data = &dataout;
		var A CD Y ;
	run;
		
	/*Create dataset with overall survival over time*/
	proc sort data=&dataout;
		by tpoint;
	run;

	proc freq data=&dataout  /*noprint*/;
		tables tpoint*CD tpoint*A /nopercent nocol;
	run;

	
	proc freq data=&dataout  /*noprint*/;
		tables tpoint*Y / out=counts nopercent nocol;
	run;

	proc transpose data=counts out=output (drop=_name_ _label_) prefix=Y_;
		by tpoint;
  		id Y;
		var COUNT;
	run;
	data output;
		set output;
		id = 1;
	proc transpose data=output out=output2 prefix=Y_0_;
		by id;
		id tpoint;
		var Y_0 ;
	proc transpose data=output out=output3 prefix=Y_1_;
		by id;
		id tpoint;
		var Y_1 ;
	data wide;
   		merge  output2(drop=_name_) output3(drop=_name_);
    		by id;
	data &datasum;
		set wide;
		s0 = 1;
		array sY(&tpoints) sY1 - sY&tpoints;
		array pY(&tpoints) pY1 - pY&tpoints;
		array Y0(&tpoints) Y_0_1 - Y_0_&tpoints;
		array Y1(&tpoints) Y_1_1 - Y_1_&tpoints;

		do i=1 to &tpoints;
			if i = 1 then do;
				pY(i) = ( Y1(i))/(Y0(i)+Y1(i));
				sY(i) = s0*(1-pY(i));
			end;
			else do;
				pY(i) = ( Y1(i))/(Y0(i)+Y1(i));
				sY(i) = sY(i-1)*(1-pY(i));
			end;
		end;
		               
    		do i=1 to  &tpoints ; 
			tpoint = i;
			p_Y  =	pY(i);
			s_tpoint  =   sY(i);
    		output;
		end;

		keep tpoint p_Y s_tpoint;     
	run;

	/*dataset for g-formula (time from 0 to 11)*/
	data &dataout;
		set &dataout;
		tpoint  = tpoint - 1;
		tpoint2 = tpoint2 - 1;	
		keep  id  tpoint tpoint2 Y  A A_l1 A_l2 A_0 CD CD_l1 CD_l2 CD_0 IT0 T maxT;
	run;


	proc datasets library=work nolist;
		delete counts output output2 output3 wide temp1;
	quit;


%mend simulate;


/*******************************************************************************************************************************/
/******************************************************************************************************************************/
/**Macro METHOD_COMPARISON: simulates data under AFT, runs g-formula on data, runs agent-based model with specified parameters*/
%macro method_comparison(
	tpoints = , 			/*number of time points for simulation*/
	subjects = , 			/*sample size for g-formula and ABMs*/
	subjectsSim = ,			/*sample size of simulated cohort*/
	nsamples = , 			/*number of bootstrap samples for g-formula*/

    lam = , cut= , contU = 0,					/*parameters for prevalence of U in simulated data*/
	psi1 = 0,						/*true effect of A on Y; psi1 = 0 for null effect*/

	
/* model: logit(L) = sgamma0 + sgamma1*U +sgamma2*L_1 + sgamma3*L_2 +sgamma4*A_1 + sgamma5*A_2*/
   sgamma0 = ,sgamma1 = , sgamma2 = , sgamma3 = , sgamma4 = , sgamma5 = , sgamma6 = ,	
	kappa = , tau = ,


/* model: logit(A) = salpha0 + salpha1*L + salpha2*L_1+ salpha3*L_2 + salpha4*A_1 + salpha5*A_2 */
   salpha0 = , salpha1 = , salpha2 = , salpha3 =  , 
	

	refint = ,
	
	pCD0 = ,

	betadata= betas,
	abm_paramdata = , 			/*name of external dataset containing parameter estimate(s) for use in ABM simulation */
	scenario = ,				/*scenario number for record keeping*/
	paramdest = 				/*output file for beta parameters*/

);

/*Record start time and date*/
   %let timenow=%sysfunc(time(), time.);
   %let datenow=%sysfunc(date(), date9.);
   %put Beginning Simulation;
   %put Start time is &datenow / &timenow ;
   %put ;

   /*set up interventions for g-formula*/
		%let interv1  =
		     intno     = 1, 
		     intlabel  = 'Always treat',
		     nintvar   = 1,
		     intvar1   = A,
		     inttype1  = 1,
		     intvalue1 = 1,
      	  	     intpr1    = 1,
		     inttimes1 = 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20;
		%let interv2  =
		     intno     = 2, 
		     intlabel  = 'Never treat',
		     nintvar   = 1,
		     intvar1   = A,
		     inttype1  = 1,
		     intvalue1 = 0,
		     intpr1    = 1,
		     inttimes1 = 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20;

	/*simulate random seeds*/
	%let lstream = %sysevalf(&subjects)*%sysevalf(&tpoints);
	%let fseed =%eval(123 + %sysevalf(&scenario));

	%seedgen(fseed=&fseed,                 /*initial seed for seed generator*/
			lstream= &lstream ,    /* number of seeds needed for one stream  */
			seeds_per_line= 10 ,   /*number of seeds per individual*/
			nseeds = 1
			);
 
	%let seeds=%sysfunc(open(work.Seeds));
	%let NObs=%sysfunc(attrn(&seeds,NOBS));
	%syscall set(seeds);

	%let rc=%sysfunc(fetchobs(&seeds, 1));
	
	/*Assign seed values for data simulation and ABM*/
		     	 %let sT0 =   %sysevalf(&seed1);
       		   	 %let sCD =  %sysevalf(&seed2);
       	 	 	 %let sA =    %sysevalf(&seed3);
       	  		 %let sY =    %sysevalf(&seed4);

			 %let sABM1 = %sysevalf(&seed5);
			 %let sABM2   = %sysevalf(&seed6);
	         	 %let sABM3   = %sysevalf(&seed7);
	         	 %let sABM4   = %sysevalf(&seed8);
	         	 %let sABM5   = %sysevalf(&seed9);
		 
			%put Simulating initial data;
			%put ;

	  /*simulate dataset, store in Simdata_all*/
	        %simulate(
				subjects =&subjectsSim,
				tpoints = &tpoints,
				sgamma0 = &sgamma0, sgamma1 = &sgamma1, sgamma2 = &sgamma2, sgamma3 = &sgamma3, sgamma4 = &sgamma4, sgamma5 = &sgamma5, sgamma6 = &sgamma6,
				salpha0 = &salpha0, salpha1 = &salpha1, salpha2 = &salpha2, salpha3 = &salpha3,  
				kappa = &kappa, tau = &tau, 
				psi1=&psi1,
				lam=&lam,
				cut=&cut, 
				contU = &contU, pCD0 = &pCD0
				);

			data _null_;
				set temp;
				call symput("pA0", pA_0);
				call symput("pCD0",pCD_0);
			run;


			data Simdata_combined;
				set temp;
				keep pA_0 pCD_0;
			run;
			data results_obs;
				set Simdata_surv;
			run;

			proc sort data= Simdata_all;
				by id;
			run;

			data first;
				set Simdata_all end = _end_;
				by id;
				if first.id then output;
			run;

			proc printto print = &outputdest;
			run;
			%if &contU = 0 %then %do;
				proc freq data=first;
					tables IT0;
					title "Prevalence of U (cut = &cut), Scenario: &scen";					
				run;
			%end;
			%else %if &contU = 1 %then %do;
				proc means data = first mean std min max p5 p95;
					var IT0;
					title "Distribution of continuous U, Scenario: &scen";
				run;
			%end;
			proc printto;
			run;

			%put Data simulation complete;
			%put ;

			
			/*Run G-formula on Simulated Dataset*/
			%put Running G-formula;
			%put ;


			/*g-formula without confounding adjustment*/
			title 'GFORMULA SAMPLE';
			%gformula(
			    data= Simdata_all,
			    id=id,
			    time=tpoint,
			    timeptype =conspl,
			    timeknots = 1 3 6 9 12,
			    timeinc = 1,
			    timepoints = &tpoints, 
				
			    outc=Y,
			    outctype = binsurv,

			    ncov=1,
			    cov1 = A,      cov1otype = 1, cov1ptype = lag1bin,

			    seed=  %eval(1234+%sysevalf(&scenario)), 

			    nsimul= ,
			    nsamples = &nsamples, 

			    numint=2,
			    intervname=interv,

				betadata= betas_conf,
			 	survdata=work.survdth,
				savelib=work,
				covmeandata=covmeandth,
				observed_surv=obsurvdth,
				refint = &refint, /*never treat as reference*/
				print_stats=0, 
				/*pringlogstats=0,
				outputs=no*/
				resultsdata=results_gform_conf
			 );


			 /*Combine g-formula results across runs*/
			data results_gform_conf;
				set results_gform_conf;
				keep int int2 pD pD_llim95 pD_ulim95 rd rd_llim95 rd_ulim95 rr rr_llim95 rr_ulim95 pD_mean pD_std intervened averinterv;
			run;	

			
			/*G-formula with confounding adjustment*/
			title 'GFORMULA SAMPLE';
			%gformula(
			    data= Simdata_all,
			    id=id,
			    time=tpoint,
			    timeptype =conspl,
			    timeknots = 1 3 6 9 12,
			    timeinc = 1,
			    timepoints = &tpoints, 
				
			    outc=Y,
			    outctype = binsurv,
			    outcinteract= 1*2,

			    ncov=2,
				cov1 = CD,     cov1otype = 1, cov1ptype = lag2bin,
				cov2 = A,      cov2otype = 1, cov2ptype = lag2bin,

			    seed=  %eval(1234+%sysevalf(&scenario)), 

			    nsimul= ,
			    nsamples = &nsamples, 

			    numint=2,
			    intervname=interv,

				betadata= betas,
			 	survdata=work.survdth,
				savelib=work,
				covmeandata=covmeandth,
				observed_surv=obsurvdth,
				refint = &refint, /*never treat as reference*/
				print_stats=0, 
				/*pringlogstats=0,
				outputs=no*/
				resultsdata=results_gform
			 );

			 /*Combine g-formula results across runs*/
			data results_gform;
				set results_gform;
				keep int int2 pD pD_llim95 pD_ulim95 rd rd_llim95 rd_ulim95 rr rr_llim95 rr_ulim95 pD_mean pD_std intervened averinterv;
				run;	


			/*output beta parameter estimates*/
			%let n = betas_;
			%let x = &n.&scenario;

				data &x;
					set betas /*(drop=_sample_)*/;
				run;
				
			proc datasets library=work nolist;
				delete betas _beta_   _paramdata_  _simuldata_ _inputd_ _ref_ fin finfin; 
			quit;
			
			proc printto print = &paramdest;
			run;
			proc print data = &x (obs = 1);
			title "G-formula parameter estimates";
			run;
			proc printto; 
			run;
	
			%put G-formula complete;
			%put ;

			
			%put Running ABM;
			%put ;


		/*Run ABM*/
	
		/*ABM type 1 only for base case scenario*/
		%if &scenario = 1 %then %do;
			

			%gformula(
			    data= Simdata_all,
			    id=id,
			    time=tpoint,
			    timeptype =conspl,
			    timeknots = 1 3 6 9 12,
			    timeinc = 1,
			    timepoints = &tpoints, 
				
			    outc=Y,
			    outctype = binsurv,
			    outcinteract= 1*2,
		    
			    ncov=2,
				cov1 = CD,     cov1otype = 1, cov1ptype = lag2bin,
				cov2 = A,      cov2otype = 1, cov2ptype = lag2bin,

			    seed=  %eval(1234+%sysevalf(&scenario)), 

			    nsimul= ,
			    nsamples = 0, 

			    numint=2,
			    intervname=interv,

				betadata= work.&x,
				usebetadata=1,
			 	survdata=work.survdth,
				savelib=work,
				covmeandata=covmeandth,
				observed_surv=obsurvdth,
				refint = &refint, 
				print_stats=0, 
				resultsdata= ABM_summary
			 );

			data results_ABM;
				set ABM_summary;
				keep int int2 pD pD_llim95 pD_ulim95 rd rd_llim95 rd_ulim95 rr rr_llim95 rr_ulim95 pD_mean pD_std intervened averinterv;

			run;
			proc datasets library=work nolist;
				delete ABM_summary; 
			quit;

			proc printto print = &paramdest;
			run;
			proc print data = &x (obs = 1);
			title "ABM1 parameter estimates";
			run;
			proc printto; 
			run;



		%end;

		/*ABM types 1,2, 3 for other scenarios*/
		%else %if &scenario > 1 %then %do;
		

			%put Running ABM - all params current pop;
			%put ;

			%gformula(
			    data= Simdata_all,
			    id=id,
			    time=tpoint,
			    timeptype =conspl,
			    timeknots = 1 3 6 9 12,
			    timeinc = 1,
			    timepoints = &tpoints, 
				
			    outc=Y,
			    outctype = binsurv,
			    outcinteract= 1*2,
			    
			    ncov=2,
				cov1 = CD,     cov1otype = 1, cov1ptype = lag2bin, 
				cov2 = A,      cov2otype = 1, cov2ptype = lag2bin, 

			    seed=  %eval(1234+%sysevalf(&scenario)), 

			    nsimul= ,
			    nsamples = 0, 

			    numint=2,
			    intervname=interv,

				betadata= work.&x,
				usebetadata=1,
			 	survdata=work.survdth,
				savelib=work,
				covmeandata=covmeandth,
				observed_surv=obsurvdth,
				refint = &refint, 
				print_stats=0, 
				resultsdata= ABM_summary_a
			 );

				data results_ABM_a;
					set ABM_summary_a;
				keep int int2 pD pD_llim95 pD_ulim95 rd rd_llim95 rd_ulim95 rr rr_llim95 rr_ulim95 pD_mean pD_std intervened averinterv;

				run;
				proc datasets library=work nolist;
					delete ABM_summary_a; 
				quit;

			proc printto print = &paramdest;
			run;
			proc print data = &x (obs = 1);
			title "ABM1 parameter estimates";
			run;
			proc printto; 
			run;
			%put Running ABM - all params base case;
			%put ;

			%gformula(
			    data= Simdata_all,
			    id=id,
			    time=tpoint,
			    timeptype =conspl,
			    timeknots = 1 3 6 9 12,
			    timeinc = 1,
			    timepoints = &tpoints, 
				
			    outc=Y,
			    outctype = binsurv,
			    outcinteract= 1*2,
	    
			    ncov=2,
				cov1 = CD,     cov1otype = 1, cov1ptype = lag2bin, 
				cov2 = A,      cov2otype = 1, cov2ptype = lag2bin, 

			    seed=  %eval(1234+%sysevalf(&scenario)), 

			    nsimul= ,
			    nsamples = 0, 

			    numint=2,
			    intervname=interv,

				betadata=  work.&&abm_paramdata,
				usebetadata=1,
			 	survdata=work.survdth,
				savelib=work,
				covmeandata=covmeandth,
				observed_surv=obsurvdth,
				refint = &refint, 
				print_stats=0, 
				resultsdata= ABM_summary_b
			 );
			
				data results_ABM_b;
					set ABM_summary_b;
				keep int int2 pD pD_llim95 pD_ulim95 rd rd_llim95 rd_ulim95 rr rr_llim95 rr_ulim95 pD_mean pD_std intervened averinterv;

				run;
				proc datasets library=work nolist;
					delete ABM_summary_b; 
				quit;

			proc printto print = &paramdest;
			run;
			proc print data = work.&&abm_paramdata (obs = 1);
			title "ABM2 parameter estimates";
			run;
			proc printto; 
			run;


			%put Running ABM - only past A params base case;
			%put ;


			data abm5_l_a;
				set work.&x;
				keep _sample_ bvar1_0-bvar1_8 bvar2_0 - bvar2_9 boutc0-boutc6 boutc8-boutc9;
			run;
			data abm5_y;	
				set work.&&abm_paramdata;
				keep _sample_ boutc7;
			run;

			data abm5;
				merge abm5_l_a abm5_y;
				by _sample_;
			run;

			%gformula(
			    data= Simdata_all,
			    id=id,
			    time=tpoint,
			    timeptype =conspl,
			    timeknots = 1 3 6 9 12,
			    timeinc = 1,
			    timepoints = &tpoints, 
				
			    outc=Y,
			    outctype = binsurv,
			    outcinteract= 1*2,
			    
			    ncov=2,
				cov1 = CD,     cov1otype = 1, cov1ptype = lag2bin,
				cov2 = A,      cov2otype = 1, cov2ptype = lag2bin,
			    seed=  %eval(1234+%sysevalf(&scenario)), 

			    nsimul= ,
			    nsamples = 0, 

			    numint=2,
			    intervname=interv,

				betadata= work.abm5,
				usebetadata=1,
			 	survdata=work.survdth,
				savelib=work,
				covmeandata=covmeandth,
				observed_surv=obsurvdth,
				refint = &refint, 
				print_stats=0, 
				resultsdata= ABM_summary_e
			 );
			
				data results_ABM_e;
					set ABM_summary_e;
				keep int int2 pD pD_llim95 pD_ulim95 rd rd_llim95 rd_ulim95 rr rr_llim95 rr_ulim95 pD_mean pD_std intervened averinterv;

				run;
				proc datasets library=work nolist;
					delete ABM_summary_e; 
				quit;

			proc printto print = &paramdest;
			run;
			proc print data = work.abm5 (obs = 1);
			title "ABM3 parameter estimates";
			run;
			proc printto; 
			run;


			%put Running ABM - only past A params cohort of interest;
			%put ;

			data abm4_l_a;
				set work.&&abm_paramdata;
				keep _sample_ bvar1_0-bvar1_8 bvar2_0 - bvar2_9 boutc0-boutc6 boutc8-boutc9;
			run;
		
			data abm4_y;	
				set work.&x;
				keep _sample_ boutc7;
			run;

			data abm4;
				merge abm4_l_a abm4_y;
				by _sample_;
			run;

			%gformula(
			    data= Simdata_all,
			    id=id,
			    time=tpoint,
			    timeptype =conspl,
			    timeknots = 1 3 6 9 12,
			    timeinc = 1,
			    timepoints = &tpoints, 
				
			    outc=Y,
			    outctype = binsurv,
			    outcinteract= 1*2,

			    ncov=2,
				cov1 = CD,     cov1otype = 1, cov1ptype = lag2bin, 
				cov2 = A,      cov2otype = 1, cov2ptype = lag2bin, 

			    seed=  %eval(1234+%sysevalf(&scenario)), 

			    nsimul= ,
			    nsamples = 0, 

			    numint=2,
			    intervname=interv,

				betadata= work.abm4,
				usebetadata=1,
			 	survdata=work.survdth,
				savelib=work,
				covmeandata=covmeandth,
				observed_surv=obsurvdth,
				refint = &refint, 
				print_stats=0, 
				resultsdata= ABM_summary_d
			 );
			
				data results_ABM_d;
					set ABM_summary_d;
				keep int int2 pD pD_llim95 pD_ulim95 rd rd_llim95 rd_ulim95 rr rr_llim95 rr_ulim95 pD_mean pD_std intervened averinterv;

				run;
				proc datasets library=work nolist;
					delete ABM_summary_d; 
				quit;


			proc printto print = &paramdest;
			run;
			proc print data =  work.abm4 (obs = 1);
			title "ABM4 parameter estimates";
			run;
			proc printto; 
			run;

			proc datasets library=work nolist;
				delete &x; 
			quit;
		%end;				
	
		%put ABM complete;
		%put ;


	%let rc2=%sysfunc(close(&seeds));

	/*Print Final output*/

	%let nametemp10 = cuminc_sim;
	%let nametemp11 = final_OBScuminc;
	%let nametemp12 =  &nametemp10.&scenario;
	%let namefinal1 =  &nametemp11.&scenario;

	%let nametemp21 = final_gformula;
	%let namefinal2 =  &nametemp21.&scenario;


	%let nametemp21b = final_gformula_conf;
	%let namefinal2b =  &nametemp21b.&scenario;

	%if &scenario = 1 %then %do;

	%let nametemp31 = final_abmcuminc;
	%let namefinal3 =  &nametemp31.&scenario;

	%end;

	%else %if &scenario > 1 %then %do;
	
	%let nametemp31a = final_abm1cuminc;
	%let namefinal3a =  &nametemp31a.&scenario;

	%let nametemp31b = final_abm2cuminc;
	%let namefinal3b =  &nametemp31b.&scenario;
	
	%let nametemp31e = final_abm5cuminc;
	%let namefinal3e =  &nametemp31e.&scenario;	

	%let nametemp31d = final_abm4cuminc;
	%let namefinal3d =  &nametemp31d.&scenario;

	%end;

		/*Simulated Data*/
		data &nametemp12;
			set results_obs;
			if tpoint = symget('tpoints') then do;
				pY = (1-s_tpoint)*100;
				St = s_tpoint*100;
			end;
		run;

		proc univariate data=&nametemp12 noprint;
			var pY St;
			output out = &namefinal1
			mean =  pY_mean St_mean;
		data &namefinal1;
			set &namefinal1;
			%labels2;
	

		/*G-formula*/
		data &namefinal2;
			set results_gform;
		run;

		data &namefinal2b;
			set results_gform_conf;
		run;


		/*ABM*/

		%if &scenario = 1 %then %do;
		data &namefinal3;
			set results_ABM ;
		run;

		%end;

		%else %if &scenario > 1 %then %do;

			data &namefinal3a;
				set results_ABM_a ;
				%labels2;
			run;

			data &namefinal3b;
				set results_ABM_b ;
				%labels2;
			run;

			data &namefinal3e;
				set results_ABM_e ;
				%labels2;
			run;

			data &namefinal3d;
				set results_ABM_d ;
				%labels2;
			run;
		%end;

		

   %let timenow2=%sysfunc(time(), time.);
   %let datenow2=%sysfunc(date(), date9.);
   %put ;
   %put End time is &datenow2 / &timenow2 ;
   %put ;

%mend  method_comparison;


/******************************************************************************/
/******************************************************************************/
/*Macro SEEDGEN: Generate random seeds*/
%macro seedgen(fseed=,lstream=,seeds_per_line =,nseeds=);
%if %bquote(&seeds_per_line) = %then %let seeds_per_line = 1 ;
data seeds( keep = bsample  %do i = 1 %to %eval(&seeds_per_line); seed&i %end;);
	 retain seed &fseed;
	 bsample = 0 ;
	 do i=1 to min((&nseeds)*&lstream,2**31-1) by &lstream;
		      do index = 1 to %eval(&seeds_per_line);
			    
			        do j=1 to &lstream;       
			            call ranuni(seed,x);             
			        end;
			        
			        if i = 1  and index = 1 then do;
			             index = 2;
			             seed1 = &fseed ;
			        end;
			        
			        %do i = 1 %to %eval(&seeds_per_line);
			            if index = &i then seed&i = seed;
			        %end;
			end;
			bsample = bsample + 1;    
		 output;
	end;
 run;

proc print data = seeds;
run;

 %mend seedgen;

 
/******************************************************************************/
/******************************************************************************/
/*Macro LABELS2: format output for simulation and agent-based model*/
%macro labels2;
   	label int       = 'Intervention';
  	label int2      = 'Description';
   
   	label tpoint    ='Time';

   	label p_Y = 'Risk (%)';
	label pY_mean = 'Risk (%)';
   	label pY_std  = 'Risk (%), SD';
   	label pY_llim95 ='Lower limit 95% CI';
   	label pY_ulim95 ='Upper limit 95% CI';

	label rd = 'Risk Difference (%)';
	label rd_mean = 'Risk Difference (%)';
   	label rd_std = 'Risk Difference (%), SD';
	label rd_ulim95 = 'Upper limit 95% CI';
   	label rd_llim95 = 'Lower limit 95% CI';
%mend labels2;


/******************************************************************************/
/******************************************************************************/
/*Macro PRINTRESULTS: formats and prints results from specified number of scenarios*/
%macro printresults(
		scenarios = , basecase = , outputdest=
		);
%let varSim = pY_mean;
%let varGform = int int2 pd pd_llim95 pd_ulim95 rd RD_llim95 RD_ulim95 rr RR_llim95 RR_ulim95 intervened averinterv;
%let varABM =  int int2 pd  rd  rr  intervened averinterv ;

%do i= 1 %to &scenarios;
	%let scenario = &i;

	%if &i =2 %then %do;
		%let namesuff = High-risk; 
	%end;
	%else %if &i =3 %then %do;
		%let namesuff = Low-risk; 
	%end;

	%let nametemp11 = final_OBScuminc;
	%let namefinal1 =  &nametemp11.&scenario;

	%let nametemp21 = final_gformula;
	%let namefinal2 =  &nametemp21.&scenario;

	%let nametemp21b = final_gformula_conf;
	%let namefinal2b =  &nametemp21b.&scenario;


	%if &scenario = 1 %then %do;

		%let nametemp31 = final_abmcuminc;
		%let namefinal3 =  &nametemp31.&scenario;

	%end;

	%else %if &scenario > 1 %then %do;
	
		%let nametemp31a = final_abm1cuminc;
		%let namefinal3a =  &nametemp31a.&scenario;

		%let nametemp31b = final_abm2cuminc;
		%let namefinal3b =  &nametemp31b.&scenario;

		%let nametemp31e = final_abm5cuminc;
		%let namefinal3e =  &nametemp31e.&scenario;	

		%let nametemp31d = final_abm4cuminc;
		%let namefinal3d =  &nametemp31d.&scenario;	

	%end;


		%if &basecase = &scenario %then %do;
			proc printto print = &outputdest;
			run;
			proc print data=&namefinal1 label;
				var &varSIm;
				title "Summary of Simulated Data: Basecase, Scenario &scenario";
			run;
			proc printto print = &outputdest;
			run;
			proc print data= &namefinal2 label;
				var &varGform;
				title "Summary of G-formula Results: Basecase, Scenario &scenario";
			run;
			proc printto print = &outputdest;
			run;
			proc print data= &namefinal2b label;
				var &varGform;
				title "Summary of G-formula Results, without confounding adjustment: Basecase, Scenario &scenario";
			run;
			proc printto print = &outputdest;
			run;
			proc print data=&namefinal3 label;
				var &varABM;
				title "Summary of ABM type 1 Results: Basecase, Scenario &scenario";
				title2 "All parameters from cohort of interest";
			run;
			proc printto print = &outputdest;
			run;
		%end;
		
		%else %if &basecase ne &scenario %then %do;
			proc printto print = &outputdest;
			run;
			proc print data=&namefinal1 label;
				var &varSim;
				title "Summary of Simulated Data: &namesuff population, Scenario &scenario ";
			run;
			proc printto print = &outputdest;
			run;
			proc print data= &namefinal2 label;
				var &varGform;
				title "Summary of G-formula Results: &namesuff population, Scenario &scenario";
			run;
			proc printto print = &outputdest;
			run;
			proc print data= &namefinal2b label;
				var &varGform;
				title "Summary of G-formula Results, without confounding adjustment: Basecase, Scenario &scenario";
			run;
			proc printto print = &outputdest;
			run;
			proc print data=&namefinal3a label;
				var &varABM;
				title "Summary of ABM TYPE 1 Results: &namesuff population, Scenario &scenario";
				title2 "All parameters from cohort of interest";
			run;
			proc print data=&namefinal3b label;
				var &varABM;
				title "Summary of ABM TYPE 2 Results: &namesuff population, Scenario &scenario";
				title2 "All parameters from base-case cohort";
			run;
			proc print data=&namefinal3e label;
				var &varABM;
				title "Summary of ABM TYPE 3 Results: &namesuff population, Scenario &scenario";
				title2 "All past A parameters from base-case cohort";
			run;

			proc print data=&namefinal3d label;
				var &varABM;
				title "Summary of ABM TYPE 4 Results: &namesuff population, Scenario &scenario";
				title2 "All past A parameters from cohort of interest";
			run;
		%end;
%end;
%mend printresults;


/************************************************/
/*Run Method Comparison under range of scenarios*/
/***********************************************/


%macro final_results(outputdest="output.rtf", logdest="log.rtf", paramdest = "params.rtf", subjects = 0, nsamples =0 , 
			psi1 = 0, int1lam = 0, int2lam=0, int3lam =0,
			int1U =1.10, g2 = 0, g3 = 0, g4 = 0, g5 =0 , g6 = 0, kappa =0 , tau = 0, 
			a0 =0 , a1 =0 , a2 =0 , a3 = 0, cut = 0, contU =0, pCD0 =  );

proc printto log=&logdest;
run;

%let timenow=%sysfunc(time(), time.);
%let datenow=%sysfunc(date(), date9.);
%put Initiating Program;
%put Start time is &datenow / &timenow ;
%put ; 

/*Specify macro variables for use in all scenarios*/
%let  scenarios = 3;
%let  refint = 2; 
%let  tpoints = 12;

%createparams( int1U = &int1U, int1lam = &int1lam, int2lam=&int2lam, int3lam = &int3lam, scen = &scenarios, kappa = &kappa, tau = &tau, 
		g2 = &g2, g3 = &g3, g4 = &g4, g5 =&g5, g6 = &g6, a0 = &a0, a1 =&a1, a2 =&a2 , a3 = &a3);

%let simulinp = %sysfunc(open(work.simulinputs));
%let NObs3 = %sysfunc(attrn(&simulinp,NOBS));
%syscall set(simulinp);

%do scen = 1 %to &scenarios;
	%put &scen;
	%let set = %sysfunc(fetchobs(&simulinp,&scen));
		%let  sgamma0 = %sysevalf(&g0);
		%let  sgamma1 = %sysevalf(&g1);		/*U->L*/
		%let  sgamma2 = %sysevalf(&g2);
		%let  sgamma3 = %sysevalf(&g3);
		%let  sgamma4 = %sysevalf(&g4);
		%let  sgamma5 = %sysevalf(&g5);
		%let  sgamma6 = %sysevalf(&g6);
		%let  kappa   = %sysevalf(&kappa);
		%let  tau     = %sysevalf(&tau);

		%let  salpha0 = %sysevalf(&a0);
		%let  salpha1 = %sysevalf(&a1);
		%let  salpha2 = %sysevalf(&a2);
		%let  salpha3 = %sysevalf(&a3);

		%let  lam = %sysevalf(&lambda);		/*prev U*/

		%let abm_paramdata = Betas_1;

		%method_comparison(
			lam = &lam, cut= &cut, 		/*parameters for prevalence of U in simulated data*/
			contU = &contU, 		/*switch for binary or continuous U in simulated data*/
      			psi1 = &psi1,				
			abm_paramdata = &abm_paramdata,	 			
			scenario = &scen,				/*scenario number for record keeping*/

			tpoints = &tpoints, 		
			subjects = &subjects, 
			subjectsSim = &subjects,   
			nsamples = &nsamples, 			

			sgamma0 = &sgamma0, sgamma1 = &sgamma1, sgamma2 = &sgamma2, 
				sgamma3 =&sgamma3, sgamma4 = &sgamma4, sgamma5 =&sgamma5, sgamma6 = &sgamma6, 
				kappa = &kappa, tau= &tau,
			
			salpha0 = &salpha0, salpha1 = &salpha1 , salpha2 =  &salpha2, salpha3 =  &salpha3,

			refint = &refint,
			paramdest =  &paramdest, 
			pCD0 = &pCD0
			);

		
	proc printto log=&logdest;
	run;
	%put Scenario &scen complete;
	%put %sysfunc(time(), time.);
	proc printto log=&logdest;
	run;
		

%end;

%let set2=%sysfunc(close(&simulinp));

%printresults(basecase = 1, scenarios = &scenarios, outputdest=&outputdest);


proc printto log=&logdest;
run;

%let timenow2=%sysfunc(time(), time.);
%let datenow2=%sysfunc(date(), date9.);
%put Program Complete;
%put End time is &datenow2 / &timenow2 ;
%put ; 	

proc printto;
run;

%mend final_results;



 

/**********************************************************************************/
/****Input parameters for creating simulated population to dataset*/

%macro createparams(int1U = , int2U = , int3U = , int1lam = , int2lam=, int3lam =, scen = , 
			g2 = 4.00, g3 = 2.50, g4 = 0.1, g5 =0.1, g6 = 1, a0 = 0.33, a1 = 5.0, a2 =3.00 , a3 = 1.50, kappa = , tau = );
data simulinputs;
  do i = 1 to &scen;
	/*constant parameters*/
	pvalue = cdf('normal', ((200-588)/300), 0, 1);
	g0 = pvalue; /*sgamma0 - logit prob of CD4<200*: p(CD4<200) = 0.09 - Althoff 2007 (CEPAC ref) */

	g1 = &int1U;/*sgamma1 - Effect of U on L*/
	
	/*constant parameters*/
	g2 = &g2; 	/*sgamma2 - log RR low CD4 given low CD4 last month*/
	g3 = &g3; 	/*sgamma3 - log RR low CD4 given low CD4 2 months ago*/
	g4 = &g4; /*sgamma4 - log RR low CD4 given on treatment last month*/  
	g5 = &g5; /*sgamma5 - log RR low CD4 given on treatment 2 months ago*/ 
	g6 = &g6;
	kappa = &kappa;
	tau = &tau;

	a0 = log(&a0/(1-&a0));     /*salpha0 - logit prob of on treatment*/
	a1 = &a1; 	/*salpha1 - log RR for on treatment when CD4 <200 in current month*/
	a2 = &a2; 	/*salpha2 - log RR for on treatment when CD4 <200 in last month*/
	a3 = &a3;  /*salpha3 - log RR for on treatment when CD4 <200 in 2 months ago*/

	/*scenario dependent parameter: lambda*/
	if i = 1 then
		lambda = &int1lam; /*Prevalence of U*/
	else if i = 2 then 
		lambda = &int2lam;
	else if i =3 then 
		lambda = &int3lam;
	output;
  end;
run;

proc printto print = &outputdest;
run;

proc print data=simulinputs;
run;

proc printto;
run;

%mend createparams;
