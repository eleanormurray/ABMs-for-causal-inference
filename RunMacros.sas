
libname abm '/ABM';

%include 'macros.sas';
%include '/gformula/June-2015/gformula.sas';
%include 'rcspline.sas';
%include 'lst.sas';

options minoperator mindelimiter=',';
options nonotes nocenter ;



/**Main analyses**/

%let start = 1;
%let numruns = 0;

%let name =null.;
%let psi1 = 0;

%let nameA =  out.rtf;
%let nameB =  log.rtf;
%let nameC  = params.rtf;

%let logdest	= &name.&nameB;
%let paramdest	= &name.&nameC;

data simparams_test;
 int1U =25;
 g2 = 0;
 g3 = 0;
 g4 = 0.675;
 g5 = 0;
 g6 = 0;
 a0 = 0.1;
 a1 = 0.3;
 a2 = 0.25;
 a3 = 0.10;
 kappa = -11.2;
 tau =0.38 ;
 pCD0 = 0.20;
run;

proc print data = simparams_test;
run;


%let newsim = %sysfunc(open(simparams_test));
%let NObs3 = %sysfunc(attrn(&newsim,NOBS));
%syscall set(newsim);

%macro factorial(startrun = &start, num = &numruns, name1 = , psi1 = , nsamples = 500);

	%let endrun = %sysevalf(&startrun + &num);

	%do blockrun = &startrun %to &endrun;

		%let set99 = %sysfunc(fetchobs(&newsim,&blockrun));

		%let int1U_new = %sysevalf(&int1U);
		%let g4_new    = %sysevalf(&g4);
		%let kappa_new = %sysevalf(&kappa);
		%let tau_new   = %sysevalf(&tau);
		%let a0_new    = %sysevalf(&a0);
		%let a1_new    = %sysevalf(&a1);
		%let a2_new    = %sysevalf(&a2);
		%let a3_new    = %sysevalf(&a3);
		
		%let pCD0_new = %sysevalf(&pCD0);

		/*%let name1 = Benef;*/
		%let name2 = &name1.&blockrun;
		%let outputdest = &name2.&nameA;

		%final_results(outputdest="&outputdest", logdest= "&logdest", paramdest = "&paramdest",
			subjects =100000, nsamples = &nsamples, psi1 = &psi1, int1lam = 0.010 , int2lam=0.1, int3lam = 0.001, 
			int1U = &int1U_new, g4 = &g4_new, kappa = &kappa_new, tau = &tau_new, 
			a0 = &a0_new, a1 = &a1_new, a2 = &a2_new, a3 = &a3_new, contU = 1, pCD0 = &pCD0_new);

	%end;

	%let set299=%sysfunc(close(&newsim));
%mend;

%factorial(startrun = &start, num = &numruns, name1 = &name, psi1 = &psi1);
