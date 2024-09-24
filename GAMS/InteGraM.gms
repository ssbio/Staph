*************************************************************
************************InteGraM****************************
*************************************************************
*************************************************************
*************************************************************

$INLINECOM /*  */
$ONEMPTY

OPTIONS

	limrow = 1000
       optCR = 0
        optCA = 0
        iterlim = 100000
        decimals = 7
        reslim = 100000
        work = 5000000;

*********Defining Sets**************************************
SETS

	i					set of metabolites

$include "metabolites.txt"	

	j					set of reactions

$include "reactions.txt"

	SMP					set for concnetration sampling	/1*100/

	cs(i)    		set of metabolites constrained by thermodynamics
$include "s_metabolites.txt"

	spec(j)
$include "s_reactions.txt"

	off(j)
$include "turn_off_rxns.txt"

	aauprxns2(j)
$include "aa_uptake_rxns2.txt"

	fix(i)
$include "conc_fix.txt"

	mem(j)
$include "membrane_reactions.txt"

	mc(j)
$include "vector.txt"
;

*************************************************************

***********Defining Parameters*******************************
PARAMETERS

	S(i,j)					stoichiometric matrix

$include "sij.txt"

	v_max(j)				maximum flux of v(j)
	
$include "upper_bound.txt"

	v_min(j)				minimum flux of v(j)

$include "lower_bound.txt"


	delta_G_o(j)				maximum flux of v(j)

$include "deltaGo_s.txt"

	Cmax(i)					highest concentration

$include "cmax_s.txt"

	Cmin(i)					Lowest concentration

$include "cmin_s.txt"

	R					Gas constant
/0.008314/

	T					Temperature
/310.15/
;
**************************************************************
*********Defining Equations***********************************
EQUATIONS

	objective				objective function
	constraint_1(j)				constranit 1 of the formulation
	upper_bound(i)				upper bound
	lower_bound(i)				lower bound
	ATPtoADP				ATP to ADP ratio fixed according to Noor et al. 2014
	NADHtoNAD				NADH to NAD ratio fixed according to Noor et al. 2014
	glu_to_ac_con
	mass_balance(i)				steady state mass balance
	lower_bound_rxn(j)				lower bounds on reactions
	upper_bound_rxn(j)				upper bounds on reactions
	constraint_5(j)
	constraint_6(j)
	constraint_56
	therm
	pFBA1(j)
	pFBA2(j)
;
**************************************************************

*********Defining Variables***********************************

SCALARS

	M		/100000/;
POSITIVE VARIABLES

	x(i)					concentration of a metabolite

FREE VARIABLES

	v(j)					reaction flux
	B					Maximum driving force
	deltaG(j)				driving force in the given condition
	Z
	dummy(j)					objective value
;
BINARY VARIABLES

	n(j);
****************************************************************
v.up(aauprxns2(j))= 2.75;

v.up(off(j))= 0;

v.lo('EX_ser__L_e_b') = 2.75 ;
v.up('EX_ser__L_e_b') = 2.75 ;

v.up('EX_ac_e_b')= 0;

v.up('EX_lac__L_e')= 0;


*turn off "PPPGO2"
v.up('PPPGO2') = 0;
v.lo('PPPGO2') = 0;

*turn off "PFLr"
v.up('PFLr') = 0;
v.lo('PFLr') = 0;

v.up('Sa_biomass_universal') = 0.06127;
v.lo('Sa_biomass_universal') = 0.06127;

v.up('EX_glc__D_e_b') = 10;
v.lo('EX_glc__D_e_b') = 10;


*v.lo('ACKr')= 9.08243;

*v.lo('PTAr')= 9.08243;


***************Defining Model***********************************
objective..			Z =e= v('EX_ac_e_f');

mass_balance(i)..		sum(j,S(i,j)*v(j)) =e= 0;

lower_bound_rxn(j)..		v_min(j) =l= v(j);

upper_bound_rxn(j)..		v(j) =l= v_max(j);

pFBA1(j)..			v(j) =l= dummy(j);

pFBA2(j)..			-v(j) =l= dummy(j);

upper_bound(i)$cs(i)..		x(i) =l= log(Cmax(i));

lower_bound(i)$cs(i)..		x(i) =g= log(Cmin(i));


ATPtoADP..			x('atp_c[c]') =e= 7 * x('adp_c[c]');

NADHtoNAD..			x('nadh_c[c]') =e= (2.51/6.42) * x('nad_c[c]');


glu_to_ac_con..		x('glc__D_e[c]') =e= (1.65/3.4) * x('ac_c[c]');


constraint_5(j)$spec(j)..		v(j) =l= n(j) * v_max(j);

constraint_6(j)$spec(j)..		B =l= -deltaG(j) + M * (1-n(j));

constraint_1(j)$spec(j)..		-deltaG(j) =e= -delta_G_o(j) - R * T * sum(i, S(i,j) * x(i));

therm..			B =e= 1.5;

constraint_56..		sum(j$mem(j),mc(j)*v(j)) =e= 61;


Model TCA_MDF/all/;
******************************************************************

********************************************
solve TCA_MDF using mip maximizing Z;
****************Output File*****************

*********Output File for values*********************
FILE RESULTS4 /InteGraM.txt/;

PUT RESULTS4;

PUT "reaction	FLUX	DeltaG"/;

LOOP(j$spec(j),
	
	PUT j.tl:0:100,"    ", v.l(j):20:5,"	",deltaG.l(j):20:5/;
		
);

PUTCLOSE;


