/*********************************************************/
/*  kpga.c                                               */
/*                                                       */
/*  ナップサック問題                                     */
/*  遺伝的アルゴリズム（GA）による求解プログラム         */
/*                                                       */
/*                                                       */
/*********************************************************/
#include "kpga.h"


int kpga(){
	
	/*  荷物の初期化   */
	Parcel pal, *ppal;
	ppal = &pal;
	
	initparcel(ppal);
	getparcel(ppal);

	/*  染色体の初期化    */
	Ga_pool ga, *pga;
	pga = &ga;
	Ga_ngpool gang, *pgang;
	pgang = &gang;
	
	initpool(pga);
	getpool(pga);
	
	initngpool(pgang);
	// getngpool(pgang);
	
	srand((unsigned)time(NULL));    /* time関数を用いた乱数初期化の定番式 */
	
	int generation=0;
	for(generation=0; generation < LASTG; ++generation){
		printf("%d世代\n", generation);
		mating(pga, pgang, ppal);
		mutation(pgang);
		selectng(pga, pgang, ppal);
		printp(pga, ppal);
	}
	
	// getpool(pga);
	
	return 0;
}


