/*********************************************************/
/*  ga.c                                                 */
/*                                                       */
/*  遺伝的アルゴリズム（GA）プログラム                   */
/*                                                       */
/*                                                       */
/*********************************************************/
#include "ga.h"

/*********************************************************/
/*  GA 関数                                              */
/*********************************************************/

/*********************************************************/
/* initparcel() : 荷物の初期化                           */
/*********************************************************/
void initparcel(Parcel *pa){
	
	int i=0;
	
	srand((unsigned)time(NULL));    /* time関数を用いた乱数初期化の定番式 */
	
	for(i=0; i < N; ++i){
		pa->parcel[i][0] = rndn(MAXVALUE);
		pa->parcel[i][1] = rndn(MAXVALUE);
	}
	
	// while( (i < N) && (scanf("%d %d", &pa->parcel[i][0], &pa->parcel[i][1]) != EOF) ){
	// 	++i;
	// }
	
}


/*********************************************************/
/* getparcel() : 荷物の値出力                            */
/*********************************************************/
void getparcel(Parcel *pa){
	
	int i=0;
	
	for(i=0; i < N; ++i){
		printf("parcel %d:\t", i);
		printf("weight = %d\t", pa->parcel[i][0]);
		printf("value = %d\n", pa->parcel[i][1]);
	}
	
}


/*********************************************************/
/* initpool() : 初期集団の生成                           */
/*********************************************************/
void initpool(Ga_pool *gapl){
	
	int i=0, j=0;
	
	srand((unsigned)time(NULL));    /* time関数を用いた乱数初期化の定番式 */
	
	for(i=0; i < POOLSIZE; ++i){
		for(j=0; j < N; ++j){
			gapl->pool[i][j] = rndn(2);
		}
	}
	
}


/*********************************************************/
/* getpool() : 集団の出力                                */
/*********************************************************/
void getpool(Ga_pool *gapl){
	
	int i=0, j=0;
	
	printf("*** pool  ***\n");
	for(i=0; i < POOLSIZE; ++i){
		printf("i=%d  ", i); 
		for(j=0; j < N; ++j){
			printf("%d", gapl->pool[i][j]);
		}
		printf("\n");
	}
	
}


/*********************************************************/
/* printp() : 結果出力                                   */
/*********************************************************/
void printp(Ga_pool *gapl, Parcel *pa){
	
	int i=0, j=0;
	int fitness;
	double totalfitness = 0;
	int elite, bestfit = 0;
	
	for(i=0; i < POOLSIZE; ++i){
		for(j=0; j < N; ++j){
			printf("%1d", gapl->pool[i][j]);
		}
		fitness = evalfit(i, gapl,pa);
		printf(":i=%d \t%d\n", i, fitness);
		if(fitness > bestfit){ /*  エリート解  */
			bestfit = fitness;
			elite = i;
		}
		totalfitness += fitness;
	}
	
	/*  エリート解の適応度を出力   */
	printf("elite=%d\tmax=%d \t", elite, bestfit);
	/*  平均の適応度を出力   */
	printf("average=%lf\n", totalfitness / POOLSIZE);
	
}


/*********************************************************/
/* initngpool() : 次世代染色体プールの初期化                */
/*********************************************************/
void initngpool(Ga_ngpool *gang){
	
	int i=0, j=0;
	
	srand((unsigned)time(NULL));    /* time関数を用いた乱数初期化の定番式 */
	
	for(i=0; i < POOLSIZE * 2; ++i){
		for(j=0; j < N; ++j){
			gang->ngpool[i][j] = rndn(2);
		}
	}
	
}

/*********************************************************/
/* getngpool() : 次世代染色体プールの出力                */
/*********************************************************/
void getngpool(Ga_ngpool *gang){
	
	int i=0, j=0;
	
	printf("*** ngpool  ***\n");
	for(i=0; i < POOLSIZE * 2; ++i){
		printf("i=%d  ", i); 
		for(j=0; j < N; ++j){
			printf("%d", gang->ngpool[i][j]);
		}
		printf("\n");
	}

}


/*********************************************************/
/* mating() : 交叉                                       */
/*********************************************************/
void mating(Ga_pool *gapl, Ga_ngpool *gang, Parcel *pa){
	
	int i=0;
	int totalfitness = 0;
	int roulette[POOLSIZE];
	int mama, papa;               /* 親の遺伝子の番号    */
	
	/*  ルーレット作成       */
	for(i=0; i < POOLSIZE; ++i){
		roulette[i] = evalfit(i, gapl, pa);
		/*  適応度の合計値を計算     */
		totalfitness += roulette[i];
	}
	
	/*  選択と交叉を繰り返す     */
	for(i=0; i < POOLSIZE; ++i){
		do{ /*  親の選択    */
			mama = selectp(roulette, totalfitness);
			papa = selectp(roulette, totalfitness);
		}while(mama == papa); /* 重複の削除    */
		
		/*  特定の2染色体の交叉    */
		crossing(mama, papa, i, gapl, gang);
	}
	
}


/*********************************************************/
/* evalfit() : 適応度の計算                              */
/*********************************************************/
int evalfit(int c, Ga_pool *gapl, Parcel *pa){
	
	int pos;           /*  遺伝子座の指定   */
	int value = 0;     /*  評価値    */
	int weight = 0;    /*  重量  */
	
	for(pos=0; pos < N; ++pos){
		weight += pa->parcel[pos][0] * gapl->pool[c][pos];
		value += pa->parcel[pos][1] * gapl->pool[c][pos];
	}
	
	if(weight >= WEIGHTLIMIT) value = 0;
	
	return value;
	
}

/*********************************************************/
/* selectng() : 親の選択                                 */
/*********************************************************/
int selectp(int roulette[POOLSIZE], int totalfitness){
	
	int i=0;
	int ball;
	int acc;
	
	// srand((unsigned)time(NULL));    /* time関数を用いた乱数初期化の定番式 */
	
	ball = rndn(totalfitness);
	for(i=0; i < POOLSIZE; ++i){
		acc += roulette[i];
		if(acc > ball) break;
	}
	
	return i;
	
}


/*********************************************************/
/* crossing() : 特定の2染色体の交叉                      */
/*********************************************************/
void crossing(int mama, int papa, int i, Ga_pool *gapl, Ga_ngpool *gang){
	
	int j=0;
	int cp;    /*  交叉する点    */
	
	// srand((unsigned)time(NULL));    /* time関数を用いた乱数初期化の定番式 */
	
	/*  交叉点の決定    */
	cp = rndn(N);
	
	/*  前半部分のコピー     */
	for(j=0; j < cp; ++j){
		gang->ngpool[i * 2][j]     = gapl->pool[mama][j];
		gang->ngpool[i * 2 + 1][j] = gapl->pool[papa][j];
	}
	
	/*  後半部分のコピー     */
	for(j=cp; j < N; ++j){
		gang->ngpool[i * 2 + 1][j] = gapl->pool[mama][j];
		gang->ngpool[i * 2][j]     = gapl->pool[papa][j];
	}
	
}


/*********************************************************/
/* mutation() : 突然変異                                 */
/*********************************************************/
void mutation(Ga_ngpool *gang){
	
	int i=0, j=0;
	
	for(i=0; i < POOLSIZE; ++i){
		for(j=0; j < N; ++j){
			if((double)rndn(100) / 100.0 <= MRATE){
				/*  反転の突然変異   */
				gang->ngpool[i][j] = notval(gang->ngpool[i][j]);
			}
		}
	}
	
}


/*********************************************************/
/* notval() : 真理値の反転                               */
/*********************************************************/
int notval(int v){
	
	if(v == YES) return NO;
	else return YES;
	
}


/*********************************************************/
/* selectng() : 選択                                     */
/*********************************************************/
void selectng(Ga_pool *gapl, Ga_ngpool *gang, Parcel *pa){
	
	int i, j, c;
	int totalfitness = 0;         /*  適応度の合計値  */
	int roulette[POOLSIZE * 2];   /*  適応度を格納   */
	int ball;                     /*  玉（選択位置の数値）  */
	int acc = 0;
	
	/* 選択を繰り返す   */
	for(i=0; i < POOLSIZE; ++i){
		/*  ルーレット作成    */
		totalfitness = 0;
		for(c=0; c < POOLSIZE * 2; ++c){
			roulette[c] = evalfitng(c, gang, pa);
			/* 適応度の合計値を計算     */
			totalfitness += roulette[c];
		}
		
		// srand((unsigned)time(NULL));    /* time関数を用いた乱数初期化の定番式 */
		
		/*  染色体を一つ選ぶ   */
		ball = rndn(totalfitness);
		acc = 0;
		for(c=0; c < POOLSIZE; ++c){
			acc += roulette[c];
			if(acc > ball) break;
		}
		
		/*  染色体のコピー    */
		for(j=0; j < N; ++j){
			gapl->pool[i][j] = gang->ngpool[c][j];
		}
		
	}
	
}


/*********************************************************/
/* evalfitng() : 適応度の計算                            */
/*********************************************************/
int evalfitng(int c, Ga_ngpool *gang, Parcel *pa){
	
	int pos;           /*  遺伝子座の指定   */
	int value = 0;     /*  評価値    */
	int weight = 0;    /*  重量  */
	
	for(pos=0; pos < N; ++pos){
		weight += pa->parcel[pos][0] * gang->ngpool[c][pos];
		value += pa->parcel[pos][1] * gang->ngpool[c][pos];
	}
	
	if(weight >= WEIGHTLIMIT) value = 0;
	
	return value;
	
}

