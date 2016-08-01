/*********************************************************/
/*  ga.h                                                 */
/*                                                       */
/*  遺伝的アルゴリズム（GA）プログラム                   */
/*                                                       */
/*                                                       */
/*********************************************************/
#ifndef _GA_H_
#define _GA_H_

#include <stdio.h>
#include <stdlib.h>
#include "rnd.h"

/*  データ構造                               */
#define MAXVALUE 100                     /* 重さと価値の最大値    */
#define N 30                             /* 荷物の個数    */
#define POOLSIZE 30                      /* プールサイズ  */
#define WEIGHTLIMIT (N * MAXVALUE / 4)   /* 重量制限  */
#define MRATE 0.01                       /* 突然変異の確率  */
#define YES 1                            /* YESに対する整数値  */
#define NO 0                             /* NOに対する整数値  */

typedef struct{
	int parcel[N][2];        /* 荷物     */
} Parcel;


typedef struct{
	int pool[POOLSIZE][N];        /* 染色体プール     */
} Ga_pool;


typedef struct{
	int ngpool[POOLSIZE * 2][N];  /* 次世代染色体プール    */
} Ga_ngpool;


/*  関数                       */
void initparcel(Parcel *parcel);   // 荷物の初期化
void getparcel(Parcel *pa);        // 荷物の値出力
void initpool(Ga_pool *gapl);      // 初期集団の生成
void getpool(Ga_pool *gapl);       // 集団の出力
void printp(Ga_pool *gapl, Parcel *pa);  // 結果出力
	
void initngpool(Ga_ngpool *gang);      // 次世代染色体プールの初期化
void getngpool(Ga_ngpool *gang);       // 次世代染色体プールの出力


/*  交叉  */
void mating(Ga_pool *gapl, Ga_ngpool *gang, Parcel *pa);
int evalfit(int c, Ga_pool *gapl, Parcel *pa);              /*  適応度の計算   */
int selectp(int roulette[POOLSIZE], int totalfitness);      /*  親の選択   */
void crossing(int mama, int papa, int i, Ga_pool *gapl, Ga_ngpool *gang);   /*  特定の2染色体の交叉   */

/*  突然変異  */
void mutation(Ga_ngpool *gang);
int notval(int v);   /*  真理値の反転  */

/*  選択  */
void selectng(Ga_pool *gapl, Ga_ngpool *gang, Parcel *pa);
int evalfitng(int c, Ga_ngpool *gang, Parcel *pa);          /*  適応度の計算   */


#endif

