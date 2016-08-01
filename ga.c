/*********************************************************/
/*  ga.c                                                 */
/*                                                       */
/*  ��`�I�A���S���Y���iGA�j�v���O����                   */
/*                                                       */
/*                                                       */
/*********************************************************/
#include "ga.h"

/*********************************************************/
/*  GA �֐�                                              */
/*********************************************************/

/*********************************************************/
/* initparcel() : �ו��̏�����                           */
/*********************************************************/
void initparcel(Parcel *pa){
	
	int i=0;
	
	srand((unsigned)time(NULL));    /* time�֐���p���������������̒�Ԏ� */
	
	for(i=0; i < N; ++i){
		pa->parcel[i][0] = rndn(MAXVALUE);
		pa->parcel[i][1] = rndn(MAXVALUE);
	}
	
	// while( (i < N) && (scanf("%d %d", &pa->parcel[i][0], &pa->parcel[i][1]) != EOF) ){
	// 	++i;
	// }
	
}


/*********************************************************/
/* getparcel() : �ו��̒l�o��                            */
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
/* initpool() : �����W�c�̐���                           */
/*********************************************************/
void initpool(Ga_pool *gapl){
	
	int i=0, j=0;
	
	srand((unsigned)time(NULL));    /* time�֐���p���������������̒�Ԏ� */
	
	for(i=0; i < POOLSIZE; ++i){
		for(j=0; j < N; ++j){
			gapl->pool[i][j] = rndn(2);
		}
	}
	
}


/*********************************************************/
/* getpool() : �W�c�̏o��                                */
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
/* printp() : ���ʏo��                                   */
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
		if(fitness > bestfit){ /*  �G���[�g��  */
			bestfit = fitness;
			elite = i;
		}
		totalfitness += fitness;
	}
	
	/*  �G���[�g���̓K���x���o��   */
	printf("elite=%d\tmax=%d \t", elite, bestfit);
	/*  ���ς̓K���x���o��   */
	printf("average=%lf\n", totalfitness / POOLSIZE);
	
}


/*********************************************************/
/* initngpool() : ��������F�̃v�[���̏�����                */
/*********************************************************/
void initngpool(Ga_ngpool *gang){
	
	int i=0, j=0;
	
	srand((unsigned)time(NULL));    /* time�֐���p���������������̒�Ԏ� */
	
	for(i=0; i < POOLSIZE * 2; ++i){
		for(j=0; j < N; ++j){
			gang->ngpool[i][j] = rndn(2);
		}
	}
	
}

/*********************************************************/
/* getngpool() : ��������F�̃v�[���̏o��                */
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
/* mating() : ����                                       */
/*********************************************************/
void mating(Ga_pool *gapl, Ga_ngpool *gang, Parcel *pa){
	
	int i=0;
	int totalfitness = 0;
	int roulette[POOLSIZE];
	int mama, papa;               /* �e�̈�`�q�̔ԍ�    */
	
	/*  ���[���b�g�쐬       */
	for(i=0; i < POOLSIZE; ++i){
		roulette[i] = evalfit(i, gapl, pa);
		/*  �K���x�̍��v�l���v�Z     */
		totalfitness += roulette[i];
	}
	
	/*  �I���ƌ������J��Ԃ�     */
	for(i=0; i < POOLSIZE; ++i){
		do{ /*  �e�̑I��    */
			mama = selectp(roulette, totalfitness);
			papa = selectp(roulette, totalfitness);
		}while(mama == papa); /* �d���̍폜    */
		
		/*  �����2���F�̂̌���    */
		crossing(mama, papa, i, gapl, gang);
	}
	
}


/*********************************************************/
/* evalfit() : �K���x�̌v�Z                              */
/*********************************************************/
int evalfit(int c, Ga_pool *gapl, Parcel *pa){
	
	int pos;           /*  ��`�q���̎w��   */
	int value = 0;     /*  �]���l    */
	int weight = 0;    /*  �d��  */
	
	for(pos=0; pos < N; ++pos){
		weight += pa->parcel[pos][0] * gapl->pool[c][pos];
		value += pa->parcel[pos][1] * gapl->pool[c][pos];
	}
	
	if(weight >= WEIGHTLIMIT) value = 0;
	
	return value;
	
}

/*********************************************************/
/* selectng() : �e�̑I��                                 */
/*********************************************************/
int selectp(int roulette[POOLSIZE], int totalfitness){
	
	int i=0;
	int ball;
	int acc;
	
	// srand((unsigned)time(NULL));    /* time�֐���p���������������̒�Ԏ� */
	
	ball = rndn(totalfitness);
	for(i=0; i < POOLSIZE; ++i){
		acc += roulette[i];
		if(acc > ball) break;
	}
	
	return i;
	
}


/*********************************************************/
/* crossing() : �����2���F�̂̌���                      */
/*********************************************************/
void crossing(int mama, int papa, int i, Ga_pool *gapl, Ga_ngpool *gang){
	
	int j=0;
	int cp;    /*  ��������_    */
	
	// srand((unsigned)time(NULL));    /* time�֐���p���������������̒�Ԏ� */
	
	/*  �����_�̌���    */
	cp = rndn(N);
	
	/*  �O�������̃R�s�[     */
	for(j=0; j < cp; ++j){
		gang->ngpool[i * 2][j]     = gapl->pool[mama][j];
		gang->ngpool[i * 2 + 1][j] = gapl->pool[papa][j];
	}
	
	/*  �㔼�����̃R�s�[     */
	for(j=cp; j < N; ++j){
		gang->ngpool[i * 2 + 1][j] = gapl->pool[mama][j];
		gang->ngpool[i * 2][j]     = gapl->pool[papa][j];
	}
	
}


/*********************************************************/
/* mutation() : �ˑR�ψ�                                 */
/*********************************************************/
void mutation(Ga_ngpool *gang){
	
	int i=0, j=0;
	
	for(i=0; i < POOLSIZE; ++i){
		for(j=0; j < N; ++j){
			if((double)rndn(100) / 100.0 <= MRATE){
				/*  ���]�̓ˑR�ψ�   */
				gang->ngpool[i][j] = notval(gang->ngpool[i][j]);
			}
		}
	}
	
}


/*********************************************************/
/* notval() : �^���l�̔��]                               */
/*********************************************************/
int notval(int v){
	
	if(v == YES) return NO;
	else return YES;
	
}


/*********************************************************/
/* selectng() : �I��                                     */
/*********************************************************/
void selectng(Ga_pool *gapl, Ga_ngpool *gang, Parcel *pa){
	
	int i, j, c;
	int totalfitness = 0;         /*  �K���x�̍��v�l  */
	int roulette[POOLSIZE * 2];   /*  �K���x���i�[   */
	int ball;                     /*  �ʁi�I���ʒu�̐��l�j  */
	int acc = 0;
	
	/* �I�����J��Ԃ�   */
	for(i=0; i < POOLSIZE; ++i){
		/*  ���[���b�g�쐬    */
		totalfitness = 0;
		for(c=0; c < POOLSIZE * 2; ++c){
			roulette[c] = evalfitng(c, gang, pa);
			/* �K���x�̍��v�l���v�Z     */
			totalfitness += roulette[c];
		}
		
		// srand((unsigned)time(NULL));    /* time�֐���p���������������̒�Ԏ� */
		
		/*  ���F�̂���I��   */
		ball = rndn(totalfitness);
		acc = 0;
		for(c=0; c < POOLSIZE; ++c){
			acc += roulette[c];
			if(acc > ball) break;
		}
		
		/*  ���F�̂̃R�s�[    */
		for(j=0; j < N; ++j){
			gapl->pool[i][j] = gang->ngpool[c][j];
		}
		
	}
	
}


/*********************************************************/
/* evalfitng() : �K���x�̌v�Z                            */
/*********************************************************/
int evalfitng(int c, Ga_ngpool *gang, Parcel *pa){
	
	int pos;           /*  ��`�q���̎w��   */
	int value = 0;     /*  �]���l    */
	int weight = 0;    /*  �d��  */
	
	for(pos=0; pos < N; ++pos){
		weight += pa->parcel[pos][0] * gang->ngpool[c][pos];
		value += pa->parcel[pos][1] * gang->ngpool[c][pos];
	}
	
	if(weight >= WEIGHTLIMIT) value = 0;
	
	return value;
	
}

