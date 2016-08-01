/*********************************************************/
/*  ga.h                                                 */
/*                                                       */
/*  ��`�I�A���S���Y���iGA�j�v���O����                   */
/*                                                       */
/*                                                       */
/*********************************************************/
#ifndef _GA_H_
#define _GA_H_

#include <stdio.h>
#include <stdlib.h>
#include "rnd.h"

/*  �f�[�^�\��                               */
#define MAXVALUE 100                     /* �d���Ɖ��l�̍ő�l    */
#define N 30                             /* �ו��̌�    */
#define POOLSIZE 30                      /* �v�[���T�C�Y  */
#define WEIGHTLIMIT (N * MAXVALUE / 4)   /* �d�ʐ���  */
#define MRATE 0.01                       /* �ˑR�ψق̊m��  */
#define YES 1                            /* YES�ɑ΂��鐮���l  */
#define NO 0                             /* NO�ɑ΂��鐮���l  */

typedef struct{
	int parcel[N][2];        /* �ו�     */
} Parcel;


typedef struct{
	int pool[POOLSIZE][N];        /* ���F�̃v�[��     */
} Ga_pool;


typedef struct{
	int ngpool[POOLSIZE * 2][N];  /* ��������F�̃v�[��    */
} Ga_ngpool;


/*  �֐�                       */
void initparcel(Parcel *parcel);   // �ו��̏�����
void getparcel(Parcel *pa);        // �ו��̒l�o��
void initpool(Ga_pool *gapl);      // �����W�c�̐���
void getpool(Ga_pool *gapl);       // �W�c�̏o��
void printp(Ga_pool *gapl, Parcel *pa);  // ���ʏo��
	
void initngpool(Ga_ngpool *gang);      // ��������F�̃v�[���̏�����
void getngpool(Ga_ngpool *gang);       // ��������F�̃v�[���̏o��


/*  ����  */
void mating(Ga_pool *gapl, Ga_ngpool *gang, Parcel *pa);
int evalfit(int c, Ga_pool *gapl, Parcel *pa);              /*  �K���x�̌v�Z   */
int selectp(int roulette[POOLSIZE], int totalfitness);      /*  �e�̑I��   */
void crossing(int mama, int papa, int i, Ga_pool *gapl, Ga_ngpool *gang);   /*  �����2���F�̂̌���   */

/*  �ˑR�ψ�  */
void mutation(Ga_ngpool *gang);
int notval(int v);   /*  �^���l�̔��]  */

/*  �I��  */
void selectng(Ga_pool *gapl, Ga_ngpool *gang, Parcel *pa);
int evalfitng(int c, Ga_ngpool *gang, Parcel *pa);          /*  �K���x�̌v�Z   */


#endif

