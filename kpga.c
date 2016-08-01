/*********************************************************/
/*  kpga.c                                               */
/*                                                       */
/*  �i�b�v�T�b�N���                                     */
/*  ��`�I�A���S���Y���iGA�j�ɂ�鋁���v���O����         */
/*                                                       */
/*                                                       */
/*********************************************************/
#include "kpga.h"


int kpga(){
	
	/*  �ו��̏�����   */
	Parcel pal, *ppal;
	ppal = &pal;
	
	initparcel(ppal);
	getparcel(ppal);

	/*  ���F�̂̏�����    */
	Ga_pool ga, *pga;
	pga = &ga;
	Ga_ngpool gang, *pgang;
	pgang = &gang;
	
	initpool(pga);
	getpool(pga);
	
	initngpool(pgang);
	// getngpool(pgang);
	
	srand((unsigned)time(NULL));    /* time�֐���p���������������̒�Ԏ� */
	
	int generation=0;
	for(generation=0; generation < LASTG; ++generation){
		printf("%d����\n", generation);
		mating(pga, pgang, ppal);
		mutation(pgang);
		selectng(pga, pgang, ppal);
		printp(pga, ppal);
	}
	
	// getpool(pga);
	
	return 0;
}


