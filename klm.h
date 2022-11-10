/*************************************************************************
    > File Name: krm.h
    > Author: Zhang Yuteng
    > Mail:fellylanma@aliyun.com
    > Created Time: 2019��05��13�� ����һ 22ʱ06��13��
 ************************************************************************/

#ifndef _KLM_H
#define _KLM_H
#include "easyMatrix.h"


#define MAX_N 5 //״̬����X��ά��//����ȡ�󲻿�С��������
#define MAX_M 4 //��������Z��ά��//����ȡ�󲻿�С��������

#define CREATE_KLM_MATRIX_STATIC(x,y,matrix,_data) \
static struct easyMatrix matrix; \
    matrix.rows = x; \
    matrix.cols = y; \
    matrix.element = _data; 

typedef struct {
    float P;//����Э���� ��ֵ������Ϊ0 ! ! ! ! ! 
    float out;//�������˲������
    float Kg;//����������
    float Q;//��������Э����
    float R;//�۲�����Э����
}KLM_SIP_T;//�������ļ򵥣�һά��

typedef struct {
    struct easyMatrix* X;    /* ״̬���� */ //���ֵ
    struct easyMatrix* Q;  /* ��������Э���� */ //������
    struct easyMatrix* R;  /* �������Э���� */ //������
    struct easyMatrix* F;  /* �û��������� f() state-transition function */
    struct easyMatrix* F_t;  /* �û��������� f() state-transition function */
    struct easyMatrix* B;  /* ���ƾ��� */
    struct easyMatrix* u;  /* �������� */
    struct easyMatrix* Z;  /* �������� */ //����������ֵ

    struct easyMatrix* H;  /* �ſɱȾ���Ĳ���ģ�� */
    struct easyMatrix* Ht; /* �ſɱȾ����ת�ò��� */

    struct easyMatrix* P;  /* Ԥ�����Э���� */ //��ֵ����Ϊ0
    struct easyMatrix* G;  /* ����������;Ҳ��ΪK */

    struct easyMatrix* I;  /*  */

    //struct easyMatrix* F;  /* �ſɱȾ���Ĺ���ģ�� */
    //struct easyMatrix* Ft; /* �ſɱȾ����ת�ù��� */
   // struct easyMatrix* Pp; /* P, ����Э���� */

}KLM_T;

#endif//_MAGRIDE_PLANNING_EASYMATRIX_H_
