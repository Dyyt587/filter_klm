#include"klm.h"


/*how to use*/
//first �������󲢸���ֵps��P��ֵ����Ϊ0
//second ���������ӵ��������ṹ����
//third ����KLM_Hander() ֮���ٶ�ȡ����
KLM_T krm;
void KLM_Init(void)
{
	//int a[2*3] = { 2 };
#define xx_N 2
#define xx_M 1


	static val_xx_G[xx_N][xx_M] = { 0 };//d�Ծ����ʼ��
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_M, xx_G, val_xx_G);
	krm.G = &xx_G;

	static val_xx_Q[xx_N][xx_N] = { 0 };//d�Ծ����ʼ��
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_N, xx_Q, val_xx_Q);
	krm.Q = &xx_Q;

	static val_xx_F[xx_M][xx_M] = { 0 };//d�Ծ����ʼ��
	CREATE_KLM_MATRIX_STATIC(xx_M, xx_M, xx_F, val_xx_F);
	krm.F = &xx_F;

	static val_xx_X[xx_N][1] = { 0 };//d�Ծ����ʼ��
	CREATE_KLM_MATRIX_STATIC(xx_N, 1, xx_X, val_xx_X); //X���ǽṹ�屾��
	krm.X = &xx_X;

	static val_xx_Z[xx_M][1] = { 0 };//d�Ծ����ʼ��
	CREATE_KLM_MATRIX_STATIC(xx_M, 1, xx_Z, val_xx_Z);
	krm.Z = &xx_Z;

	static val_xx_H[xx_M][xx_N] = { 0 };//d�Ծ����ʼ��
	CREATE_KLM_MATRIX_STATIC(xx_M, xx_N, xx_H, val_xx_H);
	krm.X = &xx_H;

	static val_xx_B[xx_N][1] = { 0 };//d�Ծ����ʼ��
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_N, xx_B, val_xx_B);
	krm.B = &xx_B;

	static val_xx_u[xx_N][1] = { 0 };//d�Ծ����ʼ��
	CREATE_KLM_MATRIX_STATIC(xx_N, 1, xx_u, val_xx_u);
	krm.u = &xx_u;

	static val_xx_P[xx_N][xx_N] = { 1 };//d�Ծ����ʼ��,��ʼֵ��Ϊ0���������趨
	for (int i = 0; i < xx_N; ++i)
	{
		for (int j = 0; i < xx_N; ++i)
		{
			val_xx_P[i][j] = 1.0f;
		}
	}
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_N, xx_P, val_xx_P);
	krm.P = &xx_P;


//���¾��󲻳�ʼ����ֵ
	static val_xx_I[xx_N][xx_N] = { 0 };//�Ծ����ʼ������ʼ��Ϊ��λ��
	for (int i = 0; i < xx_N; ++i)
		val_xx_I[i][i] = 1;
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_N, xx_I, val_xx_I);
	krm.I = &xx_I;

	static val_xx_F_t[xx_N][xx_N] = { 0 };//d�Ծ����ʼ��
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_N, xx_F_t, val_xx_F_t);
	krm.F_t = &xx_F_t;

	static val_xx_Ht[xx_N][xx_M] = { 0 };//d�Ծ����ʼ��
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_M, xx_Ht, val_xx_Ht);
	krm.Ht = &xx_Ht;

}


/**
* �������˲���
* @param 	KLM_SIP_T* kfp �������ṹ�����
* float input ��Ҫ�˲��Ĳ����Ĳ���ֵ�����������Ĳɼ�ֵ��
* @return �˲���Ĳ���������ֵ��
*/
float KalmanFilter(KLM_SIP_T* kfp, float input)
{
	//����״̬���̣�X{k-1}= F (X{k-1_1}) + B{k}u{K}
	//kfp->out = kfp->out;//�����б�д
	//Ԥ��Э����̣�kʱ��ϵͳ����Э���� = k-1ʱ�̵�ϵͳЭ���� + ��������Э����
	kfp->P = kfp->P + kfp->Q;
	//���������淽�̣����������� = kʱ��ϵͳ����Э���� / ��kʱ��ϵͳ����Э���� + �۲�����Э���
	kfp->Kg = kfp->P / (kfp->P + kfp->R);
	//��������ֵ���̣�kʱ��״̬����������ֵ = ״̬������Ԥ��ֵ + ���������� * ������ֵ - ״̬������Ԥ��ֵ��
	kfp->out = kfp->out + kfp->Kg * (input - kfp->out);//��Ϊ��һ�ε�Ԥ��ֵ������һ�ε����ֵ
	//����Э�����: ���ε�ϵͳЭ����� kfp->LastP ����һ������׼����
	kfp->P = (1 - kfp->Kg) * kfp->P;
	return kfp->out;
}

/*	tmp0��N * 1	tmp1��N * 1
	tmp2��N * N	tmp3��N * M
	tmp4��M * N	tmp5��M * M
	tmp6��M * M	temp7 N * N
	tmp8��M * N tmp9��M * 1
	tmp10��M * 1
*/
//������ʱ����ռ�ô����ռ䣬��5����ʽ�ֱ�д�ɺ��������ڱ����������ٶȻ�ռ�ÿռ���Ż�
void _formula1(KLM_T* krm)
{
	float val0[MAX_N * 1] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_N, 1, tmp0, val0);

	float val1[MAX_N * 1] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_N, 1, tmp1, val1);
	//��ʽһ ����״̬����
	/* X{k-1}= F (X{k-1_1}) + B{k}u{K} */
	//X: N * 1 ��������
	//F��N*N��״̬ת�ƾ���
	//B��N*N�Ŀ��ƾ��� 
	//u��N*1�Ŀ�������
	//tmp0��N*1 tmp1��N*1
	multiMatrix(krm->F, krm->X, &tmp0);
	multiMatrix(krm->B, krm->u, &tmp1);
	addMatrix(&tmp0, &tmp1, krm->X);
}
void _formula2(KLM_T* krm)
{
	float val2[MAX_N * MAX_N] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_N, MAX_N, tmp2, val2);
	float val7[MAX_N * MAX_N] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_N, MAX_N, tmp7, val7);
	/* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
	// F_ ��״̬ת�ƾ��� Q����������
	//F ��F_t��N*N��״̬ת�ƾ���
	//P: N*N��Э�������
	//Q : N*N�Ĺ�����������
	//tmp2��N*N�ľ���
	//CREATE_MATRIX_ONSTACK(n,n,);
	multiMatrix(krm->F, krm->P, &tmp2);
	transMatrix(krm->F, krm->F_t);
	multiMatrix(&tmp2, krm->F_t, &tmp7);
	addMatrix(&tmp7, krm->Q, krm->P);
}
void _formula3(KLM_T* krm)
{

	float val3[MAX_N * MAX_M] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_N, MAX_M, tmp3, val3);

	float val4[MAX_M * MAX_N] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_M, MAX_N, tmp4, val4);

	float val5[MAX_M * MAX_M] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_M, MAX_M, tmp5, val5);

	float val6[MAX_M * MAX_M] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_M, MAX_M, tmp6, val6);
	/* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
	//P��P_t: N*N��Э�������
	//G��N*M�Ŀ���������
	//H��ȡ����Z�����X�������Z��M*1��H��M*N
	//H_t��N*M
	//tmp3��N*M tmp4��M*N tmp5��M*M tmp6��M*M
	transMatrix(krm->H, krm->Ht);
	multiMatrix(krm->P, krm->Ht, &tmp3);

	multiMatrix(krm->H, krm->P, &tmp4);
	multiMatrix(&tmp4, krm->Ht, &tmp5);
	addMatrix(&tmp5, krm->R, &tmp6);
	if (invMatrix(&tmp6, &tmp5) == 0)while (1);//�������ô���
	multiMatrix(&tmp3, &tmp5, krm->G);
}
void _formula4(KLM_T* krm)
{
	float val0[MAX_N * 1] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_N, 1, tmp0, val0);
	float val1[MAX_N * 1] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_N, 1, tmp1, val1);
	float val9[MAX_M * 1] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_M, 1, tmp9, val9);
	float val10[MAX_M * 1] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_M, 1, tmp10, val10);
	/* \hat{X}_k = \hat{x_k} + G_k(z_k - h(\hat{X}_k)) */
	//Z��M*1�Ĳ�������
	//X��N*1��������
	//H��ȡ����Z�����X�������Z��M*1��H��M*N
	//G��N*M�Ŀ���������
	//
	multiMatrix(krm->H, krm->X, &tmp9);//temp4���õ�M*1
	subMatrix(krm->Z, &tmp9, &tmp10);
	multiMatrix(krm->G, &tmp10, &tmp0);
	addMatrix(krm->X, &tmp0, &tmp1);
	copyMatrix(&tmp1, krm->X);
}
void _formula5(KLM_T* krm)
{
	float val2[MAX_N * MAX_N] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_N, MAX_N, tmp2, val2);
	float val7[MAX_N * MAX_N] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_N, MAX_N, tmp7, val7);
	/* P_k = (I - G_k H_k) P_k */
	//I:N*N�ĵ�λ����
	//G��N*M�Ŀ���������
	//H��ȡ����Z�����X�������Z��M*1��H��M*N
	//P��P_t: N*N��Э�������
	//
	multiMatrix(krm->G, krm->H, &tmp2);
	subMatrix(krm->I, &tmp2, &tmp7);
	multiMatrix(&tmp7, krm->P, &tmp2);
	copyMatrix(&tmp2, krm->P);
}
//Z:��������������
void KLM_Hander(KLM_T* krm, struct easyMatrix* Z)
{
	/*	tmp0��N * 1	tmp1��N * 1
	tmp2��N * N	tmp3��N * M
	tmp4��M * N	tmp5��M * M
	tmp6��M * M	temp7 N * N
	tmp8��M * N tmp9��M * 1
	tmp10��M * 1
*/
//float val0[MAX_N * 1] = { 0 }; 
//CREATE_MATRIX_ONSTACK(MAX_N, 1, tmp0, val0);
//float val1[MAX_N * 1] = { 0 };
//CREATE_MATRIX_ONSTACK(MAX_N, 1, tmp1, val1);
//float val2[MAX_N * MAX_N] = { 0 };
//CREATE_MATRIX_ONSTACK(MAX_N, MAX_N, tmp2, val2);
//float val3[MAX_N * MAX_M] = { 0 };
//CREATE_MATRIX_ONSTACK(MAX_N, MAX_M, tmp3, val3);
//float val4[MAX_M * MAX_N] = { 0 };
//CREATE_MATRIX_ONSTACK(MAX_M, MAX_N, tmp4, val4);
//float val5[MAX_M * MAX_M] = { 0 };
//CREATE_MATRIX_ONSTACK(MAX_M, MAX_M, tmp5, val5);
//float val6[MAX_M * MAX_M] = { 0 };
//CREATE_MATRIX_ONSTACK(MAX_M, MAX_M, tmp6, val6);
//float val7[MAX_N * MAX_N] = { 0 };
//CREATE_MATRIX_ONSTACK(MAX_N, MAX_N, tmp7, val7);
//float val8[MAX_M * MAX_N] = { 0 };
//CREATE_MATRIX_ONSTACK(MAX_M, MAX_N, tmp8, val8);
//float val9[MAX_M * 1] = { 0 };
//CREATE_MATRIX_ONSTACK(MAX_M, 1, tmp9, val9);
//float val10[MAX_M * 1] = { 0 };
//CREATE_MATRIX_ONSTACK(MAX_M, 1, tmp10, val10);
	copyMatrix(Z, krm->Z);
	_formula1(krm);
	_formula2(krm);
	_formula3(krm);
	_formula4(krm);
	_formula5(krm);
}

void kKLM_GetResult(KLM_T* krm , struct easyMatrix* X)
{
	copyMatrix(krm->X ,X);
}
struct easyMatrix* kKLM_GetResult_p(KLM_T* krm)
{
	return krm->X;
}