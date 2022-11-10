#include"klm.h"


/*how to use*/
//first 创建矩阵并赋初值ps；P初值不可为0
//second 将矩阵链接到卡尔曼结构体上
//third 调用KLM_Hander() 之后再读取即可
KLM_T krm;
void KLM_Init(void)
{
	//int a[2*3] = { 2 };
#define xx_N 2
#define xx_M 1


	static val_xx_G[xx_N][xx_M] = { 0 };//d对矩阵初始化
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_M, xx_G, val_xx_G);
	krm.G = &xx_G;

	static val_xx_Q[xx_N][xx_N] = { 0 };//d对矩阵初始化
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_N, xx_Q, val_xx_Q);
	krm.Q = &xx_Q;

	static val_xx_F[xx_M][xx_M] = { 0 };//d对矩阵初始化
	CREATE_KLM_MATRIX_STATIC(xx_M, xx_M, xx_F, val_xx_F);
	krm.F = &xx_F;

	static val_xx_X[xx_N][1] = { 0 };//d对矩阵初始化
	CREATE_KLM_MATRIX_STATIC(xx_N, 1, xx_X, val_xx_X); //X就是结构体本身
	krm.X = &xx_X;

	static val_xx_Z[xx_M][1] = { 0 };//d对矩阵初始化
	CREATE_KLM_MATRIX_STATIC(xx_M, 1, xx_Z, val_xx_Z);
	krm.Z = &xx_Z;

	static val_xx_H[xx_M][xx_N] = { 0 };//d对矩阵初始化
	CREATE_KLM_MATRIX_STATIC(xx_M, xx_N, xx_H, val_xx_H);
	krm.X = &xx_H;

	static val_xx_B[xx_N][1] = { 0 };//d对矩阵初始化
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_N, xx_B, val_xx_B);
	krm.B = &xx_B;

	static val_xx_u[xx_N][1] = { 0 };//d对矩阵初始化
	CREATE_KLM_MATRIX_STATIC(xx_N, 1, xx_u, val_xx_u);
	krm.u = &xx_u;

	static val_xx_P[xx_N][xx_N] = { 1 };//d对矩阵初始化,初始值不为0，可自行设定
	for (int i = 0; i < xx_N; ++i)
	{
		for (int j = 0; i < xx_N; ++i)
		{
			val_xx_P[i][j] = 1.0f;
		}
	}
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_N, xx_P, val_xx_P);
	krm.P = &xx_P;


//以下矩阵不初始其他值
	static val_xx_I[xx_N][xx_N] = { 0 };//对矩阵初始化，初始化为单位阵
	for (int i = 0; i < xx_N; ++i)
		val_xx_I[i][i] = 1;
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_N, xx_I, val_xx_I);
	krm.I = &xx_I;

	static val_xx_F_t[xx_N][xx_N] = { 0 };//d对矩阵初始化
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_N, xx_F_t, val_xx_F_t);
	krm.F_t = &xx_F_t;

	static val_xx_Ht[xx_N][xx_M] = { 0 };//d对矩阵初始化
	CREATE_KLM_MATRIX_STATIC(xx_N, xx_M, xx_Ht, val_xx_Ht);
	krm.Ht = &xx_Ht;

}


/**
* 卡尔曼滤波器
* @param 	KLM_SIP_T* kfp 卡尔曼结构体参数
* float input 需要滤波的参数的测量值（即传感器的采集值）
* @return 滤波后的参数（最优值）
*/
float KalmanFilter(KLM_SIP_T* kfp, float input)
{
	//先验状态方程；X{k-1}= F (X{k-1_1}) + B{k}u{K}
	//kfp->out = kfp->out;//需自行编写
	//预测协方差方程：k时刻系统估算协方差 = k-1时刻的系统协方差 + 过程噪声协方差
	kfp->P = kfp->P + kfp->Q;
	//卡尔曼增益方程：卡尔曼增益 = k时刻系统估算协方差 / （k时刻系统估算协方差 + 观测噪声协方差）
	kfp->Kg = kfp->P / (kfp->P + kfp->R);
	//更新最优值方程：k时刻状态变量的最优值 = 状态变量的预测值 + 卡尔曼增益 * （测量值 - 状态变量的预测值）
	kfp->out = kfp->out + kfp->Kg * (input - kfp->out);//因为这一次的预测值就是上一次的输出值
	//更新协方差方程: 本次的系统协方差付给 kfp->LastP 威下一次运算准备。
	kfp->P = (1 - kfp->Kg) * kfp->P;
	return kfp->out;
}

/*	tmp0；N * 1	tmp1：N * 1
	tmp2：N * N	tmp3：N * M
	tmp4：M * N	tmp5：M * M
	tmp6：M * M	temp7 N * N
	tmp8：M * N tmp9：M * 1
	tmp10：M * 1
*/
//鉴于临时矩阵占用大量空间，将5个公式分别写成函数，便于编译器进行速度或占用空间的优化
void _formula1(KLM_T* krm)
{
	float val0[MAX_N * 1] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_N, 1, tmp0, val0);

	float val1[MAX_N * 1] = { 0 };
	CREATE_MATRIX_ONSTACK(MAX_N, 1, tmp1, val1);
	//公式一 先验状态方程
	/* X{k-1}= F (X{k-1_1}) + B{k}u{K} */
	//X: N * 1 的列向量
	//F：N*N的状态转移矩阵
	//B：N*N的控制矩阵 
	//u：N*1的控制向量
	//tmp0；N*1 tmp1：N*1
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
	// F_ ：状态转移矩阵 Q：过程噪声
	//F 、F_t：N*N的状态转移矩阵
	//P: N*N的协方差矩阵
	//Q : N*N的过程噪声矩阵
	//tmp2：N*N的矩阵
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
	//P、P_t: N*N的协方差矩阵
	//G；N*M的卡尔曼增益
	//H；取决于Z矩阵和X矩阵，如果Z：M*1，H：M*N
	//H_t；N*M
	//tmp3：N*M tmp4：M*N tmp5：M*M tmp6：M*M
	transMatrix(krm->H, krm->Ht);
	multiMatrix(krm->P, krm->Ht, &tmp3);

	multiMatrix(krm->H, krm->P, &tmp4);
	multiMatrix(&tmp4, krm->Ht, &tmp5);
	addMatrix(&tmp5, krm->R, &tmp6);
	if (invMatrix(&tmp6, &tmp5) == 0)while (1);//矩阵设置错误
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
	//Z：M*1的测量向量
	//X：N*1的列向量
	//H；取决于Z矩阵和X矩阵，如果Z：M*1，H：M*N
	//G；N*M的卡尔曼增益
	//
	multiMatrix(krm->H, krm->X, &tmp9);//temp4仅用到M*1
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
	//I:N*N的单位矩阵
	//G；N*M的卡尔曼增益
	//H；取决于Z矩阵和X矩阵，如果Z：M*1，H：M*N
	//P、P_t: N*N的协方差矩阵
	//
	multiMatrix(krm->G, krm->H, &tmp2);
	subMatrix(krm->I, &tmp2, &tmp7);
	multiMatrix(&tmp7, krm->P, &tmp2);
	copyMatrix(&tmp2, krm->P);
}
//Z:传感器测量矩阵
void KLM_Hander(KLM_T* krm, struct easyMatrix* Z)
{
	/*	tmp0；N * 1	tmp1：N * 1
	tmp2：N * N	tmp3：N * M
	tmp4：M * N	tmp5：M * M
	tmp6：M * M	temp7 N * N
	tmp8：M * N tmp9：M * 1
	tmp10：M * 1
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