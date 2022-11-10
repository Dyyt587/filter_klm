/*************************************************************************
    > File Name: krm.h
    > Author: Zhang Yuteng
    > Mail:fellylanma@aliyun.com
    > Created Time: 2019年05月13日 星期一 22时06分13秒
 ************************************************************************/

#ifndef _KLM_H
#define _KLM_H
#include "easyMatrix.h"


#define MAX_N 5 //状态向量X的维度//不可取大不可小！！！！
#define MAX_M 4 //测量向量Z的维度//不可取大不可小！！！！

#define CREATE_KLM_MATRIX_STATIC(x,y,matrix,_data) \
static struct easyMatrix matrix; \
    matrix.rows = x; \
    matrix.cols = y; \
    matrix.element = _data; 

typedef struct {
    float P;//估算协方差 初值不可以为0 ! ! ! ! ! 
    float out;//卡尔曼滤波器输出
    float Kg;//卡尔曼增益
    float Q;//过程噪声协方差
    float R;//观测噪声协方差
}KLM_SIP_T;//卡尔曼的简单（一维）

typedef struct {
    struct easyMatrix* X;    /* 状态向量 */ //输出值
    struct easyMatrix* Q;  /* 过程噪声协方差 */ //超参数
    struct easyMatrix* R;  /* 测量误差协方差 */ //超参数
    struct easyMatrix* F;  /* 用户定义的输出 f() state-transition function */
    struct easyMatrix* F_t;  /* 用户定义的输出 f() state-transition function */
    struct easyMatrix* B;  /* 控制矩阵 */
    struct easyMatrix* u;  /* 控制向量 */
    struct easyMatrix* Z;  /* 测量向量 */ //传感器测量值

    struct easyMatrix* H;  /* 雅可比矩阵的测量模型 */
    struct easyMatrix* Ht; /* 雅可比矩阵的转置测量 */

    struct easyMatrix* P;  /* 预测误差协方差 */ //初值不可为0
    struct easyMatrix* G;  /* 卡尔曼增益;也称为K */

    struct easyMatrix* I;  /*  */

    //struct easyMatrix* F;  /* 雅可比矩阵的过程模型 */
    //struct easyMatrix* Ft; /* 雅可比矩阵的转置过程 */
   // struct easyMatrix* Pp; /* P, 先验协方差 */

}KLM_T;

#endif//_MAGRIDE_PLANNING_EASYMATRIX_H_
