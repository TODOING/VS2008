#ifndef MATRIX_EIGENVALUE_H
#define MATRIX_EIGENVALUE_H

#include "Eigen/Dense"
using namespace Eigen;

// 标准矩阵，数值存在矩阵中
typedef Eigen::Matrix< double, Dynamic, Dynamic, RowMajor> MyMatrix;

// 标准列向量，数值存在矩阵中
typedef Eigen::Matrix< double, Dynamic, 1> MyVector;

class MatrixEigenvalue
{
public:
	/**
	* @brief 求实对称矩阵的特征值及特征向量
	* @param src_matrix							所求的矩阵
	* @param eigenvalues						矩阵的特征值
	* @param eigenvectors						矩阵的特征向量（按照列来存），如果特征值的序号为1，那么对应的特征向量为第1列
	* @param eigenvalues_percent				矩阵的特征值百分比
	* @param eigenvalues_accumulate_percent	    矩阵的特征值累积百分比
	* @param deps								累计误差
	* @return 返回值小于0表示超过迭代jt次仍未达到精度要求，返回值大于0表示正常返回 
	*/ 
	static bool GetMatrixEigen(MyMatrix src_matrix, MyVector &eigenvalues, MyMatrix &eigenvectors, MyVector *eigenvalues_percent = NULL, 
		MyVector *eigenvalues_accumulate_percent = NULL, double deps = 0.00000001);

private:
	/** 
	* @brief 求实对称矩阵的特征值及特征向量的雅格比法  
	* 利用雅格比(Jacobi)方法求实对称矩阵的全部特征值及特征向量  
	* @param symmetry_matrix      长度为n*n的数组，存放实对称矩阵，返回时对角线存放n个特征值  
	* @param dims				  矩阵的阶数  
	* @param eigenvectors		  长度为n*n的数组，返回特征向量(按列存储)  
	* @param eps				  控制精度要求  
	* @param jt					  整型变量，控制最大迭代次数  
	* @return 返回false表示超过迭代jt次仍未达到精度要求，返回true表示正常返回
	*/   
	static bool Jacobi(double symmetry_matrix[], int dim, double eigenvectors[], double eps, int jt);


};

#endif// MATRIX_EIGENVALUE_H