#ifndef MATRIX_EIGENVALUE_H
#define MATRIX_EIGENVALUE_H

#include "Eigen/Dense"
using namespace Eigen;

// ��׼������ֵ���ھ�����
typedef Eigen::Matrix< double, Dynamic, Dynamic, RowMajor> MyMatrix;

// ��׼����������ֵ���ھ�����
typedef Eigen::Matrix< double, Dynamic, 1> MyVector;

class MatrixEigenvalue
{
public:
	/**
	* @brief ��ʵ�Գƾ��������ֵ����������
	* @param src_matrix							����ľ���
	* @param eigenvalues						���������ֵ
	* @param eigenvectors						������������������������棩���������ֵ�����Ϊ1����ô��Ӧ����������Ϊ��1��
	* @param eigenvalues_percent				���������ֵ�ٷֱ�
	* @param eigenvalues_accumulate_percent	    ���������ֵ�ۻ��ٷֱ�
	* @param deps								�ۼ����
	* @return ����ֵС��0��ʾ��������jt����δ�ﵽ����Ҫ�󣬷���ֵ����0��ʾ�������� 
	*/ 
	static bool GetMatrixEigen(MyMatrix src_matrix, MyVector &eigenvalues, MyMatrix &eigenvectors, MyVector *eigenvalues_percent = NULL, 
		MyVector *eigenvalues_accumulate_percent = NULL, double deps = 0.00000001);

private:
	/** 
	* @brief ��ʵ�Գƾ��������ֵ�������������Ÿ�ȷ�  
	* �����Ÿ��(Jacobi)������ʵ�Գƾ����ȫ������ֵ����������  
	* @param symmetry_matrix      ����Ϊn*n�����飬���ʵ�Գƾ��󣬷���ʱ�Խ��ߴ��n������ֵ  
	* @param dims				  ����Ľ���  
	* @param eigenvectors		  ����Ϊn*n�����飬������������(���д洢)  
	* @param eps				  ���ƾ���Ҫ��  
	* @param jt					  ���ͱ�������������������  
	* @return ����false��ʾ��������jt����δ�ﵽ����Ҫ�󣬷���true��ʾ��������
	*/   
	static bool Jacobi(double symmetry_matrix[], int dim, double eigenvectors[], double eps, int jt);


};

#endif// MATRIX_EIGENVALUE_H