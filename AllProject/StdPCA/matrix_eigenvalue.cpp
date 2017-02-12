#include "matrix_eigenvalue.h"

#include <math.h>
#include <map>

using std::map;

#ifndef ABS
#  define ABS(x)        ((x<0) ? (-1*(x)) : x)
#endif

bool MatrixEigenvalue::GetMatrixEigen(MyMatrix src_matrix, MyVector &eigenvalues, MyMatrix &eigenvectors, MyVector *eigenvalues_percent /* = NULL */, MyVector *eigenvalues_accumulate_percent /* = NULL */, double deps /* = 0.00000001 */)
{
	int dim = src_matrix.rows();
	double *matrix = new double[dim*dim];
	double *matrix_eigenvectors = new double[dim*dim];

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			int index = i * dim + j;
			matrix[index] = src_matrix(i, j);
			matrix_eigenvectors[index] = 0.0;
		}
	}

	bool return_value = Jacobi(matrix, dim, matrix_eigenvectors, deps, 100);
	if (!return_value)
		return false;

	MyVector temp_eigenvalues(dim);
	MyMatrix temp_eigenvectors(dim, dim);

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			int index = i*dim + j;
			temp_eigenvectors(i, j) = matrix_eigenvectors[index];

			if (i == j)
				temp_eigenvalues[i] = matrix[index];
		}
	}

	if (matrix)
	{
		delete []matrix;
		matrix = NULL;
	}
	if (matrix_eigenvectors)
	{
		delete []matrix_eigenvectors;
		matrix_eigenvectors = NULL;
	}

	double sum_eigenvalue = 0.0;// 特征值总和
	map<double, int> map_eigenvalue;// 记录每个特征值所在的位置，即所在波段位置
	for (size_t index = 0; index < temp_eigenvalues.size(); index++)
	{
		map_eigenvalue[temp_eigenvalues[index]] = index;
		sum_eigenvalue += temp_eigenvalues[index];
	}

	// 对协方差矩阵的特征值和特征向量排序并计算累计百分比和百分比
	map<double, int>::reverse_iterator iter = map_eigenvalue.rbegin();
	for (size_t ii = 0; ii < temp_eigenvalues.size(); ii++)
	{
		if (iter == map_eigenvalue.rend())
			break;

		eigenvalues[ii] = iter->first;
		int index = iter->second;// 获取该特征值对应的特征向量对应的列号

		if (eigenvalues_percent != NULL)// 计算百分比以及累积百分比
		{
			double temp = iter->first / sum_eigenvalue;
			(*eigenvalues_percent)[ii] = temp;

			if (eigenvalues_accumulate_percent != NULL)
			{
				if (ii != 0)
					(*eigenvalues_accumulate_percent)[ii] = (*eigenvalues_accumulate_percent)[ii - 1] + temp;
				else
					(*eigenvalues_accumulate_percent)[ii] = temp;
			}
		}

		double max_value = ABS(temp_eigenvectors(0, index));
		int num = 0;
		for (int row = 0; row < dim; row++)// 获取特征向量中绝对值最大的位置
		{
			double temp = ABS(temp_eigenvectors(row, index));

			if (max_value < temp)
			{
				max_value = temp;
				num = row;
			}
		}

		bool is_positive = false;// 判断最大的特征向量中的值是否为正
		if (max_value - temp_eigenvectors(num, index) < 0.000000001)
			is_positive = true;

		for (int row = 0; row < dim; row++)// 确保每一个特征向量的绝对值最大的都为正
		{
			if (!is_positive)
				eigenvectors(row, ii) = -temp_eigenvectors(row, index);
			else
				eigenvectors(row, ii) = temp_eigenvectors(row, index);
		}

		iter++;
	}

	return true;
}

bool MatrixEigenvalue::Jacobi(double symmetry_matrix[], int dim, double eigenvectors[], double eps, int jt)
{
	int i, j, p, q, u, w, t, s, l;   
	double fm, cn, sn, omega, x, y, d; 

	l = 1;	   
	for(i = 0; i <= dim-1; i++)// 初始化特征向量矩阵使其非对角线元素全为0，对角线元素为1
	{     
		eigenvectors[i*dim + i] = 1.0;   
		for(j = 0; j <= dim-1; j++)   
		{   
			if(i != j)     
				eigenvectors[i*dim + j] = 0.0;   
		}   
	}   

	while (true)
	{
		fm = 0.0;   
		for(i = 0; i <= dim-1; i++)// 找出,矩阵symmetry_matrix(特征值),中除对角线外其他元素的最大绝对值   
		{   
			for(j = 0; j <= dim-1; j++)// 这个最大值是位于symmetry_matrix[p][q] ,等于fm   
			{   
				d = fabs(symmetry_matrix[i*dim + j]);   

				if((i != j) && (d > fm))   
				{   
					fm = d;     
					p = i;     
					q = j;   
				}   
			}   
		} 

		if(fm < eps)// 精度复合要求   
			return true;   

		if(l > jt)// 迭代次数太多   
			return false;

		l++;// 迭代计数器   
		u = p*dim + q;   
		w = p*dim + p;     
		t = q*dim + p;     
		s = q*dim + q;   
		x = -symmetry_matrix[u];   
		y = (symmetry_matrix[s]-symmetry_matrix[w]) / 2.0;// x y的求法不同   
		omega = x / sqrt(x*x + y*y);//sin2θ  

		//tan2θ=x/y = -2.0*a[u]/(a[s]-a[w])   
		if(y < 0.0)   
			omega = -omega;

		sn = 1.0 + sqrt(1.0 - omega*omega);     
		sn = omega / sqrt(2.0*sn);//sinθ   
		cn = sqrt(1.0 - sn*sn);//cosθ

		fm = symmetry_matrix[w];// 变换前的symmetry_matrix[w]   symmetry_matrix[p][p]   
		symmetry_matrix[w] = fm*cn*cn + symmetry_matrix[s]*sn*sn + symmetry_matrix[u]*omega;   
		symmetry_matrix[s] = fm*sn*sn + symmetry_matrix[s]*cn*cn - symmetry_matrix[u]*omega;   
		symmetry_matrix[u] = 0.0;   
		symmetry_matrix[t] = 0.0;

		// 以下是旋转矩阵,旋转了了p行,q行,p列,q列   
		// 但是四个特殊点没有旋转(这四个点在上述语句中发生了变化)   
		// 其他不在这些行和列的点也没变   
		// 旋转矩阵,旋转p行和q行   
		for(j = 0; j <= dim-1; j++)   
		{   
			if((j != p) && (j != q))   
			{   
				u = p*dim + j;   
				w = q*dim + j;   
				fm = symmetry_matrix[u];   
				symmetry_matrix[u] = symmetry_matrix[w]*sn + fm*cn;   
				symmetry_matrix[w] = symmetry_matrix[w]*cn - fm*sn;   
			}   
		} 
		
		// 旋转矩阵,旋转p列和q列   
		for(i = 0; i <= dim-1; i++)   
		{   
			if((i != p) && (i != q))   
			{   
				u = i*dim + p;     
				w = i*dim + q;   
				fm = symmetry_matrix[u];   
				symmetry_matrix[u] = symmetry_matrix[w]*sn + fm*cn;   
				symmetry_matrix[w] = symmetry_matrix[w]*cn - fm*sn;   
			}   
		}   

		// 记录旋转矩阵特征向量   
		for(i = 0; i <= dim-1; i++)   
		{   
			u = i*dim + p;     
			w = i*dim + q;   
			fm = eigenvectors[u];   
			eigenvectors[u] = eigenvectors[w]*sn + fm*cn;   
			eigenvectors[w] = eigenvectors[w]*cn - fm*sn;   
		}  
	}

	return true;
}