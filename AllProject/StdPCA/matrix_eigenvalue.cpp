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

	double sum_eigenvalue = 0.0;// ����ֵ�ܺ�
	map<double, int> map_eigenvalue;// ��¼ÿ������ֵ���ڵ�λ�ã������ڲ���λ��
	for (size_t index = 0; index < temp_eigenvalues.size(); index++)
	{
		map_eigenvalue[temp_eigenvalues[index]] = index;
		sum_eigenvalue += temp_eigenvalues[index];
	}

	// ��Э������������ֵ�������������򲢼����ۼưٷֱȺͰٷֱ�
	map<double, int>::reverse_iterator iter = map_eigenvalue.rbegin();
	for (size_t ii = 0; ii < temp_eigenvalues.size(); ii++)
	{
		if (iter == map_eigenvalue.rend())
			break;

		eigenvalues[ii] = iter->first;
		int index = iter->second;// ��ȡ������ֵ��Ӧ������������Ӧ���к�

		if (eigenvalues_percent != NULL)// ����ٷֱ��Լ��ۻ��ٷֱ�
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
		for (int row = 0; row < dim; row++)// ��ȡ���������о���ֵ����λ��
		{
			double temp = ABS(temp_eigenvectors(row, index));

			if (max_value < temp)
			{
				max_value = temp;
				num = row;
			}
		}

		bool is_positive = false;// �ж��������������е�ֵ�Ƿ�Ϊ��
		if (max_value - temp_eigenvectors(num, index) < 0.000000001)
			is_positive = true;

		for (int row = 0; row < dim; row++)// ȷ��ÿһ�����������ľ���ֵ���Ķ�Ϊ��
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
	for(i = 0; i <= dim-1; i++)// ��ʼ��������������ʹ��ǶԽ���Ԫ��ȫΪ0���Խ���Ԫ��Ϊ1
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
		for(i = 0; i <= dim-1; i++)// �ҳ�,����symmetry_matrix(����ֵ),�г��Խ���������Ԫ�ص�������ֵ   
		{   
			for(j = 0; j <= dim-1; j++)// ������ֵ��λ��symmetry_matrix[p][q] ,����fm   
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

		if(fm < eps)// ���ȸ���Ҫ��   
			return true;   

		if(l > jt)// ��������̫��   
			return false;

		l++;// ����������   
		u = p*dim + q;   
		w = p*dim + p;     
		t = q*dim + p;     
		s = q*dim + q;   
		x = -symmetry_matrix[u];   
		y = (symmetry_matrix[s]-symmetry_matrix[w]) / 2.0;// x y���󷨲�ͬ   
		omega = x / sqrt(x*x + y*y);//sin2��  

		//tan2��=x/y = -2.0*a[u]/(a[s]-a[w])   
		if(y < 0.0)   
			omega = -omega;

		sn = 1.0 + sqrt(1.0 - omega*omega);     
		sn = omega / sqrt(2.0*sn);//sin��   
		cn = sqrt(1.0 - sn*sn);//cos��

		fm = symmetry_matrix[w];// �任ǰ��symmetry_matrix[w]   symmetry_matrix[p][p]   
		symmetry_matrix[w] = fm*cn*cn + symmetry_matrix[s]*sn*sn + symmetry_matrix[u]*omega;   
		symmetry_matrix[s] = fm*sn*sn + symmetry_matrix[s]*cn*cn - symmetry_matrix[u]*omega;   
		symmetry_matrix[u] = 0.0;   
		symmetry_matrix[t] = 0.0;

		// ��������ת����,��ת����p��,q��,p��,q��   
		// �����ĸ������û����ת(���ĸ�������������з����˱仯)   
		// ����������Щ�к��еĵ�Ҳû��   
		// ��ת����,��תp�к�q��   
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
		
		// ��ת����,��תp�к�q��   
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

		// ��¼��ת������������   
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