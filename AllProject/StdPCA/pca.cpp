#include "pca.h"
#include "gdal_priv.h"
#include "imgalg_define.h"
#include <stdio.h>
#include <iostream>
#include "statistics_io.h"

#include "AlgProcessTime.h"

PCA::PCA(const char *src_file)
{
	m_band_mean = NULL;
	m_band_stad = NULL;
	m_covariance_or_relativity = NULL;
	m_src_dataset = NULL;
	m_statistics = NULL;
	m_dst_datatype = GDT_Float32;

	m_is_covariance = true;
	m_src_file = src_file;

	m_image_title_size = 4;
}

PCA::~PCA()
{
	if (m_src_dataset != NULL)
		GDALClose((GDALDatasetH) m_src_dataset);

	RELEASE(m_band_mean);
	RELEASE(m_band_stad);
	RELEASE(m_covariance_or_relativity);
	if (m_statistics) fclose(m_statistics);

	if (m_sta_io)
		delete m_sta_io;
}

int PCA::ExecutePCA(const char* pca_file, const char *statistics_file, int pca_band_count/* = -1*/, GDALDataType dst_type/* = 0*/, 
					bool is_covariance/* = true*/, bool is_like_envi/* = true*/, const char* format/* = "GTiff"*/)
{
	m_statistics_file = statistics_file;
	m_is_covariance = is_covariance;
	m_dst_datatype = dst_type;

	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	m_src_dataset = (GDALDataset *) GDALOpen(m_src_file, GA_ReadOnly);
	if (m_src_dataset == NULL)
		return RE_FILENOTEXIST;

	m_band_count = m_src_dataset->GetRasterCount();
	m_band_mean = new double[m_band_count];
	m_band_stad = new double[m_band_count];
	m_covariance_or_relativity = new double[m_band_count*m_band_count];	
	memset(m_band_mean, 0, sizeof(double)*m_band_count);
	memset(m_band_stad, 0, sizeof(double)*m_band_count);
	memset(m_covariance_or_relativity, 0, sizeof(double)*m_band_count*m_band_count);

	m_sta_io = new PCAStatisticsIO(m_statistics_file, m_band_count);
	if (!m_sta_io->WriteInit())
		return RE_FILENOTEXIST;

	m_statistics = fopen(m_statistics_file, "wb");
	if (m_statistics == NULL)
		return RE_FILENOTSUPPORT;

	m_src_datatype = m_src_dataset->GetRasterBand(1)->GetRasterDataType();

	bool is_float = true;
	if (m_dst_datatype == GDT_Int16 || m_dst_datatype == GDT_Int32)
		is_float = false;

	if (m_src_datatype == GDT_Byte)
	{
		if (m_dst_datatype == GDT_Int16)
			return RunPCA_New<byte, DT_16S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Int32)
			return RunPCA_New<byte, DT_32S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float32)
			return RunPCA<byte, DT_32F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float64)
			return RunPCA_New<byte, DT_64F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else
			return RE_DATATYPNOSUPPORT;
	}
	else if (m_src_datatype == GDT_UInt16)
	{
		if (m_dst_datatype == GDT_Int16)
			return RunPCA<DT_16U, DT_16S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Int32)
			return RunPCA<DT_16U, DT_32S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float32)
			return RunPCA_New<DT_16U, DT_32F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float64)
			return RunPCA<DT_16U, DT_64F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else
			return RE_DATATYPNOSUPPORT;
	}
	else if (m_src_datatype == GDT_Int16)
	{
		if (m_dst_datatype == GDT_Int16)
			return RunPCA<DT_16S, DT_16S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Int32)
			return RunPCA<DT_16S, DT_32S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float32)
			return RunPCA_New<DT_16S, DT_32F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float64)
			return RunPCA<DT_16S, DT_64F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else
			return RE_DATATYPNOSUPPORT;
	}
	else if (m_src_datatype == GDT_UInt32)
	{
		if (m_dst_datatype == GDT_Int16)
			return RunPCA<DT_32U, DT_16S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Int32)
			return RunPCA<DT_32U, DT_32S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float32)
			return RunPCA<DT_32U, DT_32F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float64)
			return RunPCA<DT_32U, DT_64F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else
			return RE_DATATYPNOSUPPORT;
	}
	else if (m_src_datatype == GDT_Int32)
	{
		if (m_dst_datatype == GDT_Int16)
			return RunPCA<DT_32S, DT_16S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Int32)
			return RunPCA<DT_32S, DT_32S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float32)
			return RunPCA<DT_32S, DT_32F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float64)
			return RunPCA<DT_32S, DT_64F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else
			return RE_DATATYPNOSUPPORT;
	}
	else if (m_src_datatype == GDT_Float32)
	{
		if (m_dst_datatype == GDT_Int16)
			return RunPCA<DT_32F, DT_16S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Int32)
			return RunPCA<DT_32F, DT_32S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float32)
			return RunPCA<DT_32F, DT_32F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float64)
			return RunPCA<DT_32F, DT_64F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else
			return RE_DATATYPNOSUPPORT;
	}
	else if (m_src_datatype == GDT_Float64)
	{
		if (m_dst_datatype == GDT_Int16)
			return RunPCA<DT_64F, DT_16S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Int32)
			return RunPCA<DT_64F, DT_32S>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float32)
			return RunPCA<DT_64F, DT_32F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else if (m_dst_datatype == GDT_Float64)
			return RunPCA<DT_64F, DT_64F>(pca_file, pca_band_count, is_float, is_like_envi, format);
		else
			return RE_DATATYPNOSUPPORT;
	}
	else
		return RE_DATATYPNOSUPPORT;
}

void PCA::CalcEigenvaluesAndEigenvectors()
{
	MyMatrix matrix = m_covar_or_relate_matrix;
	MyVector eigenvalues(m_band_count);
	MyMatrix eigenvectors(m_band_count, m_band_count);
	MyVector contribute(m_band_count);// 可输出让用户选择
	MyVector accumulate_contribute(m_band_count);// 可输出让用户选择

	//////////////////////////////////////////////////////////////////////////
	// 利用Eigen3计算特征值和特征向量
	//EigenSolver<MyMatrix> es(matrix);
	//MyMatrix D = es.pseudoEigenvalueMatrix();
	//MyMatrix V = es.pseudoEigenvectors();
	//cout << "The pseudo-eigenvalue matrix D is:" << endl << es.eigenvalues() << endl;
	//cout << "The pseudo-eigenvector matrix V is:" << endl << es.eigenvectors() << endl;
	//cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;
	//////////////////////////////////////////////////////////////////////////

	MatrixEigenvalue::GetMatrixEigen(matrix, eigenvalues, eigenvectors, &contribute, &accumulate_contribute, 0.0001);
	m_eigenvalues = eigenvalues;
	m_eigenvectors = eigenvectors;

	// 写入统计文件
	double *temp1 = new double[m_band_count*m_band_count];
	for (int i = 0; i < m_band_count; i++)
	{
		for (int j = 0; j < m_band_count; j++)
			temp1[i*m_band_count + j] = m_eigenvectors(i, j); 
	}
	m_sta_io->WriteEigenvectors(temp1);
	RELEASE(temp1);

	double *temp2 = new double[m_band_count];
	for (int i = 0; i < m_band_count; i++)
		temp2[i] = eigenvalues[i];
	m_sta_io->WriteEigenvalue(temp2);

	for (int i = 0; i < m_band_count; i++)
		temp2[i] = accumulate_contribute[i];
	
	m_sta_io->WriteAccumulateContribute(temp2);
	m_sta_io->WriteToFile();
	RELEASE(temp2);
}



int PCA::ExecuteInversePCA(const char* inverse_pca_file, const char *statistics_file, GDALDataType dst_type/* = 0*/, const char* format /* = "GTiff" */)
{
	m_sta_io = new PCAStatisticsIO(statistics_file);
	if (!m_sta_io->ReadInit())
		return RE_FILENOTEXIST;

	int dst_band_count = m_sta_io->ReadBandCount();
	m_dst_datatype = dst_type;

	// 从统计文件中获取原始影像的均值和特征向量
	double *egeinvectors = new double[dst_band_count*dst_band_count];
	m_band_mean = new double[dst_band_count];
	m_sta_io->ReadEigenvectors(egeinvectors);
	m_sta_io->ReadMean(m_band_mean);

	//Map<MyMatrix> egeinvectors_matrix(egeinvectors, dst_band_count, dst_band_count);
	MyExtMatrix egeinvectors_matrix(egeinvectors, dst_band_count, dst_band_count);
	MyMatrix inverse_egeinvectors(dst_band_count, dst_band_count);

	inverse_egeinvectors = egeinvectors_matrix.inverse();//需要注意求逆矩阵时的条件
	m_select_eigenvectors = inverse_egeinvectors;
	
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
	m_src_dataset = (GDALDataset *) GDALOpen(m_src_file, GA_ReadOnly);
	if (m_src_dataset == NULL)
		return RE_FILENOTEXIST;

	m_band_count = m_src_dataset->GetRasterCount();
	m_src_datatype = m_src_dataset->GetRasterBand(1)->GetRasterDataType();

	bool is_float = false;
	if (m_dst_datatype == GDT_Float32 || m_dst_datatype == GDT_Float64)
		is_float = true;

	if (m_src_datatype == GDT_Int16)
	{
		if (m_dst_datatype == GDT_Byte)
			return RunInversePCA<DT_16S, byte>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_UInt16)
			return RunInversePCA<DT_16S, DT_16U>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Int16)
			return RunInversePCA<DT_16S, DT_16S>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_UInt32)
			return RunInversePCA<DT_16S, DT_32U>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Int32)
			return RunInversePCA<DT_16S, DT_32S>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Float32)
			return RunInversePCA<DT_16S, DT_32F>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Float64)
			return RunInversePCA<DT_16S, DT_64F>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else
			return RE_DATATYPNOSUPPORT;
	}
	else if (m_src_datatype == GDT_Int32)
	{
		if (m_dst_datatype == GDT_Byte)
			return RunInversePCA<DT_32S, byte>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_UInt16)
			return RunInversePCA<DT_32S, DT_16U>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Int16)
			return RunInversePCA<DT_32S, DT_16S>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_UInt32)
			return RunInversePCA<DT_32S, DT_32U>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Int32)
			return RunInversePCA<DT_32S, DT_32S>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Float32)
			return RunInversePCA<DT_32S, DT_32F>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Float64)
			return RunInversePCA<DT_32S, DT_64F>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else
			return RE_DATATYPNOSUPPORT;
	}
	else if (m_src_datatype == GDT_Float32)
	{
		if (m_dst_datatype == GDT_Byte)
			return RunInversePCA<DT_32F, byte>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_UInt16)
			return RunInversePCA<DT_32F, DT_16U>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Int16)
			return RunInversePCA<DT_32F, DT_16S>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_UInt32)
			return RunInversePCA<DT_32F, DT_32U>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Int32)
			return RunInversePCA<DT_32F, DT_32S>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Float32)
			return RunInversePCA<DT_32F, DT_32F>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Float64)
			return RunInversePCA<DT_32F, DT_64F>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else
			return RE_DATATYPNOSUPPORT;
	}
	else if (m_src_datatype == GDT_Float64)
	{
		if (m_dst_datatype == GDT_Byte)
			return RunInversePCA<DT_64F, byte>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_UInt16)
			return RunInversePCA<DT_64F, DT_16U>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Int16)
			return RunInversePCA<DT_64F, DT_16S>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_UInt32)
			return RunInversePCA<DT_64F, DT_32U>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Int32)
			return RunInversePCA<DT_64F, DT_32S>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Float32)
			return RunInversePCA<DT_64F, DT_32F>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else if (m_dst_datatype == GDT_Float64)
			return RunInversePCA<DT_64F, DT_64F>(inverse_pca_file, inverse_egeinvectors, is_float, format);
		else
			return RE_DATATYPNOSUPPORT;
	}
	else
		return RE_DATATYPNOSUPPORT;
}

template <typename T1, typename T2>
int PCA::RunPCA(const char *pca_file, int pca_band_count, bool is_float, bool is_like_envi, const char* format)
{
	int w = m_src_dataset->GetRasterXSize();
	int h = m_src_dataset->GetRasterYSize();
	int band_size = w*h;

	// 分块信息
	int block_unit = w*m_band_count;									// 分块基本单位（单位：像素）
	int block_h = m_image_title_size * Mb / (block_unit*sizeof(T1));	// 每一块的高度（单位：像素）
	int block_nums = h / block_h;										// 能分成的整块数n
	int last_block_h = h % block_h;										// 剩余最后一块的高度（单位：像素）
	if (last_block_h != 0) block_nums++;								// 如果未能整分，则最后剩余的行数也分块一块，即总块数为n+1。

	int *band_map = new int[m_band_count];
	for (int i = 0; i < m_band_count; i++)
		band_map[i] = i + 1;

	CAlgProcessTime::Alg_start();

	/************************************************************************/
	/*						统计均值和标准差                                                                     */
	/************************************************************************/

	// 统计均值
	int read_h = block_h;
	T1 *block_buf_data = new T1[read_h*block_unit];

	for (int i = 0; i < block_nums; i++)
	{
		if (i == block_nums - 1)
		{
			read_h = (h - 1)%block_h + 1;
			RELEASE(block_buf_data);
			block_buf_data = new T1[read_h*block_unit];
		}

		// 按BIP的格式读取图像
		m_src_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_src_datatype,
			m_band_count, band_map, sizeof(T1)*m_band_count, sizeof(T1)*m_band_count*w, sizeof(T1));

		int band_block_size = read_h*w;
		for (int b = 0; b < m_band_count; b++)
		{			
			for (int j = 0; j < band_block_size; j++)			
				m_band_mean[b] += block_buf_data[j*m_band_count + b];
		}
	}
	RELEASE(block_buf_data);

	for (int i = 0; i < m_band_count; i++)
		m_band_mean[i] /= band_size; 

	// 统计标准差
	read_h = block_h;
	block_buf_data = new T1[read_h*block_unit];

	for (int i = 0; i < block_nums; i++)
	{
		if (i == block_nums - 1)
		{
			read_h = (h - 1)%block_h + 1;
			RELEASE(block_buf_data);
			block_buf_data = new T1[read_h*block_unit];
		}

		m_src_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_src_datatype,
			m_band_count, band_map, sizeof(T1)*m_band_count, sizeof(T1)*m_band_count*w, sizeof(T1));

		int band_block_size = read_h*w;
		for (int b = 0; b < m_band_count; b++)
		{			
			for (int j = 0; j < band_block_size; j++)			
			{
				double temp = block_buf_data[j*m_band_count + b] - m_band_mean[b];
				m_band_stad[b] += temp*temp;
			}
		}
	}
	RELEASE(block_buf_data);

	for (int i = 0; i < m_band_count; i++)
		m_band_stad[i] = sqrt(m_band_stad[i] / band_size);

	/************************************************************************/
	/*						计算协方差或相关系数矩阵                   */
	/************************************************************************/
	read_h = block_h;
	block_buf_data = new T1[read_h*block_unit];
	
	int element_index;
	for (int i = 0; i < block_nums; i++)
	{
		if (i == block_nums - 1)
		{
			read_h = (h - 1)%block_h + 1;
			RELEASE(block_buf_data);
			block_buf_data = new T1[read_h*block_unit];
		}

		m_src_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_src_datatype,
			m_band_count, band_map, sizeof(T1)*m_band_count, sizeof(T1)*m_band_count*w, sizeof(T1));

		int band_block_size = read_h*w;		
		element_index = 0;
		for (int b1 = 0; b1 < m_band_count; b1++)
		{
			for(int b2 = 0; b2 < m_band_count; b2++)
			{
				if (b2 < b1)
				{
					element_index++;
					continue;
				}

				if (b1 == b2)
				{
					if (m_is_covariance)
						m_covariance_or_relativity[element_index] = m_band_stad[b1] * m_band_stad[b1];
					else// 相关系数
						m_covariance_or_relativity[element_index] = 1.0;						

					element_index++;
					continue;
				}

				for (int k = 0; k < band_block_size; k++)
				{
					int offset = k*m_band_count;
					m_covariance_or_relativity[b1*m_band_count + b2] += 
						(block_buf_data[offset+b1] - m_band_mean[b1])*(block_buf_data[offset+b2] - m_band_mean[b2]);
				}

				element_index++;
			}
		}
	}
	RELEASE(block_buf_data);
	RELEASE(band_map);

	// 计算协方差矩阵或相关系数矩阵
	for(int r = 0; r < m_band_count; r++)
	{
		for(int c = 0; c < m_band_count; c++)
		{
			int index = r*m_band_count + c;

			if (c < r)
			{
				m_covariance_or_relativity[index] = m_covariance_or_relativity[c*m_band_count + r];
				continue;
			}

			if (r == c)
				continue;

			m_covariance_or_relativity[index] /= band_size;				
			if (!m_is_covariance)// 相关系数
				m_covariance_or_relativity[index] = m_covariance_or_relativity[index] / (m_band_stad[r]*m_band_stad[c]);
		}
	}
	m_sta_io->WriteMean(m_band_mean);

	CAlgProcessTime::Alg_end("统计均值和协方差矩阵");

	/************************************************************************/
	/*						计算特征值和特征向量                                                                     */
	/************************************************************************/
	CAlgProcessTime::Alg_start();

	//Map<MyMatrix> covariance_matrix(m_covariance_or_relativity, m_band_count, m_band_count);
	MyExtMatrix covariance_matrix(m_covariance_or_relativity, m_band_count, m_band_count);
	m_covar_or_relate_matrix = covariance_matrix;
	m_sta_io->WriteCovarianceOrCorrelation(m_covariance_or_relativity);
	CalcEigenvaluesAndEigenvectors();

	CAlgProcessTime::Alg_end("计算特征值和特征向量矩阵");

	/************************************************************************/
	/*						原始数据进行PCA变换                                                                     */
	/************************************************************************/
	CAlgProcessTime::Alg_start();

	int new_band_count = m_band_count;
	if (pca_band_count > 0)
		new_band_count = pca_band_count;

	MyMatrix select_eigenvectors(m_band_count, new_band_count);
	for (int i = 0; i < new_band_count; i++)
	{
		for (int j = 0; j < m_band_count; j++)
			select_eigenvectors(j, i) = m_eigenvectors(j, i);
	}

	m_select_eigenvectors = select_eigenvectors;
	int dst_band_count = new_band_count;	

	// 设置存储格式与原始影像一样（BSQ，BIP，BIL）
	char **dst_gdal_options = NULL;
	const char *interleave = m_src_dataset->GetMetadataItem("INTERLEAVE", "IMAGE_STRUCTURE");// 获取数据存储类型(BAND-BSQ,PIXEL-BIP,LINE-BIL)
	dst_gdal_options = CSLSetNameValue(dst_gdal_options, "INTERLEAVE", interleave);

	GDALDriver *dst_diver = (GDALDriver *)GDALGetDriverByName(format);
	GDALDataset *dst_dataset = dst_diver->Create(pca_file, w, h, dst_band_count, m_dst_datatype, dst_gdal_options);
	if ( dst_dataset == NULL )
		return RE_FILENOTSUPPORT;

	double geo_transform[6] = { 0 };
	m_src_dataset->GetGeoTransform(geo_transform);
	dst_dataset->SetGeoTransform(geo_transform);
	dst_dataset->SetProjection(m_src_dataset->GetProjectionRef());

	//////////////////////////////////////////////////////////////////////////
	read_h = block_h;
	block_buf_data = new T1[read_h*block_unit];	
	T2 *dst_buffer_data = new T2[read_h*w*dst_band_count];

	int *src_band_map = new int[m_band_count];
	int *dst_band_map = new int[dst_band_count];

	for (int i = 1; i <= m_band_count; i++)
		src_band_map[i - 1] = i;
	for (int i = 1; i <= dst_band_count; i++)
		dst_band_map[i - 1] = i;

	for (int i = 0; i < block_nums; i++)
	{
		if (i == block_nums - 1)
		{
			read_h = (h - 1)%block_h + 1;
			RELEASE(block_buf_data);
			RELEASE(dst_buffer_data);
			block_buf_data = new T1[read_h*block_unit];
			dst_buffer_data = new T2[read_h*w*dst_band_count];
		}

		m_src_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_src_datatype,
			m_band_count, src_band_map, sizeof(T1)*m_band_count, sizeof(T1)*m_band_count*w, sizeof(T1));

		int band_sample = w*read_h;		
		for (int j = 0; j < band_sample; j++)
		{
			int offset = j*m_band_count;
			int offset_out = j*dst_band_count;
			
			for (int k1 = 0; k1 < dst_band_count; k1++)
			{
				double temp = 0.0;
				
				for (int k2 = 0; k2 < m_band_count; k2++)
				{			
					temp += (block_buf_data[offset + k2] - m_band_mean[k2])* select_eigenvectors(k2, k1);			
				}

				dst_buffer_data[offset_out + k1] = is_float ? (T2)temp : (T2)floor(temp + 0.5);
			}
		}		

		dst_dataset->RasterIO(GF_Write, 0, i*block_h, w, read_h, dst_buffer_data, w, read_h, m_dst_datatype,
			dst_band_count, dst_band_map, sizeof(T2)*dst_band_count, sizeof(T2)*dst_band_count*w, sizeof(T2));
	}
	RELEASE(block_buf_data);
	RELEASE(dst_buffer_data);
	RELEASE(src_band_map);
	RELEASE(dst_band_map);

	GDALClose((GDALDatasetH) dst_dataset);

	CAlgProcessTime::Alg_end("原始数据进行PCA变换");
	//////////////////////////////////////////////////////////////////////////

	/************************************************************************/
	/*						将结果减去均值，使结果和ENVI一样                                                                     */
	/************************************************************************/
	if (!is_like_envi)
	{
		CAlgProcessTime::Alg_start();

		GDALDataset *pca_dataset = (GDALDataset *) GDALOpen(pca_file, GA_Update);
		if (pca_dataset == NULL)
			return RE_FILENOTSUPPORT;

		// 统计均值
		int band_count = dst_band_count;
		block_unit = w*band_count;
		block_h = m_image_title_size * Mb / (block_unit*sizeof(T2));
		block_nums = h / block_h;
		last_block_h = h % block_h;
		if (last_block_h != 0) block_nums++;

		// 读取的波段
		band_map = new int[band_count];
		for (int i = 0; i < band_count; i++)
			band_map[i] = i + 1;

		double *band_mean = new double[band_count];
		memset(band_mean, 0, sizeof(double)*band_count);

		read_h = block_h;
		T2 *block_buf_data = new T2[read_h*block_unit];	
		for (int i = 0; i < block_nums; i++)
		{
			if (i == block_nums - 1)
			{
				read_h = (h - 1)%block_h + 1;
				RELEASE(block_buf_data);
				block_buf_data = new T2[read_h*block_unit];
			}

			// 按BIP的格式读取图像
			pca_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_dst_datatype,
				band_count, band_map, sizeof(T2)*band_count, sizeof(T2)*band_count*w, sizeof(T2));

			int band_block_size = read_h*w;
			for (int b = 0; b < band_count; b++)
			{			
				for (int j = 0; j < band_block_size; j++)			
					band_mean[b] += block_buf_data[j*band_count + b];
			}
		}
		RELEASE(block_buf_data);

		for (int i = 0; i < band_count; i++)
			band_mean[i] /= band_size;

		// 减均值
		read_h = block_h;
		block_buf_data = new T2[read_h*block_unit];	
		for (int i = 0; i < block_nums; i++)
		{
			if (i == block_nums - 1)
			{
				read_h = (h - 1)%block_h + 1;
				RELEASE(block_buf_data);
				block_buf_data = new T2[read_h*block_unit];
			}

			pca_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_dst_datatype,
				band_count, band_map, sizeof(T2)*band_count, sizeof(T2)*band_count*w, sizeof(T2));

			int band_block_size = read_h*w;
			for (int b = 0; b < band_count; b++)
			{
				for (int j = 0; j < band_block_size; j++)			
					block_buf_data[j*band_count + b] = is_float ? (T2)(block_buf_data[j*band_count + b] - band_mean[b]) : (T2)floor(block_buf_data[j*band_count + b] - band_mean[b] + 0.5) ;
			}

			pca_dataset->RasterIO(GF_Write, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_dst_datatype,
				band_count, band_map, sizeof(T2)*band_count, sizeof(T2)*band_count*w, sizeof(T2));
		}
		RELEASE(band_map);
		RELEASE(block_buf_data);
		GDALClose((GDALDatasetH) pca_dataset);

		CAlgProcessTime::Alg_end("减去均值（ENVI）");
	}// is_like_envi

	return RE_SUCCESS;
}

template <typename T1, typename T2>
int PCA::RunInversePCA(const char *inverse_pca_file, MyMatrix &inverse_egeinvectors, bool is_float, const char* format)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
	int w = m_src_dataset->GetRasterXSize();
	int h = m_src_dataset->GetRasterYSize();
	int dst_band_count = inverse_egeinvectors.cols();

	char **dst_gdal_options = NULL;
	const char *interleave = m_src_dataset->GetMetadataItem("INTERLEAVE", "IMAGE_STRUCTURE");// 获取数据存储类型(BAND-BSQ,PIXEL-BIP,LINE-BIL)
	dst_gdal_options = CSLSetNameValue(dst_gdal_options, "INTERLEAVE", interleave);

	GDALDriver *dst_diver = (GDALDriver *)GDALGetDriverByName(format);
	GDALDataset *dst_dataset = dst_diver->Create(inverse_pca_file, w, h, dst_band_count, m_dst_datatype, dst_gdal_options);
	if ( dst_dataset == NULL )
		return RE_FILENOTSUPPORT;

	double geo_transform[6] = { 0 };
	m_src_dataset->GetGeoTransform(geo_transform);
	dst_dataset->SetGeoTransform(geo_transform);
	dst_dataset->SetProjection(m_src_dataset->GetProjectionRef());

	//////////////////////////////////////////////////////////////////////////
	int block_unit = w*m_band_count;									
	int block_h = m_image_title_size * Mb / (block_unit*sizeof(T1));	
	int block_nums = h / block_h;										
	int last_block_h = h % block_h;										
	if (last_block_h != 0) block_nums++;

	int read_h = block_h;
	T1 *block_buf_data = new T1[read_h*block_unit];	
	T2 *dst_buffer_data = new T2[read_h*w*dst_band_count];

	int *src_band_map = new int[m_band_count];
	int *dst_band_map = new int[dst_band_count];

	for (int i = 1; i <= m_band_count; i++)
		src_band_map[i - 1] = i;
	for (int i = 1; i <= dst_band_count; i++)
		dst_band_map[i - 1] = i;

	for (int i = 0; i < block_nums; i++)
	{
		if (i == block_nums - 1)
		{
			read_h = (h - 1)%block_h + 1;
			
			RELEASE(block_buf_data);
			RELEASE(dst_buffer_data);
			block_buf_data = new T1[read_h*block_unit];
			dst_buffer_data = new T2[read_h*w*dst_band_count];
		}

		// 将读取的数据按BIP格式存储在内存中
		m_src_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_src_datatype,
			m_band_count, src_band_map, sizeof(T1)*m_band_count, sizeof(T1)*m_band_count*w, sizeof(T1));

		int band_sample = w*read_h;		
		for (int j = 0; j < band_sample; j++)
		{
			int offset = j*m_band_count;
			int offset_out = j*dst_band_count;

			for (int k1 = 0; k1 < dst_band_count; k1++)
			{
				double temp = 0.0;

				for (int k2 = 0; k2 < m_band_count; k2++)			
					temp += block_buf_data[offset + k2] * inverse_egeinvectors(k2, k1);			

				dst_buffer_data[offset_out + k1] = is_float ? (T2)(temp + m_band_mean[k1]) : (T2)floor(temp + m_band_mean[k1] + 0.5);
			}
		}		

		dst_dataset->RasterIO(GF_Write, 0, i*block_h, w, read_h, dst_buffer_data, w, read_h, m_dst_datatype,
			dst_band_count, dst_band_map, sizeof(T2)*dst_band_count, sizeof(T2)*dst_band_count*w, sizeof(T2));
	}
	RELEASE(block_buf_data);
	RELEASE(dst_buffer_data);
	RELEASE(src_band_map);
	RELEASE(dst_band_map);

	GDALClose((GDALDatasetH) dst_dataset);

	return RE_SUCCESS;
}

// 实验
template <typename T1, typename T2>
int PCA::RunPCA_New(const char *pca_file, int pca_band_count, bool is_float, bool is_like_envi, const char* format)
{
	int w = m_src_dataset->GetRasterXSize();
	int h = m_src_dataset->GetRasterYSize();
	int band_size = w*h;

	// 分块信息
	int block_unit = w*m_band_count;									// 分块基本单位（单位：像素）
	int block_h = m_image_title_size * Mb / (block_unit*sizeof(T1));	// 每一块的高度（单位：像素）
	int block_nums = h / block_h;										// 能分成的整块数n
	int last_block_h = h % block_h;										// 剩余最后一块的高度（单位：像素）
	if (last_block_h != 0) block_nums++;								// 如果未能整分，则最后剩余的行数也分块一块，即总块数为n+1。

	int *band_map = new int[m_band_count];
	for (int i = 0; i < m_band_count; i++)
		band_map[i] = i + 1;

	CAlgProcessTime::Alg_start();

	/************************************************************************/
	/*						统计均值和标准差                                                                     */
	/************************************************************************/

	// 统计均值
	int read_h = block_h;
	T1 *block_buf_data = new T1[read_h*block_unit];

	for (int i = 0; i < block_nums; i++)
	{
		if (i == block_nums - 1)
		{
			read_h = (h - 1)%block_h + 1;
			RELEASE(block_buf_data);
			block_buf_data = new T1[read_h*block_unit];
		}

		// 按BIP的格式读取图像
		m_src_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_src_datatype,
			m_band_count, band_map, sizeof(T1)*m_band_count, sizeof(T1)*m_band_count*w, sizeof(T1));

		int band_block_size = read_h*w;
		for (int b = 0; b < m_band_count; b++)
		{			
			for (int j = 0; j < band_block_size; j++)		
			{
				double temp = block_buf_data[j*m_band_count + b];
				// sum of X：x1 + x2 + x3 + ... + xN
				m_band_mean[b] += temp;

				// sum of X*X：x1*x1 + x2*x2 + x3*x3 + ...
				m_band_stad[b] += temp*temp;


				// sum of X*Y：x1*y1 + x2*y2 + x3*y3 + ...
				int index1 = b*m_band_count;
				int index2 = j*m_band_count;
				for (int b1 = 0; b1 < m_band_count; b1++)
				{
					m_covariance_or_relativity[index1 + b1] += temp*block_buf_data[index2 + b1];
				}
			}			
		}
	}
	RELEASE(block_buf_data);
	RELEASE(band_map);

	// 平均值
	for (int i = 0; i < m_band_count; i++)
		m_band_mean[i] /= band_size;

	// 标准差
	for (int i = 0; i < m_band_count; i++)
	{
		m_band_stad[i] /= band_size;
		m_band_stad[i] = sqrt(m_band_stad[i] - m_band_mean[i]*m_band_mean[i]);
	}

	// 协方差矩阵
	for (int i = 0; i < m_band_count; i++)
	{
		for (int j = 0; j < m_band_count; j++)
		{
			m_covariance_or_relativity[i*m_band_count + j] /= band_size;
			m_covariance_or_relativity[i*m_band_count + j] -= m_band_mean[i]*m_band_mean[j];

			if (!m_is_covariance)// 相关系数
			{
				if (i == j) 
				{
					m_covariance_or_relativity[i*m_band_count + j] = 1.0;
					continue;
				}

				m_covariance_or_relativity[i*m_band_count + j] /= (m_band_stad[i]*m_band_stad[j]);
			}
		}
	}

	//// 统计标准差
	//read_h = block_h;
	//block_buf_data = new T1[read_h*block_unit];

	//for (int i = 0; i < block_nums; i++)
	//{
	//	if (i == block_nums - 1)
	//	{
	//		read_h = (h - 1)%block_h + 1;
	//		RELEASE(block_buf_data);
	//		block_buf_data = new T1[read_h*block_unit];
	//	}

	//	m_src_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_src_datatype,
	//		m_band_count, band_map, sizeof(T1)*m_band_count, sizeof(T1)*m_band_count*w, sizeof(T1));

	//	int band_block_size = read_h*w;
	//	for (int b = 0; b < m_band_count; b++)
	//	{			
	//		for (int j = 0; j < band_block_size; j++)			
	//		{
	//			double temp = block_buf_data[j*m_band_count + b] - m_band_mean[b];
	//			m_band_stad[b] += temp*temp;
	//		}
	//	}
	//}
	//RELEASE(block_buf_data);

	//for (int i = 0; i < m_band_count; i++)
	//	m_band_stad[i] = sqrt(m_band_stad[i] / band_size);

	///************************************************************************/
	///*						计算协方差或相关系数矩阵                   */
	///************************************************************************/
	//read_h = block_h;
	//block_buf_data = new T1[read_h*block_unit];

	//int element_index;
	//for (int i = 0; i < block_nums; i++)
	//{
	//	if (i == block_nums - 1)
	//	{
	//		read_h = (h - 1)%block_h + 1;
	//		RELEASE(block_buf_data);
	//		block_buf_data = new T1[read_h*block_unit];
	//	}

	//	m_src_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_src_datatype,
	//		m_band_count, band_map, sizeof(T1)*m_band_count, sizeof(T1)*m_band_count*w, sizeof(T1));

	//	int band_block_size = read_h*w;		
	//	element_index = 0;
	//	for (int b1 = 0; b1 < m_band_count; b1++)
	//	{
	//		for(int b2 = 0; b2 < m_band_count; b2++)
	//		{
	//			if (b2 < b1)
	//			{
	//				element_index++;
	//				continue;
	//			}

	//			if (b1 == b2)
	//			{
	//				if (m_is_covariance)
	//					m_covariance_or_relativity[element_index] = m_band_stad[b1] * m_band_stad[b1];
	//				else// 相关系数
	//					m_covariance_or_relativity[element_index] = 1.0;						

	//				element_index++;
	//				continue;
	//			}

	//			for (int k = 0; k < band_block_size; k++)
	//			{
	//				int offset = k*m_band_count;
	//				m_covariance_or_relativity[b1*m_band_count + b2] += 
	//					(block_buf_data[offset+b1] - m_band_mean[b1])*(block_buf_data[offset+b2] - m_band_mean[b2]);
	//			}

	//			element_index++;
	//		}
	//	}
	//}
	//RELEASE(block_buf_data);
	//RELEASE(band_map);

	//// 计算协方差矩阵或相关系数矩阵
	//for(int r = 0; r < m_band_count; r++)
	//{
	//	for(int c = 0; c < m_band_count; c++)
	//	{
	//		int index = r*m_band_count + c;

	//		if (c < r)
	//		{
	//			m_covariance_or_relativity[index] = m_covariance_or_relativity[c*m_band_count + r];
	//			continue;
	//		}

	//		if (r == c)
	//			continue;

	//		m_covariance_or_relativity[index] /= band_size;				
	//		if (!m_is_covariance)// 相关系数
	//			m_covariance_or_relativity[index] = m_covariance_or_relativity[index] / (m_band_stad[r]*m_band_stad[c]);
	//	}
	//}
	m_sta_io->WriteMean(m_band_mean);

	CAlgProcessTime::Alg_end("统计均值和协方差矩阵");

	/************************************************************************/
	/*						计算特征值和特征向量                                                                     */
	/************************************************************************/
	CAlgProcessTime::Alg_start();

	//Map<MyMatrix> covariance_matrix(m_covariance_or_relativity, m_band_count, m_band_count);
	MyExtMatrix covariance_matrix(m_covariance_or_relativity, m_band_count, m_band_count);
	m_covar_or_relate_matrix = covariance_matrix;
	m_sta_io->WriteCovarianceOrCorrelation(m_covariance_or_relativity);
	CalcEigenvaluesAndEigenvectors();

	CAlgProcessTime::Alg_end("计算特征值和特征向量矩阵");

	/************************************************************************/
	/*						原始数据进行PCA变换                                                                     */
	/************************************************************************/
	CAlgProcessTime::Alg_start();

	int new_band_count = m_band_count;
	if (pca_band_count > 0)
		new_band_count = pca_band_count;

	MyMatrix select_eigenvectors(m_band_count, new_band_count);
	for (int i = 0; i < new_band_count; i++)
	{
		for (int j = 0; j < m_band_count; j++)
			select_eigenvectors(j, i) = m_eigenvectors(j, i);
	}

	m_select_eigenvectors = select_eigenvectors;
	int dst_band_count = new_band_count;	

	// 设置存储格式与原始影像一样（BSQ，BIP，BIL）
	char **dst_gdal_options = NULL;
	const char *interleave = m_src_dataset->GetMetadataItem("INTERLEAVE", "IMAGE_STRUCTURE");// 获取数据存储类型(BAND-BSQ,PIXEL-BIP,LINE-BIL)
	dst_gdal_options = CSLSetNameValue(dst_gdal_options, "INTERLEAVE", interleave);

	GDALDriver *dst_diver = (GDALDriver *)GDALGetDriverByName(format);
	GDALDataset *dst_dataset = dst_diver->Create(pca_file, w, h, dst_band_count, m_dst_datatype, dst_gdal_options);
	if ( dst_dataset == NULL )
		return RE_FILENOTSUPPORT;

	double geo_transform[6] = { 0 };
	m_src_dataset->GetGeoTransform(geo_transform);
	dst_dataset->SetGeoTransform(geo_transform);
	dst_dataset->SetProjection(m_src_dataset->GetProjectionRef());

	//////////////////////////////////////////////////////////////////////////
	read_h = block_h;
	block_buf_data = new T1[read_h*block_unit];	
	T2 *dst_buffer_data = new T2[read_h*w*dst_band_count];

	int *src_band_map = new int[m_band_count];
	int *dst_band_map = new int[dst_band_count];

	for (int i = 1; i <= m_band_count; i++)
		src_band_map[i - 1] = i;
	for (int i = 1; i <= dst_band_count; i++)
		dst_band_map[i - 1] = i;

	for (int i = 0; i < block_nums; i++)
	{
		if (i == block_nums - 1)
		{
			read_h = (h - 1)%block_h + 1;
			RELEASE(block_buf_data);
			RELEASE(dst_buffer_data);
			block_buf_data = new T1[read_h*block_unit];
			dst_buffer_data = new T2[read_h*w*dst_band_count];
		}

		m_src_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_src_datatype,
			m_band_count, src_band_map, sizeof(T1)*m_band_count, sizeof(T1)*m_band_count*w, sizeof(T1));

		int band_sample = w*read_h;		
		for (int j = 0; j < band_sample; j++)
		{
			int offset = j*m_band_count;
			int offset_out = j*dst_band_count;

			for (int k1 = 0; k1 < dst_band_count; k1++)
			{
				double temp = 0.0;

				for (int k2 = 0; k2 < m_band_count; k2++)
				{			
					temp += (block_buf_data[offset + k2]  - m_band_mean[k2]) * select_eigenvectors(k2, k1);			
				}

				dst_buffer_data[offset_out + k1] = is_float ? (T2)temp : (T2)floor(temp + 0.5);
			}
		}		

		dst_dataset->RasterIO(GF_Write, 0, i*block_h, w, read_h, dst_buffer_data, w, read_h, m_dst_datatype,
			dst_band_count, dst_band_map, sizeof(T2)*dst_band_count, sizeof(T2)*dst_band_count*w, sizeof(T2));
	}
	RELEASE(block_buf_data);
	RELEASE(dst_buffer_data);
	RELEASE(src_band_map);
	RELEASE(dst_band_map);

	GDALClose((GDALDatasetH) dst_dataset);

	CAlgProcessTime::Alg_end("原始数据进行PCA变换");
	//////////////////////////////////////////////////////////////////////////

	/************************************************************************/
	/*						将结果减去均值，使结果和ENVI一样                                                                     */
	/************************************************************************/
	if (!is_like_envi)
	{
		CAlgProcessTime::Alg_start();

		GDALDataset *pca_dataset = (GDALDataset *) GDALOpen(pca_file, GA_Update);
		if (pca_dataset == NULL)
			return RE_FILENOTSUPPORT;

		// 统计均值
		int band_count = dst_band_count;
		block_unit = w*band_count;
		block_h = m_image_title_size * Mb / (block_unit*sizeof(T2));
		block_nums = h / block_h;
		last_block_h = h % block_h;
		if (last_block_h != 0) block_nums++;

		// 读取的波段
		band_map = new int[band_count];
		for (int i = 0; i < band_count; i++)
			band_map[i] = i + 1;

		double *band_mean = new double[band_count];
		memset(band_mean, 0, sizeof(double)*band_count);

		read_h = block_h;
		T2 *block_buf_data = new T2[read_h*block_unit];	
		for (int i = 0; i < block_nums; i++)
		{
			if (i == block_nums - 1)
			{
				read_h = (h - 1)%block_h + 1;
				RELEASE(block_buf_data);
				block_buf_data = new T2[read_h*block_unit];
			}

			// 按BIP的格式读取图像
			pca_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_dst_datatype,
				band_count, band_map, sizeof(T2)*band_count, sizeof(T2)*band_count*w, sizeof(T2));

			int band_block_size = read_h*w;
			for (int b = 0; b < band_count; b++)
			{			
				for (int j = 0; j < band_block_size; j++)			
					band_mean[b] += block_buf_data[j*band_count + b];
			}
		}
		RELEASE(block_buf_data);

		for (int i = 0; i < band_count; i++)
			band_mean[i] /= band_size;

		// 减均值
		read_h = block_h;
		block_buf_data = new T2[read_h*block_unit];	
		for (int i = 0; i < block_nums; i++)
		{
			if (i == block_nums - 1)
			{
				read_h = (h - 1)%block_h + 1;
				RELEASE(block_buf_data);
				block_buf_data = new T2[read_h*block_unit];
			}

			pca_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_dst_datatype,
				band_count, band_map, sizeof(T2)*band_count, sizeof(T2)*band_count*w, sizeof(T2));

			int band_block_size = read_h*w;
			for (int b = 0; b < band_count; b++)
			{
				for (int j = 0; j < band_block_size; j++)			
					block_buf_data[j*band_count + b] = is_float ? (T2)(block_buf_data[j*band_count + b] - band_mean[b]) : (T2)floor(block_buf_data[j*band_count + b] - band_mean[b] + 0.5) ;
			}

			pca_dataset->RasterIO(GF_Write, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_dst_datatype,
				band_count, band_map, sizeof(T2)*band_count, sizeof(T2)*band_count*w, sizeof(T2));
		}
		RELEASE(band_map);
		RELEASE(block_buf_data);
		GDALClose((GDALDatasetH) pca_dataset);

		CAlgProcessTime::Alg_end("减去均值（ENVI）");
	}// is_like_envi

	return RE_SUCCESS;
}