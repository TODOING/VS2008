#include "pca.h"
#include "gdal_priv.h"
#include "imgalg_define.h"
#include <stdio.h>
#include <iostream>
#include "statistics_io.h"

PCA::PCA(const char *src_file)
{
	m_band_mean = NULL;
	m_band_stad = NULL;
	m_relativity = NULL;
	m_src_dataset = NULL;
	m_statistics = NULL;

	m_is_covariance = true;
	m_src_file = src_file;
}

PCA::~PCA()
{
	if (m_src_dataset != NULL)
		GDALClose((GDALDatasetH) m_src_dataset);

	RELEASE(m_band_mean);
	RELEASE(m_band_stad);
	RELEASE(m_relativity);
	if (m_statistics) fclose(m_statistics);

	if (m_sta_io)
		delete m_sta_io;
}

int PCA::ExecutePCA(const char* pca_file, int pca_band_count /* = -1 */, bool is_covariance /* = true */, bool is_like_envi /* = true */, const char* format /* = "GTiff" */)
{
	m_is_covariance = is_covariance;

	// 第一步，数据预处理，统计每一波段的均值和标准差
	int return_value = PreProcessData();
	if (return_value != RE_SUCCESS)
		return return_value;

	// 第二步，计算协方差矩阵和相关系数矩阵R，以及内部求特征值和特征向量
	return_value = CalcCovarianceMartix();
	if (return_value != RE_SUCCESS)
		return return_value;

	// 计算主成分得分，并写入到文件中
	return_value = CreatePCAFile(pca_file, pca_band_count, format);
	if (return_value != RE_SUCCESS)
		return return_value;

	if (is_like_envi)
	{
		return_value = CalcSubAvg(pca_file);
		if (return_value != RE_SUCCESS)
			return return_value;
	}

	return RE_SUCCESS;
}

int PCA::ExecutePCA(const char* pca_file, const char *statistics_file, int pca_band_count/* = -1*/, bool is_covariance/* = true*/, 
			   bool is_like_envi/* = true*/, const char* format/* = "GTiff"*/)
{
	m_statistics_file = statistics_file;
	return ExecutePCA(pca_file, pca_band_count, is_covariance, is_like_envi, format);
}

int PCA::PreProcessData()
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	m_src_dataset = (GDALDataset *) GDALOpen(m_src_file, GA_ReadOnly);
	if (m_src_dataset == NULL)
		return RE_FILENOTEXIST;

	m_band_count = m_src_dataset->GetRasterCount();
	m_band_mean = new double[m_band_count];
	m_band_stad = new double[m_band_count];

	m_sta_io = new PCAStatisticsIO(m_statistics_file, m_band_count);
	m_sta_io->WriteInit();

	m_statistics = fopen(m_statistics_file, "w");
	if (m_statistics == NULL)
		return RE_FILENOTSUPPORT;

	double max_value, min_value;
	//fprintf(m_statistics, "band			min			max			mean		stddev\n");
	for (int i = 1; i <= m_band_count; i++)
	{
		m_src_dataset->GetRasterBand(i)->ComputeStatistics(FALSE, &min_value, &max_value, 
			m_band_mean + (i - 1), m_band_stad + (i - 1), NULL, NULL);
		
		//fprintf(m_statistics, "%-6d %-16.6f %-16.6f %-16.6f %-16.6f\n", i, min_value, max_value, m_band_mean[i-1], m_band_stad[i-1]);
	}
	m_sta_io->WriteMean(m_band_mean);
	//fprintf(m_statistics, "\n\n");

	return RE_SUCCESS;
}

int PCA::CalcCovarianceMartix()
{
	if (m_src_dataset == NULL)
		return RE_FILENOTEXIST;

	int width = m_src_dataset->GetRasterXSize();
	int height = m_src_dataset->GetRasterYSize();
	int image_dims = width*height;

	int element_num = m_band_count*m_band_count;
	int element_index = 0;
	
	m_relativity = new double[element_num];
	Map<MyMatrix> covariance_matrix(m_relativity, m_band_count, m_band_count);
	
	// 求协方差或相关系数矩阵
	for (int i1 = 1; i1 <= m_band_count; i1++)
	{
		for (int i2 = 1; i2 <= m_band_count; i2++)
		{
			if (i2 < i1)
			{
				// TODO：有待验证
				m_relativity[element_index] = m_relativity[(i2 - 1)*m_band_count + (i1 - 1)];
				
				element_index++;
				continue;
			}

			if (i1 == i2)
			{
				if (!m_is_covariance)
					m_relativity[element_index] = 1.0;
				else
					m_relativity[element_index] = m_band_stad[i1 - 1] * m_band_stad[i1 -1];
				
				element_index++;
				continue;
			}

			//////////////////////////////////////////////////////////////////////////
			// 分块处理
			GDALRasterBand *band1 = m_src_dataset->GetRasterBand(i1);
			GDALRasterBand *band2 = m_src_dataset->GetRasterBand(i2);

			DT_64F *buffer_data1 = new DT_64F[width];
			DT_64F *buffer_data2 = new DT_64F[width];

			double temp = 0.0;
			for (int j = 0; j < height; j++)
			{
				band1->RasterIO(GF_Read, 0, j, width, 1, buffer_data1, width, 1, GDT_Float64, 0, 0);
				band2->RasterIO(GF_Read, 0, j, width, 1, buffer_data2, width, 1, GDT_Float64, 0, 0);

				for (int i = 0; i < width; i++)
					temp += (buffer_data1[i] - m_band_mean[i1 - 1]) * (buffer_data2[i] - m_band_mean[i2 - 1]);
			}

			RELEASE(buffer_data1);
			RELEASE(buffer_data2);
			//////////////////////////////////////////////////////////////////////////

			m_relativity[element_index] = temp / image_dims;

			if (!m_is_covariance)
				m_relativity[element_index] = m_relativity[element_index] / (m_band_stad[i1 - 1]*m_band_stad[i2 - 1]);

			element_index++;
		}
	}

	m_relate_matrix = covariance_matrix;

	m_sta_io->WriteCovarianceOrCorrelation(m_relativity);
	
	//fprintf(m_statistics, "Covariance\n");
	//for (int i = 0; i < m_band_count; i++)
	//{
	//	for (int j = 0; j < m_band_count; j++)
	//	{
	//		fprintf(m_statistics, "%-16.6f", m_relate_matrix(i, j));
	//	}
	//	fprintf(m_statistics, "\n");
	//}
	//fprintf(m_statistics, "\n\n");

	CalcEigenvaluesAndEigenvectors();
	return RE_SUCCESS;
}

void PCA::CalcEigenvaluesAndEigenvectors()
{
	MyMatrix matrix;
	matrix = m_relate_matrix;

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

	double *temp1 = new double[m_band_count*m_band_count];
	for (int i = 0; i < m_band_count; i++)
	{
		for (int j = 0; j < m_band_count; j++)
			temp1[i*m_band_count + j] = m_eigenvectors(i, j); 
	}
	m_sta_io->WriteEigenvectors(temp1);

	if (temp1)
	{
		delete []temp1;
		temp1 = NULL;
	}

	double *temp2 = new double[m_band_count];
	int i;
	for (i = 0; i < m_band_count; i++)
		temp2[i] = eigenvalues[i];
	m_sta_io->WriteEigenvalue(temp2);

	for (i = 0; i < m_band_count; i++)
		temp2[i] = accumulate_contribute[i];
	m_sta_io->WriteAccumulateContribute(temp2);
	m_sta_io->WriteToFile();

	if (temp2)
	{
		delete []temp2;
		temp2 = NULL;
	}

	//fprintf(m_statistics, "Eigenvectors\n");
	//for (int i = 0; i < m_band_count; i++)
	//{
	//	for (int j = 0; j < m_band_count; j++)
	//		fprintf(m_statistics, "%-16.6f", m_eigenvectors(j, i));

	//	fprintf(m_statistics, "\n");
	//}
	//fprintf(m_statistics, "\n\n");

	//fprintf(m_statistics, "Eigenvalues\n");
	//for (int i = 0; i < m_band_count; i++)
	//{
	//	fprintf(m_statistics, "%-16.6f\n", m_eigenvalues[i]);
	//}
	//fprintf(m_statistics, "\n\n");

	//fprintf(m_statistics, "contribute\n");
	//for (int i = 0; i < m_band_count; i++)
	//{
	//	fprintf(m_statistics, "%-16.6f\n", contribute[i]);
	//}
	//fprintf(m_statistics, "\n\n");

	//fprintf(m_statistics, "accumulate_contribute\n");
	//for (int i = 0; i < m_band_count; i++)
	//{
	//	fprintf(m_statistics, "%-16.6f\n", accumulate_contribute[i]);
	//}
	//fprintf(m_statistics, "\n\n");
}

int PCA::CreatePCAFile(const char* pca_file, int pca_band_count, const char* format)
{
	if (m_src_dataset == NULL)
		return RE_FILENOTEXIST;

	int width = m_src_dataset->GetRasterXSize();
	int height = m_src_dataset->GetRasterYSize();

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

	return LinearCombination(pca_file, m_select_eigenvectors, NULL, format);// Combination，组合
}

int PCA::LinearCombination(const char *pca_file, MyMatrix &select_eigenvectors,  double *mean, const char *format)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	int width = m_src_dataset->GetRasterXSize();
	int height = m_src_dataset->GetRasterYSize();
	int band_count = select_eigenvectors.cols();

	GDALDriver *pca_driver = (GDALDriver *)GDALGetDriverByName(format);
	GDALDataset *pca_dataset = pca_driver->Create(pca_file, width, height, band_count, GDT_Float32, NULL);
	if ( pca_dataset == NULL )
		return RE_FILENOTSUPPORT;

	double geo_transform[6] = { 0 };
	m_src_dataset->GetGeoTransform(geo_transform);
	pca_dataset->SetGeoTransform(geo_transform);
	pca_dataset->SetProjection(m_src_dataset->GetProjectionRef());

	DT_64F *src_buffer_data = new DT_64F[width*m_band_count];
	DT_32F *pca_buffer_data = new DT_32F[width*band_count];
	int *src_band_map = new int[m_band_count];
	int *pca_band_map = new int[band_count];

	for (int i = 1; i <= m_band_count; i++)
		src_band_map[i - 1] = i;
	for (int i = 1; i <= band_count; i++)
		pca_band_map[i - 1] = i;

	//////////////////////////////////////////////////////////////////////////
	// 分块处理
	for (int h = 0; h < height; h++)
	{
		m_src_dataset->RasterIO(GF_Read, 0, h, width, 1, src_buffer_data, width, 1, GDT_Float64, m_band_count, src_band_map, 0, 0, 0);

		for (int j = 0; j < width; j++)
		{
			for (int k1 = 0; k1 < band_count; k1++)
			{
				double temp = 0.0;
				for (int k2 = 0; k2 < m_band_count; k2++)
				{			
					temp += src_buffer_data[k2*width + j] * select_eigenvectors(k2, k1);			
				}

				if (mean == NULL)
					pca_buffer_data[k1*width + j] = (float)temp;
				else
					pca_buffer_data[k1*width + j] = (float)temp + (float)mean[k1];
			}
		}

		pca_dataset->RasterIO(GF_Write, 0, h, width, 1, pca_buffer_data, width, 1, GDT_Float32, band_count, pca_band_map, 0, 0, 0);
	}

	RELEASE(src_buffer_data);
	RELEASE(pca_buffer_data);
	RELEASE(src_band_map);
	RELEASE(pca_band_map);
	//////////////////////////////////////////////////////////////////////////

	GDALClose((GDALDatasetH) pca_dataset);

	return RE_SUCCESS;
}

int PCA::CalcSubAvg(const char* pca_file)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDataset *pca_dataset = (GDALDataset *) GDALOpen(pca_file, GA_Update);
	if (pca_dataset == NULL)
		return RE_FILENOTSUPPORT;

	int width = pca_dataset->GetRasterXSize();
	int height = pca_dataset->GetRasterYSize();
	int band_count = pca_dataset->GetRasterCount();

	DT_32F *buffer_data = new DT_32F[width];
	for (int b = 1; b <= band_count; b++)
	{
		double min_value = 0.0;
		double max_value = 0.0;
		double mean_value = 0.0;
		double stddev_value = 0.0;

		GDALRasterBand *band = pca_dataset->GetRasterBand(b);
		band->ComputeStatistics(FALSE, &min_value, &max_value, &mean_value, &stddev_value, NULL, NULL);
		for (int i = 0; i < height; i++)
		{
			band->RasterIO(GF_Read, 0, i, width, 1, buffer_data, width, 1, GDT_Float32, 0, 0);
			
			for (int j = 0; j < width; j++)
				buffer_data[j] = buffer_data[j] - mean_value;

			band->RasterIO(GF_Write, 0, i, width, 1, buffer_data, width, 1, GDT_Float32, 0, 0);
		}
	}

	GDALClose((GDALDatasetH) pca_dataset);
	RELEASE(buffer_data);

	string temp_file = string(pca_file) + ".aux.xml";
	remove(temp_file.c_str());

	return RE_SUCCESS;
}

int PCA::ExecuteInversePCA(const char* inverse_pca_file, const char *statistics_file, const char* format /* = "GTiff" */)
{
	m_sta_io = new PCAStatisticsIO(statistics_file);
	m_sta_io->ReadInit();
	int band_count = m_sta_io->ReadBandCount();

	double *egeinvectors = new double[band_count*band_count];
	m_band_mean = new double[band_count];
	m_sta_io->ReadEigenvectors(egeinvectors);
	m_sta_io->ReadMean(m_band_mean);

	Map<MyMatrix> egeinvectors_matrix(egeinvectors, band_count, band_count);
	MyMatrix inverse_egeinvectors(band_count, band_count);

	inverse_egeinvectors = egeinvectors_matrix.inverse();
	
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
	m_src_dataset = (GDALDataset *) GDALOpen(m_src_file, GA_ReadOnly);
	if (m_src_dataset == NULL)
		return RE_FILENOTEXIST;

	m_band_count = m_src_dataset->GetRasterCount();
	return LinearCombination(inverse_pca_file, inverse_egeinvectors, m_band_mean, format);
}



