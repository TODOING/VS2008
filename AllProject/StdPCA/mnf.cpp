#include "mnf.h"
#include "gdal_priv.h"
#include "imgalg_define.h"
#include <stdio.h>
#include <iostream>
#include "statistics_io.h"

#include "AlgProcessTime.h"

MNF::MNF(const char *src_file)
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

MNF::~MNF()
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

int MNF::ExecuteMNF(const char* mnf_file, const char *statistics_file, int mnf_band_count/* = -1*/, GDALDataType dst_type/* = 0*/, 
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

	m_statistics = fopen(m_statistics_file, "wb");
	if (m_statistics == NULL)
		return RE_FILENOTSUPPORT;

	m_src_datatype = m_src_dataset->GetRasterBand(1)->GetRasterDataType();

	int w = m_src_dataset->GetRasterXSize();
	int h = m_src_dataset->GetRasterYSize();
	int band_size = w*h;

	// 分块信息
	int block_unit = w*m_band_count;									// 分块基本单位（单位：像素）
	int block_h = m_image_title_size * Mb / (block_unit*sizeof(double));	// 每一块的高度（单位：像素）
	int block_nums = h / block_h;										// 能分成的整块数n
	int last_block_h = h % block_h;										// 剩余最后一块的高度（单位：像素）
	if (last_block_h != 0) block_nums++;								// 如果未能整分，则最后剩余的行数也分块一块，即总块数为n+1。

	int *band_map = new int[m_band_count];
	for (int i = 0; i < m_band_count; i++)
		band_map[i] = i + 1;

	// 统计均值
	int read_h = block_h;
	double *block_buf_data;// = new double[read_h*block_unit];
	double *block_buf_extdata;// = new double[(read_h+1)*(w+1)*m_band_count];
	
	for (int i = 0; i < block_nums; i++)
	{
		if (i == block_nums - 1)
		{
			read_h = (h - 1)%block_h + 1;

			RELEASE(block_buf_data);
			RELEASE(block_buf_extdata);
			//block_buf_data = new double[read_h*block_unit];
		}

		// 读取数据
		if (i == 0)
		{
			block_buf_data = new double[read_h*block_unit];
			block_buf_extdata = new double[(read_h+1)*(w+1)*m_band_count];

			// 按BIP的格式读取图像,按波段存
			m_src_dataset->RasterIO(GF_Read, 0, 0, w, read_h, block_buf_data, w, read_h, GDT_Float64,
				m_band_count, band_map, 0, 0, 0);

			for (int i = 0; i<m_band_count; i++)
			{
				double *pSrcTemp = block_buf_data + w*read_h*i;
				double *pNewTemp = block_buf_extdata + (read_h+1)*(w + 1)*i;// 与原始影像左上对齐
				int srcLineBytes = sizeof(double)*w;

				// 将原始影像数据复制到新数组中相应位置
				memcpy(pNewTemp, pSrcTemp, srcLineBytes);
				for(int y=0; y<read_h; y++)
					memcpy(pNewTemp+w+1+y*(w+1), pSrcTemp+y*w, srcLineBytes);

				// 扩展新数组左边界和右边界
				for (int y=0; y<read_h+1; y++)
				{
					// 右边界
					*(pNewTemp + y*(w+1) + w) = *(pNewTemp + y*(w+1) + w-1);
				}

			}
			
		}
		else
		{
			block_buf_data = new double[(read_h+1)*w*m_band_count];
			block_buf_extdata = new double[(read_h+1)*(w+1)*m_band_count];

			// 按BIP的格式读取图像,按波段存
			m_src_dataset->RasterIO(GF_Read, 0, i*block_h - 1, w, read_h+1, block_buf_data, w, read_h+1, GDT_Float64,
				m_band_count, band_map, 0, 0, 0);


			for (int i = 0; i<m_band_count; i++)
			{
				double *pSrcTemp = block_buf_data + w*(read_h+1)*i;
				double *pNewTemp = block_buf_extdata + (read_h+1)*(w + 1)*i;// 与原始影像左上对齐
				int srcLineBytes = sizeof(double)*w;

				// 将原始影像数据复制到新数组中相应位置
				for(int y=0; y<read_h+1; y++)
					memcpy(pNewTemp+y*(w+1), pSrcTemp+y*w, srcLineBytes);

				// 扩展新数组左边界和右边界
				for (int y=0; y<read_h+1; y++)
				{
					*(pNewTemp + y*(w+1) + w) = *(pSrcTemp + y*w+w-1);
				}
			}
		}

		for (int i = 0; i < read_h; i++)
		{
			for (int j = 0; j < w; j++)
			{
				double *noise = new double[m_band_count];
 				for (int b1 = 0; b1 < m_band_count; b1++)
				{
					double *start = block_buf_extdata + w+1 + b1*(read_h+1)*(w+1);// 每个波段的起始位置

					noise[b1] = start[i*(w+1)+j] - (start[-(w+1)] + start[i*(w+1)+j+1])/2;//取一个光谱
					m_band_mean[b1] += noise[b1];	

		/*			int index1 = b1*m_band_count;
					for (int b2 = 0; b2 < m_band_count; b2++)
					{
						m_covariance_or_relativity[index1 + b2] += noise[b1]*noise[b2];
					}*/
				}		

				for(int ir=0;ir<m_band_count;ir++)
				{
					for(int il=0;il<m_band_count;il++)
					{
						m_covariance_or_relativity[ir*m_band_count+il] += noise[ir]*noise[il];
					}
				}
				RELEASE(noise);
			}
		}

		RELEASE(block_buf_data);
		RELEASE(block_buf_extdata);

	}
	RELEASE(block_buf_data);
	RELEASE(band_map);

	// 平均值
	for (int i = 0; i < m_band_count; i++)
		m_band_mean[i] /= band_size;

	// 协方差矩阵
	for (int i = 0; i < m_band_count; i++)
	{
		for (int j = 0; j < m_band_count; j++)
		{
			m_covariance_or_relativity[i*m_band_count + j] /= band_size;
			m_covariance_or_relativity[i*m_band_count + j] -= m_band_mean[i]*m_band_mean[j];
		}
	}
	

	return RE_SUCCESS;
}


// 实验
template <typename T1, typename T2>
int MNF::RunMNF_New(const char *mnf_file, int mnf_band_count, bool is_float, bool is_like_envi, const char* format)
{
	//int w = m_src_dataset->GetRasterXSize();
	//int h = m_src_dataset->GetRasterYSize();
	//int band_size = w*h;

	//// 分块信息
	//int block_unit = w*m_band_count;									// 分块基本单位（单位：像素）
	//int block_h = m_image_title_size * Mb / (block_unit*sizeof(T1));	// 每一块的高度（单位：像素）
	//int block_nums = h / block_h;										// 能分成的整块数n
	//int last_block_h = h % block_h;										// 剩余最后一块的高度（单位：像素）
	//if (last_block_h != 0) block_nums++;								// 如果未能整分，则最后剩余的行数也分块一块，即总块数为n+1。

	//int *band_map = new int[m_band_count];
	//for (int i = 0; i < m_band_count; i++)
	//	band_map[i] = i + 1;

	//CAlgProcessTime::Alg_start();

	///************************************************************************/
	///*						统计均值，标准差，协方差或相关系数矩阵                                                                     */
	///************************************************************************/

	//// 统计均值
	//int read_h = block_h;
	//T1 *block_buf_data = new T1[read_h*block_unit];

	//for (int i = 0; i < block_nums; i++)
	//{
	//	if (i == block_nums - 1)
	//	{
	//		read_h = (h - 1)%block_h + 1;
	//		RELEASE(block_buf_data);
	//		block_buf_data = new T1[read_h*block_unit];
	//	}

	//	// 按BIP的格式读取图像
	//	m_src_dataset->RasterIO(GF_Read, 0, i*block_h, w, read_h, block_buf_data, w, read_h, m_src_datatype,
	//		m_band_count, band_map, sizeof(T1)*m_band_count, sizeof(T1)*m_band_count*w, sizeof(T1));

	//	int band_block_size = read_h*w;
	//	for (int b = 0; b < m_band_count; b++)
	//	{			
	//		for (int j = 0; j < band_block_size; j++)		
	//		{
	//			double temp = block_buf_data[j*m_band_count + b];
	//			// sum of X：x1 + x2 + x3 + ... + xN
	//			m_band_mean[b] += temp;

	//			// sum of X*X：x1*x1 + x2*x2 + x3*x3 + ...
	//			m_band_stad[b] += temp*temp;


	//			// sum of X*Y：x1*y1 + x2*y2 + x3*y3 + ...
	//			int index1 = b*m_band_count;
	//			int index2 = j*m_band_count;
	//			for (int b1 = 0; b1 < m_band_count; b1++)
	//			{
	//				m_covariance_or_relativity[index1 + b1] += temp*block_buf_data[index2 + b1];
	//			}
	//		}			
	//	}
	//}
	//RELEASE(block_buf_data);
	//RELEASE(band_map);

	//// 平均值
	//for (int i = 0; i < m_band_count; i++)
	//	m_band_mean[i] /= band_size;

	//// 标准差
	//for (int i = 0; i < m_band_count; i++)
	//{
	//	m_band_stad[i] /= band_size;
	//	m_band_stad[i] = sqrt(m_band_stad[i] - m_band_mean[i]*m_band_mean[i]);
	//}

	//// 协方差矩阵
	//for (int i = 0; i < m_band_count; i++)
	//{
	//	for (int j = 0; j < m_band_count; j++)
	//	{
	//		m_covariance_or_relativity[i*m_band_count + j] /= band_size;
	//		m_covariance_or_relativity[i*m_band_count + j] -= m_band_mean[i]*m_band_mean[j];

	//		if (!m_is_covariance)// 相关系数
	//		{
	//			if (i == j) 
	//			{
	//				m_covariance_or_relativity[i*m_band_count + j] = 1.0;
	//				continue;
	//			}

	//			m_covariance_or_relativity[i*m_band_count + j] /= (m_band_stad[i]*m_band_stad[j]);
	//		}
	//	}
	//}
	//m_sta_io->WriteMean(m_band_mean);

	//CAlgProcessTime::Alg_end("统计均值和协方差矩阵");

	///************************************************************************/
	///*						计算特征值和特征向量                                                                     */
	///************************************************************************/
	//CAlgProcessTime::Alg_start();

	//MyExtMatrix covariance_matrix(m_covariance_or_relativity, m_band_count, m_band_count);
	//m_covar_or_relate_matrix = covariance_matrix;
	//m_sta_io->WriteCovarianceOrCorrelation(m_covariance_or_relativity);
	//CalcEigenvaluesAndEigenvectors();

	return RE_SUCCESS;
}

