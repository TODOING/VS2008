#include "ssdc.h"
#include <stdio.h>
#include <iostream>
#include "gdal_priv.h"
#include "imgalg_define.h"
#include "statistics_io.h"
#include "AlgProcessTime.h"

using std::vector;

SSDC::SSDC(const char *src_file)
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

	m_block_w = 16;
	m_block_h = 16;
}

SSDC::~SSDC()
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

int SSDC::ExecuteSSDC(const char* ssdc_file, const char *statistics_file, int ssdc_band_count/* = -1*/, GDALDataType dst_type/* = 0*/, 
					bool is_covariance/* = true*/, bool is_like_envi/* = true*/, const char* format/* = "GTiff"*/)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	m_src_dataset = (GDALDataset *) GDALOpen(m_src_file, GA_ReadOnly);
	if (m_src_dataset == NULL)
		return RE_FILENOTEXIST;

	m_band_count = m_src_dataset->GetRasterCount();
	int width = m_src_dataset->GetRasterXSize();
	int height = m_src_dataset->GetRasterYSize();
	int band_size = width*height;

	// 分块信息
	//int block_unit = w*m_band_count;									// 分块基本单位（单位：像素）
	//int block_h = m_image_title_size * Mb / (block_unit*sizeof(T1));	// 每一块的高度（单位：像素）
	//int block_nums = h / block_h;										// 能分成的整块数n
	//int last_block_h = h % block_h;										// 剩余最后一块的高度（单位：像素）
	//if (last_block_h != 0) block_nums++;								// 如果未能整分，则最后剩余的行数也分块一块，即总块数为n+1。

	//int *band_map = new int[m_band_count];
	//for (int i = 0; i < m_band_count; i++)
	//	band_map[i] = i + 1;


	// 分块处理，将影像分成很多16*16大小的块，通过循环对每一块进行处理
	int x_block_num = width / m_block_w;  // 计算列方向上块数
	int y_block_num = height / m_block_h; // 计算行方向块数
	int block_size = m_block_w*m_block_h;
	int band_block_num = x_block_num*x_block_num;

	double *block_buf_data= new double[block_size*m_band_count];

	int coef_rows = m_block_w*m_block_h - 1;// i = j 不参与计算
	int coef_cols = 4;// a,b,c,d

	vector<double> band_noise_vec;// 存放所有波段的噪声估计
	int start_block_index = band_block_num * 0.1;
	int end_block_index = band_block_num * 0.9;
	int calc_block_num = end_block_index - start_block_index;
	for (int xb = 0; xb < x_block_num; xb++)
	{
		for (int yb = 0; yb < y_block_num; yb++)
		{
			// 将读取的数据按BSQ方式存在内存中
			m_src_dataset->RasterIO(GF_Read, xb*m_block_w, yb*m_block_h, m_block_w,
				m_block_h, block_buf_data, m_block_w, m_block_h, GDT_Float64, m_band_count, 0, 0, 0, 0);

			// 系数矩阵：coef = (X'X)^-1X'Y
			

			// 循环处理每一波段中的每一地块
			
			for (int b = 0; b < m_band_count; b++)
			{
				// b == 0
				vector<double> block_noise_vec;
				if (b == 0)
				{
					coef_cols = 3;
					double *X= new double[coef_rows*coef_cols];
					double *Y = new double[coef_rows];

					MyExtMatrix X_mat(X, coef_rows, coef_cols);
					MyExtVector Y_vec(Y, coef_rows);
					MyVector coef_vec(coef_cols);

					MyVector YY_vec(coef_rows);// 保存估计结果
					
					// 循环处理每一块，初始化求解回归系数的矩阵数据
					for (int i = 0; i < m_block_h; i++)
					{
						for (int j = 0; j < m_block_w; j++)
						{
							if (i > 0)
							{
								Y[i*m_block_w + j - 1] = block_buf_data[b*block_size + i*m_block_w + j];
								X[(i*m_block_w + j - 1)*3] = block_buf_data[(b+1)*block_size + i*m_block_w + j];// a
								X[(i*m_block_w + j - 1)*3 + 1] = block_buf_data[b*block_size + i*m_block_w + j - 1];// c
								X[(i*m_block_w + j - 1)*3 + 2] = 1;// d
							}

							if (i == 0 && j > 0)
							{
								Y[i*m_block_w + j - 1] = block_buf_data[b*block_size + i*m_block_w + j];
								X[(i*m_block_w + j - 1)*3] = block_buf_data[(b+1)*block_size + i*m_block_w + j];// a
								X[(i*m_block_w + j - 1)*3 + 1] = block_buf_data[b*block_size + i*m_block_w + j - 1];// c
								X[(i*m_block_w + j - 1)*3 + 2] = 1;// d
							}

						}// end for w of block
					}// end for h of block
					
					// 求解系数向量：coef = (X'X)^-1X'Y
					MyMatrix temp_mat = X_mat.transpose();
					coef_vec = (temp_mat*X_mat).inverse()*temp_mat*Y_vec;// 回归系数向量

					YY_vec = X_mat*coef_vec;// 估计值
					YY_vec = Y_vec - YY_vec;// 残差
					double SS = YY_vec.dot(YY_vec);// 点乘求r^2和
					double unbiased_esti_vari = SS / (coef_rows - 3);
					block_noise_vec.push_back(unbiased_esti_vari);

					printf("a = %.6lf\n", coef_vec[0]);
					printf("c = %.6lf\n", coef_vec[1]);
					printf("d = %.6lf\n", coef_vec[2]);
					printf("SS = %.6lf\n", SS);
					printf("unbiased_esti_vari = %.6lf\n", unbiased_esti_vari);
					
					continue;
				}


				// TODO
				// 0 < b < (m_band_count - 1)
				if (b > 0 && b < m_band_count - 1)
				{

					continue;
				}

				// TODO
				// b == m_band_count - 1
				if (b == m_band_count - 1)
				{

					continue;
				}

				// 先进行排序，然后在按照
				block_noise_vec.sort();// 从小到大排序				
				double sum = 0.0;
				for (int i = start_block_index; i < end_block_index; i++)
					sum += block_noise_vec[i];

				band_noise_vec.push_back(sum / calc_block_num);

			}// end for band

		}// end for yb of block

	}// end for xb of block




	


	return RE_SUCCESS;
}
