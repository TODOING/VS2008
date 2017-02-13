#ifndef PCA_H
#define PCA_H

#include "matrix_eigenvalue.h"
#include "gdal.h"
class GDALDataset;
class PCAStatisticsIO;

class PCA
{
public:
	PCA(const char *src_file);
	~PCA();

	/**
	* @brief PC变换
	* @param pszPCAFile			输出主成分变换后文件的路径
	* @param iBandCount			主成分变换后文件的的波段个数（默认为全部-1）
	* @param bIsCovariance		采用相关系数还是方差-协方差矩阵来计算，默认为协方差矩阵
	* @param bIsLikeEnvi		计算结果是否按照ENVI方式输出，即将所有的数据减去均值，使得每个波段的均值为0
	* @param pszFormat			输出文件格式，默认为GeoTiff格式
	* @return 返回代码
	*/
	int ExecutePCA(const char* pca_file, int pca_band_count = -1, bool is_covariance = true, 
		bool is_like_envi = true, const char* format = "GTiff");

	int ExecutePCA(const char* pca_file, const char *statistics_file, int pca_band_count = -1, GDALDataType dst_type = GDT_Float64, bool is_covariance = true, 
		bool is_like_envi = true, const char* format = "GTiff");

	/**
	* @brief PC逆变换
	* @param pszPCAFile			输出主成分变换后文件的路径
	* @param pmMatrix			主成分变换的特征向量矩阵
	* @param pvMeanVector		原始图像的均值向量，可以为NULL
	* @param pszFormat			输出文件格式，默认为GeoTiff格式
	* @return 返回代码
	*/
	int ExecuteInversePCA(const char* inverse_pca_file, const char *statistics_file,  GDALDataType dst_type = GDT_Float64, const char* format = "GTiff");

private:
	/**
	* @brief 数据预处理，进行图像信息的统计等
	* @return 返回代码
	*/
	int PreProcessData();

	/**
	* @brief 计算协方差矩阵和相关系数矩阵R，第二步
	* @return 返回代码
	*/
	int CalcCovarianceMartix();

	/**
	* @brief 计算特征值和特征向量，第三步和计算贡献率以及累积贡献率，第四步
	*/
	void CalcEigenvaluesAndEigenvectors();

	/**
	* @brief 计算主成分得分，并写入到文件中，第五和第六步
	* @param pszPCAFile		输出主成分变换后文件的路径
	* @param iBandCount		主成分变换后文件的的波段个数（默认为全部-1）
	* @param pszFormat		输出文件格式
	* @return 返回代码
	*/
	int CreatePCAFile(const char* pca_file, int pca_band_count, const char* format);

	/**
	* @brief 计算主成分得分，并写入到文件中，第五和第六步
	* @param pszPCAFile 输出主成分变换后文件的路径
	* @return 返回代码
	*/
	int CalcSubAvg(const char* pca_file);

	template <typename T>
	int CalcSubAvg(const char* pca_file);

	// TODO：写成模板
	int LinearCombination(const char *dst_file, MyMatrix &select_eigenvectors, double *mean, const char *format);


	template <typename T>
	int LinearCombination(const char *dst_file, MyMatrix &select_eigenvectors,  double *mean, const char *format);

	template <typename T>
	void CoreCalcInt16OrInt32(double *src_buffer_data, T *dst_buffer_data, double *mean, int dst_band_count, int block_sample);
	template <typename T>
	void CoreCalcInt16OrInt32(double *src_buffer_data, T *dst_buffer_data, int dst_band_count, int block_sample);
	
	template <typename T>
	void CoreCalcFloat32OrFloat64(double *src_buffer_data, T *dst_buffer_data, double *mean, int dst_band_count, int block_sample);
	template <typename T>
	void CoreCalcFloat32OrFloat64(double *src_buffer_data, T *dst_buffer_data, int dst_band_count, int block_sample);

private:	
	const char *m_src_file;		/*<! 要变换的文件路径 */
	const char *m_statistics_file; // 统计文件
	FILE *m_statistics;
	bool m_is_covariance;		/*<! PCA变换方式，true为协方差，false为相关系数 */

	GDALDataset *m_src_dataset;	/*<! 要变换的文件指针 */

	int m_band_count;			/*<! 波段个数 */	
	double *m_band_mean;			/*<! 波段均值 */	
	double *m_band_stad;			/*<! 波段标准差 */

	double *m_relativity;			/*<! 相关系数矩阵(协方差矩阵)中的元素 */
	MyMatrix m_relate_matrix;		/*<! 相关系数矩阵(协方差矩阵) */
	MyVector m_eigenvalues;			/*<! 相关系数矩阵(协方差矩阵)的特征值 */
	MyMatrix m_eigenvectors;		/*<! 相关系数矩阵(协方差矩阵)的特征向量 */
	MyMatrix m_select_eigenvectors; /*<! 构建选择后的特征向量矩阵 */
	PCAStatisticsIO *m_sta_io;
	
	GDALDataType m_dst_type;	// 目标数据类型，参照gdal.h中的GDALDataType, // 只允许保存（GDT_Int16, GDT_Int32, GDT_Float32, GDT_Float64等4种类型）
	int m_cache_size;
	int m_image_title_size;
};

#endif// PCA_H