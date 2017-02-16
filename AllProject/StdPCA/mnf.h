#ifndef MNF_H
#define MNF_H

#include "matrix_eigenvalue.h"
#include "gdal.h"
class GDALDataset;
class MNFStatisticsIO;

class MNF
{
public:
	MNF(const char *src_file);
	~MNF();

	/**
	* @brief PC变换
	* @param pszMNFFile			输出主成分变换后文件的路径
	* @param iBandCount			主成分变换后文件的的波段个数（默认为全部-1）
	* @param bIsCovariance		采用相关系数还是方差-协方差矩阵来计算，默认为协方差矩阵
	* @param bIsLikeEnvi		计算结果是否按照ENVI方式输出，即将所有的数据减去均值，使得每个波段的均值为0
	* @param pszFormat			输出文件格式，默认为GeoTiff格式
	* @return 返回代码
	*/
	int ExecuteMNF(const char* mnf_file, const char *statistics_file, int mnf_band_count = -1, GDALDataType dst_type = GDT_Float32, bool is_covariance = true, 
		bool is_like_envi = true, const char* format = "GTiff");

	/**
	* @brief PC逆变换
	* @param pszMNFFile			输出主成分变换后文件的路径
	* @param pmMatrix			主成分变换的特征向量矩阵
	* @param pvMeanVector		原始图像的均值向量，可以为NULL
	* @param pszFormat			输出文件格式，默认为GeoTiff格式
	* @return 返回代码
	*/
	int ExecuteInverseMNF(const char* inverse_mnf_file, const char *statistics_file,  GDALDataType dst_type = GDT_Float32, const char* format = "GTiff");

private:

	/**
	* @brief 计算特征值和特征向量，第三步和计算贡献率以及累积贡献率，第四步
	*/
	void CalcEigenvaluesAndEigenvectors();

	template <typename T1, typename T2>
	int RunMNF(const char *mnf_file, int mnf_band_count, bool is_float, bool is_like_envi, const char* format);

	template <typename T1, typename T2>
	int RunMNF_New(const char *mnf_file, int mnf_band_count, bool is_float, bool is_like_envi, const char* format);

	template <typename T1, typename T2>
	int RunInverseMNF(const char *inverse_mnf_file, MyMatrix &inverse_egeinvectors, bool is_float, const char* format);

private:	
	const char *m_src_file;		/*<! 要变换的文件路径 */
	GDALDataset *m_src_dataset;	/*<! 要变换的文件指针 */
	int m_band_count;			/*<! 波段个数 */	
	double *m_band_mean;			/*<! 波段均值 */	
	double *m_band_stad;			/*<! 波段标准差 */

	bool m_is_covariance;		/*<! MNF变换方式，true为协方差，false为相关系数 */
	
	const char *m_statistics_file; // 统计文件	
	FILE *m_statistics;
	MNFStatisticsIO *m_sta_io;

	double *m_covariance_or_relativity;		/*<! 相关系数矩阵(协方差矩阵)中的元素 */
	MyMatrix m_covar_or_relate_matrix;		/*<! 相关系数矩阵(协方差矩阵) */
	MyVector m_eigenvalues;					/*<! 相关系数矩阵(协方差矩阵)的特征值 */
	MyMatrix m_eigenvectors;				/*<! 相关系数矩阵(协方差矩阵)的特征向量 */
	MyMatrix m_select_eigenvectors;			/*<! 构建选择后的特征向量矩阵 */
	
	GDALDataType m_src_datatype;	
	GDALDataType m_dst_datatype;	// 目标数据类型，参照gdal.h中的GDALDataType, // 只允许保存（GDT_Int16, GDT_Int32, GDT_Float32, GDT_Float64等4种类型）

	int m_image_title_size;			// 分块大小,以MB为单位
};


#endif