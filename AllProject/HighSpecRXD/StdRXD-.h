#ifndef STDRXD_H
#define STDRXD_H

#include "Eigen/Dense"

#include <vector>

using namespace Eigen;

typedef Eigen::Matrix<double, Dynamic, Dynamic, RowMajor> RMatrix;
typedef Eigen::Matrix< double , 1, Dynamic > RVector;
typedef Eigen::Matrix< double , Dynamic, 1 > CVector;

class GDALDataset;

class CHighSpecRXD
{
public: 
	CHighSpecRXD(const char *pszSrcFile);
	~CHighSpecRXD();

	int ExecuteStdRXD(const char* pszRXDFile, int iBandCount = -1, const char* pszFormat = "GTiff");
private:
	/**
	* @brief 数据预处理，读取影像数据至内存，进行图像信息的统计等(包括波段数，波段均值，标准差)
	* @return 返回代码
	*/
	int PreProcessData();

	/**
	* @brief 计算协方差矩阵的逆矩阵
	* @return 返回代码
	*/
	int CalcInvCovarianceMartix();

	/**
	* @brief 计算主成分得分，并写入到文件中，第五和第六步
	* @param pszRXDFile		输出主成分变换后文件的路径
	* @param iBandCount		主成分变换后文件的的波段个数（默认为全部-1）
	* @param pszFormat		输出文件格式
	* @return 返回代码
	*/
	int CreateRXDFile(const char* pszRXDFile, int iBandCount, const char* pszFormat = "GTiff");

private: 

	const char* m_pszSrcFile;
	GDALDataset *m_pSrcDS;
	double *m_pBandStad;			/*<! 波段标准差 */
	double *m_pSrcDataBuf;			/*<! 存储整幅影像的数据，速度会很快 */

	int m_iBandCount;			 	/*<! 波段个数 */	
	int m_imgWidth;					/*<! 影像宽 */
	int m_imgHeight;				/*<! 影像高 */
	
	RMatrix m_InvCovarianceMatrix;	/*<~! 协方差矩阵 */
	CVector m_MeanVector;			/*<~! 样本均值列向量 */

	std::vector<RVector> m_rowVector;
	std::vector<CVector> m_colVector;


};


#endif // STDRXD_H