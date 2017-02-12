#ifndef FILTER_H
#define FILTER_H

#include <string>
#include "gdal_priv.h"
using std::string;

// 平滑（Box模板(只考虑领域影响)，高斯(Gauss)低通滤波(考虑领域和位置影响)，中值滤波）

/**
* @brief 滤波工具类
*/
class CFilter
{
public:
	/**
	* @brief 低通滤波（平滑），可进行m*n模板滤波(m,n为奇数)
	* @param strSrcFile 影像文件
	* @param strDstFile 输出的滤波文件
	* @param m 模板宽
	* @param n 模板高
	* @param pCoefArray 卷积数组
	* @param coef 卷积系数
	* @param return 操作成功返回true，否则返回false
	*/
	static bool BoxLowPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n, float *pCoefArray, float coef=1);
	static bool GaussLowPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n, float *pCoefArray, float coef=1);
	static bool MedianLowPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n);

	static bool LaplacianHighPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n, float *pCoefArray, float coef=1);


private:
	/**
	* @brief 低通滤波函数模板（Box和Gauss），可进行不同数据类型的低通滤波
	* @param pSrcDs 输入影像数据集
	* @param pDstDs 输出滤波后的影像数据集
	* @param eDataType 影像的数据类型
	* @param m 模板宽
	* @param n 模板高
	* @param pCoefArray 卷积数组
	* @param coef 卷积系数
	* @param nBandCount 需要处理的波段数
	* @param panBandMap 需要处理的波段序号数组(从1开始)
	* @param return
	*/
	template<typename T>
	static void LowPassMxNGenericT(GDALDataset *pSrcDs, GDALDataset *pDstDs, GDALDataType eDataType, int m, int n, float *pCoefArray, 
		float coef=1, int nBandCount=0, int *panBandMap=NULL);

	/**
	* @brief 低通滤波的核心函数，进行模板与窗口运算
	* @param pCenterDataBuf 窗口中心位置(重建后的数组里)
	* @param pCoefArray 卷积数组
	* @param rr 模板行半径
	* @param rc 模板列半径
	* @param newDataBufX 重建后的数组的宽
	* @param return 返回模板与窗口运算的结果(相乘)
	*/
	template<typename T>
	static double GetCenterValueT(T *pCenterDataBuf, float *pCoefArray, int rr, int rc, int newDataBufX);
	
	/**
	* @brief 扩展重建后的数组的外边界，复制影像最外圈的值（扩展方式：重复）
	* @param pSrcDB 影像数据缓冲区
	* @param pNewDB 重建后的数组的缓冲区
	* @param srcSizeX 影像宽
	* @param srcSizeY 影像高
	* @param newSizeX 重建后的数组的宽
	* @param rr 模板行半径
	* @param rc 模板列半径
	* @param return
	*/
	template<typename T>
	static void BorderExtendT(T *pSrcDB, T *pNewDB, int srcSizeX, int srcSizeY, int newSizeX, int rr, int rc);


	/**
	* @brief 获取中值滤波中的中值
	* @param pCenterDataBuf 窗口中心位置(重建后的数组里)
	* @param rr 模板行半径
	* @param rc 模板列半径
	* @param newDataBufX 重建后的数组的宽
	* @param return 返回窗口中的中值
	*/
	template<typename T>
	static T GetMedianCenterValueT(T *pCenterDataBuf, int rr, int rc, int newDataBufX);

	/**
	* @brief 低通滤波函数模板(中值滤波)，可进行不同数据类型的低通滤波
	* @param pSrcDs 输入影像数据集
	* @param pDstDs 输出滤波后的影像数据集
	* @param eDataType 影像的数据类型
	* @param m 模板宽
	* @param n 模板高
	* @param nBandCount 需要处理的波段数
	* @param panBandMap 需要处理的波段序号数组(从1开始)
	* @param return
	*/
	template<typename T>
	static void MedianLowPassMxNGenericT(GDALDataset *pSrcDs, GDALDataset *pDstDs, GDALDataType eDataType, int m, int n, int nBandCount=0, int *panBandMap=NULL);

	template<typename T>
	static void LaplacianHighPassMxNGenericT(GDALDataset *pSrcDs, GDALDataset *pDstDs, GDALDataType eDataType, int m, int n, float *pCoefArray, 
		float coef=1, int nBandCount=0, int *panBandMap=NULL);

};

// 交换两个数
template<typename T>
inline void exch (T *array, int i, int j);

// 在数组中查找中位数
template<typename T>
T select_middle (T *array, int beg, int end, int n);

#endif // FILTER_H