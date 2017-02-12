#ifndef CONVOLUTIONS_H
#define CONVOLUTIONS_H

#include <string>
#include "gdal_priv.h"
using std::string;
/**
* @brief 低通滤波
*/

enum BORDER_EXT
{
	BORDER_REPLOCATE,/*aaaaaa|abcdefgh|hhhhhhh，重复(ENVI这么处理)*/
	BORDER_CONSTANT, /*iiiiiiii|abcdefgh|iiiiiiiii，常量i填充，一般为0*/
	BORDER_REFLECT, /*fedcba|abcdefgh|hgfedcb，反射*/
};

bool LowPass(string strSrcFile, string strDstFile, int kernelSize, BORDER_EXT eBorderType);


bool LowPassMxN(string strSrcFile, string strDstFile, int m, int n, float *fCoefArray);
bool LowPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n, float *pCoefArray, float coef=1);

/**
* @brief 获取低通滤波后窗口中心的值
* @param pCenterDataBuf 待处理的中心值的地址(在重建后的数组里)
* @param pCoefArray 模板数组地址
* @param rRow 模板窗口行半径
* @param rCol 模板窗口列半径
* @param newDataBufX 重建后的数组的列大小
* @param return 返回处理后的中心值
*/
float GetCenterValue(int *pCenterDataBuf, float *pCoefArray, int rRow, int rCol, int newDataBufX);

/************************************************************************/
/*							以下为函数模板                         */
/************************************************************************/

template<typename T>
void LowPassMxNGenericT(GDALDataset *pSrcDs, GDALDataset *pDstDs, GDALDataType eDataType, int m, int n, float *pCoefArray, 
					   float coef=1, int nBandCount=0, int *panBandMap=NULL);

// 返回模板和窗口相乘的结果
template<typename T>
double GetCenterValueT(T *pCenterDataBuf, float *pCoefArray, int rRow, int rCol, int newDataBufX);
template<typename T>
void BorderExtendT(T *pSrcDB, T *pNewDB, int srcSizeX, int srcSizeY, int newSizeX, int rr, int rc);

#endif // CONVOLUTIONS_H
