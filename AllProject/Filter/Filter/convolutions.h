#ifndef CONVOLUTIONS_H
#define CONVOLUTIONS_H

#include <string>
#include "gdal_priv.h"
using std::string;
/**
* @brief ��ͨ�˲�
*/

enum BORDER_EXT
{
	BORDER_REPLOCATE,/*aaaaaa|abcdefgh|hhhhhhh���ظ�(ENVI��ô����)*/
	BORDER_CONSTANT, /*iiiiiiii|abcdefgh|iiiiiiiii������i��䣬һ��Ϊ0*/
	BORDER_REFLECT, /*fedcba|abcdefgh|hgfedcb������*/
};

bool LowPass(string strSrcFile, string strDstFile, int kernelSize, BORDER_EXT eBorderType);


bool LowPassMxN(string strSrcFile, string strDstFile, int m, int n, float *fCoefArray);
bool LowPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n, float *pCoefArray, float coef=1);

/**
* @brief ��ȡ��ͨ�˲��󴰿����ĵ�ֵ
* @param pCenterDataBuf �����������ֵ�ĵ�ַ(���ؽ����������)
* @param pCoefArray ģ�������ַ
* @param rRow ģ�崰���а뾶
* @param rCol ģ�崰���а뾶
* @param newDataBufX �ؽ����������д�С
* @param return ���ش���������ֵ
*/
float GetCenterValue(int *pCenterDataBuf, float *pCoefArray, int rRow, int rCol, int newDataBufX);

/************************************************************************/
/*							����Ϊ����ģ��                         */
/************************************************************************/

template<typename T>
void LowPassMxNGenericT(GDALDataset *pSrcDs, GDALDataset *pDstDs, GDALDataType eDataType, int m, int n, float *pCoefArray, 
					   float coef=1, int nBandCount=0, int *panBandMap=NULL);

// ����ģ��ʹ�����˵Ľ��
template<typename T>
double GetCenterValueT(T *pCenterDataBuf, float *pCoefArray, int rRow, int rCol, int newDataBufX);
template<typename T>
void BorderExtendT(T *pSrcDB, T *pNewDB, int srcSizeX, int srcSizeY, int newSizeX, int rr, int rc);

#endif // CONVOLUTIONS_H
