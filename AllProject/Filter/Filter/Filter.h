#ifndef FILTER_H
#define FILTER_H

#include <string>
#include "gdal_priv.h"
using std::string;

// ƽ����Boxģ��(ֻ��������Ӱ��)����˹(Gauss)��ͨ�˲�(���������λ��Ӱ��)����ֵ�˲���

/**
* @brief �˲�������
*/
class CFilter
{
public:
	/**
	* @brief ��ͨ�˲���ƽ�������ɽ���m*nģ���˲�(m,nΪ����)
	* @param strSrcFile Ӱ���ļ�
	* @param strDstFile ������˲��ļ�
	* @param m ģ���
	* @param n ģ���
	* @param pCoefArray �������
	* @param coef ���ϵ��
	* @param return �����ɹ�����true�����򷵻�false
	*/
	static bool BoxLowPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n, float *pCoefArray, float coef=1);
	static bool GaussLowPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n, float *pCoefArray, float coef=1);
	static bool MedianLowPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n);

	static bool LaplacianHighPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n, float *pCoefArray, float coef=1);


private:
	/**
	* @brief ��ͨ�˲�����ģ�壨Box��Gauss�����ɽ��в�ͬ�������͵ĵ�ͨ�˲�
	* @param pSrcDs ����Ӱ�����ݼ�
	* @param pDstDs ����˲����Ӱ�����ݼ�
	* @param eDataType Ӱ�����������
	* @param m ģ���
	* @param n ģ���
	* @param pCoefArray �������
	* @param coef ���ϵ��
	* @param nBandCount ��Ҫ����Ĳ�����
	* @param panBandMap ��Ҫ����Ĳ����������(��1��ʼ)
	* @param return
	*/
	template<typename T>
	static void LowPassMxNGenericT(GDALDataset *pSrcDs, GDALDataset *pDstDs, GDALDataType eDataType, int m, int n, float *pCoefArray, 
		float coef=1, int nBandCount=0, int *panBandMap=NULL);

	/**
	* @brief ��ͨ�˲��ĺ��ĺ���������ģ���봰������
	* @param pCenterDataBuf ��������λ��(�ؽ����������)
	* @param pCoefArray �������
	* @param rr ģ���а뾶
	* @param rc ģ���а뾶
	* @param newDataBufX �ؽ��������Ŀ�
	* @param return ����ģ���봰������Ľ��(���)
	*/
	template<typename T>
	static double GetCenterValueT(T *pCenterDataBuf, float *pCoefArray, int rr, int rc, int newDataBufX);
	
	/**
	* @brief ��չ�ؽ�����������߽磬����Ӱ������Ȧ��ֵ����չ��ʽ���ظ���
	* @param pSrcDB Ӱ�����ݻ�����
	* @param pNewDB �ؽ��������Ļ�����
	* @param srcSizeX Ӱ���
	* @param srcSizeY Ӱ���
	* @param newSizeX �ؽ��������Ŀ�
	* @param rr ģ���а뾶
	* @param rc ģ���а뾶
	* @param return
	*/
	template<typename T>
	static void BorderExtendT(T *pSrcDB, T *pNewDB, int srcSizeX, int srcSizeY, int newSizeX, int rr, int rc);


	/**
	* @brief ��ȡ��ֵ�˲��е���ֵ
	* @param pCenterDataBuf ��������λ��(�ؽ����������)
	* @param rr ģ���а뾶
	* @param rc ģ���а뾶
	* @param newDataBufX �ؽ��������Ŀ�
	* @param return ���ش����е���ֵ
	*/
	template<typename T>
	static T GetMedianCenterValueT(T *pCenterDataBuf, int rr, int rc, int newDataBufX);

	/**
	* @brief ��ͨ�˲�����ģ��(��ֵ�˲�)���ɽ��в�ͬ�������͵ĵ�ͨ�˲�
	* @param pSrcDs ����Ӱ�����ݼ�
	* @param pDstDs ����˲����Ӱ�����ݼ�
	* @param eDataType Ӱ�����������
	* @param m ģ���
	* @param n ģ���
	* @param nBandCount ��Ҫ����Ĳ�����
	* @param panBandMap ��Ҫ����Ĳ����������(��1��ʼ)
	* @param return
	*/
	template<typename T>
	static void MedianLowPassMxNGenericT(GDALDataset *pSrcDs, GDALDataset *pDstDs, GDALDataType eDataType, int m, int n, int nBandCount=0, int *panBandMap=NULL);

	template<typename T>
	static void LaplacianHighPassMxNGenericT(GDALDataset *pSrcDs, GDALDataset *pDstDs, GDALDataType eDataType, int m, int n, float *pCoefArray, 
		float coef=1, int nBandCount=0, int *panBandMap=NULL);

};

// ����������
template<typename T>
inline void exch (T *array, int i, int j);

// �������в�����λ��
template<typename T>
T select_middle (T *array, int beg, int end, int n);

#endif // FILTER_H