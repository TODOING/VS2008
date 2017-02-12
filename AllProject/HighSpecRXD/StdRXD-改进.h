#ifndef STDRXD_H
#define STDRXD_H

#include "Eigen/Dense"

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

	void GetRXDMatrix(RMatrix &mEigenVectors, CVector &vMeanValues);
private:
	/**
	* @brief ����Ԥ��������ͼ����Ϣ��ͳ�Ƶ�(���������������ξ�ֵ����׼��)
	* @return ���ش���
	*/
	int PreProcessData();

	/**
	* @brief ����Э����������R���ڶ���
	* @return ���ش���
	*/
	int CalcCovarianceMartix();

	/**
	* @brief �������ɷֵ÷֣���д�뵽�ļ��У�����͵�����
	* @param pszPCAFile		������ɷֱ任���ļ���·��
	* @param iBandCount		���ɷֱ任���ļ��ĵĲ��θ�����Ĭ��Ϊȫ��-1��
	* @param pszFormat		����ļ���ʽ
	* @return ���ش���
	*/
	int CreateRXDFile(const char* pszRXDFile, int iBandCount, const char* pszFormat = "GTiff");

private: 

	const char* m_pszSrcFile;

	GDALDataset *m_pSrcDS;

	int m_iBandCount;			 	/*<! ���θ��� */	
	double *m_pBandMean;			/*<! ���ξ�ֵ */	
	double *m_pBandStad;			/*<! ���α�׼�� */
	double *m_pCovarianceElem;		/*<~! Э��������е�Ԫ�� */
	RMatrix m_InvCovarianceMatrix;		/*<~! Э������� */
	CVector m_PixelVector;          /*<~! ���������� */
	CVector m_MeanVector;			/*<~! ������ֵ������ */

	double *m_pSrcDataBuf;
	int m_imgWidth;
	int m_imgHeight;

};


#endif // STDRXD_H