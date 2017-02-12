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
	* @brief ����Ԥ������ȡӰ���������ڴ棬����ͼ����Ϣ��ͳ�Ƶ�(���������������ξ�ֵ����׼��)
	* @return ���ش���
	*/
	int PreProcessData();

	/**
	* @brief ����Э�������������
	* @return ���ش���
	*/
	int CalcInvCovarianceMartix();

	/**
	* @brief �������ɷֵ÷֣���д�뵽�ļ��У�����͵�����
	* @param pszRXDFile		������ɷֱ任���ļ���·��
	* @param iBandCount		���ɷֱ任���ļ��ĵĲ��θ�����Ĭ��Ϊȫ��-1��
	* @param pszFormat		����ļ���ʽ
	* @return ���ش���
	*/
	int CreateRXDFile(const char* pszRXDFile, int iBandCount, const char* pszFormat = "GTiff");

private: 

	const char* m_pszSrcFile;
	GDALDataset *m_pSrcDS;
	double *m_pBandStad;			/*<! ���α�׼�� */
	double *m_pSrcDataBuf;			/*<! �洢����Ӱ������ݣ��ٶȻ�ܿ� */

	int m_iBandCount;			 	/*<! ���θ��� */	
	int m_imgWidth;					/*<! Ӱ��� */
	int m_imgHeight;				/*<! Ӱ��� */
	
	RMatrix m_InvCovarianceMatrix;	/*<~! Э������� */
	CVector m_MeanVector;			/*<~! ������ֵ������ */

	std::vector<RVector> m_rowVector;
	std::vector<CVector> m_colVector;


};


#endif // STDRXD_H