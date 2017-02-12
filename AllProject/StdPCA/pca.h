#ifndef PCA_H
#define PCA_H

#include "matrix_eigenvalue.h"
class GDALDataset;
class PCAStatisticsIO;

class PCA
{
public:
	PCA(const char *src_file);
	~PCA();

	/**
	* @brief PC�任
	* @param pszPCAFile			������ɷֱ任���ļ���·��
	* @param iBandCount			���ɷֱ任���ļ��ĵĲ��θ�����Ĭ��Ϊȫ��-1��
	* @param bIsCovariance		�������ϵ�����Ƿ���-Э������������㣬Ĭ��ΪЭ�������
	* @param bIsLikeEnvi		�������Ƿ���ENVI��ʽ������������е����ݼ�ȥ��ֵ��ʹ��ÿ�����εľ�ֵΪ0
	* @param pszFormat			����ļ���ʽ��Ĭ��ΪGeoTiff��ʽ
	* @return ���ش���
	*/
	int ExecutePCA(const char* pca_file, int pca_band_count = -1, bool is_covariance = true, 
		bool is_like_envi = true, const char* format = "GTiff");

	int ExecutePCA(const char* pca_file, const char *statistics_file, int pca_band_count = -1, bool is_covariance = true, 
		bool is_like_envi = true, const char* format = "GTiff");

	/**
	* @brief PC��任
	* @param pszPCAFile			������ɷֱ任���ļ���·��
	* @param pmMatrix			���ɷֱ任��������������
	* @param pvMeanVector		ԭʼͼ��ľ�ֵ����������ΪNULL
	* @param pszFormat			����ļ���ʽ��Ĭ��ΪGeoTiff��ʽ
	* @return ���ش���
	*/
	int ExecuteInversePCA(const char* inverse_pca_file, const char *statistics_file,  const char* format = "GTiff");

private:
	/**
	* @brief ����Ԥ��������ͼ����Ϣ��ͳ�Ƶ�
	* @return ���ش���
	*/
	int PreProcessData();

	/**
	* @brief ����Э�����������ϵ������R���ڶ���
	* @return ���ش���
	*/
	int CalcCovarianceMartix();

	/**
	* @brief ��������ֵ�������������������ͼ��㹱�����Լ��ۻ������ʣ����Ĳ�
	*/
	void CalcEigenvaluesAndEigenvectors();

	/**
	* @brief �������ɷֵ÷֣���д�뵽�ļ��У�����͵�����
	* @param pszPCAFile		������ɷֱ任���ļ���·��
	* @param iBandCount		���ɷֱ任���ļ��ĵĲ��θ�����Ĭ��Ϊȫ��-1��
	* @param pszFormat		����ļ���ʽ
	* @return ���ش���
	*/
	int CreatePCAFile(const char* pca_file, int pca_band_count, const char* format);

	/**
	* @brief �������ɷֵ÷֣���д�뵽�ļ��У�����͵�����
	* @param pszPCAFile ������ɷֱ任���ļ���·��
	* @return ���ش���
	*/
	int CalcSubAvg(const char* pca_file);

	int LinearCombination(const char *pca_file, MyMatrix &select_eigenvectors, double *mean, const char *format);

private:	
	const char *m_src_file;		/*<! Ҫ�任���ļ�·�� */
	const char *m_statistics_file; // ͳ���ļ�
	FILE *m_statistics;
	bool m_is_covariance;		/*<! PCA�任��ʽ��trueΪЭ���falseΪ���ϵ�� */

	GDALDataset *m_src_dataset;	/*<! Ҫ�任���ļ�ָ�� */

	int m_band_count;			/*<! ���θ��� */	
	double *m_band_mean;			/*<! ���ξ�ֵ */	
	double *m_band_stad;			/*<! ���α�׼�� */

	double *m_relativity;			/*<! ���ϵ������(Э�������)�е�Ԫ�� */
	MyMatrix m_relate_matrix;		/*<! ���ϵ������(Э�������) */
	MyVector m_eigenvalues;			/*<! ���ϵ������(Э�������)������ֵ */
	MyMatrix m_eigenvectors;		/*<! ���ϵ������(Э�������)���������� */
	MyMatrix m_select_eigenvectors; /*<! ����ѡ���������������� */
	PCAStatisticsIO *m_sta_io;
};

#endif// PCA_H