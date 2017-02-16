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
	* @brief PC�任
	* @param pszMNFFile			������ɷֱ任���ļ���·��
	* @param iBandCount			���ɷֱ任���ļ��ĵĲ��θ�����Ĭ��Ϊȫ��-1��
	* @param bIsCovariance		�������ϵ�����Ƿ���-Э������������㣬Ĭ��ΪЭ�������
	* @param bIsLikeEnvi		�������Ƿ���ENVI��ʽ������������е����ݼ�ȥ��ֵ��ʹ��ÿ�����εľ�ֵΪ0
	* @param pszFormat			����ļ���ʽ��Ĭ��ΪGeoTiff��ʽ
	* @return ���ش���
	*/
	int ExecuteMNF(const char* mnf_file, const char *statistics_file, int mnf_band_count = -1, GDALDataType dst_type = GDT_Float32, bool is_covariance = true, 
		bool is_like_envi = true, const char* format = "GTiff");

	/**
	* @brief PC��任
	* @param pszMNFFile			������ɷֱ任���ļ���·��
	* @param pmMatrix			���ɷֱ任��������������
	* @param pvMeanVector		ԭʼͼ��ľ�ֵ����������ΪNULL
	* @param pszFormat			����ļ���ʽ��Ĭ��ΪGeoTiff��ʽ
	* @return ���ش���
	*/
	int ExecuteInverseMNF(const char* inverse_mnf_file, const char *statistics_file,  GDALDataType dst_type = GDT_Float32, const char* format = "GTiff");

private:

	/**
	* @brief ��������ֵ�������������������ͼ��㹱�����Լ��ۻ������ʣ����Ĳ�
	*/
	void CalcEigenvaluesAndEigenvectors();

	template <typename T1, typename T2>
	int RunMNF(const char *mnf_file, int mnf_band_count, bool is_float, bool is_like_envi, const char* format);

	template <typename T1, typename T2>
	int RunMNF_New(const char *mnf_file, int mnf_band_count, bool is_float, bool is_like_envi, const char* format);

	template <typename T1, typename T2>
	int RunInverseMNF(const char *inverse_mnf_file, MyMatrix &inverse_egeinvectors, bool is_float, const char* format);

private:	
	const char *m_src_file;		/*<! Ҫ�任���ļ�·�� */
	GDALDataset *m_src_dataset;	/*<! Ҫ�任���ļ�ָ�� */
	int m_band_count;			/*<! ���θ��� */	
	double *m_band_mean;			/*<! ���ξ�ֵ */	
	double *m_band_stad;			/*<! ���α�׼�� */

	bool m_is_covariance;		/*<! MNF�任��ʽ��trueΪЭ���falseΪ���ϵ�� */
	
	const char *m_statistics_file; // ͳ���ļ�	
	FILE *m_statistics;
	MNFStatisticsIO *m_sta_io;

	double *m_covariance_or_relativity;		/*<! ���ϵ������(Э�������)�е�Ԫ�� */
	MyMatrix m_covar_or_relate_matrix;		/*<! ���ϵ������(Э�������) */
	MyVector m_eigenvalues;					/*<! ���ϵ������(Э�������)������ֵ */
	MyMatrix m_eigenvectors;				/*<! ���ϵ������(Э�������)���������� */
	MyMatrix m_select_eigenvectors;			/*<! ����ѡ���������������� */
	
	GDALDataType m_src_datatype;	
	GDALDataType m_dst_datatype;	// Ŀ���������ͣ�����gdal.h�е�GDALDataType, // ֻ�����棨GDT_Int16, GDT_Int32, GDT_Float32, GDT_Float64��4�����ͣ�

	int m_image_title_size;			// �ֿ��С,��MBΪ��λ
};


#endif