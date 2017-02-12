#include "stdafx.h"
#include "Filter.h"

#include "gdal_priv.h"
#include <vector>
#include <algorithm> // ʹ��std::sort()
using std::vector;

template<typename T>
inline void exch (T *array, int i, int j)  
{  
	if (i!=j)  
	{
		T temp = array[i];
		array[i] = array[j];
		array[j] = temp;
	}  
}  

template<typename T>
T select_middle (T *array, int beg, int end, int n)  
{  
	if (n == 1)  
		return array[0];  

	int i = beg, j;  
	for (j=i+1; j<=end; ++j)  
		if (array[j]<=array[beg])  
		{  
			++i;  
			exch<T>(array, i, j);  
		}  
		exch<T>(array, beg, i);  

		if (i < n/2)  
			return select_middle (array, i+1, end, n);  
		else if (i > n/2)  
			return select_middle (array, beg, i-1, n);  
		else  
		{  
			if (n%2)  
				return array[i];  
			else  
			{  
				T m = array[0];  
				for (int j=1; j<i; ++j)  
					if (array[j] > m)  
						m=array[j];  
				return static_cast<T>((array[i]+m)/2);  
			}  
		}  
}



bool CFilter::BoxLowPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n, float *pCoefArray, float coef/* =1 */)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	GDALDataset *pSrcDs = (GDALDataset *)GDALOpen(strSrcFile.c_str(), GA_ReadOnly);
	if (pSrcDs == NULL)
	{
		printf("���ļ�%s���ݼ�ʧ��!\n", strSrcFile.c_str());
		return false;
	}

	GDALDataType eDataType = pSrcDs->GetRasterBand(1)->GetRasterDataType();
	int sizeX = pSrcDs->GetRasterXSize();
	int sizeY = pSrcDs->GetRasterYSize();
	double adfGeoTransform[6] = {0};
	pSrcDs->GetGeoTransform(adfGeoTransform);

	// ��Ҫ����Ĳ���
	int nBandCount = pSrcDs->GetRasterCount();
	int *panBandMap = new int[nBandCount];
	for (int i=0; i<nBandCount; i++)
		panBandMap[i] = i+1;

	GDALDriver *pDriver = (GDALDriver *) GDALGetDriverByName("GTiff");
	GDALDataset *pDstDs = NULL;
	switch (eDataType)
	{
	case GDT_Byte:
		//LowPassMxNGenericT<unsigned char>(pSrcDs, pDstDs, eDataType, m, n, pCoefArray, 
		//	coef, nBandCount, panBandMap);
		//break;
	case GDT_UInt16:
		//LowPassMxNGenericT<unsigned short int>(pSrcDs, pDstDs, eDataType, m, n, pCoefArray, 
		//	coef, nBandCount, panBandMap);
		//break;
	case GDT_Int16:
		//LowPassMxNGenericT<short int>(pSrcDs, pDstDs, eDataType, m, n, pCoefArray, 
		//	coef, nBandCount, panBandMap);
		//break;
	case GDT_UInt32:
		//LowPassMxNGenericT<unsigned int>(pSrcDs, pDstDs, eDataType, m, n, pCoefArray, 
		//	coef, nBandCount, panBandMap);
		//break;
	case GDT_Int32:
		
		// ��������ļ�
		pDstDs = pDriver->Create(strDstFile.c_str(), sizeX, sizeY, nBandCount, GDT_Int32, NULL);
		if (pDstDs == NULL)
		{
			printf("�����ļ�%s���ݼ�ʧ��!\n", strDstFile.c_str());
			return false;
		}		
		pDstDs->SetGeoTransform(adfGeoTransform);
		pDstDs->SetProjection(pSrcDs->GetProjectionRef());

		LowPassMxNGenericT<int>(pSrcDs, pDstDs, GDT_Int32, m, n, pCoefArray, 
			coef, nBandCount, panBandMap);
		break;
	case GDT_Float32:
		//LowPassMxNGenericT<float>(pSrcDs, pDstDs, eDataType, m, n, pCoefArray, 
		//	coef, nBandCount, panBandMap);
		break;
	case GDT_Float64:
		
		// ��������ļ�
		pDstDs = pDriver->Create(strDstFile.c_str(), sizeX, sizeY, nBandCount, GDT_Float64, NULL);
		if (pDstDs == NULL)
		{
			printf("�����ļ�%s���ݼ�ʧ��!\n", strDstFile.c_str());
			return false;
		}
		pDstDs->SetGeoTransform(adfGeoTransform);
		pDstDs->SetProjection(pSrcDs->GetProjectionRef());

		LowPassMxNGenericT<double>(pSrcDs, pDstDs, GDT_Float64, m, n, pCoefArray, 
			coef, nBandCount, panBandMap);
		break;
	}

	GDALClose((GDALDatasetH)pSrcDs);
	GDALClose((GDALDatasetH)pDstDs);

	return true;
}

template<typename T>
void CFilter::LowPassMxNGenericT(GDALDataset *pSrcDs, GDALDataset *pDstDs, GDALDataType eDataType, int m, int n, float *pCoefArray, 
						float coef/* = 1*/, int nBandCount/* = 0*/, int *panBandMap/* = NULL*/)
{
	int srcSizeX = pSrcDs->GetRasterXSize();
	int srcSizeY = pSrcDs->GetRasterYSize();
	int srcDims = srcSizeX*srcSizeY;

	int rr = (n-1)/2;// ģ���а뾶
	int rc = (m-1)/2;// ģ���а뾶
	int newSizeX = srcSizeX + 2*rc;
	int newSizeY = srcSizeY + 2*rr;
	int newDims = newSizeX*newSizeY;

	for (int i=0; i<nBandCount; i++) // ѭ��������
	{
		GDALRasterBand *pSrcBandi = pSrcDs->GetRasterBand(panBandMap[i]);

		T *pSrcDB = new T[srcDims];
		pSrcBandi->RasterIO(GF_Read, 0, 0, srcSizeX, srcSizeY, pSrcDB, srcSizeX, srcSizeY, eDataType, 0, 0);

		// �ؽ�����
		T *pNewDB = new T[newDims];
		memset(pNewDB, 0, sizeof(T)*newDims);

		// ���������չ�߽�--��ԭʼӰ������һȦ��ֵ������չ�ı߽磨�ظ�ģʽ��
		BorderExtendT<T>(pSrcDB, pNewDB, srcSizeX, srcSizeY, newSizeX, rr, rc);

		double tempNum = 0;
		T *pSrcTemp = NULL;
		T *pNewTemp = NULL;
		for (int y=rr; y<newSizeY-rr; y++)
		{
			for (int x=rc; x<newSizeX-rc; x++)
			{
				pNewTemp = pNewDB + y*newSizeX + x;// ��������λ��
				pSrcTemp = pSrcDB + (y-rr)*srcSizeX + (x-rc);
				tempNum = GetCenterValueT<T>(pNewTemp, pCoefArray, rr, rc, newSizeX);
				tempNum *= coef;
				*pSrcTemp = static_cast<T>(tempNum);
			}
		}

		// д���ļ�
		GDALRasterBand *pDstBandi = pDstDs->GetRasterBand(panBandMap[i]);
		pDstBandi->RasterIO(GF_Write, 0, 0, srcSizeX, srcSizeY, pSrcDB, srcSizeX, srcSizeY, eDataType, 0, 0);

		delete []pSrcDB;
		delete []pNewDB;
		pSrcDB = NULL;
		pNewDB = NULL;

	}// for band
}

template<typename T>
double CFilter::GetCenterValueT(T *pCenterDataBuf, float *pCoefArray, int rr, int rc, int newDataBufX)
{
	double resValue=0;
	for (int y=rr; y>=-rr; y--)// ��
	{
		for (int x=rc; x>=-rc; x--)// ��
			resValue += (*(pCenterDataBuf-y*newDataBufX-x))*(*(pCoefArray++));// �������
	}

	return resValue;
}

template<typename T>
void CFilter::BorderExtendT(T *pSrcDB, T *pNewDB, int srcSizeX, int srcSizeY, int newSizeX, int rr, int rc)
{
	T *pSrcTemp = pSrcDB;
	T *pNewTemp = pNewDB + rr*newSizeX + rc;// ��ԭʼӰ�����϶���
	int srcLineBytes = sizeof(T)*srcSizeX;

	// ��ԭʼӰ�����ݸ��Ƶ�����������Ӧλ��
	for(int y=0; y<srcSizeY; y++)
		memcpy(pNewTemp+y*newSizeX, pSrcTemp+y*srcSizeX, srcLineBytes);

	// ��չ�������ϱ߽�
	for (int i=1; i<=rr; i++)
	{
		memcpy(pNewTemp-i*newSizeX, pSrcTemp, srcLineBytes);

		// �ϱ߽����β��ֵ
		for (int j=1; j<=rc; j++)
		{
			// �׸�ֵ
			*(pNewTemp-i*newSizeX-j) = *pSrcTemp;
			// β��ֵ
			*(pNewTemp-i*newSizeX+srcSizeX-1+j) = *(pSrcTemp + srcSizeX-1);
		}
	}


	// ��չ�������±߽�
	for (int i=1; i<=rr; i++)
	{
		memcpy(pNewTemp+(srcSizeY-1+i)*newSizeX, pSrcTemp+(srcSizeY-1)*srcSizeX, srcLineBytes);

		// �±߽����β��ֵ
		for (int j=1; j<=rc; j++)
		{
			// �׸�ֵ
			*(pNewTemp+(srcSizeY-1+i)*newSizeX-j) = *(pSrcTemp+(srcSizeY-1)*srcSizeX);
			// β��ֵ
			*(pNewTemp+(srcSizeY-1+i)*newSizeX+srcSizeX-1+j) = *(pSrcTemp+(srcSizeY-1)*srcSizeX + srcSizeX-1);
		}
	}

	// ��չ��������߽���ұ߽�
	for (int y=0; y<srcSizeY; y++)
	{
		for (int i=1; i<=rc; i++)
		{
			// ��߽�
			*(pNewTemp + y*newSizeX-i) = *(pSrcTemp + y*srcSizeX);
			// �ұ߽�
			*(pNewTemp + y*newSizeX+srcSizeX-1+i) = *(pSrcTemp + y*srcSizeX+srcSizeX-1);
		}
	}
}

bool CFilter::GaussLowPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n, float *pCoefArray, float coef/* =1 */)
{
	return BoxLowPassMxNGeneric(strSrcFile, strDstFile, m, n, pCoefArray, coef);
}

template<typename T>
T CFilter::GetMedianCenterValueT(T *pCenterDataBuf, int rr, int rc, int newDataBufX)
{
	int winDims = (2*rr+1)*(2*rc+1);
	T *pTemp = new T[winDims];
	T *pCenterTemp = pTemp + rr*(2*rc+1) + rc;
	for (int y=rr; y>=-rr; y--)// ��
		for (int x=rc; x>=-rc; x--)// ��		
			*(pCenterTemp-y*(2*rc+1)-x) = *(pCenterDataBuf-y*newDataBufX-x);

	return select_middle<T>(pTemp, 0, winDims-1, winDims);
}

bool CFilter::MedianLowPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	GDALDataset *pSrcDs = (GDALDataset *)GDALOpen(strSrcFile.c_str(), GA_ReadOnly);
	if (pSrcDs == NULL)
	{
		printf("���ļ�%s���ݼ�ʧ��!\n", strSrcFile.c_str());
		return false;
	}

	GDALDataType eDataType = pSrcDs->GetRasterBand(1)->GetRasterDataType();
	int sizeX = pSrcDs->GetRasterXSize();
	int sizeY = pSrcDs->GetRasterYSize();

	// ��Ҫ����Ĳ���
	int nBandCount = pSrcDs->GetRasterCount();
	int *panBandMap = new int[nBandCount];
	for (int i=0; i<nBandCount; i++)
		panBandMap[i] = i+1;

	// ��������ļ�
	GDALDriver *pDriver = (GDALDriver *) GDALGetDriverByName("GTiff");
	GDALDataset *pDstDs = pDriver->Create(strDstFile.c_str(), sizeX, sizeY, nBandCount, eDataType, NULL);
	if (pDstDs == NULL)
	{
		printf("�����ļ�%s���ݼ�ʧ��!\n", strDstFile.c_str());
		return false;
	}
	double adfGeoTransform[6] = {0};
	pSrcDs->GetGeoTransform(adfGeoTransform);
	pDstDs->SetGeoTransform(adfGeoTransform);
	pDstDs->SetProjection(pSrcDs->GetProjectionRef());

	switch (eDataType)
	{
	case GDT_Byte:
		MedianLowPassMxNGenericT<unsigned char>(pSrcDs, pDstDs, eDataType, m, n, nBandCount, panBandMap);
		break;
	case GDT_UInt16:
		MedianLowPassMxNGenericT<unsigned short int>(pSrcDs, pDstDs, eDataType, m, n, nBandCount, panBandMap);
		break;
	case GDT_Int16:
		MedianLowPassMxNGenericT<short int>(pSrcDs, pDstDs, eDataType, m, n, nBandCount, panBandMap);
		break;
	case GDT_UInt32:
		MedianLowPassMxNGenericT<unsigned int>(pSrcDs, pDstDs, eDataType, m, n, nBandCount, panBandMap);
		break;
	case GDT_Int32:
		MedianLowPassMxNGenericT<int>(pSrcDs, pDstDs, eDataType, m, n, nBandCount, panBandMap);
		break;
	case GDT_Float32:
		MedianLowPassMxNGenericT<float>(pSrcDs, pDstDs, eDataType, m, n, nBandCount, panBandMap);
		break;
	case GDT_Float64:
		MedianLowPassMxNGenericT<double>(pSrcDs, pDstDs, eDataType, m, n, nBandCount, panBandMap);
		break;
	}

	GDALClose((GDALDatasetH)pSrcDs);
	GDALClose((GDALDatasetH)pDstDs);

	return true;
}

template<typename T>
void CFilter::MedianLowPassMxNGenericT(GDALDataset *pSrcDs, GDALDataset *pDstDs, GDALDataType eDataType, int m, int n, int nBandCount/* =0 */, int *panBandMap/* =NULL */)
{
	int srcSizeX = pSrcDs->GetRasterXSize();
	int srcSizeY = pSrcDs->GetRasterYSize();
	int srcDims = srcSizeX*srcSizeY;

	int rr = (n-1)/2;// ģ���а뾶
	int rc = (m-1)/2;// ģ���а뾶
	int newSizeX = srcSizeX ;//+ 2*rc;
	int newSizeY = srcSizeY ;//+ 2*rr;
	int newDims = newSizeX*newSizeY;

	T *pSrcDB = new T[srcDims];
	T *pNewDB = new T[newDims];
	for (int i=0; i<nBandCount; i++) // ѭ��������
	{
		GDALRasterBand *pSrcBandi = pSrcDs->GetRasterBand(panBandMap[i]);
		pSrcBandi->RasterIO(GF_Read, 0, 0, srcSizeX, srcSizeY, pSrcDB, srcSizeX, srcSizeY, eDataType, 0, 0);

		T *pSrcTemp = pSrcDB;
		T *pNewTemp = pNewDB;// ��ԭʼӰ�����϶���
		memcpy(pNewTemp, pSrcTemp, sizeof(T)*srcDims);

		pSrcTemp = NULL;
		pNewTemp = NULL;
		for (int y=rr; y<newSizeY-rr; y++)
		{
			for (int x=rc; x<newSizeX-rc; x++)
			{
				pNewTemp = pNewDB + y*newSizeX + x;// ��������λ��
				pSrcTemp = pSrcDB + y*srcSizeX + x;
				*pSrcTemp = GetMedianCenterValueT(pNewTemp, rr, rc, newSizeX);
			}
		}

		// д���ļ�
		GDALRasterBand *pDstBandi = pDstDs->GetRasterBand(panBandMap[i]);
		pDstBandi->RasterIO(GF_Write, 0, 0, srcSizeX, srcSizeY, pSrcDB, srcSizeX, srcSizeY, eDataType, 0, 0);
	
	}// for band
	delete []pSrcDB;
	delete []pNewDB;
	pSrcDB = NULL;
	pNewDB = NULL;
}

bool CFilter::LaplacianHighPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n, float *pCoefArray, float coef/* =1 */)
{
	return BoxLowPassMxNGeneric(strSrcFile, strDstFile, m, n, pCoefArray, coef);
}