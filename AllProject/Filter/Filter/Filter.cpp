#include "stdafx.h"
#include "Filter.h"

#include "gdal_priv.h"
#include <vector>
#include <algorithm> // 使用std::sort()
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
		printf("打开文件%s数据集失败!\n", strSrcFile.c_str());
		return false;
	}

	GDALDataType eDataType = pSrcDs->GetRasterBand(1)->GetRasterDataType();
	int sizeX = pSrcDs->GetRasterXSize();
	int sizeY = pSrcDs->GetRasterYSize();
	double adfGeoTransform[6] = {0};
	pSrcDs->GetGeoTransform(adfGeoTransform);

	// 需要处理的波段
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
		
		// 创建输出文件
		pDstDs = pDriver->Create(strDstFile.c_str(), sizeX, sizeY, nBandCount, GDT_Int32, NULL);
		if (pDstDs == NULL)
		{
			printf("创建文件%s数据集失败!\n", strDstFile.c_str());
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
		
		// 创建输出文件
		pDstDs = pDriver->Create(strDstFile.c_str(), sizeX, sizeY, nBandCount, GDT_Float64, NULL);
		if (pDstDs == NULL)
		{
			printf("创建文件%s数据集失败!\n", strDstFile.c_str());
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

	int rr = (n-1)/2;// 模板行半径
	int rc = (m-1)/2;// 模板列半径
	int newSizeX = srcSizeX + 2*rc;
	int newSizeY = srcSizeY + 2*rr;
	int newDims = newSizeX*newSizeY;

	for (int i=0; i<nBandCount; i++) // 循环处理波段
	{
		GDALRasterBand *pSrcBandi = pSrcDs->GetRasterBand(panBandMap[i]);

		T *pSrcDB = new T[srcDims];
		pSrcBandi->RasterIO(GF_Read, 0, 0, srcSizeX, srcSizeY, pSrcDB, srcSizeX, srcSizeY, eDataType, 0, 0);

		// 重建数组
		T *pNewDB = new T[newDims];
		memset(pNewDB, 0, sizeof(T)*newDims);

		// 填充数组扩展边界--用原始影像最外一圈的值赋给扩展的边界（重复模式）
		BorderExtendT<T>(pSrcDB, pNewDB, srcSizeX, srcSizeY, newSizeX, rr, rc);

		double tempNum = 0;
		T *pSrcTemp = NULL;
		T *pNewTemp = NULL;
		for (int y=rr; y<newSizeY-rr; y++)
		{
			for (int x=rc; x<newSizeX-rc; x++)
			{
				pNewTemp = pNewDB + y*newSizeX + x;// 窗口中心位置
				pSrcTemp = pSrcDB + (y-rr)*srcSizeX + (x-rc);
				tempNum = GetCenterValueT<T>(pNewTemp, pCoefArray, rr, rc, newSizeX);
				tempNum *= coef;
				*pSrcTemp = static_cast<T>(tempNum);
			}
		}

		// 写入文件
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
	for (int y=rr; y>=-rr; y--)// 行
	{
		for (int x=rc; x>=-rc; x--)// 列
			resValue += (*(pCenterDataBuf-y*newDataBufX-x))*(*(pCoefArray++));// 窗口相乘
	}

	return resValue;
}

template<typename T>
void CFilter::BorderExtendT(T *pSrcDB, T *pNewDB, int srcSizeX, int srcSizeY, int newSizeX, int rr, int rc)
{
	T *pSrcTemp = pSrcDB;
	T *pNewTemp = pNewDB + rr*newSizeX + rc;// 与原始影像左上对齐
	int srcLineBytes = sizeof(T)*srcSizeX;

	// 将原始影像数据复制到新数组中相应位置
	for(int y=0; y<srcSizeY; y++)
		memcpy(pNewTemp+y*newSizeX, pSrcTemp+y*srcSizeX, srcLineBytes);

	// 扩展新数组上边界
	for (int i=1; i<=rr; i++)
	{
		memcpy(pNewTemp-i*newSizeX, pSrcTemp, srcLineBytes);

		// 上边界的首尾赋值
		for (int j=1; j<=rc; j++)
		{
			// 首赋值
			*(pNewTemp-i*newSizeX-j) = *pSrcTemp;
			// 尾赋值
			*(pNewTemp-i*newSizeX+srcSizeX-1+j) = *(pSrcTemp + srcSizeX-1);
		}
	}


	// 扩展新数组下边界
	for (int i=1; i<=rr; i++)
	{
		memcpy(pNewTemp+(srcSizeY-1+i)*newSizeX, pSrcTemp+(srcSizeY-1)*srcSizeX, srcLineBytes);

		// 下边界的首尾赋值
		for (int j=1; j<=rc; j++)
		{
			// 首赋值
			*(pNewTemp+(srcSizeY-1+i)*newSizeX-j) = *(pSrcTemp+(srcSizeY-1)*srcSizeX);
			// 尾赋值
			*(pNewTemp+(srcSizeY-1+i)*newSizeX+srcSizeX-1+j) = *(pSrcTemp+(srcSizeY-1)*srcSizeX + srcSizeX-1);
		}
	}

	// 扩展新数组左边界和右边界
	for (int y=0; y<srcSizeY; y++)
	{
		for (int i=1; i<=rc; i++)
		{
			// 左边界
			*(pNewTemp + y*newSizeX-i) = *(pSrcTemp + y*srcSizeX);
			// 右边界
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
	for (int y=rr; y>=-rr; y--)// 行
		for (int x=rc; x>=-rc; x--)// 列		
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
		printf("打开文件%s数据集失败!\n", strSrcFile.c_str());
		return false;
	}

	GDALDataType eDataType = pSrcDs->GetRasterBand(1)->GetRasterDataType();
	int sizeX = pSrcDs->GetRasterXSize();
	int sizeY = pSrcDs->GetRasterYSize();

	// 需要处理的波段
	int nBandCount = pSrcDs->GetRasterCount();
	int *panBandMap = new int[nBandCount];
	for (int i=0; i<nBandCount; i++)
		panBandMap[i] = i+1;

	// 创建输出文件
	GDALDriver *pDriver = (GDALDriver *) GDALGetDriverByName("GTiff");
	GDALDataset *pDstDs = pDriver->Create(strDstFile.c_str(), sizeX, sizeY, nBandCount, eDataType, NULL);
	if (pDstDs == NULL)
	{
		printf("创建文件%s数据集失败!\n", strDstFile.c_str());
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

	int rr = (n-1)/2;// 模板行半径
	int rc = (m-1)/2;// 模板列半径
	int newSizeX = srcSizeX ;//+ 2*rc;
	int newSizeY = srcSizeY ;//+ 2*rr;
	int newDims = newSizeX*newSizeY;

	T *pSrcDB = new T[srcDims];
	T *pNewDB = new T[newDims];
	for (int i=0; i<nBandCount; i++) // 循环处理波段
	{
		GDALRasterBand *pSrcBandi = pSrcDs->GetRasterBand(panBandMap[i]);
		pSrcBandi->RasterIO(GF_Read, 0, 0, srcSizeX, srcSizeY, pSrcDB, srcSizeX, srcSizeY, eDataType, 0, 0);

		T *pSrcTemp = pSrcDB;
		T *pNewTemp = pNewDB;// 与原始影像左上对齐
		memcpy(pNewTemp, pSrcTemp, sizeof(T)*srcDims);

		pSrcTemp = NULL;
		pNewTemp = NULL;
		for (int y=rr; y<newSizeY-rr; y++)
		{
			for (int x=rc; x<newSizeX-rc; x++)
			{
				pNewTemp = pNewDB + y*newSizeX + x;// 窗口中心位置
				pSrcTemp = pSrcDB + y*srcSizeX + x;
				*pSrcTemp = GetMedianCenterValueT(pNewTemp, rr, rc, newSizeX);
			}
		}

		// 写入文件
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