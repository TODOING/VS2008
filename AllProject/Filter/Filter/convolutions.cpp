#include "stdafx.h"

#include "convolutions.h"
#include "gdal_priv.h"
#include <vector>
using std::vector;

// ptr(row, col) = ptr + row * sizeCol + col;
bool LowPass(string strSrcFile, string strDstFile, int iKernelSize, BORDER_EXT eBorderType)
{
	GDALAllRegister();
	
	GDALDataset *pSrcDs = (GDALDataset *)GDALOpen(strSrcFile.c_str(), GA_ReadOnly);
	if (pSrcDs == NULL)
	{
		printf("打开文件%s失败!\n", strSrcFile.c_str());
		return false;
	}

	int iBbandCount = pSrcDs->GetRasterCount();
	GDALRasterBand *pSrcBand = pSrcDs->GetRasterBand(1);
	if (pSrcBand == NULL)
	{
		printf("打开波段1失败!\n");
		return false;
	}

	int sizeX = pSrcBand->GetXSize();
	int sizeY = pSrcBand->GetYSize();
	int imgDims = sizeX * sizeY;

	int *pSrcDB = new int[imgDims];
	pSrcBand->RasterIO(GF_Read, 0, 0, sizeX, sizeY, pSrcDB, sizeX, sizeY, GDT_Int32, 0, 0);

	int r = (iKernelSize - 1) / 2;// 相当于模板窗口半径
	int iNewX = sizeX + 2*r;
	int iNewY = sizeY + 2*r;
	int iNewDims = iNewX*iNewY;
	int *pNewDB = new int[iNewDims];
	memset(pNewDB, 0, sizeof(int)*iNewDims);

	int *pSrc = pSrcDB;
	int *pNew = pNewDB + r*iNewX + r;// 与原始影像左上对齐
	int iBufSize = sizeof(int)*sizeX;


	// 将原始影像数据复制到新数组中相应位置
	for(int y=0; y<sizeY; y++)
		memcpy(pNew+y*iNewX, pSrc+y*sizeX, iBufSize);

	// 扩展新数组上边界
	for (int i=1; i<=r; i++)
	{
		memcpy(pNew-i*iNewX, pSrc, iBufSize);
		
		// 上边界的首尾赋值
		for (int j=1; j<=r; j++)
		{
			*(pNew-i*iNewX-j) = *pSrc;
			*(pNew-i*iNewX+sizeX-1+j) = *(pSrc + sizeX-1);
		}
	}


	// 扩展新数组下边界
	for (int i=1; i<=r; i++)
	{
		memcpy(pNew+sizeY*iNewX, pSrc+(sizeY-1)*sizeX, iBufSize);

		// 下边界的首尾赋值
		for (int j=1; j<=r; j++)
		{
			*(pNew+sizeY*iNewX-j) = *(pSrc+(sizeY-1)*sizeX);
			*(pNew+sizeY*iNewX+sizeX-1+j) = *(pSrc+(sizeY-1)*sizeX + sizeX-1);
		}
	}

	// 扩展新数组左边界和右边界
	for (int y=0; y<sizeY; y++)
	{
		for (int i=1; i<=r; i++)
		{
			// 左边界
			*(pNew + y*iNewX-i) = *(pSrc + y*sizeX);
			// 右边界
			*(pNew + y*iNewX+sizeX-1+i) = *(pSrc + y*sizeX+sizeX-1);
		}
	}

	float CoefArray[9] = { 1, 1, 1,   
						   1, 1, 1,   
						   1, 1, 1 };
	float coef = 1.0/9;
	float TempNum = 0;
	pSrc = NULL;
	pNew = NULL;
	for (int y=r; y<iNewY-r; y++)// r为窗口半径
	{
		for (int x=r; x<iNewX-r; x++)
		{
			pNew = pNewDB + y*iNewX + x;
			pSrc = pSrcDB + (y-r)*sizeX + (x-r);

			TempNum = (*(pNew-iNewX-r))*CoefArray[0];
			TempNum += (*(pNew-iNewX))*CoefArray[1];
			TempNum += (*(pNew-iNewX+r))*CoefArray[2];
			TempNum += (*(pNew-r))*CoefArray[3];
			TempNum += (*(pNew))*CoefArray[4];
			TempNum += (*(pNew+r))*CoefArray[5];
			TempNum += (*(pNew+iNewX-r))*CoefArray[6];
			TempNum += (*(pNew+iNewX))*CoefArray[7];
			TempNum += (*(pNew+iNewX+r))*CoefArray[8];

			TempNum *= coef;			*pSrc = static_cast<int>(TempNum);
		}
	}

	GDALDriver *pDriver = (GDALDriver *) GDALGetDriverByName("GTiff");
	GDALDataset *pDstDs = pDriver->Create(strDstFile.c_str(), sizeX, sizeY, 1, GDT_Int32, NULL);
	double adfGeoTransform[6] = {0};
	pSrcDs->GetGeoTransform(adfGeoTransform);
	pDstDs->SetGeoTransform(adfGeoTransform);
	pDstDs->SetProjection(pSrcDs->GetProjectionRef());

	int panMap[1] = {1};
	pDstDs->RasterIO(GF_Write, 0, 0, sizeX, sizeY, pSrcDB, sizeX, sizeY, GDT_Int32, 1, panMap, 0, 0, 0);

	delete []pSrcDB;
	delete []pNewDB;
	pSrcDB = NULL;
	pNewDB = NULL;

	GDALClose((GDALDatasetH)pSrcDs);
	GDALClose((GDALDatasetH)pDstDs);

	return true;
}

bool LowPassMxN(string strSrcFile, string strDstFile, int m, int n, float *fCoefArray)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	GDALDataset *pSrcDs = (GDALDataset *)GDALOpen(strSrcFile.c_str(), GA_ReadOnly);
	if (pSrcDs == NULL)
	{
		printf("打开文件%s失败!\n", strSrcFile.c_str());
		return false;
	}

	int iBbandCount = pSrcDs->GetRasterCount();
	GDALRasterBand *pSrcBand = pSrcDs->GetRasterBand(1);
	if (pSrcBand == NULL)
	{
		printf("打开波段1失败!\n");
		return false;
	}

	int sizeX = pSrcBand->GetXSize();
	int sizeY = pSrcBand->GetYSize();
	int imgDims = sizeX * sizeY;

	int *pSrcDB = new int[imgDims];
	pSrcBand->RasterIO(GF_Read, 0, 0, sizeX, sizeY, pSrcDB, sizeX, sizeY, GDT_Int32, 0, 0);

	int rr = (n-1)/2; // 模板行半径y, m*n模板窗口
	int rc = (m-1)/2; // 模板列半径x

	int iNewX = sizeX + 2*rc;
	int iNewY = sizeY + 2*rr;
	int iNewDims = iNewX*iNewY;
	int *pNewDB = new int[iNewDims];
	memset(pNewDB, 0, sizeof(int)*iNewDims);

	int *pSrc = pSrcDB;
	int *pNew = pNewDB + rr*iNewX + rc;// 与原始影像左上对齐
	int iBufSize = sizeof(int)*sizeX;

	// 将原始影像数据复制到新数组中相应位置
	for(int y=0; y<sizeY; y++)
		memcpy(pNew+y*iNewX, pSrc+y*sizeX, iBufSize);

	// 扩展新数组上边界
	for (int i=1; i<=rr; i++)
	{
		memcpy(pNew-i*iNewX, pSrc, iBufSize);

		// 上边界的首尾赋值
		for (int j=1; j<=rc; j++)
		{
			// 首赋值
			*(pNew-i*iNewX-j) = *pSrc;
			// 尾赋值
			*(pNew-i*iNewX+sizeX-1+j) = *(pSrc + sizeX-1);
		}
	}


	// 扩展新数组下边界
	for (int i=1; i<=rr; i++)
	{
		memcpy(pNew+(sizeY-1+i)*iNewX, pSrc+(sizeY-1)*sizeX, iBufSize);

		// 下边界的首尾赋值
		for (int j=1; j<=rc; j++)
		{
			// 首赋值
			*(pNew+(sizeY-1+i)*iNewX-j) = *(pSrc+(sizeY-1)*sizeX);
			// 尾赋值
			*(pNew+(sizeY-1+i)*iNewX+sizeX-1+j) = *(pSrc+(sizeY-1)*sizeX + sizeX-1);
		}
	}

	// 扩展新数组左边界和右边界
	for (int y=0; y<sizeY; y++)
	{
		for (int i=1; i<=rc; i++)
		{
			// 左边界
			*(pNew + y*iNewX-i) = *(pSrc + y*sizeX);
			// 右边界
			*(pNew + y*iNewX+sizeX-1+i) = *(pSrc + y*sizeX+sizeX-1);
		}
	}


	float TempNum = 0;
	float coef = 1.0/25;
	pSrc = NULL;
	pNew = NULL;
	for (int y=rr; y<iNewY-rr; y++)
	{
		for (int x=rc; x<iNewX-rc; x++)
		{
			pNew = pNewDB + y*iNewX + x;
			pSrc = pSrcDB + (y-rr)*sizeX + (x-rc);
			//TempNum = GetCenterValue(pNew, fCoefArray, rr, rc, iNewX);
			TempNum *= coef;
			*pSrc = static_cast<int>(TempNum);
		}
	}

	GDALDriver *pDriver = (GDALDriver *) GDALGetDriverByName("GTiff");
	GDALDataset *pDstDs = pDriver->Create(strDstFile.c_str(), sizeX, sizeY, 1, GDT_Int32, NULL);
	double adfGeoTransform[6] = {0};
	pSrcDs->GetGeoTransform(adfGeoTransform);
	pDstDs->SetGeoTransform(adfGeoTransform);
	pDstDs->SetProjection(pSrcDs->GetProjectionRef());

	int panMap[1] = {1};
	pDstDs->RasterIO(GF_Write, 0, 0, sizeX, sizeY, pSrcDB, sizeX, sizeY, GDT_Int32, 1, panMap, 0, 0, 0);

	delete []pSrcDB;
	delete []pNewDB;
	pSrcDB = NULL;
	pNewDB = NULL;

	GDALClose((GDALDatasetH)pSrcDs);
	GDALClose((GDALDatasetH)pDstDs);

	return true;
}


template<typename T>
void BorderExtendT(T *pSrcDB, T *pNewDB, int srcSizeX, int srcSizeY, int newSizeX, 
				  int rr, int rc)
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


template<typename T>
double GetCenterValueT(T *pCenterDataBuf, float *pCoefArray, int rRow, int rCol, int newDataBufX)
{
	double resValue=0;
	for (int y=rRow; y>=-rRow; y--)// 行
	{
		for (int x=rCol; x>=-rCol; x--)// 列
			resValue += (*(pCenterDataBuf-y*newDataBufX-x))*(*(pCoefArray++));// 窗口相乘
	}

	return resValue;
}

template<typename T>
void LowPassMxNGenericT(GDALDataset *pSrcDs, GDALDataset *pDstDs, GDALDataType eDataType, int m, int n, float *pCoefArray, 
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

		// 填充数组扩展边界--用原始影像最外一圈的值赋给扩展的边界
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


bool LowPassMxNGeneric(string strSrcFile, string strDstFile, int m, int n, float *pCoefArray, float coef/*=1*/)
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
		LowPassMxNGenericT<unsigned char>(pSrcDs, pDstDs, eDataType, m, n, pCoefArray, 
			coef, nBandCount, panBandMap);
		break;
	case GDT_UInt16:
		LowPassMxNGenericT<unsigned short int>(pSrcDs, pDstDs, eDataType, m, n, pCoefArray, 
			coef, nBandCount, panBandMap);
		break;
	case GDT_Int16:
		LowPassMxNGenericT<short int>(pSrcDs, pDstDs, eDataType, m, n, pCoefArray, 
			coef, nBandCount, panBandMap);
		break;
	case GDT_UInt32:
		LowPassMxNGenericT<unsigned int>(pSrcDs, pDstDs, eDataType, m, n, pCoefArray, 
			coef, nBandCount, panBandMap);
		break;
	case GDT_Int32:
		LowPassMxNGenericT<int>(pSrcDs, pDstDs, eDataType, m, n, pCoefArray, 
			coef, nBandCount, panBandMap);
		break;
	case GDT_Float32:
		LowPassMxNGenericT<float>(pSrcDs, pDstDs, eDataType, m, n, pCoefArray, 
			coef, nBandCount, panBandMap);
		break;
	case GDT_Float64:
		LowPassMxNGenericT<double>(pSrcDs, pDstDs, eDataType, m, n, pCoefArray, 
			coef, nBandCount, panBandMap);
		break;
	}

	GDALClose((GDALDatasetH)pSrcDs);
	GDALClose((GDALDatasetH)pDstDs);

	return true;
}






































































































//float GetCenterValue(int *pCenterDataBuf, float *pCoefArray, int rRow, int rCol, int newDataBufX)
//{
//	float resValue=0;
//	for (int y=rRow; y>=-rRow; y--)// 行
//	{
//		for (int x=rCol; x>=-rCol; x--)// 列
//			resValue += (*(pCenterDataBuf-y*newDataBufX-x))*(*(pCoefArray++));// 窗口相乘
//	}
//
//	return resValue;
//}

//int iBbandCount = pSrcDs->GetRasterCount();
//GDALRasterBand *pSrcBand = pSrcDs->GetRasterBand(1);
//if (pSrcBand == NULL)
//{
//	printf("打开波段1失败!\n");
//	return false;
//}
//
//int sizeX = pSrcBand->GetXSize();
//int sizeY = pSrcBand->GetYSize();
//int imgDims = sizeX * sizeY;
//
//int *pSrcDB = new int[imgDims];
//pSrcBand->RasterIO(GF_Read, 0, 0, sizeX, sizeY, pSrcDB, sizeX, sizeY, GDT_Int32, 0, 0);
//
//int rr = (n-1)/2; // 模板行半径y, m*n模板窗口
//int rc = (m-1)/2; // 模板列半径x
//
//int iNewX = sizeX + 2*rc;
//int iNewY = sizeY + 2*rr;
//int iNewDims = iNewX*iNewY;
//int *pNewDB = new int[iNewDims];
//memset(pNewDB, 0, sizeof(int)*iNewDims);
//
//int *pSrc = pSrcDB;
//int *pNew = pNewDB + rr*iNewX + rc;// 与原始影像左上对齐
//int iBufSize = sizeof(int)*sizeX;
//
//// 将原始影像数据复制到新数组中相应位置
//for(int y=0; y<sizeY; y++)
//memcpy(pNew+y*iNewX, pSrc+y*sizeX, iBufSize);
//
//// 扩展新数组上边界
//for (int i=1; i<=rr; i++)
//{
//	memcpy(pNew-i*iNewX, pSrc, iBufSize);
//
//	// 上边界的首尾赋值
//	for (int j=1; j<=rc; j++)
//	{
//		// 首赋值
//		*(pNew-i*iNewX-j) = *pSrc;
//		// 尾赋值
//		*(pNew-i*iNewX+sizeX-1+j) = *(pSrc + sizeX-1);
//	}
//}
//
//
//// 扩展新数组下边界
//for (int i=1; i<=rr; i++)
//{
//	memcpy(pNew+(sizeY-1+i)*iNewX, pSrc+(sizeY-1)*sizeX, iBufSize);
//
//	// 下边界的首尾赋值
//	for (int j=1; j<=rc; j++)
//	{
//		// 首赋值
//		*(pNew+(sizeY-1+i)*iNewX-j) = *(pSrc+(sizeY-1)*sizeX);
//		// 尾赋值
//		*(pNew+(sizeY-1+i)*iNewX+sizeX-1+j) = *(pSrc+(sizeY-1)*sizeX + sizeX-1);
//	}
//}
//
//// 扩展新数组左边界和右边界
//for (int y=0; y<sizeY; y++)
//{
//	for (int i=1; i<=rc; i++)
//	{
//		// 左边界
//		*(pNew + y*iNewX-i) = *(pSrc + y*sizeX);
//		// 右边界
//		*(pNew + y*iNewX+sizeX-1+i) = *(pSrc + y*sizeX+sizeX-1);
//	}
//}
//
//
//float TempNum = 0;
//float coef = 1.0/25;
//pSrc = NULL;
//pNew = NULL;
//for (int y=rr; y<iNewY-rr; y++)
//{
//	for (int x=rc; x<iNewX-rc; x++)
//	{
//		pNew = pNewDB + y*iNewX + x;
//		pSrc = pSrcDB + (y-rr)*sizeX + (x-rc);
//		//TempNum = GetCenterValue(pNew, fCoefArray, rr, rc, iNewX);
//		TempNum *= coef;
//		*pSrc = static_cast<int>(TempNum);
//	}
//}















































/*
for(int i = 0; i < r; i++)
{
	int sCurRow = std::max(0,i-rhSe);// rhSe--模板窗口行半径；s--start，e-end
	int eCurRow = std::min(r-1,i+rhSe);
	float *ptrPSI = psiImg->ptr<float>(i);// 取行指针
	
	for(int j = 0;j< c;j++)
	{
		int sCurCol = std::max(0,j-chSe);// chSe--模板窗口列半径
		int eCurCol = std::min(c-1,j+chSe);
		//获取以i,j为中心点的矩形感兴区域存储为imgROI
		int width = eCurCol - sCurCol + 1;
		int height = eCurRow - sCurRow + 1;
		// 加入判断保持一致
		if(width>=se.cols)// se--模板
			width = se.cols;
		if(height>=se.rows)
			height = se.rows;

		int rtSeSY =(se.rows - height)/2;
		int rtSeSX = (se.cols - width)/2;
		cv::Rect rtSe(rtSeSX,rtSeSY,width,height);
		se = se(rtSe).clone();

		cv::Rect rt(sCurCol,sCurRow,width,height);
		cv::Mat imgROI;
		imgROI = img(rt).clone();
		//计算该方向同质区域的延伸长度
		std::vector<cv::Mat> vecROI;
		cv::split(imgROI,vecROI);
		cv::Mat tmp = cv::Mat(height,width,CV_32F,cv::Scalar(0));
		int k = 0;
		while(k<bands)
		{ 
			tmp = tmp + vecROI[k].mul(se);
			tmp = tmp - img.at<float>(i,j);
			tmp = cv::abs(tmp);
			++k;
		}// while
		int tmpLen = getLength(tmp,thresd1*bands);
		*(ptrPSI+j) += float(tmpLen);
	}// for j
}// for i
*/


// 扩展新数组上边界
//for (int i=1; i<=r; i++)
//{
//	memcpy(pNew-i*iNewX, pSrc, iBufSize);
//
//	// 上边界的首尾赋值
//	for (int j=1; j<=r; j++)
//	{
//		*(pNew-i*iNewX-j) = *pSrc;
//		*(pNew-i*iNewX+sizeX-1+j) = *(pSrc + sizeX-1);
//		//pNew[-i*iNewX-j] = pSrc[0];
//		//pNew[-i*iNewX+sizeX-1+j] = pSrc[sizeX-1];
//	}
//}
//
//
//// 扩展新数组下边界
//for (int i=1; i<=r; i++)
//{
//	memcpy(pNew+sizeY*iNewX, pSrc+(sizeY-1)*sizeX, iBufSize);
//
//	// 下边界的首尾赋值
//	for (int j=1; j<=r; j++)
//	{
//		*(pNew+sizeY*iNewX-j) = *(pSrc+(sizeY-1)*sizeX);
//		*(pNew+sizeY*iNewX+sizeX-1+j) = *(pSrc+(sizeY-1)*sizeX + sizeX-1);
//		//pNew[sizeY*iNewX-j] = pSrc[(sizeY-1)*sizeX];
//		//pNew[sizeY*iNewX+sizeX-1+j] = pSrc[(sizeY-1)*sizeX + sizeX-1];
//	}
//}
//
//// 扩展新数组左边界和右边界
//for (int y=0; y<sizeY; y++)
//{
//	for (int i=1; i<=r; i++)
//	{
//		// 左边界
//		*(pNew + y*iNewX-i) = *(pSrc + y*sizeX);
//		//pNew[y*iNewX-i] = pSrc[y*sizeX];
//		// 右边界
//		*(pNew + y*iNewX+sizeX-1+i) = *(pSrc + y*sizeX+sizeX-1);
//		//pNew[y*iNewX+sizeX-1+i] = pSrc[y*sizeX+sizeX-1];
//	}
//}

			
			// 3*3窗口
			/*TempNum = pNew[-iNewX-r]*CoefArray[0];
			TempNum += pNew[-iNewX]*CoefArray[1];
			TempNum += pNew[-iNewX+r]*CoefArray[2];
			TempNum += pNew[-r]*CoefArray[3];
			TempNum += pNew[0]*CoefArray[4];
			TempNum += pNew[r]*CoefArray[5];
			TempNum += pNew[iNewX-r]*CoefArray[6];
			TempNum += pNew[iNewX]*CoefArray[7];
			TempNum += pNew[iNewX+r]*CoefArray[8];
			TempNum *= coef;			*pSrc = static_cast<int>(TempNum);*//*是否需要单独写成一个函数？？？*/
//T *pSrcTemp = pSrcDB;
//T *pNewTemp = pNewDB + rr*newSizeX + rc;// 与原始影像左上对齐
//int srcLineBytes = sizeof(T)*srcSizeX;
//
//// 将原始影像数据复制到新数组中相应位置
//for(int y=0; y<srcSizeY; y++)
//memcpy(pNewTemp+y*newSizeX, pSrcTemp+y*srcSizeX, srcLineBytes);
//
//// 扩展新数组上边界
//for (int i=1; i<=rr; i++)
//{
//	memcpy(pNewTemp-i*newSizeX, pSrcTemp, srcLineBytes);
//
//	// 上边界的首尾赋值
//	for (int j=1; j<=rc; j++)
//	{
//		// 首赋值
//		*(pNewTemp-i*newSizeX-j) = *pSrcTemp;
//		// 尾赋值
//		*(pNewTemp-i*newSizeX+srcSizeX-1+j) = *(pSrcTemp + srcSizeX-1);
//	}
//}
//
//
//// 扩展新数组下边界
//for (int i=1; i<=rr; i++)
//{
//	memcpy(pNewTemp+(srcSizeY-1+i)*newSizeX, pSrcTemp+(srcSizeY-1)*srcSizeX, srcLineBytes);
//
//	// 下边界的首尾赋值
//	for (int j=1; j<=rc; j++)
//	{
//		// 首赋值
//		*(pNewTemp+(srcSizeY-1+i)*newSizeX-j) = *(pSrcTemp+(srcSizeY-1)*srcSizeX);
//		// 尾赋值
//		*(pNewTemp+(srcSizeY-1+i)*newSizeX+srcSizeX-1+j) = *(pSrcTemp+(srcSizeY-1)*srcSizeX + srcSizeX-1);
//	}
//}
//
//// 扩展新数组左边界和右边界
//for (int y=0; y<srcSizeY; y++)
//{
//	for (int i=1; i<=rc; i++)
//	{
//		// 左边界
//		*(pNewTemp + y*newSizeX-i) = *(pSrcTemp + y*srcSizeX);
//		// 右边界
//		*(pNewTemp + y*newSizeX+srcSizeX-1+i) = *(pSrcTemp + y*srcSizeX+srcSizeX-1);
//	}
//}