#include "StdRXD.h"
#include "common.h"

#include "gdal_priv.h"
#include "AlgProcessTime.h"

using namespace std;

CHighSpecRXD::CHighSpecRXD(const char *pszSrcFile)
{
	m_pSrcDS		= NULL;
	m_pSrcDataBuf = NULL;
	m_pBandStad = NULL;
	m_pszSrcFile = pszSrcFile;
	m_iBandCount = 0;
	m_imgHeight = 0;
	m_imgWidth = 0;		
}

CHighSpecRXD::~CHighSpecRXD()
{
	if (m_pSrcDS != NULL)
		GDALClose( (GDALDatasetH) m_pSrcDS );

	RELEASE(m_pBandStad);
	RELEASE(m_pSrcDataBuf);
}

int CHighSpecRXD::PreProcessData()
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	m_pSrcDS = (GDALDataset*) GDALOpen(m_pszSrcFile, GA_ReadOnly);
	if (m_pSrcDS == NULL)
	{
		return RE_FILENOTEXIST;
	}

	m_iBandCount = m_pSrcDS->GetRasterCount();
	m_imgWidth = m_pSrcDS->GetRasterXSize();
	m_imgHeight = m_pSrcDS->GetRasterYSize();
	int iImageSize = m_imgWidth*m_imgHeight;

	m_pBandStad = new double [m_iBandCount];
	m_pSrcDataBuf = new double[m_iBandCount*m_imgWidth*m_imgHeight];

	int *panBandMap = new int[m_iBandCount];
	for (int i=1; i<=m_iBandCount; ++i)
		panBandMap[i-1] = i;

	CAlgProcessTime::Alg_start();
	m_pSrcDS->RasterIO(GF_Read, 0, 0, m_imgWidth, m_imgHeight, m_pSrcDataBuf, 
		m_imgWidth, m_imgHeight, GDT_Float64, m_iBandCount, panBandMap, 0, 0, 0);
	CAlgProcessTime::Alg_end();
	printf("读取数据时间：%lf s\n", CAlgProcessTime::GetAlgProcessTime());

	delete []panBandMap;
	panBandMap = NULL;

	double *pBuf = NULL;
	CVector meanVector(m_iBandCount);

	CAlgProcessTime::Alg_start();
	for (int i=1; i<=m_iBandCount; i++)
	{
		double dMaxValue, dMinValue;
		//pBuf = m_pSrcDataBuf+iImageSize*(i-1);// 定位到每一个波段缓冲区

		/*GDALRasterBand *pBand = m_pSrcDS->GetRasterBand(i);
		pBand->ComputeStatistics(TRUE, &dMinValue, &dMaxValue, 
			&meanVector[i-1], m_pBandStad+(i-1), NULL, NULL);*/

		m_pSrcDS->GetRasterBand(i)->ComputeStatistics(TRUE, &dMinValue, &dMaxValue, 
			&meanVector[i-1], m_pBandStad+(i-1), NULL, NULL);
		
		// 读取全部数据至内存
		//pBand->RasterIO(GF_Read, 0, 0, m_imgWidth, m_imgHeight, pBuf, m_imgWidth, m_imgHeight, GDT_Float64, 0, 0);
	}
	CAlgProcessTime::Alg_end();
	printf("统计时间：%lf s\n", CAlgProcessTime::GetAlgProcessTime());
	m_MeanVector = meanVector;

	return RE_SUCCESS;
}

int CHighSpecRXD::CalcInvCovarianceMartix()
{
	if (m_pSrcDS == NULL)
	{
		return RE_FILENOTEXIST;
	}

	int iImageSize = m_imgWidth * m_imgHeight;
	RMatrix covMatrix(m_iBandCount, m_iBandCount);	//协方差矩阵
	for (int i=0; i<m_iBandCount; ++i)
	{
		for (int j=0; j<m_iBandCount; ++j)
		{
			covMatrix(i, j) = 0;
		}
	}

	RVector rowVector(m_iBandCount);
	CVector colVector(m_iBandCount);

	CAlgProcessTime::Alg_start();
	for (int i=0; i<iImageSize; ++i)
	{
		for (int j=0; j<m_iBandCount; ++j)
		{
			double *pBuff1 = m_pSrcDataBuf + iImageSize*j;

			//rowVector[j] = pBuff1[i] - m_MeanVector[j];
			colVector[j] = pBuff1[i] - m_MeanVector[j];
			
		}
		m_colVector.push_back(colVector);

		//covMatrix = covMatrix + colVector*rowVector;
		covMatrix = covMatrix + colVector*colVector.transpose();
	}
	covMatrix /= iImageSize;
	m_InvCovarianceMatrix = covMatrix.inverse();
	CAlgProcessTime::Alg_end();
	printf("计算协方差时间：%lf s\n", CAlgProcessTime::GetAlgProcessTime());

	/*CAlgProcessTime::Alg_start();
	for (int i=0; i<m_imgHeight; ++i)
	{
		for (int j=0; j<m_imgWidth; ++j)
		{
			for (int b=0; b<m_iBandCount; ++b)
			{
				double *pBuff1 = m_pSrcDataBuf + iImageSize*b;

				rowVector[b] = pBuff1[i*m_imgHeight + m_imgWidth] - m_MeanVector[b];
				colVector[b] = pBuff1[i*m_imgHeight + m_imgWidth] - m_MeanVector[b];
			}

			covMatrix = covMatrix + colVector*rowVector;
		}
	}
	covMatrix /= iImageSize;
	m_InvCovarianceMatrix = covMatrix.inverse();
	CAlgProcessTime::Alg_end();
	printf("计算协方差时间：%lf s\n", CAlgProcessTime::GetAlgProcessTime());*/


	//CAlgProcessTime::Alg_start();
	//for (int i1=1; i1<=m_iBandCount; i1++)
	//{
	//	for (int i2=1; i2<=m_iBandCount; i2++)
	//	{
	//		if (i2<i1)
	//			continue;

	//		if (i1 == i2)	//若是同一个波段，不用计算了
	//		{
	//			covMatrix(i1-1, i2-1) = m_pBandStad[i1-1]*m_pBandStad[i1-1];
	//			continue;
	//		}

	//		double *pBuff1 = m_pSrcDataBuf + iImageSize*(i1-1);
	//		double *pBuff2 = m_pSrcDataBuf + iImageSize*(i2-1);

	//		double dTemp = 0.0;
	//		for (int i=0; i<iImageSize; i++)
	//			dTemp += ((pBuff1[i] - m_MeanVector[i1-1]) * (pBuff2[i] - m_MeanVector[i2-1]));

	//		covMatrix(i1-1, i2-1) = dTemp / iImageSize;
	//		covMatrix(i2-1, i1-1) = dTemp / iImageSize;
	//	}
	//}
	/*CAlgProcessTime::Alg_end();
	printf("计算协方差时间：%lf s\n", CAlgProcessTime::GetAlgProcessTime());
	m_InvCovarianceMatrix = covMatrix.inverse();*/

	
	return RE_SUCCESS;
}

int CHighSpecRXD::ExecuteStdRXD(const char* pszRXDFile, int iBandCount /* = -1 */, const char* pszFormat /* = "GTiff" */)
{
	// 数据预处理，统计均值和协方差，以及读取影像数据至内存
	int iRev = PreProcessData();
	if(iRev != RE_SUCCESS)
		return iRev;

	// 计算协防差矩阵的逆矩阵
	iRev = CalcInvCovarianceMartix();
	if(iRev != RE_SUCCESS)
		return iRev;

	CVector ruVector(m_iBandCount);// (r-u)
	RVector rutVector(m_iBandCount);// (r-u)T

	double *pOutData = new double[m_imgWidth*m_imgHeight];
	int iImageSize = m_imgWidth * m_imgHeight;
	RMatrix resRXD(1,1);
	double *pBuf = NULL;

	// 计算RXD
	CAlgProcessTime::Alg_start();
	for (int i=0; i<iImageSize; ++i)
	{
		//for (int j=0; j<m_iBandCount; ++j)
		//{
		//	pBuf = m_pSrcDataBuf + iImageSize*j;// 定位到每个波段的起始位置

		//	ruVector[j] = pBuf[i] - m_MeanVector[j];// (r-u)
		//	rutVector[j] = pBuf[i] - m_MeanVector[j];// (r-u)T

		//	
		//}

		//resRXD = rutVector*m_InvCovarianceMatrix*ruVector;// 按公式计算
		//pOutData[i] = resRXD(0, 0);
		pOutData[i] = m_colVector[i].transpose()*m_InvCovarianceMatrix*m_colVector[i];
 	}
	CAlgProcessTime::Alg_end();
	printf("RXD计算时间：%lf s\n", CAlgProcessTime::GetAlgProcessTime());

	// 创建RXD文件
	GDALDriver *pDriver = (GDALDriver *) GDALGetDriverByName("GTiff");
	GDALDataset *pDstDs = pDriver->Create(pszRXDFile, m_imgWidth, m_imgHeight, 1, GDT_Float64, NULL);
	if (pDstDs == NULL)
	{
		printf("创建文件%s数据集失败!\n", pszRXDFile);
		return 0;
	}		
	double adfGeoTransform[6] = {0};
	m_pSrcDS->GetGeoTransform(adfGeoTransform);
	pDstDs->SetGeoTransform(adfGeoTransform);
	pDstDs->SetProjection(m_pSrcDS->GetProjectionRef());

	GDALRasterBand *pDstBand = pDstDs->GetRasterBand(1);
	pDstBand->RasterIO(GF_Write, 0, 0, m_imgWidth, m_imgHeight, pOutData, m_imgWidth, m_imgHeight, GDT_Float64, 0, 0);

	// 释放资源
	delete []pOutData;
	pOutData = NULL;
	GDALClose((GDALDatasetH)pDstDs);

	return RE_SUCCESS;
}
