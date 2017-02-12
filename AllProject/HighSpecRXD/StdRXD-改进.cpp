#include "StdRXD.h"
#include "common.h"

#include "gdal_priv.h"

CHighSpecRXD::CHighSpecRXD(const char *pszSrcFile)
{
	m_pBandMean		= NULL;
	m_pCovarianceElem	= NULL;
	m_pSrcDS		= NULL;

	m_pszSrcFile = pszSrcFile;

	m_pSrcDataBuf = NULL;
	m_imgHeight = 0;
	m_imgWidth = 0;
}

CHighSpecRXD::~CHighSpecRXD()
{
	if (m_pSrcDS != NULL)
		GDALClose( (GDALDatasetH) m_pSrcDS );

	//RELEASE(m_pBandMean);
	RELEASE(m_pBandStad);
	//RELEASE(m_pCovarianceElem);

	RELEASE(m_pSrcDataBuf);
}

int CHighSpecRXD::PreProcessData()
{
	GDALAllRegister();
	m_pSrcDS = (GDALDataset*) GDALOpen(m_pszSrcFile, GA_ReadOnly);
	if (m_pSrcDS == NULL)
	{
		/*if (m_pProcess != NULL)
			m_pProcess->SetMessage("输入文件不能打开！");*/

		return RE_FILENOTEXIST;
	}

	m_iBandCount = m_pSrcDS->GetRasterCount();
	m_imgWidth = m_pSrcDS->GetRasterXSize();
	m_imgHeight = m_pSrcDS->GetRasterYSize();
	int iImageSize = m_imgWidth*m_imgHeight;

	//m_pBandMean = new double [m_iBandCount];
	m_pBandStad = new double [m_iBandCount];
	m_pSrcDataBuf = new double[m_iBandCount*m_imgWidth*m_imgHeight];

	double *pBuf = NULL;
	CVector meanVector(m_iBandCount);
	for (int i=1; i<=m_iBandCount; i++)	//获取每个波段的均值 //和标准差
	{
		double dMaxValue, dMinValue;
		pBuf = m_pSrcDataBuf+iImageSize*(i-1);// 定位到每一个波段缓冲区
		GDALRasterBand *pBand = m_pSrcDS->GetRasterBand(i);
		/*pBand->ComputeStatistics(FALSE, &dMinValue, &dMaxValue, 
			m_pBandMean+(i-1), m_pBandStad+(i-1), NULL, NULL);*/

		pBand->ComputeStatistics(FALSE, &dMinValue, &dMaxValue, 
			&meanVector[i-1], m_pBandStad+(i-1), NULL, NULL);
		
		// 读取全部数据至内存
		pBand->RasterIO(GF_Read, 0, 0, m_imgWidth, m_imgHeight, pBuf, m_imgWidth, m_imgHeight, GDT_Float64, 0, 0);
		
		//meanVector[i-1] = *(m_pBandMean+(i-1));
	}
	m_MeanVector = meanVector;

	return RE_SUCCESS;
}

int CHighSpecRXD::CalcCovarianceMartix()
{
	if (m_pSrcDS == NULL)
	{
		return RE_FILENOTEXIST;
	}

	int iImageSize = m_imgWidth * m_imgHeight;

	RMatrix covMatrix(m_iBandCount, m_iBandCount);	//协方差矩阵

	for (int i1=1; i1<=m_iBandCount; i1++)
	{
		for (int i2=1; i2<=m_iBandCount; i2++)
		{
			if (i2<i1)
				continue;

			if (i1 == i2)	//若是同一个波段，不用计算了
			{
				covMatrix(i1-1, i2-1) = m_pBandStad[i1-1]*m_pBandStad[i1-1];
				continue;
			}

			double *pBuff1 = m_pSrcDataBuf + iImageSize*(i1-1);
			double *pBuff2 = m_pSrcDataBuf + iImageSize*(i2-1);

			double dTemp = 0.0;
			for (int i=0; i<iImageSize; i++)
				dTemp += ((pBuff1[i] - m_MeanVector[i1-1]) * (pBuff2[i] - m_MeanVector[i2-1]));

			covMatrix(i1-1, i2-1) = dTemp / iImageSize;
			covMatrix(i2-1, i1-1) = dTemp / iImageSize;
		}
	}

	m_InvCovarianceMatrix = covMatrix.inverse();
	return RE_SUCCESS;
}

int CHighSpecRXD::ExecuteStdRXD(const char* pszRXDFile, int iBandCount /* = -1 */, const char* pszFormat /* = "GTiff" */)
{
	int iRev = PreProcessData();
	if(iRev != RE_SUCCESS)
		return iRev;

	iRev = CalcCovarianceMartix();
	if(iRev != RE_SUCCESS)
		return iRev;

	CVector ruVector(m_iBandCount);// (r-u)
	RVector rutVector(m_iBandCount);// (r-u)T

	double *pOutData = new double[m_imgWidth*m_imgHeight];
	int iImageSize = m_imgWidth * m_imgHeight;
	RMatrix resRXD(1,1);
	double *pBuf = NULL;

	for (int i=0; i<iImageSize; ++i)
	{
		for (int j=0; j<m_iBandCount; ++j)
		{
			pBuf = m_pSrcDataBuf + iImageSize*j;// 定位到每个波段的起始位置

			ruVector[j] = pBuf[i] - m_MeanVector[j];// (r-u)
			rutVector[j] = pBuf[i] - m_MeanVector[j];// (r-u)T
		}

		resRXD = rutVector*m_InvCovarianceMatrix*ruVector;// 按公式计算
		pOutData[i] = resRXD(0, 0);
	}

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
