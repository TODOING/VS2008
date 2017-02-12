#include "StdRXD.h"
#include "common.h"

#include "gdal_priv.h"

CHighSpecRXD::CHighSpecRXD(const char *pszSrcFile)
{
	m_pBandMean		= NULL;
	m_pCovarianceElem	= NULL;
	m_pSrcDS		= NULL;

	m_pszSrcFile = pszSrcFile;
}

CHighSpecRXD::~CHighSpecRXD()
{
	if (m_pSrcDS != NULL)
		GDALClose( (GDALDatasetH) m_pSrcDS );

	RELEASE(m_pBandMean);
	RELEASE(m_pBandStad);
	RELEASE(m_pCovarianceElem);
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
	m_pBandMean = new double [m_iBandCount];
	m_pBandStad = new double [m_iBandCount];
	CVector meanVector(m_iBandCount);

	for (int i=1; i<=m_iBandCount; i++)	//获取每个波段的均值 //和标准差
	{
		double dMaxValue, dMinValue;
		m_pSrcDS->GetRasterBand(i)->ComputeStatistics(FALSE, &dMinValue, &dMaxValue, 
			m_pBandMean+(i-1), m_pBandStad+(i-1), NULL, NULL);
		meanVector[i-1] = *(m_pBandMean+(i-1));
	}
	m_MeanVector = meanVector;

	return RE_SUCCESS;
}

int CHighSpecRXD::CalcCovarianceMartix()
{
	if (m_pSrcDS == NULL)
	{
		/*if (m_pProcess != NULL)
			m_pProcess->SetMessage("输入文件不能打开！");*/

		return RE_FILENOTEXIST;
	}

	int iWidth = m_pSrcDS->GetRasterXSize();
	int iHeight = m_pSrcDS->GetRasterYSize();

	int iImageSize = iWidth * iHeight;

	//int iElementNum = m_iBandCount*m_iBandCount;	//相关系数矩阵中的个数，只保存对角阵，因为是实对称阵
	//int iElementIndex = 0;								//用于遍历相关系数矩阵中元素的索引
	m_pCovarianceElem = new double [m_iBandCount*m_iBandCount];			//分配相关系数矩阵的大小

	RMatrix covMatrix(m_iBandCount, m_iBandCount);	//相关系数矩阵

	for (int i1=1; i1<=m_iBandCount; i1++)
	{
		for (int i2=1; i2<=m_iBandCount; i2++)
		{
			if (i2<i1)
				continue;

			if (i1 == i2)	//若是同一个波段，不用计算了
			{
				//if (!m_bIsCovariance)
				//	m_pRelativity[iElementIndex] = 1.0;		//相关系数矩阵
				//else
					//m_pRelativity[iElementIndex] = m_pBandStad[i1-1]*m_pBandStad[i1-1];	//方差-协方差矩阵
					covMatrix(i1-1, i2-1) = m_pBandStad[i1-1]*m_pBandStad[i1-1];

				//iElementIndex++;
				continue;
			}

			GDALRasterBand *ptrBandI1 = m_pSrcDS->GetRasterBand(i1);
			GDALRasterBand *ptrBandI2 = m_pSrcDS->GetRasterBand(i2);

			double *pBuff1 = new double[iWidth];
			double *pBuff2 = new double[iWidth];

			double dTemp = 0.0;
			for(int j=0; j<iHeight; j++)//行
			{
				ptrBandI1->RasterIO(GF_Read, 0, j, iWidth, 1, pBuff1, iWidth, 1, GDT_Float64, 0, 0);		//读取第一波段数据块
				ptrBandI2->RasterIO(GF_Read, 0, j, iWidth, 1, pBuff2, iWidth, 1, GDT_Float64, 0, 0);		//读取第二波段数据块

				for (int i=0; i<iWidth; i++)
					dTemp += ((pBuff1[i] - m_pBandMean[i1-1]) * (pBuff2[i] - m_pBandMean[i2-1]));
			}

			RELEASE(pBuff1);
			RELEASE(pBuff2);

			//m_pRelativity[iElementIndex] = dTemp / iImageSize;	//方差-协方差矩阵
			covMatrix(i1-1, i2-1) = dTemp / iImageSize;
			covMatrix(i2-1, i1-1) = dTemp / iImageSize;

			//if (!m_bIsCovariance)
				//m_pRelativity[iElementIndex] = m_pRelativity[iElementIndex] / (m_pBandStad[i1-1]*m_pBandStad[i2-1]); //相关系数矩阵

			//iElementIndex++;
		}
	}

	m_CovarianceMatrix = covMatrix;
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

	int imgWidth = m_pSrcDS->GetRasterXSize();
	int imgHeight = m_pSrcDS->GetRasterYSize();

	double** pInputData = new double*[m_iBandCount];
	for(int i = 0; i < m_iBandCount; i++)
	{
		pInputData[i] = new double[imgWidth*imgHeight];
		assert(pInputData[i]);
	} 

	for (int m=0; m<m_iBandCount; m++)//读取采样后的数据？？？？？？？为什么要重采样？？
	{
		GDALRasterBand *pBand = m_pSrcDS->GetRasterBand(m+1);
		CPLErr result = pBand->RasterIO(GF_Read, 0, 0, imgWidth, imgHeight, pInputData[m], imgWidth, imgHeight, GDT_Float64, 0, 0);

		if(CE_Failure == result)
		{
			for(int i = 0; i < m_iBandCount; i++)
			{
				delete[]pInputData[i];
				pInputData[i] = NULL;
			}

			return CE_Failure;
		}
	}

	CVector ruVector(m_iBandCount);// (r-u)
	RVector rutVector(m_iBandCount);// (r-u)T

	double *pOutData = new double[imgWidth*imgHeight];
	RMatrix revCovarianceMatrix = m_CovarianceMatrix.inverse();

	int offset = 0;
	for (int row=0; row<imgHeight; row++)
	{
		for (int col=0; col<imgWidth; col++)
		{
			//double *pBuff1 = new double[imgWidth*m_iBandCount];

			for (int band=0; band<m_iBandCount; band++)
			{
				//float meanss = KMeans(b,0);
				ruVector[band] = pInputData[band][offset + col] - m_MeanVector(band);// (r - u)
				rutVector[band] = pInputData[band][offset + col] - m_MeanVector(band);// (r-u)T
				//matrixline(0,b) = pInputData[band][offset + col] - KMeans(b,0);// (r - u)T，相当于转置
			}

			RMatrix resRXD(1,1);
			resRXD = rutVector*revCovarianceMatrix*ruVector;

			pOutData[offset + col] = resRXD(0, 0);

			//cbmatix = cbmatix + (matrixcol*matrixline);//bands*bands // 1. 先求和
		}
		offset += imgWidth;
	}

	GDALDriver *pDriver = (GDALDriver *) GDALGetDriverByName("GTiff");
	GDALDataset *pDstDs = pDriver->Create(pszRXDFile, imgWidth, imgHeight, 1, GDT_Float64, NULL);
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
	pDstBand->RasterIO(GF_Write, 0, 0, imgWidth, imgHeight, pOutData, imgWidth, imgHeight, GDT_Float64, 0, 0);


	for(int i = 0; i < m_iBandCount; i++)
	{
		delete[]pInputData[i];
		pInputData[i] = NULL;
	}
	delete[] pInputData;
	delete []pOutData;
	pInputData = NULL;
	pOutData = NULL;

	GDALClose((GDALDatasetH)pDstDs);
	return RE_SUCCESS;
}
