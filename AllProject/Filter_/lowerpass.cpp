#include "gdal_priv.h"

int TemplateOperation(const char *pszSrcFile, const char *pszDstFile, int a)
{
	int sizeConv = 2*a + 1;

	GDALAllRegister();
	
	GDALDataset *pSrcDS = (GDALDataset *)GDALOpen(pszSrcFile, GA_ReadOnly);
	if (pSrcDS == NULL)
	{
		printf("文件%s打开失败！\n", pszSrcFile);
		return -1;
	}

	GDALRasterBand *pBand1 = pSrcDS->GetRasterBand(1);
	int sizeX = pBand1->GetXSize();
	int sizeY = pBand1->GetYSize();
	int tempX = sizeX + 2*a;
	int tempY = sizeY + 2*a;

	int *pSrcDataBuffer = new int[sizeX*sizeY];
	int *pTempDataBuffer = new int[tempX*tempY];
	float *pCoefArray = new float[sizeConv*sizeConv];

	memset(pTempDataBuffer, 0, tempX*tempY*sizeof(int));
	for (int i=0; i<sizeConv*sizeConv; i++)
		pCoefArray[i] = 1.0;
	float coef = 1.0 / 9;

	pBand1->RasterIO(GF_Read, 0, 0, sizeX, sizeY, pSrcDataBuffer, sizeX, sizeY, GDT_Int32, 0, 0);

	int *pSrc = NULL;
	int *pTemp = NULL;
	for (int y=1; y<tempY-1; y++)
	{
		for (int x=1; x<tempX-1; x++)
		{
			pTemp = pTempDataBuffer + y*tempX + x;
			pSrc = pSrcDataBuffer + (y-1)*sizeX + (x-1);

			*pTemp = *pSrc;
		}
	}

	pSrc = NULL;
	pTemp = NULL;
	float TempNum = 0;
	for (int y=1; y<tempY-1; y++)// 1 相当于a
	{
		for (int x=1; x<tempX-1; x++)
		{
			pTemp = pTempDataBuffer + y*tempX + x;
			pSrc = pSrcDataBuffer + (y-1)*sizeX + (x-1);

			TempNum=(float)(*(pTemp-tempX-1))*pCoefArray[0];
			TempNum+=(float)(*(pTemp-tempX))*pCoefArray[1];
			TempNum+=(float)(*(pTemp-tempX+1))*pCoefArray[2];
			TempNum+=(float)(*(pTemp-1))*pCoefArray[3];
			TempNum+=(float)(*pTemp)*pCoefArray[4];
			TempNum+=(float)(*(pTemp+1))*pCoefArray[5];
			TempNum+=(float)(*(pTemp+tempX-1))*pCoefArray[6];
			TempNum+=(float)(*(pTemp+tempX))*pCoefArray[7];
			TempNum+=(float)(*(pTemp+tempX+1))*pCoefArray[8];
			TempNum*=coef;
			*pSrc = static_cast<int>(TempNum);
		}
	}

	// 创建输出文件并设置空间参考和坐标信息
	GDALDriver *pDriver = (GDALDriver *) GDALGetDriverByName("GTiff");
	GDALDataset *pDstDS = pDriver->Create(pszDstFile, sizeX, sizeY, 1, GDT_Int32, NULL);
	double adfGeoTransform[6] = {0};
	pSrcDS->GetGeoTransform(adfGeoTransform);
	pDstDS->SetGeoTransform(adfGeoTransform);
	pDstDS->SetProjection(pSrcDS->GetProjectionRef());

	int panMap[1] = {1};
	pDstDS->RasterIO(GF_Write, 0, 0, sizeX, sizeY, pSrcDataBuffer, sizeX, sizeY, GDT_Int32, 1, panMap, 0, 0, 0);

	delete []pTempDataBuffer;
	delete []pSrcDataBuffer;
	delete []pCoefArray;

	GDALClose((GDALDatasetH)pSrcDS);
	GDALClose((GDALDatasetH)pDstDS);

	return 0;
}


