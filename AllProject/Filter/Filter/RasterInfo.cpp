#include "stdafx.h"
#include "RasterInfo.h"

CRasterInfo::CRasterInfo()
{

}
CRasterInfo::~CRasterInfo()
{

}
void CRasterInfo::Init(const std::vector<GDALDataType>& vecDataType, int imgWidth, int imgHeight, 
					   const string& strWkt, double* pGeoTransform)
{
	m_vecBandType = vecDataType;
	m_imgWidth = imgWidth;
	m_imgHeight = imgHeight;
	m_strWkt = strWkt;
	m_pGeoTrans = pGeoTransform;

}