#ifndef RASTERINFO_H
#define RASTERINFO_H

#include <string>
#include "gdal_priv.h"

using std::string;
using std::vector;

enum enStoreType
{
	StoreBIP = 0,
	StoreBSQ = 1,
	StoreBIL = 2
};

// 定义一个栅格信息类
class CRasterInfo
{
public:
	CRasterInfo();
	CRasterInfo(const CRasterInfo& otherRasterInfo);
	~CRasterInfo();

	void Init(const std::vector<GDALDataType>& vecDataType,	int imgWidth, int imgHeight,
		const string& strWkt, double* pGeoTransform);

	CRasterInfo& operator=(const CRasterInfo& otherRasterInfo);

public: 
	int m_imgWidth;
	int m_imgHeight;
	string m_strWkt;
	double* m_pGeoTrans;

	bool GetGeoTransform(double* pGeoTransform) const;
	void SetGeoTransform(double* pGeoTransform);

	// 获取横向分辨率
	double GetHorResolution();
	// 获取纵向分辨率
	double GetVerResolution();

	
	enStoreType  m_StoreType; // 图像的存储类型：bsq,bip或bil,默认为bsq

	//std::vector<unsigned short> m_BitCount;
	std::vector<GDALDataType> m_vecBandType;  // 波段类型
	std::vector<double> m_vecInvalidValue; // 无效值
	std::vector<std::pair<short int, string> > m_vecBandDes; // 波段描述


};

// 获取波段的存储类型
//const char *pszItem = pSrcDs->GetMetadataItem("INTERLEAVE", "IMAGE_STRUCTURE");
//if (pszItem != NULL)
//{
//	if (EQUAL(pszItem, "BAND"))// bsq
//	{
//		//m_rasterInfo.m_StoreType = SysDataSource::StoreBSQ;
//
//	}
//	if (EQUAL(pszItem, "PIXEL"))// bip
//	{
//
//	}
//	if (EQUAL(pszItem, "LINE"))// bil
//	{
//
//	}
//}

#endif // RASTERINFO_H