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

// ����һ��դ����Ϣ��
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

	// ��ȡ����ֱ���
	double GetHorResolution();
	// ��ȡ����ֱ���
	double GetVerResolution();

	
	enStoreType  m_StoreType; // ͼ��Ĵ洢���ͣ�bsq,bip��bil,Ĭ��Ϊbsq

	//std::vector<unsigned short> m_BitCount;
	std::vector<GDALDataType> m_vecBandType;  // ��������
	std::vector<double> m_vecInvalidValue; // ��Чֵ
	std::vector<std::pair<short int, string> > m_vecBandDes; // ��������


};

// ��ȡ���εĴ洢����
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