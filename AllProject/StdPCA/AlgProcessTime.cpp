//#include "stdafx.h"
#include "AlgProcessTime.h"
#include <stdio.h>

LARGE_INTEGER CAlgProcessTime::m_nFreq = { 0 };
LARGE_INTEGER CAlgProcessTime::m_nBeginTime = { 0 };
LARGE_INTEGER CAlgProcessTime::m_nEndTime = { 0 };

void CAlgProcessTime::Alg_start()
{
	QueryPerformanceCounter(&m_nBeginTime);
}

void CAlgProcessTime::Alg_end()
{
	QueryPerformanceCounter(&m_nEndTime);
}

void CAlgProcessTime::Alg_end(const char *info)
{
	QueryPerformanceCounter(&m_nEndTime);
	printf(info);
	printf(": [");
	printf("ºÄÊ± %lf Ãë]\n", GetAlgProcessTime());
}

double CAlgProcessTime::GetAlgProcessTime()
{
	QueryPerformanceFrequency(&m_nFreq);
	return (double)(m_nEndTime.QuadPart - m_nBeginTime.QuadPart)/(double)m_nFreq.QuadPart;
}