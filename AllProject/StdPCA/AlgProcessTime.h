#ifndef ALGPROCESSTIME_H
#define ALGPROCESSTIME_H

#include <windows.h>

class CAlgProcessTime
{
public: 
	static void Alg_start();
	static void Alg_end();
	static void Alg_end(const char *info);
	static double GetAlgProcessTime();// 以s为单位
private:
	static LARGE_INTEGER m_nFreq;
	static LARGE_INTEGER m_nBeginTime;
	static LARGE_INTEGER m_nEndTime;
};

#endif // ALGPROCESSTIME_H