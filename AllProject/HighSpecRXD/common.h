/***************************************************************************
*
* Time: 2009-09-21
* Project: ң��ƽ̨
* Purpose: ���Ŀ��ļ�
* Author:  ����¼
* Copyright (c) 2009, liminlu0314@gmail.com
* Describe:�ṩ���õ��������Ͷ����
*
****************************************************************************/

#ifndef IMGALG_H
#define IMGALG_H

/**
* \file ImgCore.h
* @brief �������Ͷ���
*
* �����ӿڣ�ʹ��C���Է�ʽ�����������Ͷ���
*/

/**
* ������MS Widowsƽ̨�ϵľ��� lml 2010-10-19
* warning C4100: ��*��: δ���õ��β�
* warning C4190: ��identifier1����ָ���� C ���ӣ����������� C �����ݵ� UDT��identifier2��
* warning C4251: �ࡰtype����Ҫ���ࡰtype2���Ŀͻ���ʹ�� dll �ӿ�
* warning C4275: �� DLL �ӿ������identifier����Ϊ DLL �ӿ������identifier���Ļ�ʹ��
* warning C4305: �ӡ�type1������type2���ض�
* warning C4309: �ضϳ���ֵ
* warning C4819: ���ļ����������ڵ�ǰ����ҳ(936)�б�ʾ���ַ����뽫���ļ�����Ϊ Unicode ��ʽ�Է�ֹ���ݶ�ʧ
* warning C4996: ʹ���˷Ǳ�׼��չ: �޶�����ʹ����ö��
*/
#if defined(_MSC_VER) && (_MSC_VER >= 1400)
#pragma warning(disable: 4100 4190 4251 4275 4305 4309 4819 4996 )
#endif

/**
* @brief �Ƿ�ʹ��Vld�ڴ�й¶��⹤��
*/
#if _USEVLD
#if _DEBUG	//��debugģʽ�¼���ڴ�й¶
#include "vld.h"
#endif
#endif

/**
* @brief �Ƿ�ʹ��LOG���߽���д��־
*/
#if _USELOG
#define USE_LOG4CPP
#endif

#include <float.h>
#include <algorithm>
#include <deque>
#include <fstream>
#include <limits>
#include <map>
#include <stack>
#include <string>
#include <vector>
using namespace std;

/**
* @brief �������Ŷ���
*/
#ifdef IMGALG_EXPORTS
#define IMGALG_API __declspec(dllexport)
#else
#define IMGALG_API __declspec(dllimport)
#endif

/**
* @brief ����NULL
*/
#ifndef NULL
#  define NULL  0
#endif

/**
* @brief ����FALSE
*/
#ifndef FALSE
#  define FALSE 0
#endif

/**
* @brief ����TRUE
*/
#ifndef TRUE
#  define TRUE  1
#endif

#ifndef MAX
/*! �����ֵ */
#  define MIN(a, b)      ((a<b) ? a : b)
/*! ����Сֵ */
#  define MAX(a, b)      ((a>b) ? a : b)
#endif

/**
* @brief ����ABS�������ֵ
*/
#ifndef ABS
#  define ABS(x)        ((x<0) ? (-1*(x)) : x)
#endif

/**
* @brief ����PI=3.141592653...�Լ��Ⱥͻ���ת��
*/
#ifndef M_PI
/*! ����Բ����PI */
# define M_PI  3.1415926535897932384626433832795
/*! ����ת�� */
# define DEG_PER_RAD      ((double)(180.0/M_PI))
/*! ��ת���� */
# define RAD_PER_DEG      ((double)(M_PI/180.0))
#endif

/**
* @brief ����ƽ��
*/
#ifndef M_SQUARE
# define M_SQUARE(x)  (x)*(x)
#endif

/**
* @brief ��������
*/
#ifndef M_CUBE
# define M_CUBE(x)  (x)*(x)*(x)
#endif

/*! �жϸ������Ƿ�NaNֵ */
inline bool isnan(const float& v)  { return _isnan(v) ? true : false; }
/*! �ж�double���Ƿ�NaNֵ */
inline bool isnan(const double& v) { return _isnan(v) ? true : false; }
/*! ��ȡdouble��NaNֵ */
inline double nan() { return numeric_limits<double>::quiet_NaN(); }

/**
* @brief float���͵ļ�ֵ
*/
#ifndef FLT_EQUALS
/*! �������Ƿ���� */
#define FLT_EQUALS(x, y)  (fabs((double)x-y)<FLT_EPSILON)
/*! �������Ƿ����(ָ���Ƚ���ֵ) */
#define FLT_EQUALS_N(x, y, z)  (fabs((double)x-y)<z)
#endif

#ifndef FLT_ZERO
/*! �������Ƿ�Ϊ0 */
#define FLT_ZERO(x)  (fabs(x)<FLT_EPSILON)
#endif

/**
* @brief �ͷ�����
*/
#define RELEASE(x)	if(x!=NULL) {delete []x; x = NULL;}

/**
* @brief ��ȫ��������Ϊ����ϵͳĬ������
*/
#define SET_LOCAL	{ locale::global(locale("")); setlocale(LC_ALL,"Chinese-simplified"); }
/**
* @brief ��ԭȫ�������趨
*/
#define REVERT_LOCAL	locale::global(locale("C"))


#ifndef EQUAL
#if defined(WIN32) || defined(WIN32CE)
/*! �Ƚ��ַ����Ƿ���� */
#  define EQUALN(a, b, n)           (_strnicmp(a, b, n) == 0)
/*! �Ƚ��ַ����Ƿ���� */
#  define EQUAL(a, b)              (_stricmp(a, b) == 0)
#else
/*! �Ƚ��ַ����Ƿ���� */
#  define EQUALN(a, b, n)           (strncasecmp(a, b, n) == 0)
/*! �Ƚ��ַ����Ƿ���� */
#  define EQUAL(a, b)              (strcasecmp(a, b) == 0)
#endif
#endif

/*! byte */
typedef unsigned char	byte;
/*! 8U */
typedef unsigned char	DT_8U;
/*! 16U */
typedef unsigned short	DT_16U;
/*! 16S */
typedef short			DT_16S;
/*! 32U */
typedef unsigned int	DT_32U;
/*! 32S */
typedef int				DT_32S;
/*! 32F */
typedef float			DT_32F;
/*! 64F */
typedef double			DT_64F;

/*! �ɹ�ִ�� */
const int RE_SUCCESS		= 0;
/*! �ļ������� */
const int RE_FILENOTEXIST	= 1;
/*! �ļ���ʽ����֧�� */
const int RE_FILENOTSUPPORT	= 2;
/*! ͼ���������Ͳ���ȷ */
const int RE_FILETYPEERROR	= 3;
/*! ����ͼ��ʧ�� */
const int RE_CREATEFAILED	= 4;
/*! ����������� */
const int RE_PARAMERROR		= 5;
/*! �������� */
const int RE_FAILED			= 6;
/*! ͼ�񲻴��ڹ������� */
const int RE_NOSAMEEXTENT	= 7;
/*! �û�ȡ������ */
const int RE_USERCANCEL		= 8;
/*! �ļ��Ѿ���ʹ�� */
const int RE_FILEISUESED	= 9;
/*! ��֧�ֵ�������� */
const int RE_DEPTHNOTSUPPORT	= 10;
/*! ��������������Ҫ�� */
const int RE_BANDCOUNTERROR		= 11;
/*! �ļ�������ͶӰ */
const int RE_NOPROJECTION		= 12;
/*! ͶӰ��һ�� */
const int RE_PROJECTIONDIFF		= 13;

#endif