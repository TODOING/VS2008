/***************************************************************************
*
* Time: 2009-09-21
* Project: 遥感平台
* Purpose: 核心库文件
* Author:  李民录
* Copyright (c) 2009, liminlu0314@gmail.com
* Describe:提供常用的数据类型定义等
*
****************************************************************************/

#ifndef IMGALG_DEFINE_H
#define IMGALG_DEFINE_H

/**
* \file ImgCore.h
* @brief 核心类型定义
*
* 导出接口（使用C语言方式），核心类型定义
*/

/**
* 忽略在MS Widows平台上的警告 lml 2010-10-19
* warning C4100: “*”: 未引用的形参
* warning C4190: “identifier1”有指定的 C 链接，但返回了与 C 不兼容的 UDT“identifier2”
* warning C4251: 类“type”需要由类“type2”的客户端使用 dll 接口
* warning C4275: 非 DLL 接口类键“identifier”作为 DLL 接口类键“identifier”的基使用
* warning C4305: 从“type1”到“type2”截断
* warning C4309: 截断常数值
* warning C4819: 该文件包含不能在当前代码页(936)中表示的字符。请将该文件保存为 Unicode 格式以防止数据丢失
* warning C4996: 使用了非标准扩展: 限定名中使用了枚举
*/
#if defined(_MSC_VER) && (_MSC_VER >= 1400)
#pragma warning(disable: 4100 4190 4251 4275 4305 4309 4819 4996 )
#endif

/**
* @brief 是否使用Vld内存泄露监测工具
*/
#if _USEVLD
#if _DEBUG	//在debug模式下检测内存泄露
#include "vld.h"
#endif
#endif

/**
* @brief 是否使用LOG工具进行写日志
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
* @brief 导出符号定义
*/
#ifdef IMGALG_EXPORTS
#define IMGALG_API __declspec(dllexport)
#else
#define IMGALG_API __declspec(dllimport)
#endif

/**
* @brief 定义NULL
*/
#ifndef NULL
#  define NULL  0
#endif

/**
* @brief 定义FALSE
*/
#ifndef FALSE
#  define FALSE 0
#endif

/**
* @brief 定义TRUE
*/
#ifndef TRUE
#  define TRUE  1
#endif

#ifndef MAX
/*! 求最大值 */
#  define MIN(a, b)      ((a<b) ? a : b)
/*! 求最小值 */
#  define MAX(a, b)      ((a>b) ? a : b)
#endif

/**
* @brief 定义ABS，求绝对值
*/
#ifndef ABS
#  define ABS(x)        ((x<0) ? (-1*(x)) : x)
#endif

/**
* @brief 定义PI=3.141592653...以及度和弧度转换
*/
#ifndef M_PI
/*! 定义圆周率PI */
# define M_PI  3.1415926535897932384626433832795
/*! 弧度转度 */
# define DEG_PER_RAD      ((double)(180.0/M_PI))
/*! 度转弧度 */
# define RAD_PER_DEG      ((double)(M_PI/180.0))
#endif

/**
* @brief 定义平方
*/
#ifndef M_SQUARE
# define M_SQUARE(x)  (x)*(x)
#endif

/**
* @brief 定义立方
*/
#ifndef M_CUBE
# define M_CUBE(x)  (x)*(x)*(x)
#endif

/*! 判断浮点数是否NaN值 */
inline bool isnan(const float& v)  { return _isnan(v) ? true : false; }
/*! 判断double数是否NaN值 */
inline bool isnan(const double& v) { return _isnan(v) ? true : false; }
/*! 获取double的NaN值 */
inline double nan() { return numeric_limits<double>::quiet_NaN(); }

/**
* @brief float类型的极值
*/
#ifndef FLT_EQUALS
/*! 浮点数是否相等 */
#define FLT_EQUALS(x, y)  (fabs((double)x-y)<FLT_EPSILON)
/*! 浮点数是否相等(指定比较阈值) */
#define FLT_EQUALS_N(x, y, z)  (fabs((double)x-y)<z)
#endif

#ifndef FLT_ZERO
/*! 浮点数是否为0 */
#define FLT_ZERO(x)  (fabs(x)<FLT_EPSILON)
#endif

/**
* @brief 释放数组
*/
#define RELEASE(x)	if(x!=NULL) {delete []x; x = NULL;}

/**
* @brief 将全局区域设为操作系统默认区域
*/
#define SET_LOCAL	{ locale::global(locale("")); setlocale(LC_ALL,"Chinese-simplified"); }
/**
* @brief 还原全局区域设定
*/
#define REVERT_LOCAL	locale::global(locale("C"))


#ifndef EQUAL
#if defined(WIN32) || defined(WIN32CE)
/*! 比较字符串是否相等 */
#  define EQUALN(a, b, n)           (_strnicmp(a, b, n) == 0)
/*! 比较字符串是否相等 */
#  define EQUAL(a, b)              (_stricmp(a, b) == 0)
#else
/*! 比较字符串是否相等 */
#  define EQUALN(a, b, n)           (strncasecmp(a, b, n) == 0)
/*! 比较字符串是否相等 */
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

/*! 成功执行 */
const int RE_SUCCESS		= 0;
/*! 文件不存在 */
const int RE_FILENOTEXIST	= 1;
/*! 文件格式不被支持 */
const int RE_FILENOTSUPPORT	= 2;
/*! 图像数据类型不正确 */
const int RE_FILETYPEERROR	= 3;
/*! 创建图像失败 */
const int RE_CREATEFAILED	= 4;
/*! 输入参数错误 */
const int RE_PARAMERROR		= 5;
/*! 其他错误 */
const int RE_FAILED			= 6;
/*! 图像不存在公共区域 */
const int RE_NOSAMEEXTENT	= 7;
/*! 用户取消操作 */
const int RE_USERCANCEL		= 8;
/*! 文件已经被使用 */
const int RE_FILEISUESED	= 9;
/*! 不支持的像素深度 */
const int RE_DEPTHNOTSUPPORT	= 10;
/*! 波段数量不符合要求 */
const int RE_BANDCOUNTERROR		= 11;
/*! 文件不存在投影 */
const int RE_NOPROJECTION		= 12;
/*! 投影不一致 */
const int RE_PROJECTIONDIFF		= 13;

#endif// IMGALG_DEFINE_H