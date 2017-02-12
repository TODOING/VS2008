#ifndef _CUL_LOG_H_
#define _CUL_LOG_H_
#include <stdio.h>
#include <errno.h>
/**
* 说明:
* snprintf()函数的格式跟printf差不多一样，是在c里面用的函数，包含在 #include <stdio.h>头文件中。
* 但snprintf()函数并不是标准c/c++中规定的函数，所以在许多编译器中，厂商提供了其相应的实现的版本。
* 在gcc中，该函数名称就是snprintf()，而在VS中称为_snprintf。 所以在需要使用snprintf()时改成_snprintf就可以了，或则在预编译处加入：
*/
#if _MSC_VER
#define snprintf _snprintf
#endif

// 调试时使用控制台输出，Release时使用文件输出
#ifdef _DEBUG
#define USE_CONSOLE_LOG 1
#else
#define USE_CONSOLE_LOG 0
#endif

#define  USE_LOG4CPLUS
#ifdef USE_LOG4CPLUS
#include <log4cplus/configurator.h>
#include <string>

static log4cplus::Logger logger= log4cplus::Logger::getInstance("Log");
static void init_log(const std::string & path)
{
	log4cplus::PropertyConfigurator::doConfigure(path);   
}

// CL -- CL的缩写
#define CL_ALL        log4cplus::TRACE_LOG_LEVEL
#define CL_TRACE        log4cplus::TRACE_LOG_LEVEL
#define CL_DEBUG        log4cplus::DEBUG_LOG_LEVEL
#define CL_INFO            log4cplus::INFO_LOG_LEVEL
#define CL_NOTICE        log4cplus::INFO_LOG_LEVEL
#define CL_WARNING    log4cplus::WARN_LOG_LEVEL
#define CL_ERROR        log4cplus::ERROR_LOG_LEVEL
#define CL_CRITICAL    log4cplus::FATAL_LOG_LEVEL

#define CULLog(l, ...) \
	do { \
		if(logger.isEnabledFor(l)) { \
			char __buf__internal__[2046]={0}; \
			snprintf(__buf__internal__,2045,__FUNCTION__"; -"##__VA_ARGS__); \
			logger.forcedLog(l, __buf__internal__, __FILE__, __LINE__); \
		} \
	} while(0);
#elif define USE_ACE_LOG
#include "ace/Log_Msg.h"
#include "ace/Trace.h"
#define CL_ALL        LM_TRACE
#define CL_TRACE        LM_TRACE
#define CL_DEBUG        LM_DEBUG
#define CL_INFO             LM_INFO
#define CL_NOTICE        LM_NOTICE
#define CL_WARNING    LM_WARNING
#define CL_ERROR        LM_ERROR
#define CL_CRITICAL        LM_CRITICAL
#define CULLog(l,) do{ \
	ACE_DEBUG((l,"[%T|%t] %s-%s:%

d",__FILE__,__FUNCTION__,__LINE__)); \
 ACE_DEBUG((l,__VA_ARGS__)); \
 }while(0)
#else
#include <pthread.h>
#include <time.h>
#include <sys/time.h>
#define CL_LEVEL 0xFF
typedef enum{
	CL_CRITICAL=0x1,
	CL_ERROR=0x2,
	CL_WARNING=0x4,
	CL_NOTICE=0x8,
	CL_INFO=0x10,
	CL_DEBUG=0x20,
	CL_TRACE=0x40,
	CL_ALL=0x80
}CLLevel;
#define CULLog(l, ...) do{ \
	if(CL_LEVEL&l){ \
struct timeval now;\
	gettimeofday(&now,0); \
struct tm *ptm=localtime(&(now.tv_sec)); \
	printf("[%d|%d:%d:%d.%d] [%s/%s/%d] 

",pthread_self(),ptm->tm_hour,ptm->tm_min,ptm-

>tm_sec,now.tv_usec,__FILE__,__FUNCTION__,__LINE__); \
printf( __VA_ARGS__); \
	} \
	}while(0)
#endif // #ifdef USE_LOG4CPLUS 

#define die(str) { CULLog(CL_WARNING,str); return; }
#define die_0(str) { CULLog(CL_WARNING,str); return 0; }
#define die_1(str) { CULLog(CL_WARNING,str); return -1; }
#define die_ns(str) { CULLog(CL_WARNING,str); return ""; }

/*safe func return empty,0,-1*/
#define SAFE_FUNC(func) if((func)<0) \
	{ \
		CULLog(CL_ERROR,"[%s:%d]error!error[%d:%s]   \n",\
			__FILE__,__LINE__,errno,strerror(errno)); \
		exit(-1); \
	}

/*safe func but 1 err return empty,0,-1*/
#define SAFE_FUNC_BUT_ERR1(func,ERR1) do \
{ \
	if((func)<0){ \
		CULLog(CL_ERROR,"[%s:%d]error!error[%d:%s]\n",\
			__FILE__,__LINE__,errno,strerror(errno)); \
		if(errno!=ERR1) exit(-1); \
	} \
	else break; \
}while(1)

/*safe func but 2 err return empty,0,-1*/
#define SAFE_FUNC_BUT_ERR2(func,ERR1,ERR2) do \
{ \
	if((func)<0){ \
		CULLog(CL_ERROR,"[%s:%d]error!error[%d:%s] \n",\
			__FILE__,__LINE__,errno,strerror(errno)); \
		if(errno!=ERR1&&errno!=ERR2)  exit(-1); \
	} \
	else break; \
}while(1)

#endif // #ifndef _CUL_LOG_H_