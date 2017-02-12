#ifndef _UTILITYLOG_H_
#define _UTILITYLOG_H_
#include <stdio.h>
#include <errno.h>

/** 
* @brief 导出符号定义 
*/  
#ifdef UTILITYLOG_EXPORTS 
#define UTILITYLOG_API __declspec(dllexport)  
#else  
#define UTILITYLOG_API __declspec(dllimport)  
#endif

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
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#include

#endif // _UTILITYLOG_H_