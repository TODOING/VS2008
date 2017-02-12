#ifndef STATISTICS_IO_H
#define STATISTICS_IO_H
#include <stdio.h>

class PCAStatisticsIO
{
public:
	PCAStatisticsIO(const char *file, int band_count = 0);
	~PCAStatisticsIO();
	
	bool WriteInit();
	int WriteMean(double *mean);
	int WriteCovarianceOrCorrelation(double *convariance_or_correlation);
	int WriteEigenvectors(double *eigenvectors);
	int WriteEigenvalue(double *eigenvalue);
	int WriteAccumulateContribute(double *accumulate_contribute);
	int WriteToFile();

	bool ReadInit();
	int ReadBandCount();
	int ReadMean(double *mean);
	int ReadCovarianceOrCorrelation(double *convariance_or_correlation);
	int ReadEigenvectors(double *eigenvectors);
	int ReadEigenvalue(double *eigenvalue);
	int ReadAccumulateContribute(double *accumulate_contribute);

private:
	FILE *m_file;
	double *m_buffer;
	const char *m_statistics_file;
	int m_band_count;
};

#endif// STATISTICS_IO_H