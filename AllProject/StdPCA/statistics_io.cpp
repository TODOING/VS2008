#include "statistics_io.h"
#include <string.h>

PCAStatisticsIO::PCAStatisticsIO(const char *file, int band_count/* = 0*/)
{
	m_file = NULL;
	m_buffer = NULL;
	m_statistics_file = file;
	m_band_count = band_count;
}

PCAStatisticsIO::~PCAStatisticsIO()
{
	if (m_buffer)
	{
		delete []m_buffer;
		m_buffer = NULL;
	}

	if (m_file)
		fclose(m_file);
}

bool PCAStatisticsIO::WriteInit()
{
	m_file = fopen(m_statistics_file, "wb");
	if (m_file == NULL)
		return false;
	
	int total_count = 1 + m_band_count*(3 + 2*m_band_count); 
	m_buffer = new double[total_count];
	
	m_buffer[0] = m_band_count;

	return true;
}

int PCAStatisticsIO::WriteMean(double *mean)
{
	memcpy(m_buffer + 1, mean, m_band_count*sizeof(double));
	return 0;
}

int PCAStatisticsIO::WriteCovarianceOrCorrelation(double *convariance_or_correlation)
{
	memcpy(m_buffer + 1 + m_band_count, convariance_or_correlation, m_band_count*m_band_count*sizeof(double));
	return 0;
}

int PCAStatisticsIO::WriteEigenvectors(double *eigenvectors)
{
	memcpy(m_buffer + 1 + m_band_count + m_band_count*m_band_count, eigenvectors, m_band_count*m_band_count*sizeof(double));
	return 0;
}

int PCAStatisticsIO:: WriteEigenvalue(double *eigenvalue)
{
	memcpy(m_buffer + 1 + m_band_count + 2*m_band_count*m_band_count, eigenvalue, m_band_count*sizeof(double));
	return 0;
}

int PCAStatisticsIO::WriteAccumulateContribute(double *accumulate_contribute)
{
	memcpy(m_buffer + 1 + 2*m_band_count + 2*m_band_count*m_band_count, accumulate_contribute, m_band_count*sizeof(double));
	return 0;
}

int PCAStatisticsIO::WriteToFile()
{
	if (m_file && m_buffer)
		return fwrite(m_buffer, sizeof(double), 1 + m_band_count*(3 + 2*m_band_count), m_file);
}

bool PCAStatisticsIO::ReadInit()
{
	m_file = fopen(m_statistics_file, "rb");
	if (m_file == NULL)
		return false;

	double band_count;
	fread(&band_count, sizeof(double), 1, m_file);
	m_band_count = band_count;
	return true;
}

int PCAStatisticsIO::ReadMean(double *mean)
{
	if (!m_file)
		return -1;

	fseek(m_file, sizeof(double)*1, SEEK_SET);
	fread(mean, sizeof(double), m_band_count, m_file);

	return 0;
}

int PCAStatisticsIO::ReadCovarianceOrCorrelation(double *convariance_or_correlation)
{
	if (!m_file)
		return -1;

	fseek(m_file, sizeof(double)*(1 + m_band_count), SEEK_SET);
	fread(convariance_or_correlation, sizeof(double), m_band_count*m_band_count, m_file);

	return 0;
}

int PCAStatisticsIO::ReadEigenvectors(double *eigenvectors)
{
	if (!m_file)
		return -1;

	fseek(m_file, sizeof(double)*(1 + m_band_count + m_band_count*m_band_count), SEEK_SET);
	fread(eigenvectors, sizeof(double), m_band_count*m_band_count, m_file);

	return 0;
}

int PCAStatisticsIO::ReadEigenvalue(double *eigenvalue)
{
	if (!m_file)
		return -1;

	fseek(m_file, sizeof(double)*(1 + m_band_count + 2*m_band_count*m_band_count), SEEK_SET);
	fread(eigenvalue, sizeof(double), m_band_count, m_file);

	return 0;
}

int PCAStatisticsIO::ReadAccumulateContribute(double *accumulate_contribute)
{
	if (!m_file)
		return -1;

	fseek(m_file, sizeof(double)*(1 + 2*m_band_count + 2*m_band_count*m_band_count), SEEK_SET);
	fread(accumulate_contribute, sizeof(double), m_band_count, m_file);

	return 0;
}

int PCAStatisticsIO::ReadBandCount()
{
	return m_band_count;
}

