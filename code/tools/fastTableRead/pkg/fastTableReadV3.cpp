/*
 * fastTableRead.cpp
 *
 *  Created on: Aug 5th, 2014
 *      Author: mhofree
 *
 */
#include "mex.h"
#include <iterator>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sys/stat.h>
#include <algorithm>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
//#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/detail/ios.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>


#include <cassert>

// TODO: fix empty delim bug!
#define xMM(A,m,n,D) (A[m+(D*n)])
using namespace std;


uint32_T countLines(const char * fname)
{

	ifstream myfile(fname);
	size_t numLines = 0;

	if (myfile.is_open() && myfile.good())
	{
		// new lines will be skipped unless we stop it from happening:
		myfile.unsetf(std::ios_base::skipws);

		// count the newlines with an algorithm specialized for counting:
		numLines = std::count(
				std::istream_iterator<char>(myfile),
				std::istream_iterator<char>(),'\n');
	}
    return numLines;
}

bool isGzipSuffix(char* fileName) 
{

    std::string::size_type idx;
    idx = string(fileName).rfind('.');
    bool isGzip = false;

    if(idx != std::string::npos)
    {
        std::string extension = string(fileName).substr(idx+1);
        // mexPrintf("Ext:(%s)\n",extension.c_str());
        isGzip = strncmp(extension.c_str(),"gz",2) == 0;
    }

    return isGzip;
}

int readFileData(char* fileName,int &rowLen, char* delimChars, int skipLines, char* commentChar, vector<string> &errorLines,vector<vector<string> > &dataTable)
{
    // Read gzipped files not yet supported
    uint32_T cntLines = 0;    
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inStreamFilter;

    bool isGzip = isGzipSuffix(fileName);

    std::ios_base::openmode fopenMode = std::ios_base::in;
    if (isGzip)
    {
        fopenMode = fopenMode | std::ios_base::binary;
    }

    boost::iostreams::file_source inF(fileName, fopenMode);
    if (!inF.is_open())
    {
        return -2;
    }
    if (isGzip) {
        // mexPrintf("Pre-gzip decomp\n");
        inStreamFilter.push(boost::iostreams::gzip_decompressor());
        // mexPrintf("gzip decomp\n");
    }
    inStreamFilter.push(inF);

    std::istream inStream(&inStreamFilter);
    string line;
    while (std::getline(inStream, line)) {

        vector<string> tokens;
        boost::split(tokens, line, boost::is_any_of(delimChars));

        // Skip header lines
        if (skipLines > 0)
        {
            errorLines.push_back(line);
            --skipLines;
            continue;
        }

        bool foundComment = 0;
        if (commentChar)
        {
            // Simple first in line comment masker
            if (line[0] == commentChar[0])
            {
                errorLines.push_back(line);
                continue;
            }
        }

        if ( rowLen == 0 )
        {
            rowLen = tokens.size();
        };

        if (rowLen == tokens.size())
        {
            dataTable.push_back(tokens);
            ++cntLines;
        }
        else
        {
            errorLines.push_back(line);
        }
    }
    return cntLines;
}

int readFileDataSparse(char* fileName,int &rowLen, mwSize &nzmax, char* delimChars, int skipLines, char* commentChar,int rowHeaderN,int colHeaderN, vector<string> &errorLines,vector<vector<pair<mwIndex,double> > > &dataTable,vector<string> &rowHeaderV,vector<string> &colHeaderV)
{
    // Read gzipped files not yet supported
    uint32_T cntLines = 0;    
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inStreamFilter;
    double nanV = mxGetNaN();

    bool isGzip = isGzipSuffix(fileName);

    std::ios_base::openmode fopenMode = std::ios_base::in;
    if (isGzip)
    {
        fopenMode = fopenMode | std::ios_base::binary;
    }

    boost::iostreams::file_source inF(fileName, fopenMode);
    if (!inF.is_open())
    {
        return -2;
    }
    if (isGzip) {
        // mexPrintf("Pre-gzip decomp\n");
        inStreamFilter.push(boost::iostreams::gzip_decompressor());
        // mexPrintf("gzip decomp\n");
    }
    inStreamFilter.push(inF);

    std::istream inStream(&inStreamFilter);
    string line;
    while (std::getline(inStream, line)) {

        vector<string> tokens;
        boost::split(tokens, line, boost::is_any_of(delimChars));

        // Skip header lines
        if (skipLines > 0)
        {
            errorLines.push_back(line);
            --skipLines;
            continue;
        }
        
        if (commentChar)
        {
            // Simple first in line comment masker
            if (line[0] == commentChar[0])
            {
                errorLines.push_back(line);
                continue;
            }
        }

        if (rowLen > 0 && cntLines == 0 && ((rowLen+1) == tokens.size()))
        {
            // There is no row header column -- this occurs often when writing R data.frames to table.
            ++rowLen;
        }

        if ( rowLen == 0 )
        {
            rowLen = tokens.size();
            colHeaderV = tokens;
        }
        else if (rowLen == tokens.size())
        {

            // Iterate on row
            vector<pair<mwIndex,double> > sparseRowV; 
            vector<string>::iterator dataIt = tokens.begin();
            if (rowHeaderN > 0)
            {
                rowHeaderV.push_back(dataIt->c_str());
                ++dataIt;
            }

            mwIndex spIndex = 0;
            double val = nanV;
            for ( ;dataIt != tokens.end();++dataIt)
            {
                try {
                    val = boost::lexical_cast<double>(*dataIt);
                }
                catch (const boost::bad_lexical_cast &)
                {
                    val = nanV;
                }

                if (val != 0)
                {
                    sparseRowV.push_back({spIndex,val});
                    ++nzmax;
                }
                ++spIndex;
            }

            dataTable.push_back(sparseRowV);
            ++cntLines;
        }
        else
        {
            errorLines.push_back(line);
        }
    }

    if (colHeaderV.size() == rowLen)
    {
        colHeaderV.erase(colHeaderV.begin());
    }

    return cntLines;
}

//int ends_with(const char* name, const char* extension, size_t length)
//{
//    const char* ldot = strrchr(name, '.');
//    if (ldot != NULL)
//    {
//        if (length == 0)
//            length = strlen(extension);
//        return strncmp(ldot + 1, extension, length) == 0;
//    }
//    return 0;
//}


/*******************************************************************************
 **entry point for the mex call
 **nlhs - number of outputs
 **plhs - pointer to array of outputs
 **nrhs - number of inputs
 **prhs - pointer to array of inputs
 *******************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*****************************************************************************
	 * A very simple reader function for reading a text table from file
	 * Input:
	 * flag: (1) Read text output will be a simple cell array
	 *       (2) (TODO) Read data frame
	 *       (3)
	 * filename: name of file to write
	 * colHeaderLines:
	 * rowHeadersLines:
	 * delim:
	 *
	 * Output:
	 * outCellData
	 * outErrorLines
	 *
	 * TODO:
	 * - add simple error checking.
	 *****************************************************************************/

	const mxArray* mxFlag;
	const mxArray* mxFileName;
	const mxArray* mxData = NULL;
	const mxArray* mxDataRow;
	const mxArray* mxDataCol;
	const mxArray* mxTxtRow;
	const mxArray* mxTxtCol;
	const mxArray* mxRowNum = NULL;
	const mxArray* mxDelim = NULL;
    const mxArray* mxCommentChar = NULL;


    int skipLines = 0;
    int rowHeaderN = 0;
    int colHeaderN = 0;
    char* commentChar = "";

	switch (nrhs)
	{
		case 2:
		{
			mxFlag = prhs[0];
			mxFileName = prhs[1];
			/*printf("Must be H(X), calculateProbability(X), merge(X), normaliseArray(X)\n");*/
			break;
		}
		case 3:
		{
			mxFlag = prhs[0];
			mxFileName = prhs[1];
			mxRowNum = prhs[2];
			/*printf("Must be H(X), calculateProbability(X), merge(X), normaliseArray(X)\n");*/
			break;
		}
		case 4:
		{
			mxFlag = prhs[0];
			mxFileName = prhs[1];
			mxRowNum = prhs[2];
			mxDelim = prhs[3];
			break;
		}
		case 5:
		{
			mxFlag = prhs[0];
			mxFileName = prhs[1];
			mxRowNum = prhs[2];
			mxDelim = prhs[3];
            if (prhs[4]) skipLines = *mxGetPr(prhs[4]);
			break;
		}
		case 7:
		{
            mxFlag = prhs[0];
            mxFileName = prhs[1];
            mxRowNum = prhs[2];
            mxDelim = prhs[3];
            if (prhs[4]) skipLines = *mxGetPr(prhs[4]);
            if (prhs[5]) rowHeaderN = *mxGetPr(prhs[5]);
            if (prhs[6]) colHeaderN = *mxGetPr(prhs[6]);
            break;
		}
        case 8:
        {
            mxFlag = prhs[0];
            mxFileName = prhs[1];
            mxRowNum = prhs[2];
            mxDelim = prhs[3];
            if (prhs[4]) skipLines = *mxGetPr(prhs[4]);
            if (prhs[5]) rowHeaderN = *mxGetPr(prhs[5]);
            if (prhs[6]) colHeaderN = *mxGetPr(prhs[6]);
            if (prhs[7]) commentChar = mxArrayToString(prhs[7]);
            break;
        }
		default:
		{
			mexPrintf("Error incorrect number of arguments, format is SimpleByteReadWrite(fcn,data)\n");
			return;
		}
	}

	int flag = *mxGetPr(mxFlag);

	switch (flag)
	{
		case 1: /* read simple text data */
		{

			char* fileName = mxArrayToString(mxFileName);
            int rowLen = 0;
            vector<string> errorLines;
            vector<vector<string> > dataTable;

            char * delimChars;
            if (mxDelim != NULL)
            {
                delimChars = mxArrayToString(mxDelim);
                // mexPrintf("Delim:(%s)\n",delimChars);
            }

            int errorState = readFileData(fileName,rowLen,delimChars,skipLines,commentChar,errorLines,dataTable);
            if ( errorState < 0 )
            {
                mexPrintf("Error no data read\n");
            }
            mxArray* outCellData = mxCreateCellMatrix(dataTable.size(),rowLen);
            mxArray* outErrorLines = mxCreateCellMatrix(errorLines.size(),1);
            mwIndex cellPos[2] = {0,0};

            mwIndex cntCopy = 0;
            for (vector<vector<string> >::iterator dataItRow = dataTable.begin(); dataItRow != dataTable.end();++dataItRow)
            {
                cellPos[1]=0;
                for (vector<string>::iterator dataItCol = dataItRow->begin(); dataItCol != dataItRow->end();++dataItCol)
                {
                    const char * zstr = dataItCol->c_str();
                    cntCopy = mxCalcSingleSubscript(outCellData, rowLen, cellPos);

                    mxSetCell(outCellData, cntCopy, mxCreateString(zstr));
                    ++cellPos[1];
                }
                ++cellPos[0];
            }

            cntCopy = 0;
            for (vector<string>::iterator dataItCol = errorLines.begin(); dataItCol != errorLines.end();++dataItCol)
            {
                const char * zstr = dataItCol->c_str();
                mxSetCell(outErrorLines, cntCopy++, mxCreateString(zstr));
            }


            plhs[0] = outCellData;
            plhs[1] = outErrorLines;
            break;
		}
		case 2: /* read simple numeric only data with top and bottom headers */
		{
            char* fileName = mxArrayToString(mxFileName);
            int dataDim = 0;
            vector<string> errorLines;
            vector<vector<string> > dataTable;

            char * delimChars;
            if (mxDelim != NULL)
            {
                delimChars = mxArrayToString(mxDelim);
                // mexPrintf("Delim:(%s)\n",delimChars);
            }

            int errorState = readFileData(fileName,dataDim,delimChars,skipLines,commentChar,errorLines,dataTable);

            mwIndex nDim = dataDim - rowHeaderN;
            mwIndex nRow = dataTable.size() - colHeaderN;

            // if isSparse
            // {
            
            // }
            // else {
            //     mxArray* outData = mxCreateDoubleMatrix(nRow,nDim,mxREAL);
            // }
            mxArray* outData = mxCreateDoubleMatrix(nRow,nDim,mxREAL);

            mxArray* outErrorLines = mxCreateCellMatrix(errorLines.size(),1);
            mxArray* outColNames;
            mxArray* outRowNames;
            if (colHeaderN>0)
            {
                outColNames = mxCreateCellMatrix(nDim,colHeaderN);
            }
            if (rowHeaderN>0)
            {
                outRowNames = mxCreateCellMatrix(nRow,rowHeaderN);
            }


            double nanV = mxGetNaN();
            double* outDataPr = mxGetPr(outData);

            mwIndex cellPos[2] = {0,0};
            mwIndex dataPos[2] = {0,0};
            mwIndex nColHeader = 0;
            mwIndex nRowHeader = 0;
            mwIndex cntCopy = 0;

            for (vector<vector<string> >::iterator dataItRow = dataTable.begin(); dataItRow != dataTable.end();++dataItRow)
            {
                cellPos[1]=0;
                dataPos[1]=0;

                for (vector<string>::iterator dataItCol = dataItRow->begin(); dataItCol != dataItRow->end();++dataItCol)
                {

                    if (rowHeaderN > 0 && cellPos[1] < rowHeaderN)
                    {
                        // Row header
                        if (cellPos[0] >= colHeaderN ) {
                            const char *zstr = dataItCol->c_str();
                            mxSetCell(outRowNames, nRowHeader, mxCreateString(zstr));
                            assert(nRowHeader < nRow);
                            ++nRowHeader;
                        }
                    }
                    else if (colHeaderN > 0 && cellPos[0] < colHeaderN)
                    {
                        // Col header
                        if (cellPos[1] >= rowHeaderN) {
                            const char *zstr = dataItCol->c_str();
                            mxSetCell(outColNames, nColHeader, mxCreateString(zstr));
                            assert(nColHeader < nDim);
                            ++nColHeader;
                        }
                    }
                    else {
                        cntCopy = mxCalcSingleSubscript(outData, 2, dataPos);

                        try {
                            *(outDataPr + cntCopy) = boost::lexical_cast<double>(*dataItCol);
                            ++dataPos[1];
                        }
                        catch (const boost::bad_lexical_cast &)
                        {
                            *(outDataPr + cntCopy) = nanV;
                            ++dataPos[1];
                        }
                    }
                    ++cellPos[1];
                }

                if (cellPos[0] >= colHeaderN) ++dataPos[0];
                ++cellPos[0];
            }

            cntCopy = 0;
            for (vector<string>::iterator dataItCol = errorLines.begin(); dataItCol != errorLines.end();++dataItCol)
            {
                const char * zstr = dataItCol->c_str();
                mxSetCell(outErrorLines, cntCopy++, mxCreateString(zstr));
            }


            plhs[0] = outData;
            plhs[1] = outErrorLines;
            if (rowHeaderN>0) plhs[2] = outRowNames;
            if (colHeaderN>0) plhs[3] = outColNames;

            break;
		}
		case 3: /* read sparse numeric only data with top and bottom headers */
		{
            char* fileName = mxArrayToString(mxFileName);
            int dataDim = 0;
            mwSize nzmax = 0;

            vector<string> rowHeaderV;
            vector<string> colHeaderV;
            vector<string> errorLines;

            vector<vector<pair<mwIndex,double> > > dataTable;

            char * delimChars;
            if (mxDelim != NULL)
            {
                delimChars = mxArrayToString(mxDelim);
            }

            int errorState = readFileDataSparse(fileName,dataDim,nzmax,delimChars,skipLines,commentChar,rowHeaderN,colHeaderN,errorLines,dataTable,rowHeaderV,colHeaderV);

            mwSize nDim = dataDim-1;
            mwSize nRow = dataTable.size();
    
            mexPrintf("Sparse table:(%d,%d,nnz: %d)\n",nRow,nDim,nzmax);

            // double percent_sparse = 0.2;
            // nzmax=(mwSize)ceil((double)nRow*(double)nDim*percent_sparse);

            // Dimensions are reversed due to how data is organized by rows. Will transpose upstairs
            mxArray* outData = mxCreateSparse(nDim,nRow,nzmax,mxREAL);

            mxArray* outErrorLines = mxCreateCellMatrix(errorLines.size(),1);
            mxArray* outColNames;
            mxArray* outRowNames;

            if (colHeaderN>0)
            {
                outColNames = mxCreateCellMatrix(nDim,colHeaderN);
            }
            if (rowHeaderN>0)
            {
                outRowNames = mxCreateCellMatrix(nRow,rowHeaderN);
            }


            double nan = mxGetNaN();

            // double* outDataPr = mxGetPr(outData);
            // mwIndex cellPos[2] = {0,0};
            // mwIndex dataPos[2] = {0,0};           
            // mwIndex nColHeader = 0;
            // mwIndex nRowHeader = 0;
            // mwIndex cntCopy = 0;            
            // mexPrintf("Transfer\n");

            double* sr = mxGetPr(outData);
            mwIndex* irs = mxGetIr(outData);
            mwIndex* jcs = mxGetJc(outData);

            mwIndex j = 0;
            mwIndex k = 0;
            for (vector<vector<pair<mwIndex,double> > >::iterator dataItRow = dataTable.begin(); dataItRow != dataTable.end();++dataItRow)
            {
                jcs[j] = k;
                for (vector<pair<mwIndex,double> >::iterator dataItCol = dataItRow->begin(); dataItCol != dataItRow->end();++dataItCol)
                {
                    irs[k] = dataItCol->first;
                    sr[k] = dataItCol->second;
                    // mexPrintf("%d %d %f\n",irs[k],j,sr[k]);
                    ++k;
                }
                ++j;
            }
            jcs[j] = k;

            mwIndex cntCopy = 0;
            for (vector<string>::iterator dataItCol = errorLines.begin(); dataItCol != errorLines.end();++dataItCol)
            {
                const char * zstr = dataItCol->c_str();
                mxSetCell(outErrorLines, cntCopy++, mxCreateString(zstr));
            }
            
            cntCopy = 0;
            for (vector<string>::iterator dataItCol = rowHeaderV.begin(); dataItCol != rowHeaderV.end();++dataItCol)
            {
                const char * zstr = dataItCol->c_str();
                mxSetCell(outRowNames, cntCopy++, mxCreateString(zstr));
            }

            cntCopy = 0;
            // vector<string>::iterator dataItCol = colHeaderV.begin();
            // if (dataItCol != colHeaderV.end()) ++dataItCol;
            for (vector<string>::iterator dataItCol = colHeaderV.begin(); dataItCol != colHeaderV.end();++dataItCol)
            {            
                const char * zstr = dataItCol->c_str();
                mxSetCell(outColNames, cntCopy++, mxCreateString(zstr));
            }


            plhs[0] = outData;
            plhs[1] = outErrorLines;
            if (rowHeaderN>0) plhs[2] = outRowNames;
            if (colHeaderN>0) plhs[3] = outColNames;

            break;
		}
		default:
		{
			printf("Unrecognised flag %d\n",flag);
			break;
		}/*default*/
	}

	return;
}/*mexFunction()*/




