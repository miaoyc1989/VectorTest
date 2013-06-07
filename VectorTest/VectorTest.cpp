//VectorTest.cpp : 定义控制台应用程序的入口点。
 
#include "stdafx.h"
#include "DspFilters/Dsp.h"

#include "iostream"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex> 
#include "vector.h"
#include "matrix.h"
#include <bitset>

#include "time.h"

#include "windows.h"
#include "cstring"

#include "StdString.h"
#include <map>


#define BOUNDS_CHECK
#define IDS_STRING1                     101
#define IDS_STRING2                     102


using namespace std;
using namespace splab;


const int Vector_Size = 1024000;
const int Matrix_Row = 102400;
const int Matrix_Col = 50;


#define WRITE_LOG(str) WriteLog("Log.txt", str);     

void WriteLog(char* filename,char* str)   
{   
	time_t t;   
	time(&t);   
	struct tm* tp= localtime(&t);   

	char now_str[100];   
	strftime(now_str, 100, "%Y-%m-%d %H:%M:%S", tp);   

	FILE *fo;   
	fo = fopen(filename, "a");   
	if (fo == 0) {     
		return;   
	}   

	fprintf(fo, "%s %s\r",now_str, str);   
	fclose(fo);   
}  

void Display(const double *p, int length )
{
	for( int i=0; i<length; ++i )
	{
		cout << p[i] << "\t" ;
	}
	cout << endl;
}

void TestVector()
{
	int t_Start, t_Stop, t_result; 
	int t_init, t_getsize, t_setsize, t_scale, t_set, t_resize,t_append, t_minmax, t_abs, t_sum, t_subset, t_equal;

	//
	// 向量初始化
	//
	t_Start = GetTickCount();
	Vector<double> v1(Vector_Size);				
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "init cost time:" 
		 << t_result 
		 << endl;
	cout << "element numbers :" 
		 << v1.GetSize() 
		 << endl << endl;
	t_init = t_result;

	//
	//Getsize
	//
	t_Start = GetTickCount();
	v1.GetSize();						
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "getsize cost time:" 
		 << t_result 
		 << endl;
	cout << "element numbers :" 
		 << v1.GetSize() 
		 << endl << endl;
	t_getsize = t_result;

	//
	//Setsize
	//
	t_Start = GetTickCount();
	v1.SetSize(Vector_Size + 1);						
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "Setsize cost time:" 
		 << t_result 
		 << endl;
	cout << "element numbers :" 
		 << v1.GetSize() 
		 << endl << endl;
	t_setsize = t_result;

	//
	//Scale
	//
	t_Start = GetTickCount();
	v1.Scale(30,20);				
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "scale cost time:" 
		 << t_result 
		 << endl;
	cout << "element numbers :" 
		 << v1.GetSize() 
		 << endl << endl;
	t_scale = t_result;

	//
	//Set
	//
	t_Start = GetTickCount();
	v1.Set(2);						
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "Set cost time:" 
		 << t_result 
		 << endl;
	cout << "element numbers :" 
		 << v1.GetSize() 
		 << endl << endl;
	t_set = t_result;

	//
	//Resize
	//
	t_Start = GetTickCount();
	v1.Resize(Vector_Size - 2);						
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "Resize cost time:"
		 << t_result 
		 << endl;
	cout << "element numbers :" 
		 << v1.GetSize() 
		 << endl << endl;
	t_resize = t_result;

	//
	//Append
	//

	Vector<double> v2(Vector_Size);					// 向量初始化

	t_Start = GetTickCount();
	v1.Append(v2);					
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "append cost time:" 
		 << t_result 
		 << endl;
	cout << "element numbers :" 
		 << v1.GetSize() 
		 << endl << endl;
	t_append = t_result;

	//
	//MinMax
	//
	v1.Resize(Vector_Size);
	double min = 0;
	unsigned int minIndex = 0;
	double max = 0;
	unsigned int maxIndex = 0;
	t_Start = GetTickCount();
	v1.MinMax (min, minIndex, max, maxIndex);
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "MinMax cost time:" 
		 << t_result 
		 << endl;
	cout << "element numbers :" 
		 << v1.GetSize() 
		 << endl << endl;
	t_minmax = t_result;

	//
	//Abs
	//
	v1.Set(-3);
	t_Start = GetTickCount();
	v1.Abs();					
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "Abs cost time:" 
		 << t_result 
		 << endl;
	cout << "element numbers :" 
		 << v1.GetSize() 
		 << endl << endl;
	t_abs = t_result;

	//
	//Sum
	//
	t_Start = GetTickCount();
	v1.Sum();					
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "sum cost time:" 
	  	 << t_result 
		 << endl;
	cout << "element numbers :" 
		 << v1.GetSize() 
		 << endl << endl;
	t_sum = t_result;


	//
	//SubSet
	//

	v1.Resize(30);
	Vector<double> v3(10);
	t_Start = GetTickCount();

	v3 = v1.Subset(0,10);
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "subset cost time:" 
		 << t_result 
		 << endl;
	cout << "element numbers :" 
		 << v1.GetSize() 
		 << endl << endl;
	t_subset = t_result;

	//
	// "=" 
	//
 	v3.Resize(Vector_Size);
 	Vector<double> v4(10);
	t_Start = GetTickCount();
	v4 = v3;
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout <<  " = cost time:" 
		<< t_result 
		<< endl;
	cout << "element numbers :" 
		<< v1.GetSize() 
		<< endl << endl;
	t_equal = t_result;

	//
	//打印日志
	//


	char s_dest[100] = ",";
	char s_source[5];

	//t_init
	itoa(t_init, s_source, 10);
	strcat(s_dest, s_source);

	//t_getsize
	itoa(t_getsize, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_setsize
	itoa(t_setsize, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_scale
	itoa(t_scale, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_set
	itoa(t_set, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_resize
	itoa(t_resize, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_append
	itoa(t_append, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source); 

	//t_minmax
	itoa(t_minmax, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_abs
	itoa(t_abs, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_sum
	itoa(t_sum, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_subset
	itoa(t_subset, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_equal
	itoa(t_equal, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	WRITE_LOG(s_dest);
}

void TestMatrix()
{
	int t_Start, t_Stop, t_result;   
	int t_init, t_setsize, t_resize, t_scale, t_getrows, t_getcols, t_copyrow, t_copycolumn, t_assignrow, t_minmax, t_equal;

	//
	// Matrix初始化
	//
	t_Start = GetTickCount();
	Matrix<double> m1(Matrix_Row, Matrix_Col);
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "init cost time:" 
		 << t_result 
    	 << endl;
	cout << "Matrix size :"  << m1.GetRows() << " * "
							 << m1.GetCols() << endl << endl;
	t_init = t_result;

	//
	// setsize
	//
	t_Start = GetTickCount();
	m1.SetSize(Matrix_Row - 1, Matrix_Col - 1);
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "setsize cost time:" 
		 << t_result 
		 << endl;
	cout << "Matrix size :"  << m1.GetRows() << " * "
							 << m1.GetCols() << endl << endl;
	t_setsize = t_result;

	//
	// Resize
	//
	t_Start = GetTickCount();
	m1.Resize(Matrix_Row + 2, Matrix_Col + 2);
	t_Stop = GetTickCount(); 
	t_result = t_Stop - t_Start;
	cout << "resize cost time:" 
		 << t_result 
		 << endl;
	cout << "Matrix size :"  << m1.GetRows() << " * "
		                     << m1.GetCols() << endl << endl;
	t_resize = t_result;

	//
	// Scale
	//
	t_Start = GetTickCount();
	m1.Scale(3.14,1);
	t_Stop = GetTickCount(); 
	t_result = t_Stop - t_Start;
	cout << "scale cost time:" 
		 << t_result 
		 << endl;
	cout << "Matrix size :"  << m1.GetRows() << " * "
		                     << m1.GetCols() << endl << endl;
	t_scale = t_result;

	//
	// getRows
	//
	t_Start = GetTickCount();
	m1.GetRows();
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "getRows cost time:" 
		 << t_result 
		 << endl;
	cout << "Matrix size :"  << m1.GetRows() << " * "
							 << m1.GetCols() << endl << endl;
	t_getrows = t_result;

	//
	// getCols
	//
	t_Start = GetTickCount();
	m1.GetCols();
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "getCols cost time:" 
	 	 << t_result 
		 << endl;
	cout << "Matrix size :"  << m1.GetRows() << " * "
							 << m1.GetCols() << endl << endl;
	t_getcols = t_result;

	//
	// copy row
	//
	Vector<double> v1(Matrix_Row);
	t_Start = GetTickCount();
	m1.CopyRow(1, v1);
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "copyRow cost time:" 
		 << t_result 
		 << endl;
	cout << "Matrix size :"  << m1.GetRows() << " * "
		                     << m1.GetCols() << endl << endl;
	t_copyrow = t_result;
	
	//
	// copy column
	//
	t_Start = GetTickCount();
	m1.CopyColumn(1, v1);
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "copyColumn cost time:" 
		 << t_result 
		 << endl;
	cout << "Matrix size :"  << m1.GetRows() << " * "
		                     << m1.GetCols() << endl << endl;
	t_copycolumn = t_result;

	//
	// assign row
	//
	Vector<double> v2(Matrix_Row, 123);
	t_Start = GetTickCount();
	m1.AssignRow(1, v2);
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "assign row cost time:" 
		 << t_result 
		 << endl;
	cout << "Matrix size :"  << m1.GetRows() << " * "
		                     << m1.GetCols() << endl << endl;
	t_assignrow = t_result;

	//
	// minmax
	//
	double min,max;
	unsigned int minRowIndex, minColIndex, maxRowIndex, maxColIndex;
	t_Start = GetTickCount();
	m1.MinMax(min, minRowIndex, minColIndex, max, maxRowIndex, maxColIndex);
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << "minmax row cost time:" 
		 << t_result 
		 << endl;
	cout << "Matrix size :"  << m1.GetRows() << " * "
		                     << m1.GetCols() << endl << endl;
	t_minmax = t_result;


	//
	// “=”
	//
	Matrix<double> m2(Matrix_Col + 1, Matrix_Row + 1);
	t_Start = GetTickCount();
	m2 = m1;
	t_Stop = GetTickCount();
	t_result = t_Stop - t_Start;
	cout << " =  cost time:" 
		 << t_result 
		 << endl;
	cout << "Matrix size :"  << m1.GetRows() << " * "
		                     << m1.GetCols() << endl << endl;
	t_equal = t_result;


	//
	//记录日志
	//

	char s_dest[100] = ",";
	char s_source[5];

	//t_init
	itoa(t_init, s_source, 10);
	strcat(s_dest, s_source);

	//t_setsize
	itoa(t_setsize, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_getsize
	itoa(t_resize, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_scale
	itoa(t_scale, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_getrows
	itoa(t_getrows, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_getcols
	itoa(t_getcols, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_copyrow
	itoa(t_copyrow, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_copycolumn
	itoa(t_copycolumn, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_assignrow
	itoa(t_assignrow, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_minmax
	itoa(t_minmax, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	//t_equal
	itoa(t_equal, s_source, 10);
	strcat(s_dest,",");
	strcat(s_dest, s_source);

	WRITE_LOG(s_dest);
}



	/// 重采样
	/** 如果期望采样频率小于实际采样频率，则可以进行重采样
	*	\param[in/out] nivfData_ 波形数据
	*   \param[in] nSampRate_ 实际采样率
	*	\param[in] nExpectedSampRate_ 期望采样率
	*	\param[out] nActualSampRate_ 重采样后得到的采样率
	*   \param[in] nExpectedSampNum_ 期望采样点数：<=0表示最大限度的重采样
	*	\param[in] nOrder_ 抗混滤波器的阶次
	*	\return 无
	*/
template <class Type>
bool Resample( Vector<Type> &nivfData_, int nSampRate_, int nExpectedSampRate_, int &nActualSampRate_, int nExpectedSampNum_ /*= 0*/ )
{
	// 输入校验
	if (nivfData_.GetSize() < 1)
	{
		return false; 
	}
	if (nSampRate_ < 1)
	{
		return false;
	}
	if (nExpectedSampRate_ < 1 || nExpectedSampRate_ > nSampRate_)
	{
		return false;
	}
	int nInterval = nSampRate_ / nExpectedSampRate_;          // 抽样间隔
	nActualSampRate_ = nSampRate_ / nInterval;                   // 实际将得到的采样率
	float fCutoffFreq = nActualSampRate_ / 2.56;                   // 截止频率
	float fAvg = nivfData_.Sum() / nivfData_.GetSize();            // 平均值
	int nActualSampNum = nivfData_.GetSize() / nInterval;     // 实际得到的采样点数
	if (nivfData_.GetSize() % nInterval != 0)
	{
		++nActualSampNum;
	}
	if (nExpectedSampRate_ < nSampRate_)
	{
		Vector<Type> nivdTmp = nivfData_;
		nivdTmp.Scale(1, -fAvg);   // 去直流
		// 抗混滤波
		EllipticHighPass(nivdTmp, nSampRate_, fCutoffFreq);
		
		int j = 0;
		for (int i = 0; i < nivdTmp.GetSize(); i += nInterval, ++j)
		{
			nivfData_[j] = nivdTmp[i];
		}
		nivfData_.Resize(nActualSampNum);
		nivfData_.Scale(1, fAvg);    // 加直流
	}
	if (nExpectedSampNum_>0 && nExpectedSampNum_<nivfData_.GetSize())
	{
		nivfData_.Resize(nExpectedSampNum_);
	}
	return true;
}

template <class Type>
void EllipticHighPass (Vector<Type>& nivdTmp, int nSampRate, float fCutoffFreq)
{
	int numSamples = nivdTmp.GetSize();
	float* data[1];
	data[0] = new float[numSamples];

	for (int i = 0; i < numSamples; i++)
 	{
 		data[0][i] = nivdTmp[i];
 	}

	Dsp::Filter* f = new Dsp::SmoothedFilterDesign <Dsp::RBJ::Design::LowPass, 1> (numSamples);
	Dsp::Params params;
	params[0] = 44100;	// sample rate
	params[1] = 400;	// cutoff frequency
	params[2] = 1.25;   // Q
	f->setParams (params);
	f->process (numSamples, data);

}

void TestResample()
{
	Vector<float> nivfData(Vector_Size);
	int nSampRate = 102400;
	int nExpectedSampRate = 2560;
	int nActualSampRate = 0;
	int nExpectedSampNum = 0;

//	Resample(nivfData, nSampRate, nExpectedSampRate, nActualSampRate, nExpectedSampNum);
}

unsigned char LRC(unsigned char *auchMsg, unsigned short usDataLen)
{
	unsigned char uchLRC = 0 ; /* LRC char initialized */
	while (usDataLen--) /* pass through message buffer */
	{
		uchLRC += *auchMsg++ ; /* add buffer byte without carry */
	}
	return ((unsigned char)(-((char)uchLRC))) ; /* return twos complement */
}

int _tmain(int argc, _TCHAR* argv[])
{
	// You can add just about any form of string to a CStdString
	// with operator+()
	CStdString strVal1(_T("THIS IS A STRING   "));
	OutputDebugString(strVal1 + _T("\n"));
	strVal1 += _bstr_t(" plus a BSTR string");
	strVal1 += '.' ;

	// Some conversion functions can be chained together

	strVal1.ToLower().TrimRight();

	// Case INsensitive comparison via Equals() or CompareNoCase()

	strVal1 = _T("THIS IS A STRING");
	CStdString strVal2(_T("thIs Is a sTRing"));
	_ASSERTE(strVal1 != strVal2);
	_ASSERTE(strVal1.Equals(strVal2));
	_ASSERTE(strVal1.CompareNoCase(strVal2) == 0);

	// Format() works just like CString's

	strVal1.Format(_T("This %s a string named strVal%d"), _T("IS"), 1);
	OutputDebugString(strVal1 + _T("\n"));

	// Declare an STL map class which maps strings to integers.  The
	// keys are case insensitive, so an integer stored under the key
	// _T("MYKEY") could be retrieved with the value _T("mykey")

	typedef std::map<CStdString, int, StdStringLessNoCase> CMyMap;
	CMyMap myMap;
	myMap[_T("MYKEY")] = 7;
	_ASSERTE(myMap.find(_T("mykey")) != myMap.end());


	// If we were calling some windows function that fills out a
	// buffer for us we can use the GetBuffer() function.  .

	CStdString strPath;
	::GetTempPath( MAX_PATH, strPath.GetBuffer(MAX_PATH+1));
	strPath.ReleaseBuffer();

	// You can set the resource handle for loading string resources
	// and then load them via either the constructor or the Load()
	// function.

	CStdString::SetResourceHandle(::GetModuleHandle(NULL));
	CStdString strString(MAKEINTRESOURCE(IDS_STRING1));
	_ASSERTE(_T("All your base are belong to us!") == strString);

	strString.Load(IDS_STRING2);
	_ASSERTE(_T("I see dead people") == strString);

	return 0;
} 
