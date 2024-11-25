#pragma once
#include <DECODE.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <iostream>
#include <numeric>

//#define SCORE_MAX 140 //
#define PEPTIDE_LENGTH 12
__constant__ float g_AAMap[21 * 21];

int SCORE_MAX = 10;
char AAMap[21 * 21];

size_t task = 0;
size_t task_target = 0;
size_t samplesize = 0;

unsigned int Cores;
string scoretableFP = {};
//#define FloatSHIFT 100000.0
#define FloatSHIFT 1000.0
#define FloatSHIFT_INT 1000.0
#define FloatSHIFT_Ratio 1000.0

#define CUDA_SAFE_CALL(func) \
do { \
     cudaError_t err = (func); \
     if (err != cudaSuccess) { \
         fprintf(stderr, "[Error] %s (error code: %d) at %s line %d\n", cudaGetErrorString(err), err, __FILE__, __LINE__); \
         exit(err); \
     } \
} while(0)

double P_TH[12] = { 1,0.05,0.01,0.005,0.001,0.0005,0.0001, 0.00005, 0.00001, 0.000005, 0.000001, 0.0000001 };


struct Option {
	string FPr = "./";
	string FPw = "./";
	string P_map = "./";
	string AAtable;
	string Target;
	string Start = "0";
	int threads = 1;
	int PepLimit = 0;
	int Q_Export = 0;
	int P_Export = 0;
	string Plistexport;
	int Calc = 0;               // 0: read of peptide, 1;kind of peptide
	int Distance_function = 0;  // 0:KL, 1:L1, 2:L2, 3:
	double TH = 0.05*FloatSHIFT;
} OPs;

struct PV {
	// vector<vector<float>> scorelist;
	string *filename;
	float *scorelist;  // P value
	char *Seq;         // Protein seq
	int *Seqlength;
	int scoresize = 0;
	int targetsize = 0;
	size_t totalAA = 0;
	size_t totalread = 0;
};

struct Pep_data {
	char *pep;
	float *read;
	int totalPeps = 0;  // kinds of pep
	int peplength = 0;  // kinds of pep x peplength
};

struct mems {
	//host
	size_t GPU_mem;
	size_t Host_mem;
	float *h_distance;
	float *h_AAratio;
	float *h_Scoreratio;
	float *h_read_TH;
	float *h_result;
	double *h_Read_Ratio;

	//device
	char *d_pep, *d_pro, *d_AAmap;
	float *d_AAratio;
	float *d_P_list;
	float *d_Scoreratio;
	float *d_distance;
	float *d_read;
	float *d_read_TH;
};


vector<vector<vector<string>>> queue;
vector<bool> queflag;
inline char S2C(char src) {
	switch (src) {
	case 'A':
		return 0;
	case 'R':
		return 1;
	case 'N':
		return 2;
	case 'D':
		return 3;
		break;
	case 'C':
		return 4;
		break;
	case 'Q':
		return 5;
		break;
	case 'E':
		return 6;
		break;
	case 'G':
		return 7;
		break;
	case 'H':
		return 8;
		break;
	case 'I':
		return 9;
		break;
	case 'L':
		return 10;
		break;
	case 'K':
		return 11;
		break;
	case 'M':
		return 12;
		break;
	case 'F':
		return 13;
		break;
	case 'P':
		return 14;
		break;
	case 'S':
		return 15;
		break;
	case 'T':
		return 16;
		break;
	case 'W':
		return 17;
		break;
	case 'Y':
		return 18;
		break;
	case 'V':
		return 19;
		break;
	case '*':
		return 20;
		break;
	default:
		return 20;
		break;
	}
	return 20;
}
inline char *AA2Chr2(string &src) {
	int length = src.length();
	char *dst = new char[length + 1];
	for (auto i = 0; i < length; ++i) {
		switch (src[i]) {
		case 'A':
			dst[i] = 0;
			break;
		case 'R':
			dst[i] = 1;
			break;
		case 'N':
			dst[i] = 2;
			break;
		case 'D':
			dst[i] = 3;
			break;
		case 'C':
			dst[i] = 4;
			break;
		case 'Q':
			dst[i] = 5;
			break;
		case 'E':
			dst[i] = 6;
			break;
		case 'G':
			dst[i] = 7;
			break;
		case 'H':
			dst[i] = 8;
			break;
		case 'I':
			dst[i] = 9;
			break;
		case 'L':
			dst[i] = 10;
			break;
		case 'K':
			dst[i] = 11;
			break;
		case 'M':
			dst[i] = 12;
			break;
		case 'F':
			dst[i] = 13;
			break;
		case 'P':
			dst[i] = 14;
			break;
		case 'S':
			dst[i] = 15;
			break;
		case 'T':
			dst[i] = 16;
			break;
		case 'W':
			dst[i] = 17;
			break;
		case 'Y':
			dst[i] = 18;
			break;
		case 'V':
			dst[i] = 19;
			break;
		case '*':
			dst[i] = 20;
			break;
		default:
			dst[i] = 20;
			break;
		}
	}
	dst[length] = '\0';
	return dst;
}
inline char *AA2Chr_pep2(char src[PEPTIDE_LENGTH + 1]) {
	char *dst = new char[PEPTIDE_LENGTH + 1];
	for (int i = 0; i < PEPTIDE_LENGTH + 1; ++i) {
		dst[i] = '\0';
	}
	for (auto i = 0; i < PEPTIDE_LENGTH; ++i) {
		switch (src[i]) {
		case 'A':
			dst[i] = 0;
			break;
		case 'R':
			dst[i] = 1;
			break;
		case 'N':
			dst[i] = 2;
			break;
		case 'D':
			dst[i] = 3;
			break;
		case 'C':
			dst[i] = 4;
			break;
		case 'Q':
			dst[i] = 5;
			break;
		case 'E':
			dst[i] = 6;
			break;
		case 'G':
			dst[i] = 7;
			break;
		case 'H':
			dst[i] = 8;
			break;
		case 'I':
			dst[i] = 9;
			break;
		case 'L':
			dst[i] = 10;
			break;
		case 'K':
			dst[i] = 11;
			break;
		case 'M':
			dst[i] = 12;
			break;
		case 'F':
			dst[i] = 13;
			break;
		case 'P':
			dst[i] = 14;
			break;
		case 'S':
			dst[i] = 15;
			break;
		case 'T':
			dst[i] = 16;
			break;
		case 'W':
			dst[i] = 17;
			break;
		case 'Y':
			dst[i] = 18;
			break;
		case 'V':
			dst[i] = 19;
			break;
		case '*':
			dst[i] = 20;
			break;
		default:
			dst[i] = 20;
			break;
		}
	}
	dst[12] = '\0';
	return dst;
}
Pep_data Import_Peptide(const string &FP) {
	cout << "Import_Peptide:" << FP << endl;
	char buf[256];
	FILE *stream;
	fopen_s(&stream, FP.c_str(), "r");
	vector<vector<string>> v;
	while (fgets(buf, 255, stream) != NULL) {
		vector<string> vbuf = split_comma(buf);
		if (stoi(vbuf[1]) < OPs.PepLimit)
			break;
		else {
			if (vbuf[0].length() == PEPTIDE_LENGTH)
				v.push_back(vbuf);
		}
	}
	fclose(stream);
	Pep_data data;
	char *pep = new char[v.size()*PEPTIDE_LENGTH];
	float *read = new float[v.size()];
	char *pP = pep;
	float *pR = read;
	size_t maxread = 1;
	maxread = stoi(v[0][1]);
	for (int i = 0; i < v.size(); ++i) {
		for (int j = 0; j < v[i][0].length(); ++j, ++pP) {
			*pP = S2C(v[i][0][j]);
		}
		if (stof(v[i][1]) <= maxread && stof(v[i][1]) > 0)
			read[i] = stof(v[i][1]);
		else
			read[i] = 1;
	}
	data.pep = pep;
	data.peplength = v.size()*PEPTIDE_LENGTH;
	data.read = read;
	data.totalPeps = v.size();
	return data;
}

inline void CalcDistance_KL(float *Q, const PV &P, float *d) {
	// Distance calculation by Kullback–Leibler divergence

	int slen = P.scoresize;
	int targetsize = P.targetsize;
	int s = 0;
	int i = 0;
	int j = 0;
	float *pQ = Q;
	float *pP = P.scorelist;
	float *pd = d;
	//float Plimit = 1 / P.totalread;
	double Plimit = (1.0 / (double)P.totalread) / FloatSHIFT;

	for (; i < targetsize; ++i) {
		for (j = 0; j < P.Seqlength[i]; ++j, ++pd) {
			for (s = 0; s < slen; ++s, ++pQ, ++pP) {
				if (j < P.Seqlength[i] - 12) {
					*pd += (double)*pP / FloatSHIFT * log((double)*pP / FloatSHIFT / ((double)*pQ / FloatSHIFT));
					/*
					if (*pQ / FloatSHIFT >= Plimit) {
						if (*pP / FloatSHIFT < Plimit && *pQ>0) {
							*pd += abs(Plimit * log(Plimit / (double)*pQ / FloatSHIFT));
						}
						else {
							// cout << abs(log(*fPs / *fQ)) << endl;
							*pd += abs((double)*pP / FloatSHIFT * log((double)*pP / FloatSHIFT / ((double)*pQ / FloatSHIFT)));
						}
					}
					*/
				}
			}
		}
	}
}
inline void CalcDistance_reverseKL(float *Q, const PV &P, float *d) {
	// Distance calculation by reverse Kullback–Leibler divergence

	int slen = P.scoresize;
	int targetsize = P.targetsize;
	int s = 0;
	int i = 0;
	int j = 0;
	float *pQ = Q;
	float *pP = P.scorelist;
	float *pd = d;
	//float Plimit = 1 / P.totalread;
	double Plimit = (float)(1.0 / (double)P.totalread) / FloatSHIFT;

	for (; i < targetsize; ++i) {
		for (j = 0; j < P.Seqlength[i]; ++j, ++pd) {
			for (s = 0; s < slen; ++s, ++pQ, ++pP) {
				if (j < P.Seqlength[i] - 12) {
					if (*pQ / FloatSHIFT >= Plimit) {
						if (*pP / FloatSHIFT < Plimit && *pQ>0) {
							*pd += abs((double)*pQ / FloatSHIFT * log((double)*pQ / FloatSHIFT / Plimit));
						}
						else {
							// cout << abs(log(*fPs / *fQ)) << endl;
							*pd += abs((double)*pQ / FloatSHIFT * log((double)*pQ / FloatSHIFT / ((double)*pP / FloatSHIFT)));
						}
					}
				}
			}
		}
	}
}
inline void CalcDistance_Ratio(float *Q, const PV &P, float *d) {
	// Distance calculation by log(P/Q) 
	int slen = P.scoresize;
	int targetsize = P.targetsize;
	int s = 0;
	int i = 0;
	int j = 0;
	float *pQ = Q;
	float *pP = P.scorelist;
	float *pd = d;
	//float Plimit = 1 / P.totalread;
	double Plimit = (float)(1.0 / (double)P.totalread) / FloatSHIFT;

	for (; i < targetsize; ++i) {
		for (j = 0; j < P.Seqlength[i]; ++j, ++pd) {
			for (s = 0; s < slen; ++s, ++pQ, ++pP) {
				if (j < P.Seqlength[i] - 12) {
					if (*pQ / FloatSHIFT >= Plimit) {
						if (*pP / FloatSHIFT < Plimit && *pQ>0) {
							*pd += (double)*pQ / FloatSHIFT / Plimit;
						}
						else {
							*pd += ((double)*pQ / FloatSHIFT) / ((double)*pP / FloatSHIFT);
						}
					}
				}
			}
		}
	}
}
inline void CalcDistance_Ratio_log(float *Q, const PV &P, float *d) {
	// Distance calculation by log(P/Q) 
	int slen = P.scoresize;
	int targetsize = P.targetsize;
	int s = 0;
	int i = 0;
	int j = 0;
	float *pQ = Q;
	float *pP = P.scorelist;
	float *pd = d;
	//float Plimit = 1 / P.totalread;
	double Plimit = (float)(1.0 / (double)P.totalread) / FloatSHIFT;

	for (; i < targetsize; ++i) {
		for (j = 0; j < P.Seqlength[i]; ++j, ++pd) {
			for (s = 0; s < slen; ++s, ++pQ, ++pP) {
				if (j < P.Seqlength[i] - 12) {
					if (*pQ / FloatSHIFT >= Plimit) {
						if (*pP / FloatSHIFT < Plimit && *pQ>0) {
							*pd += abs(log((double)*pQ / FloatSHIFT / Plimit));
						}
						else {
							*pd += abs(log((double)*pQ / FloatSHIFT / ((double)*pP / FloatSHIFT)));
						}
					}
				}
			}
		}
	}
}
inline void CalcDistance_L1(float *Q, const PV &P, float *d) {
	// Distance calculation by log(P/Q) 
	int slen = P.scoresize;
	int targetsize = P.targetsize;
	int s = 0;
	int i = 0;
	int j = 0;
	int count = 0;
	float *pQ = Q;
	float *pP = P.scorelist;
	float *pd = d;
	//float Plimit = 1 / P.totalread;
	double Plimit = (float)(1.0 / (double)P.totalread) / FloatSHIFT;

	for (; i < targetsize; ++i) {
		for (j = 0; j < P.Seqlength[i]; ++j, ++pd) {
			for (count = 0, s = 0; s < slen; ++s, ++pQ, ++pP) {
				if (j < P.Seqlength[i] - 12) {
					if (*pQ / FloatSHIFT >= Plimit) {
						if (*pP / FloatSHIFT < Plimit && *pQ>0) {
							*pd += abs((double)*pQ / FloatSHIFT - Plimit);
						}
						else {
							*pd += abs((double)*pQ / FloatSHIFT - (double)*pP / FloatSHIFT);
						}
						count++;
					}

				}
			}
			if (j < P.Seqlength[i] - 12) {
				*pd = *pd / count;
			}
		}
	}

}
inline void CalcDistance_L2(float *Q, const PV &P, float *d) {
	// Distance calculation by log(P/Q) 
	int slen = P.scoresize;
	int targetsize = P.targetsize;
	int s = 0;
	int i = 0;
	int j = 0;
	int count = 0;
	float *pQ = Q;
	float *pP = P.scorelist;
	float *pd = d;
	//float Plimit = 1 / P.totalread;
	double Plimit = (float)(1.0 / (double)P.totalread) / FloatSHIFT;

	for (; i < targetsize; ++i) {
		for (j = 0; j < P.Seqlength[i]; ++j, ++pd) {
			for (count = 0, s = 0; s < slen; ++s, ++pQ, ++pP) {
				if (j < P.Seqlength[i] - 12) {
					if (*pQ / FloatSHIFT >= Plimit) {
						if (*pP / FloatSHIFT < Plimit && *pQ>0) {
							*pd += pow(((double)*pQ / FloatSHIFT - Plimit), 2);
						}
						else {
							*pd += pow(((double)*pQ / FloatSHIFT - (double)*pP / FloatSHIFT), 2);
						}
						count++;
					}
				}
			}
			if (j < P.Seqlength[i] - 12) {
				*pd = *pd / count;
			}
		}
	}
}

inline void CalcDistance_Pearson(float *Q, const PV &P, float *d) {

	// Distance calculation by Pearson X2
	int slen = P.scoresize;
	int targetsize = P.targetsize;
	int s = 0;
	int i = 0;
	int j = 0;
	float *pQ = Q;
	float *pP = P.scorelist;
	float *pd = d;
	//float Plimit = 1 / P.totalread;
	double Plimit = (float)(1.0 / (double)P.totalread) / FloatSHIFT;

	for (; i < targetsize; ++i) {
		for (j = 0; j < P.Seqlength[i]; ++j, ++pd) {
			*pd = 0;
			for (s = 0; s < slen; ++s, ++pQ, ++pP) {
				if (j < P.Seqlength[i] - 12) {
					if (*pQ / FloatSHIFT >= Plimit) {

						if (*pP / FloatSHIFT < Plimit && *pQ>0) {
							*pd += pow(((double)*pQ / FloatSHIFT - Plimit), 2) / Plimit;
						}
						else {
							if (*pP > 0)
								*pd += pow(((double)*pQ / FloatSHIFT - (double)*pP / FloatSHIFT), 2) / ((double)*pP / FloatSHIFT);
						}
					}
				}
			}
		}
	}
}
inline void CalcDistance_PearsonX3(float *Q, const PV &P, float *d) {

	// Distance calculation by Pearson X2
	int slen = P.scoresize;
	int targetsize = P.targetsize;
	int s = 0;
	int i = 0;
	int j = 0;
	float *pQ = Q;
	float *pP = P.scorelist;
	float *pd = d;
	//float Plimit = 1 / P.totalread;
	double Plimit = (float)(1.0 / (double)P.totalread) / FloatSHIFT;

	for (; i < targetsize; ++i) {
		for (j = 0; j < P.Seqlength[i]; ++j, ++pd) {
			*pd = 0;
			for (s = 0; s < slen; ++s, ++pQ, ++pP) {
				if (j < P.Seqlength[i] - 12) {
					if (*pQ / FloatSHIFT >= Plimit) {

						if (*pP / FloatSHIFT < Plimit && *pQ>0) {
							*pd += pow(((double)*pQ / FloatSHIFT - Plimit), 3) / Plimit;
						}
						else {
							if (*pP > 0)
								*pd += pow(((double)*pQ / FloatSHIFT - (double)*pP / FloatSHIFT), 3) / ((double)*pP / FloatSHIFT);
						}
					}
				}
			}
		}
	}
}
inline void CalcDistance_SquaredHelinger(float *Q, const PV &P, float *d) {

	// Distance calculation by SquaredHelinger
	int slen = P.scoresize;
	int targetsize = P.targetsize;
	int s = 0;
	int i = 0;
	int j = 0;
	float *pQ = Q;
	float *pP = P.scorelist;
	float *pd = d;
	double Plimit = (float)(1.0 / (double)P.totalread) / FloatSHIFT;

	for (; i < targetsize; ++i) {
		for (j = 0; j < P.Seqlength[i]; ++j, ++pd) {
			for (s = 0; s < slen; ++s, ++pQ, ++pP) {
				if (j < P.Seqlength[i] - 12) {
					if (*pQ / FloatSHIFT >= Plimit) {
						if (*pP / FloatSHIFT < Plimit && *pQ>0) {
							*pd += pow((sqrt(Plimit) - sqrt((double)*pQ / FloatSHIFT)), 2);
						}
						else {
							*pd += pow((sqrt((double)*pP / FloatSHIFT) - sqrt((double)*pQ / FloatSHIFT)), 2);
						}
					}
				}
			}
		}
	}
}
inline void CalcDistance_JensonShannon(float *Q, const PV &P, float *d) {

	// Distance calculation by Jenson-Shannon
	int slen = P.scoresize;
	int targetsize = P.targetsize;
	int s = 0;
	int i = 0;
	int j = 0;
	float *pQ = Q;
	float *pP = P.scorelist;
	float *pd = d;
	double Plimit = (float)(1.0 / (double)P.totalread) / FloatSHIFT;

	for (; i < targetsize; ++i) {
		for (j = 0; j < P.Seqlength[i]; ++j, ++pd) {
			for (s = 0; s < slen; ++s, ++pQ, ++pP) {
				if (j < P.Seqlength[i] - 12) {
					if (*pQ >= Plimit) {
						if (*pP < Plimit && *pQ>0) {
							*pd += Plimit * (log((2 * Plimit) / ((Plimit) * ((double)*pQ / FloatSHIFT)))) + (double)*pQ / FloatSHIFT * (log(2 * ((double)*pQ / FloatSHIFT) / (((Plimit) * ((double)*pQ) / FloatSHIFT))));
						}
						else {
							*pd += (double)*pP / FloatSHIFT * (log((2 * (double)*pP / FloatSHIFT) / (((double)*pP / FloatSHIFT) * ((double)*pQ / FloatSHIFT)))) + (double)*pQ / FloatSHIFT * (log(2 * ((double)*pQ / FloatSHIFT) / ((((double)*pP / FloatSHIFT) * ((double)*pQ) / FloatSHIFT))));
						}
					}
				}
			}
			if (j < P.Seqlength[i] - 12) {
				*pd = *pd / 2;
			}
		}
	}
}

inline void CalcDistance(float *Q, const PV &P, float *d) {
	if (OPs.Distance_function == 0)
		CalcDistance_KL(Q, P, d);
	else if (OPs.Distance_function == 1)
		CalcDistance_Ratio(Q, P, d);
	else if (OPs.Distance_function == 2)
		CalcDistance_reverseKL(Q, P, d);
	else if (OPs.Distance_function == 3)
		CalcDistance_Pearson(Q, P, d);
	else if (OPs.Distance_function == 4)
		CalcDistance_SquaredHelinger(Q, P, d);
	else if (OPs.Distance_function == 5)
		CalcDistance_JensonShannon(Q, P, d);
	else if (OPs.Distance_function == 6)
		CalcDistance_Ratio_log(Q, P, d);
	else if (OPs.Distance_function == 7)
		CalcDistance_L1(Q, P, d);
	else if (OPs.Distance_function == 8)
		CalcDistance_L2(Q, P, d);
	else if (OPs.Distance_function == 30)
		CalcDistance_PearsonX3(Q, P, d);
}

inline void CalcReadRatio(float *Q, const PV &P, double *d) {
	int slen = P.scoresize;
	int targetsize = P.targetsize;
	int s = 0;
	int i = 0;
	int j = 0;
	float *pQ = Q;
	float *pP = P.scorelist;
	double *pd = d;


	for (; i < targetsize; ++i) {
		for (j = 0; j < P.Seqlength[i]; ++j, pd += 12) {
			*pd = 0;
			for (s = 0; s < slen; ++s, ++pQ, ++pP) {
				if (j < P.Seqlength[i] - 12) {
					for (int p = 0; p < 12; ++p) {
						if (*pP / FloatSHIFT >= P_TH[p]) {
							*(pd + p) = *pQ;
						}
					}
				}
			}
		}
	}


}
void Export_Q2(float *Q, const PV &P, const string &dir, const string &barcode) {
	// cout << "\rExport Q: " << dir + "\\Qlist\\" << barcode;
	Directory_check(dir + "\\Qlist2\\" + barcode);
	int targetsize = P.targetsize;
	int slen = P.scoresize;
	float *pQ = Q;
	for (int i = 0; i < targetsize; ++i) {
		// cout << "\rExport P : " << dir + "/Plist/" + P[i].filename + ".csv"
		FILE *stream;
		fopen_s(&stream, (dir + "\\Qlist2\\" + barcode + "\\" + P.filename[i] + ".csv").c_str(), "w");
		for (int j = 0; j < P.Seqlength[i]; ++j) {
			for (int k = 0; k < slen; ++k, ++pQ) {
				fprintf(stream, "%f,", *pQ);
			}
			fprintf(stream, "\n");
		}
		fclose(stream);
	}
}

void Export_Q(float *Q, const PV &P, const string &dir, const string &barcode) {
	// cout << "\rExport Q: " << dir + "\\Qlist\\" << barcode;
	Directory_check(dir + "\\Qlist\\" + barcode);
	int targetsize = P.targetsize;
	int slen = P.scoresize;
	float *pQ = Q;
	for (int i = 0; i < targetsize; ++i) {
		// cout << "\rExport P : " << dir + "/Plist/" + P[i].filename + ".csv"
		FILE *stream;
		fopen_s(&stream, (dir + "\\Qlist\\" + barcode + "\\" + P.filename[i] + ".csv").c_str(), "w");
		for (int j = 0; j < P.Seqlength[i]; ++j) {
			for (int k = 0; k < slen; ++k, ++pQ) {
				fprintf(stream, "%f,", *pQ);
			}
			fprintf(stream, "\n");
		}
		fclose(stream);
	}
}
void Export_P(const PV &P, const string &dir) {

	Directory_check(dir + "\\Plist");
	int targetsize = P.targetsize;
	int slen = P.scoresize;
	float *pP = P.scorelist;
	for (int i = 0; i < targetsize; ++i) {
		// cout << "\rExport P : " << dir + "/Plist/" + P[i].filename + ".csv"
		FILE *stream;
		fopen_s(&stream, (dir + "\\Plist\\" + P.filename[i] + ".csv").c_str(), "w");
		for (int j = 0; j < P.Seqlength[i]; ++j) {
			for (int k = 0; k < slen; ++k, ++pP) {
				fprintf(stream, "%f,", *pP);
			}
			fprintf(stream, "\n");
		}
		fclose(stream);
	}
}
inline void Export_result_header(const string &export_d, const string &fname) {
	Directory_check(export_d);
	string FP = export_d + "/" + fname + ".csv";
	FILE *stream;
	fopen_s(&stream, (FP).c_str(), "w");
	if (stream == NULL) {
		printf("%s file not open!\n", (FP).c_str());
	}
	else {
		// make header line
		// printf("Export: %s\n", (FP).c_str());
		fprintf(stream, ",pos,AA,read,distance");
		for (int n = 0; n < 21; ++n) {
			fprintf(stream, ",%c", AAlist[n]);
		}
		fprintf(stream, "\n");
		fclose(stream);
	}
}
inline void Export_result_bin(float *D, float *AA, u_int *read, const string &export_d, const string barcode, const PV &P) {
	Directory_check(export_d);

	string FP = export_d + "\\" + barcode + "_" + to_string(P.totalAA) + ".bin";

	FILE *stream;

	fopen_s(&stream, (FP).c_str(), "a");
	if (stream == NULL) {
		printf("%s file not open!\n", (FP).c_str());
	}
	else {


		char *buf = new char[1024 * 1024];
		setvbuf(stream, buf, _IOFBF, 512 * 512);
		int len = P.totalAA;
		u_int tlen = 0;
		u_int tsize = P.targetsize;
		string ID = "";
		int pos = 0;
		char *pAA = P.Seq;
		float *pD = D;
		u_int *pR = read;
		float* pRatio = AA;
		int i = 0; int j = 0; int k = 0;
		for (i = 0; i < tsize; ++i) {
			tlen = P.Seqlength[i];
			ID = P.filename[i].substr(1, P.filename[i].length() - 4);
			for (j = 0; j < tlen; ++j, ++pAA, ++pD, ++pR, ++pos) {
				fwrite(ID.c_str(), 1, ID.length(), stream);
				fwrite(&j, sizeof(int), 1, stream);
				fwrite(&AAlist[P.Seq[pos]], 1, 1, stream);
				fwrite(pR, sizeof(float), 1, stream);
				fwrite(pD, sizeof(float), 1, stream);

				for (k = 0; k < 21; ++k, ++pRatio) {
					fwrite(pRatio, sizeof(float), 1, stream);
				}
			}
		}

		fclose(stream);
		delete[] buf;
	}
	/*
	ofstream ofs(FP);

	int len = P.totalAA;
	u_int tlen = 0;
	u_int tsize = P.targetsize;
	string ID = "";
	int pos = 0;
	char *pAA = P.Seq;
	float *pD = D;
	u_int *pR = read;
	float* pRatio = AA;
	int i = 0; int j = 0; int k = 0;
	for (i = 0; i < tsize; ++i) {
		tlen = P.Seqlength[i];
		ID = P.filename[i].substr(1, P.filename[i].length() - 4);
		for (j = 0; j < tlen; ++j, ++pAA, ++pD, ++pR, ++pos) {
			ofs << ID.c_str() << "," << j << "," << AAlist[P.Seq[pos]] << "," << *pR << "," << *pD << ",";
			for (k = 0; k < 21; ++k, ++pRatio) {
				ofs << *pRatio;
			}
			ofs << "\n";
		}

	}
	ofs.close();
	*/

}
inline void checknun(float *d, size_t len) {
	for (int i = 0; i < len; ++i) {
		if (isnan(d[i]))
			d[i] = 0;
	}
}

inline void Export_result(float *D, float *AA, float *read, const string &export_d, const string barcode, const PV &P) {
	Directory_check(export_d);

	string FP = export_d + "\\" + barcode + ".csv";

	FILE *stream;

	fopen_s(&stream, (FP).c_str(), "a");
	if (stream == NULL) {
		printf("%s file not open!\n", (FP).c_str());
	}
	else {
		char *buf = new char[1024 * 1024 * 1024];
		setvbuf(stream, buf, _IOFBF, 1024 * 1024 * 1024);
		int len = P.totalAA;
		u_int tlen = 0;
		u_int tsize = P.targetsize;
		string ID = "";
		int pos = 0;
		char *pAA = P.Seq;
		float *pD = D;
		float *pR = read;
		float* pRatio = AA;
		int i = 0; int j = 0; int k = 0;
		for (i = 0; i < tsize; ++i) {
			tlen = P.Seqlength[i];
			ID = P.filename[i].substr(1, P.filename[i].length() - 4);
			for (j = 0; j < tlen; ++j, ++pAA, ++pD, ++pR, ++pos) {
				//fprintf(stream, "%s,%d,%c,%f,%f", ID.c_str(), j, AAlist[P.Seq[pos]], *pR, *pD);
				fprintf(stream, "%s,%d,%c,%f,%f", ID.c_str(), j, AAlist[P.Seq[pos]], *pR, *pD);
				for (k = 0; k < 21; ++k, ++pRatio) {
					fprintf(stream, ",%f", *pRatio);
				}
				fprintf(stream, "\n");
			}

		}

		fclose(stream);
		delete[] buf;
	}
}
inline void Export_result2(float *D, float *AA, float *read, const string &export_d, const string barcode, const PV &P, double *RR) {
	Directory_check(export_d);

	string FP = export_d + "\\" + barcode + ".csv";

	FILE *stream;

	fopen_s(&stream, (FP).c_str(), "a");
	if (stream == NULL) {
		printf("%s file not open!\n", (FP).c_str());
	}
	else {
		char *buf = new char[1024 * 1024 * 1024];
		setvbuf(stream, buf, _IOFBF, 1024 * 1024 * 1024);
		int len = P.totalAA;
		u_int tlen = 0;
		u_int tsize = P.targetsize;
		string ID = "";
		int pos = 0;
		char *pAA = P.Seq;
		float *pD = D;
		float *pR = read;
		double *pRR = RR;
		float* pRatio = AA;
		int i = 0; int j = 0; int k = 0;
		for (i = 0; i < tsize; ++i) {
			tlen = P.Seqlength[i];
			ID = P.filename[i].substr(1, P.filename[i].length() - 4);
			for (j = 0; j < tlen; ++j, ++pAA, ++pD, ++pR, ++pos) {
				//fprintf(stream, "%s,%d,%c,%f,%f", ID.c_str(), j, AAlist[P.Seq[pos]], *pR, *pD);
				fprintf(stream, "%s,%d,%c,%f,%f", ID.c_str(), j, AAlist[P.Seq[pos]], *pR, *pD);
				for (k = 0; k < 12; ++k, ++pRR) {
					fprintf(stream, ",%lf", *pRR);
				}
				for (k = 0; k < 21; ++k, ++pRatio) {
					fprintf(stream, ",%f", *pRatio);
				}
				fprintf(stream, "\n");
			}

		}

		fclose(stream);
		delete[] buf;
	}
	/*
	ofstream ofs(FP);

	int len = P.totalAA;
	u_int tlen = 0;
	u_int tsize = P.targetsize;
	string ID = "";
	int pos = 0;
	char *pAA = P.Seq;
	float *pD = D;
	u_int *pR = read;
	float* pRatio = AA;
	int i = 0; int j = 0; int k = 0;
	for (i = 0; i < tsize; ++i) {
		tlen = P.Seqlength[i];
		ID = P.filename[i].substr(1, P.filename[i].length() - 4);
		for (j = 0; j < tlen; ++j, ++pAA, ++pD, ++pR, ++pos) {
			ofs << ID.c_str() << "," << j << "," << AAlist[P.Seq[pos]] << "," << *pR << "," << *pD << ",";
			for (k = 0; k < 21; ++k, ++pRatio) {
				ofs << *pRatio;
			}
			ofs << "\n";
		}

	}
	ofs.close();
	*/

}
inline void Export_result_Q(const string &export_d, const string &fname, int length) {
	Directory_check(export_d);
	string FP = export_d + "/" + fname + ".csv";
	FILE *stream;
	// mtxf.lock();
	// cout << export_d + "/" + fname + ".csv" << endl;
	fopen_s(&stream, (FP).c_str(), "w");
	if (stream == NULL) {
		printf("%s file not open!\n", (FP).c_str());
	}
	else {
		// char buf[512 * 512];
		// setvbuf(stream, buf, _IOFBF, sizeof(buf));
		int pro = 0;
		for (const auto &i : queue) {
			if (queflag[pro])
				for (const auto &j : i) {
					for (const string &k : j) {
						fprintf(stream, "%s,", k);
					}
					fprintf(stream, "\n");
				}
			++pro;
		}

		fclose(stream);
	}
	// mtxf.unlock();
	// delete[] buf;
}
inline char* import_target(const string &FP, int len) {
	FILE *stream;
	char* seq = new char[len];
	cout << "\rimport seq :" << FP << "\t\t";
	//cout << "length: " << len << endl;
	int ret = fopen_s(&stream, (FP).c_str(), "rb");
	if (ret == 0) {
		// cout << result.scoresize <<endl;
		int iret = fread(seq, sizeof(char), len, stream);
		fclose(stream);
		if (iret != len) cout << "fread error" << iret << endl;
	}
	else {
		cout << "\rimport seq :" << FP << "error" << endl;
	}
	return seq;
}
inline void Check_P_func(const string &FP, const vector<string> &fname, vector<int> &len, vector<string> &fn, PV &P) {

	FILE *stream;
	int i = 0;
	while (1) {
		mtx.lock();
		i = task;
		task++;
		if (i >= fname.size())
		{
			mtx.unlock();
			break;
		}
		cout << i << "\rCheck : " << FP << "\\" << fname[i] << "\t\t";
		mtx.unlock();
		//cout << "\rCheck : " << FP << "\\" << fname[i] << "\t\t";
		fopen_s(&stream, (FP + "\\" + fname[i]).c_str(), "rb");
		int buf[2] = {};
		int iret = fread(buf, sizeof(int), 2, stream);
		if (iret != 0) {
			len[i] = buf[0];
			fn[i] = fname[i];
			if (i == 0)
				P.scoresize = (int)buf[1];
			fclose(stream);
		}
		//	cout << i << "\rCheck end : " << FP << "\\" << fname[i] << "\n\t";
		//	Sleep(1000);
	}
}
inline void import_Seq(const string &FP, vector<int> &len, vector<string> &fn, char *seq) {
	FILE *stream;
	int i = 0;
	while (1) {
		mtx.lock();
		i = task;
		task++;
		if (i >= fn.size())
		{
			mtx.unlock();
			break;
		}
		cout << "\rimport seq :" << fn[i].substr(1) << "\t\t";
		mtx.unlock();
		size_t pos = 0;
		for (int n = 0; n < i; ++n) {
			//cout << len[n] << endl;
			pos += len[n];
		}
		FILE *stream;
		char* s = new char[len[i]];

		int ret = fopen_s(&stream, (FP + "\\Seq" + fn[i].substr(1)).c_str(), "rb");
		if (ret == 0) {
			// cout << result.scoresize <<endl;
			int iret = fread(s, sizeof(char), len[i], stream);
			fclose(stream);
			if (iret != len[i]) cout << "fread error" << iret << endl;
		}
		else {
			cout << "\rimport seq :" << FP << "error" << endl;
		}
		char *pseq = &seq[pos];
		for (int n = 0; n < len[i]; ++n, ++s, ++pseq) {
			*pseq = *s;
			//cout << AAlist[*s] << endl;
		}
		//mtx.unlock();
	}

}
inline void Import_P_func(const string &FP, vector<int> &len, vector<string> &fn, PV &P, size_t totallength, float *buf) {
	FILE *stream;
	int i = 0;
	while (1) {
		mtx.lock();
		i = task;
		task++;
		if (i >= fn.size())
		{
			mtx.unlock();
			break;
		}
		cout << "\rImport : " << FP << "\\" << fn[i] << "\t\t";
		mtx.unlock();
		//cout << len[i] << endl;
		double* b = new double[len[i] * P.scoresize];
		fopen_s(&stream, (FP + "\\" + fn[i]).c_str(), "rb");
		_fseeki64(stream, 2 * sizeof(int), SEEK_SET);
		size_t pos = 0;
		for (int n = 0; n < i; ++n) {
			pos += len[n];
		}
		//size_t pos = std::accumulate(len.begin(), len.begin() + i - 1, 0);
		int ret = fread(b, sizeof(double), len[i] * P.scoresize, stream);
		fclose(stream);
		if (ret == len[i] * P.scoresize) {
			for (int j = 0; j < len[i] * P.scoresize; ++j) {
				buf[pos*P.scoresize + j] = (float)(b[j] * FloatSHIFT);//shift 
			}
		}
		else {
			cout << "Import P error: " << (FP + "\\" + fn[i]) << endl;
		}
		delete[] b;
		//cout << "\rImport end : " << FP << "\\" << fn[i] << "\t\t";
		//mtx.unlock();
	}

}

PV Import_P_mt(const string &FP, const vector<string> &fname) {
	PV P;
	vector<thread> threads(10);
	vector<int> len(fname.size());
	vector<string> fn(fname.size());
	task = 0;
	for (int g = 0; g < threads.size(); ++g) {
		threads[g] = thread([&] {Check_P_func(FP, fname, len, fn, P); });
	}
	for (auto &t : threads) {
		t.join();
	}
	cout << "\nCheck_P Finished" << endl;
	size_t totallength = accumulate(len.begin(), len.end(), 0);
	cout << "total length :" << totallength << " : Score size " << P.scoresize << endl;
	char *seq = new char[totallength];
	float *AA = new float[totallength * 21];
	float *Plist = new float[totallength * P.scoresize];
	string *filename = new string[fname.size()];
	int *seqlen = new int[fname.size()];
	task = 0;
	cout << "\nImport_P start" << endl;
	//COUT_VEC_String(fn);
	cout << "Plist size: " << totallength * P.scoresize << endl;
	for (int g = 0; g < threads.size(); ++g) {
		threads[g] = thread([&] {Import_P_func(FP, len, fn, P, totallength, Plist); });
	}
	for (auto &t : threads) {
		t.join();
	}

	cout << "\nImport_P Finished" << endl;
	task = 0;
	for (int g = 0; g < threads.size(); ++g) {
		threads[g] = thread([&] {import_Seq(FP, len, fn, seq); });
	}
	for (auto &t : threads) {
		t.join();
	}

	for (int i = 0; i < fn.size(); ++i) {
		filename[i] = fn[i];
		seqlen[i] = len[i];
	}


	P.filename = filename;
	P.Seq = seq;
	P.scorelist = Plist;
	P.Seqlength = seqlen;
	P.targetsize = fn.size();
	P.totalAA = totallength;
	return P;
}

void ReadP_listcondition(vector<int> &v, const string &FP,
	const vector<string> &fname, int &size) {
	FILE *stream;
	int i = 0;
	int s = v.size();
	while (i < s) {
		mtxt.lock();
		i = task;
		task++;
		if (i < s) {
			cout << "\r" << fname[i];
			mtxt.unlock();
		}
		else {
			mtxt.unlock();
			break;
		}

		fopen_s(&stream, (FP + "/" + fname[i]).c_str(), "rb");
		int buf[2] = {};
		int iret = fread(buf, sizeof(int), 2, stream);
		v[i] = buf[0];
		size = buf[1];
		fclose(stream);
	}
}

void Calc_P_LessThan_Score(PV &P, int size) {
	double sum = 0;
	float *pS = P.scorelist;
	float a = 0;
	int slen = P.scoresize;
	//	cout << "\rP convert : " << t << "/" << P.size();
	for (int i = 0; i < P.targetsize; ++i) {
		int len = P.Seqlength[i];
		for (int j = 0; j < len; ++j) {
			sum = 0;
			for (int k = 0; k < slen; ++k, pS++) {
				a = *pS;

				*pS = (float)(FloatSHIFT - sum);
				if (*pS < 0)
					*pS = 0.0;
				sum += a;
			}
		}
	}
}


__global__ void Convertratio(float *d_Scoreratio, float *d_AAratio, int slen, size_t size) {

	size_t idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx > size) return;
	int i = 0;
	int j = 0;
	float a = 1.0;
	float b = FloatSHIFT;
	float sum = 0.0;
	for (int i = 0; i < slen; ++i) {
		sum += d_Scoreratio[idx * slen + i];
	}
	if (sum == 0.0) {
		for (int i = 0; i < slen; ++i)
			d_Scoreratio[idx * slen + i] = -1;
	}
	else {
		for (int i = 0; i < slen; ++i) {
			a = d_Scoreratio[idx * slen + i];
			d_Scoreratio[idx * slen + i] = b;
			b = b - (a / sum)*FloatSHIFT;
			//printf("%f, ",b);
		}
	}

	sum = 0.0;
	for (int i = 0; i < 21; ++i) {
		sum += d_AAratio[idx * 21 + i];
	}
	if (sum == 0) {
		for (int i = 0; i < 21; ++i)
			d_AAratio[idx * 21 + i] = 0;
	}
	else {
		for (int i = 0; i < 21; ++i) {
			//printf("%f,%f\n", d_AAratio[idx * 21 + i], sum);
			d_AAratio[idx * 21 + i] = d_AAratio[idx * 21 + i] / sum;
		}
	}
}
__global__ void reset_AA(float *d_Scoreratio, float *d_AAratio, float *d_read_TH, int slen, size_t threads) {
	size_t idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx > threads) return;
	for (int i = 0; i < 21; ++i)
		d_AAratio[idx * 21 + i] = 0;
	for (int i = 0; i < slen; ++i)
		d_Scoreratio[idx * slen + i] = 0;
	d_read_TH[idx] = 0;
}
__global__ void match(char* d_AAMap, char *d_pep, char *d_pro, float *d_Scoreratio, float *d_AAratio, float *P,
	float *d_read_TH, float *d_read, float TH, int slen, size_t pepsize, size_t threads, int plength) {
	size_t idx = blockDim.x * blockIdx.x + threadIdx.x;
	if (idx > threads) return;
	int i = 0;
	int j = 0;
	size_t it = 0;
	int score = 0;
	for (int n = 0; n < plength; ++n) {
		score = 0;
		d_pep[idx*PEPTIDE_LENGTH];
		for (j = 0; j < PEPTIDE_LENGTH; ++j) {
			score += d_AAMap[d_pep[idx*PEPTIDE_LENGTH + j] * 21 + d_pro[n + j]];
		}
		atomicAdd(&d_Scoreratio[n * slen + score], (d_read[idx]) / FloatSHIFT_INT);
		//atomicAdd(&d_Scoreratio[n * slen + score], 1.0);
		//d_Scoreratio[n * slen + score] += ((float)d_read[idx]) / FloatSHIFT_INT;

		if (P[n*slen + score] < TH) {
			for (int j = 0; j < PEPTIDE_LENGTH; ++j) {
				d_AAratio[(n + j) * 21 + d_pep[idx*PEPTIDE_LENGTH + j]] += (d_read[idx]) / FloatSHIFT_INT;
			}
			atomicAdd(&d_read_TH[n], d_read[idx]);
			//d_read_TH[n] += d_read[idx];
			//printf("%d. %f, %u\n", score, P[idx*slen + score], d_read[i]);
		}
		//printf("%d\n", score);
	}

}
__global__  void ChackSeq(char* s, int size) {
	printf("ChackSeq\n");
	for (int i = 0; i < size; ++i) {
		printf("%d", s[i]);
	}
	printf("\n");
}
__global__  void ChackRead(float* s, int size) {
	printf("ChackRead\n");
	for (int i = 0; i < size; ++i) {
		printf("%f\n", s[i]);
	}
	printf("\n");
}

void Matching(const string &dir_name, const vector<string> &sample_list,
	const string &export_d, const PV P_list, string fname, const Pep_data &PepList, mems &m) {
	cout << "Start Matching " << endl;
	chrono::system_clock::time_point start, end, end1, end2, end3, end4, end5;
	start = std::chrono::system_clock::now();
	size_t maxblockx = 1024;
	size_t maxblocky = 1024;
	size_t maxgridx = 2147483647;
	size_t maxgridy = 65535;
	size_t maxgridz = 65535;
	int blocksize = 32;
	// 10496 core
	// max texture dimension size:(131072)(131072,65536)(16384,16384,16384)
	// max dimension size of  a thread block: (1024x1024x64)
	// max dimension size of a grid size : (2147483647,65535,65535)
	size_t threads = PepList.totalPeps;
	dim3 block(blocksize, 1, 1);
	dim3 grid(threads / block.x + 1, 1, 1);
	size_t i = 0;
	double duration = 0;
	double totaltime = 0;
	start = std::chrono::system_clock::now();
	//cout << fname << endl;
	//if(fname=="134")
	//ChackRead << <1,1 >> > ( m.d_read, PepList.totalPeps);
	match << <grid, block >> > (m.d_AAmap, m.d_pep, m.d_pro, m.d_Scoreratio, m.d_AAratio, m.d_P_list, m.d_read_TH, m.d_read, (float)OPs.TH, P_list.scoresize, PepList.totalPeps, threads, P_list.totalAA);
	cudaDeviceSynchronize();
	//	cout << "Calc matching finished" << endl;
	//if (fname == "134")
		//ChackRead << <1, 1 >> > (m.d_read_TH, P_list.totalAA);
	size_t threads_AA = P_list.totalAA;
	dim3 block_AA(1024, 1, 1);
	dim3 grid_AA(threads_AA / block_AA.x + 1, 1, 1);
	Convertratio << <grid_AA, 1024 >> > (m.d_Scoreratio, m.d_AAratio, P_list.scoresize, P_list.totalAA);
	cudaDeviceSynchronize();
	cout << "end Convertratio " << endl;
	cudaMemcpy(m.h_result, m.d_Scoreratio, P_list.totalAA * P_list.scoresize * sizeof(float), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	cudaMemcpy(m.h_AAratio, m.d_AAratio, P_list.totalAA * 21 * sizeof(float), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	cudaMemcpy(m.h_read_TH, m.d_read_TH, P_list.totalAA * sizeof(float), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	cout << "end cudaMemcpy " << endl;
	end1 = std::chrono::system_clock::now();
	//if (OPs.Q_Export)
	//	Export_Q2(m.h_result, P_list, export_d, fname);
	double matchingtime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end1 - start).count());
	cout << "matchingtime:  " << matchingtime / 1000000 << " s" << endl;
	//cout << "CalcDistance" << j << endl;
	//CalcReadRatio(m.h_result, P_list, m.h_Read_Ratio);
	CalcDistance(m.h_result, P_list, m.h_distance);
	end2 = std::chrono::system_clock::now();
	double distancetime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end2 - end1).count());
	cout << "distancetime:  " << distancetime / 1000000 << " s" << endl;
	checknun(m.h_distance, P_list.totalAA);
	Export_result(m.h_distance, m.h_AAratio, m.h_read_TH, export_d, fname, P_list);
	end3 = std::chrono::system_clock::now();
	if (OPs.Q_Export)
		Export_Q(m.h_result, P_list, export_d, fname);
	double exporttime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end3 - end2).count());
	cout << "exporttime:  " << exporttime / 1000000 << " s" << endl;

	reset_AA << <grid_AA, 1024 >> > (m.d_Scoreratio, m.d_AAratio, m.d_read_TH, P_list.scoresize, P_list.totalAA);
	//cudaDeviceSynchronize();
}
void MemFree1(mems &m) {
	cudaFreeHost(m.h_distance);
	cudaFreeHost(m.h_AAratio);
	cudaFreeHost(m.h_Scoreratio);
	cudaFreeHost(m.h_read_TH);
	cudaFreeHost(m.h_result);
	cudaFreeHost(m.h_Read_Ratio);
	cudaFree(m.d_pro);
	cudaFree(m.d_AAmap);
	cudaFree(m.d_P_list);
	cudaFree(m.d_read_TH);
	cudaFree(m.d_AAratio);
	cudaFree(m.d_Scoreratio);

}
void MemFree(mems &m) {
	cout << "MemFree" << endl;
	cudaFree(m.d_pep);
	cudaFree(m.d_read);

}
void MemRefresh(mems &m, const PV &P_list) {
	cout << "MemRefresh" << endl;
	for (int i = 0; i < P_list.totalAA; ++i) {
		m.h_distance[i] = 0;
		m.h_Read_Ratio[i] = 0;
		m.h_read_TH[i] = 0;
		m.h_result[i] = 0;
		for (int j = 0; j < 21; ++j)
			m.h_AAratio[i * 21 + j] = 0;
		for (int j = 0; j < P_list.scoresize; ++j)
			m.h_Scoreratio[i * P_list.scoresize + j] = 0;
	}
}
void MemPrep1(mems &m, const PV &P_list, size_t &GPU_mem, size_t &Host_mem) {
	cout << " Total Protein length: " << P_list.totalAA << " aa" << endl;
	cudaError_t err;
	GPU_mem = P_list.totalAA + 21 * 21 + sizeof(float) * P_list.totalAA * P_list.scoresize +
		sizeof(float) * P_list.totalAA + sizeof(float) * P_list.totalAA * 21 +
		sizeof(float) * P_list.totalAA * P_list.scoresize;

	Host_mem = sizeof(float) *P_list.totalAA + sizeof(float) *P_list.totalAA * 21 + sizeof(float) *P_list.totalAA * P_list.scoresize +
		sizeof(float) *P_list.totalAA +
		sizeof(float) *P_list.totalAA * P_list.scoresize;

	// on host

	CUDA_SAFE_CALL(err = cudaMallocHost(&m.h_distance, sizeof(float) * P_list.totalAA));
	CUDA_SAFE_CALL(err = cudaMallocHost(&m.h_AAratio, sizeof(float) * P_list.totalAA * 21));
	CUDA_SAFE_CALL(err = cudaMallocHost(&m.h_Scoreratio, sizeof(float) * P_list.totalAA * P_list.scoresize));
	CUDA_SAFE_CALL(err = cudaMallocHost(&m.h_read_TH, sizeof(float) * P_list.totalAA));
	CUDA_SAFE_CALL(err = cudaMallocHost(&m.h_result, sizeof(float) * P_list.totalAA * P_list.scoresize));
	CUDA_SAFE_CALL(err = cudaMallocHost(&m.h_Read_Ratio, sizeof(float) * P_list.totalAA * 12));

	MemRefresh(m, P_list);
	// on device
	CUDA_SAFE_CALL(err = cudaMalloc(&m.d_pro, P_list.totalAA));
	CUDA_SAFE_CALL(err = cudaMalloc(&m.d_AAmap, 21 * 21));
	CUDA_SAFE_CALL(err = cudaMalloc(&m.d_read_TH, sizeof(float) * P_list.totalAA));
	CUDA_SAFE_CALL(err = cudaMalloc(&m.d_AAratio, sizeof(float) * P_list.totalAA * 21));
	CUDA_SAFE_CALL(err = cudaMalloc(&m.d_Scoreratio, sizeof(float) * P_list.totalAA * P_list.scoresize));
	CUDA_SAFE_CALL(err = cudaMalloc(&m.d_P_list, sizeof(float) * P_list.totalAA * P_list.scoresize));

	CUDA_SAFE_CALL(err = cudaMemcpy(m.d_AAmap, AAMap, 21 * 21, cudaMemcpyHostToDevice));  // AAmap (21 * 21)
	CUDA_SAFE_CALL(err = cudaMemcpy(m.d_AAratio, m.h_AAratio, P_list.totalAA * 21 * sizeof(float), cudaMemcpyHostToDevice));  //
	CUDA_SAFE_CALL(err = cudaMemcpy(m.d_Scoreratio, m.h_Scoreratio, P_list.totalAA * P_list.scoresize * sizeof(float), cudaMemcpyHostToDevice));  //
	CUDA_SAFE_CALL(err = cudaMemcpy(m.d_P_list, P_list.scorelist, P_list.totalAA * P_list.scoresize * sizeof(float), cudaMemcpyHostToDevice));  //
	CUDA_SAFE_CALL(err = cudaMemcpy(m.d_pro, P_list.Seq, P_list.totalAA, cudaMemcpyHostToDevice));  //



}

void MemPrep(mems &m, Pep_data &PepList, size_t &GPU_mem, size_t &Host_mem) {
	cout << "Total peptide num: " << PepList.totalPeps << endl;
	cudaError_t err;
	GPU_mem += PepList.peplength + sizeof(float) * PepList.totalPeps + sizeof(float) * PepList.totalPeps;
	cout << "GPU mem: " << GPU_mem / 1000000 << " MB" << endl;
	cout << "Host mem: " << Host_mem / 1000000 << " MB" << endl;

	// on device
	cout << "cudaMalloc" << endl;
	char *mem1;
	float *mem2;
	m.d_pep = mem1;
	m.d_read = mem2;
	cout << "d_pep" << endl;
	cout << "malloc: " << PepList.peplength << " byte" << endl;
	int count = 0;
	err = cudaMalloc((void **)&m.d_pep, PepList.peplength);
	while (err != cudaSuccess) {
		Sleep(1000);
		cout << "malloc: " << PepList.peplength << " byte" << endl;
		err = cudaMalloc((void **)&m.d_pep, PepList.peplength);
		count++;
		if (count > 4) {
			exit(err);
		}
	}
	cout << "d_read" << endl;
	err = cudaMalloc((void **)&m.d_read, sizeof(float) * PepList.totalPeps);
	while (err != cudaSuccess) {
		Sleep(1000);
		cout << "d_read" << endl;
		err = cudaMalloc((void **)&m.d_read, sizeof(float) * PepList.totalPeps);
		count++;
		if (count > 4) {
			exit(err);
		}
	}
	cout << "cudaMemcpy" << endl;
	CUDA_SAFE_CALL(err = cudaMemcpy(m.d_pep, PepList.pep, PepList.peplength, cudaMemcpyHostToDevice));  //
	CUDA_SAFE_CALL(err = cudaMemcpy(m.d_read, PepList.read, PepList.totalPeps * sizeof(float), cudaMemcpyHostToDevice));  //
	//ChackSeq << <1, 1 >> > (m.d_pep, 100);
}
void delPep(Pep_data &PepList) {
	delete[] PepList.read;
	delete[] PepList.pep;
}
void GetGPUProf() {
	cudaDeviceProp dev;
	cudaGetDeviceProperties(&dev, 0);
	cout << "--------------------------------------------------------------" << endl;
	cout << "GPU information" << endl;
	cout << "---------------------------------------------------------------" << endl;
	printf("device %d\n", 0);
	printf(" device name : %s\n", dev.name);
	printf(" total global memory : %d (MB)\n", dev.totalGlobalMem / 1024 / 1024);
	printf(" shared memory / block : %d (KB)\n", dev.sharedMemPerBlock / 1024);
	printf(" register / block : %d\n", dev.regsPerBlock);
	printf(" warp size : %d\n", dev.warpSize);
	printf(" max pitch : %d (B)\n", dev.memPitch);
	printf(" max threads / block : %d\n", dev.maxThreadsPerBlock);
	printf(" max size of each dim. of block : (%d, %d, %d)\n", dev.maxThreadsDim[0], dev.maxThreadsDim[1], dev.maxThreadsDim[2]);
	printf(" max size of each dim. of grid  : (%d, %d, %d)\n", dev.maxGridSize[0], dev.maxGridSize[1], dev.maxGridSize[2]);
	printf(" clock rate : %d (MHz)\n", dev.clockRate / 1000);
	printf(" total constant memory : %d (KB)\n", dev.totalConstMem / 1024);
	printf(" compute capability : %d.%d\n", dev.major, dev.minor);
	printf(" alignment requirement for texture : %d\n", dev.textureAlignment);
	printf(" device overlap : %s\n", (dev.deviceOverlap ? "ok" : "not"));
	printf(" num. of multiprocessors : %d\n", dev.multiProcessorCount);
	printf(" kernel execution timeout : %s\n", (dev.kernelExecTimeoutEnabled ? "on" : "off"));
	printf(" integrated : %s\n", (dev.integrated ? "on" : "off"));
	printf(" host memory mapping : %s\n", (dev.canMapHostMemory ? "on" : "off"));

	printf(" compute mode : ");
	if (dev.computeMode == cudaComputeModeDefault) printf("default mode (multiple threads can use) \n");
	else if (dev.computeMode == cudaComputeModeExclusive) printf("exclusive mode (only one thread will be able to use)\n");
	else if (dev.computeMode == cudaComputeModeProhibited) printf("prohibited mode (no threads can use)\n");
	cout << "---------------------------------------------------------------" << endl;
}

int main(int argc, char *argv[]) {
	clock_t start = clock();
	OPs.threads = thread::hardware_concurrency() - 1;
	try {
		for (int i = 0; i < argc; ++i) {
			cout << i << ": " << argv[i] << " : ";
			if (strstr(argv[i], "--pep") != NULL) {
				++i;
				OPs.FPr = argv[i];
				cout << argv[i];
			}
			if (strstr(argv[i], "--thread") != NULL) {
				++i;
				OPs.threads = stoi(argv[i]);
				cout << argv[i];
			}
			if (strstr(argv[i], "--export") != NULL) {
				++i;
				OPs.FPw = argv[i];
				cout << argv[i];
			}
			if (strstr(argv[i], "--P_map") != NULL) {
				++i;
				OPs.P_map = argv[i];
				cout << argv[i];
			}
			if (strstr(argv[i], "--target") != NULL) {
				++i;
				OPs.Target = argv[i];
				cout << argv[i];
			}
			if (strstr(argv[i], "--Plistexport") != NULL) {
				++i;
				OPs.Plistexport = argv[i];
				cout << argv[i];
			}
			if (strstr(argv[i], "--readlimit") != NULL) {
				++i;
				OPs.PepLimit = stoi(argv[i]);
				std::cout << argv[i];
			}

			if (strstr(argv[i], "--Calc") != NULL) {
				++i;
				std::cout << argv[i];
				if (strstr(argv[i], "read") != NULL) {
					OPs.Calc = 0;
					cout << "mode read";
				}
				else {
					OPs.Calc = 1;
					cout << "mode num";
				}
			}
			if (strstr(argv[i], "--QExport") != NULL) {
				OPs.Q_Export = 1;

			}
			if (strstr(argv[i], "--PExport") != NULL) {
				OPs.P_Export = 1;

			}
			if (strstr(argv[i], "--TH") != NULL) {
				++i;
				OPs.TH = stod(argv[i])*FloatSHIFT;
				cout << OPs.TH;
			}
			if (strstr(argv[i], "--Distance_function") != NULL) {
				++i;
				OPs.Distance_function = stod(argv[i]);
				cout << OPs.Distance_function;
			}
			if (strstr(argv[i], "--Start") != NULL) {
				++i;
				OPs.Start = argv[i];
				cout << OPs.Start;
			}
			cout << endl;
		}
	}
	catch (exception &e) {
		cout << "--pep <*> \t: import directory/n";
		cout << "--export <*> \t: export directory\n";
		cout << "--P_map <*>\t: P map file path\n";
		cout << "-- sim_table <*>\t : AA similarity table\n";
		cout << "--target <*>\t: target protein direcroty\n";
		cout << "--thread <*>\t : thread size(default cours - 2)\n";
		cout << "--Plistexport <*>\t : yes(export P list, default <no> )\n";
		cout << "--readlimit <*>\t : default 0\n";
		cout << "--TH <*>\t : Threshold of cluster (default 0.001)\n";
		cout << "--Calc <read or num>\t : Use NGS read or kinds of peptide "
			"(default read)\n";
		cout << "--Distance_function <*>\t : 1: KL, 2:KL_like, 3:KL_plus\n";
		cout << "--QExport <*>\t : export Q list\n";


		return 0;
	}
	Directory_check(OPs.FPw);
	Copyfile((OPs.P_map + "\\Protein_ID_list.txt").c_str(), (OPs.FPw + "\\Protein_ID_list.txt").c_str());
	vector<string> sample_list = Get_File_path_in_dir(OPs.FPr, "csv");
	samplesize = sample_list.size();
	Load_ScoreTable_Limited2(OPs.P_map + "\\sim_table.csv", AAMap, 0, 11);
	cout << OPs.P_map << endl;
	vector<string> file_list = Get_File_path_in_dir(OPs.P_map, "bin");
	vector<string> P_file_list;
	vector<string> S_file_list;
	GetGPUProf();
	for (int n = 0; n < file_list.size(); ++n) {
		if (file_list[n][0] == 'P')
			P_file_list.push_back(file_list[n]);
		else if (file_list[n].substr(0, 3) == "Seq")
			S_file_list.push_back(file_list[n]);
	}
	sort(P_file_list.begin(), P_file_list.end());
	sort(S_file_list.begin(), S_file_list.end());
	cout << "Import P_file\t\t\t" << endl;
	int targetsize = P_file_list.size();
	int length = 0;

	PV P_list = Import_P_mt(OPs.P_map, P_file_list);

	Calc_P_LessThan_Score(P_list, targetsize);
	cout << "\nImport P_file Finsihed\t\t\t" << endl;
	/*for (int n = 0; n < 1000; ++n) {
		for (int m = 0; m < P_list.scoresize; ++m) {
			printf("%f, ", P_list.scorelist[n * P_list.scoresize + m]);
		}
		cout << endl;
	}*/
	//for (int i = 0; i < P_list.totalAA; ++i)
	//	cout << AAlist[P_list.Seq[i]];



	
	for (int i=stoi(OPs.Start); i < samplesize; ++i) {
		cout << "**********************************************" << endl;
		mems memory;
		size_t GPUmem, Hostmem;
		MemPrep1(memory, P_list, GPUmem, Hostmem);
		// queue = new Qdata[targetsize];
		Pep_data PepList = Import_Peptide(OPs.FPr + "\\" + sample_list[i]);
		size_t tRead = 0;
		for (int n = 0; n < PepList.totalPeps; ++n) {
			tRead += PepList.read[n];
		}
		P_list.totalread = tRead;
		string fname = ExtractPathWithoutExt(sample_list[i]);

		cout << "total read: " << P_list.totalread << "read " << endl;
		cout << "total peptide " << PepList.totalPeps << " peps " << endl;

		if (PepList.totalPeps < 2) {
			delPep(PepList);
			continue;
		}
		MemPrep(memory, PepList, GPUmem, Hostmem);
		MemRefresh(memory, P_list);

		Export_result_header(OPs.FPw, fname);

		Matching(OPs.FPr, sample_list, OPs.FPw, P_list, fname, PepList, memory);
		//cudaDeviceSynchronize();
		delPep(PepList);
		MemFree(memory);
		MemFree1(memory);
	}
	if (OPs.P_Export == 1)
		Export_P(P_list, OPs.FPw);
	//cudaDeviceReset();
	cout << "***************************************" << endl << endl;
	cout << " Finished " << endl << endl;
	clock_t end = clock();
	cout << "Total: " << difftime(end, start) / 60 / 1000 << "min\n";
	cout << "************************" << endl << endl;
}
