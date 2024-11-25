#pragma once
#include <stdafx.h>

#include <CheckTheFolder.h>
#include <Windows.h>
#include <algorithm>
#include <cassert>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <io.h>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <process.h>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <thread>
#include <time.h>
#include <unordered_map>
#include <vector>

#define READ_LENGTH_LIMIT 100
#define CDS_LENGTH 12

using namespace std;

mutex mtx;
mutex mtx2;
mutex mtx3;
mutex mtx4;
mutex mtx5;
mutex mtxt;
mutex mtxf;
const char *startseq = "ATACATATG";
const char *rstartseq = "CATATGTAT";
const char *startbarcode = "GTTAACTTTAAGAAGGAGATA";
const char *rstartbarcode = "CGCTGCCGCTGCCGCA";
const char *AAlist = "ARNDCQEGHILKMFPSTWYV*";
// chrono::system_clock::time_point timer;

struct Seq_read {
	int read;
	string seq;

	// char *seq[CDS_LENGTH + 1];
	bool operator<(const Seq_read &left) const {
		return read > left.read; //
	};
};
struct Seq_read_wbarcode {
	int read;
	string barcode;
	string seq;

	// char *seq[CDS_LENGTH + 1];
	bool operator<(const Seq_read &left) const {
		return read > left.read; //
	};
};
struct Fasta {
	string name;
	string seq;
};
struct Fasta_i {
	int ID;
	char *seq;
};
struct Pep_read {
	char *pep;
	u_int read = 0;
	int score = 0;
	char frag = 0;
	bool operator<(const Pep_read &right) const { return score > right.score; };
};
struct Protein_list {
	char *Peotein;
	char *Pro_seq;
};
struct Score_set {
	int pos;
	double score;
	int read;
};

inline chrono::system_clock::time_point Timer_Start() {
	chrono::system_clock::time_point start = chrono::system_clock::now();
	return start;
}
inline void Timer_End(chrono::system_clock::time_point start) {
	chrono::system_clock::time_point end = chrono::system_clock::now();
	double elapsed =
		std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
		.count(); //�����ɗv�������Ԃ��~���b�ɕϊ�
	cout << "elapsed :" << elapsed << " ms\n";
}

void Copyfile(const char *from_file_name, const char *to_file_name) {
	// streambuf_iteratorとcopy()を利用した方法

	ifstream is(from_file_name, ios::in | ios::binary);
	ofstream os(to_file_name, ios::out | ios::binary);

	// ファイルコピー
	istreambuf_iterator<char> iit(is);
	istreambuf_iterator<char> end;
	ostreambuf_iterator<char> oit(os);
	copy(iit, end, oit);
}
inline size_t GetFileSize(const char *file) {
	int fd = 0;
	struct _stat64 stbuf;
	_sopen_s(&fd, file, O_RDONLY, _SH_DENYWR, _S_IREAD);
	if (_fstat64(fd, &stbuf) == -1)
		cout << "cant get file state.\n";
	size_t size = stbuf.st_size;
	_close(fd);
	return size;
}
bool CheckFileExistence(const string &str) {
	std::ifstream ifs(str);
	return ifs.is_open();
}
bool CheckExistenceOfFolder(const string folder_name) {
	if (_mkdir(folder_name.c_str()) == 0) {
		return true;
	}
	else {
		return false;
	}
}
bool Directory_check(string dir_name) {
	string dir_s;
	int pos1 = 0;
	int pos2 = 4;
	while (1) {
		if (dir_name.find('\\', pos2) != std::string::npos) {
			pos2 = dir_name.find('\\', pos2);
			dir_s = dir_name.substr(0, pos2);
			CheckExistenceOfFolder(dir_s);
			pos2++;
			// cout << "Created directory: " << dir_s << endl;
		}
		else {
			CheckExistenceOfFolder(dir_name);
			break;
		}
	}
	return true;
}
void COUT_VEC_String(const vector<string> &vec) {
	for (auto itr = vec.begin(); itr != vec.end(); ++itr) {
		cout << *itr << "\n";
	}
}
void COUT_VEC_Int(const vector<int> &vec) {
	for (auto itr = vec.begin(); itr != vec.end(); ++itr) {
		cout << *itr << "\n";
	}
}
void COUT_VEC_Char(const vector<char *> &vec) {
	for (auto itr = vec.begin(); itr != vec.end(); ++itr) {
		cout << *itr << "\n";
	}
}
inline vector<string> split_comma(const string &src, const char *delim = ",") {
	vector<string> vec;
	string::size_type len = src.length();
	for (string::size_type i = 0, n; i < len; i = n + 1) {
		n = src.find_first_of(delim, i);
		if (n == string::npos) {
			n = len;
		}
		vec.push_back(src.substr(i, n - i));
	}
	return vec;
}
inline vector<string> split_tab(const string &src, const char *delim = "tab") {
	vector<string> vec;
	string::size_type len = src.length();
	for (string::size_type i = 0, n; i < len; i = n + 1) {
		n = src.find_first_of(delim, i);
		if (n == string::npos) {
			n = len;
		}
		vec.push_back(src.substr(i, n - i));
	}
	return vec;
}
inline vector<int> split_comma_i(const string &src, const char *delim = ",") {
	vector<int> vec;
	string::size_type len = src.length();
	for (string::size_type i = 0, n; i < len; i = n + 1) {
		n = src.find_first_of(delim, i);
		if (n == string::npos) {
			n = len;
		}
		vec.push_back(stoi(src.substr(i, n - i)));
	}
	return vec;
}
inline vector<string> split_line_string(const string &src,
	const char *delim = "\n") {
	vector<string> vec;
	string::size_type len = src.length();

	for (string::size_type i = 0, n; i < len; i = n + 1) {
		n = src.find_first_of(delim, i);
		if (n == string::npos) {
			n = len;
		}
		vec.push_back(src.substr(i, n - i));
	}
	return vec;
}
unordered_map<string, char> Import_Codon(const char *Codon_csv) {
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	cout << "Import_Codon\n";
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

	unordered_map<string, char> codon_map;
	//const char *Codon_csv = "Y:\\codon.csv";
	FILE *codon;
	char buf[16];
	int i = 0;
	char AA = '\0';
	fopen_s(&codon, Codon_csv, "r");
	while (fgets(buf, 16, codon) != NULL) {
		//printf("%s\n", buf);
		char DNA[4] = {};
		AA = buf[0];
		DNA[0] = buf[2];
		DNA[1] = buf[3];
		DNA[2] = buf[4];
		DNA[3] = '\0';
		codon_map[DNA] = AA;
	}
	/*
	for (auto itr = codon_map.begin(); itr != codon_map.end(); ++itr) {
			cout << itr->first << '\t' << itr->second << '\n';
	}*/
	// cout << "AAA" << codon_map["AAA"] << endl;
	return codon_map;
}
float *Import_Codon_Q(const string &FP) {
	float codon_Q[21];
	FILE *stream;
	ifstream ifs(FP);
	string line;
	getline(ifs, line);
	getline(ifs, line);
	getline(ifs, line);
	vector<string> strvec = split_comma(line);
	for (int i = 1; i < strvec.size(); i++) {
		codon_Q[i - 1] = stod(strvec.at(i));
		cout << AAlist[i - 1] << "/ " << codon_Q[i - 1] << endl;
	}
	return codon_Q;
}
void DelFile(const string &fileName) { remove(fileName.c_str()); }
void DelDir(const string &dirName) { remove(dirName.c_str()); }
void Load_BLOSUM62(int list[24][24]) {
	string filepath = "H:\\BLOSUM62.csv";
	FILE *blosum62;
	char buf[256];
	int line = 0;
	vector<vector<int>> scorelist;

	fopen_s(&blosum62, filepath.c_str(), "rb");
	scorelist = vector<vector<int>>(25, vector<int>(25, '\0'));
	while (fgets(buf, 256, blosum62) != NULL) {
		scorelist[line] = split_comma_i(buf);
		for (int i = 0; i < strlen(AAlist); ++i) {
			list[line][i] = scorelist[line][i];
			//			cout << i << ";" << scorelist[line][i] << endl;
		}
		line++;
		if (line >= 25)
			break;
	}
	fclose(blosum62);
}
void Load_BLOSUM62_Limited(int list[21][21], int limit_l, int limit_h) {
	string filepath = "H:\\BLOSUM62_2.csv";
	FILE *blosum62;
	char buf[256];
	int line = 0;
	vector<vector<int>> scorelist;
	fopen_s(&blosum62, filepath.c_str(), "rb");
	scorelist = vector<vector<int>>(25, vector<int>(25, '\0'));
	cout << "-----------------------------------------------------\n\n";
	cout << "Modified BLOSUM62\n\n";
	cout << "-----------------------------------------------------\n\n";
	while (fgets(buf, 256, blosum62) != NULL) {
		scorelist[line] = split_comma_i(buf);
		for (int i = 0; i < strlen(AAlist); ++i) {
			if (scorelist[line][i] < limit_l)
				list[line][i] = limit_l;
			else if (scorelist[line][i] > limit_h)
				list[line][i] = limit_h;
			else
				list[line][i] = scorelist[line][i];
			cout << list[line][i] << ", ";
		}
		cout << endl;
		line++;
		if (line >= 25)
			break;
	}
	fclose(blosum62);
	cout << "-----------------------------------------------------\n\n";
}
void Load_ScoreTable_Limited(string filepath, int list[21][21], int limit_l,
	int limit_h) {
	FILE *blosum62;
	char buf[256];
	int line = 0;
	vector<vector<int>> scorelist;
	fopen_s(&blosum62, filepath.c_str(), "rb");
	scorelist = vector<vector<int>>(21, vector<int>(21, '\0'));
	cout << "---------------------------------------------------------------\n\n";
	cout << filepath << "\n";
	cout << "---------------------------------------------------------------\n\n";
	while (fgets(buf, 256, blosum62) != NULL) {
		scorelist[line] = split_comma_i(buf);
		for (int i = 0; i < strlen(AAlist); ++i) {
			if (i >= 20) {
				list[line][i] = limit_l;
			}
			else if (scorelist[line][i] < limit_l)
				list[line][i] = limit_l;
			else if (scorelist[line][i] > limit_h)
				list[line][i] = limit_h;
			else
				list[line][i] = scorelist[line][i];
			cout << list[line][i] << ", ";
		}
		cout << endl;
		line++;
		if (line >= 20)
			break;
	}
	for (int i = 0; i < strlen(AAlist); ++i) {
		list[20][i] = limit_l;
		cout << list[line][i] << ", ";
	}
	cout << endl;
	fclose(blosum62);
	cout << "---------------------------------------------------------------\n\n";
}
void Load_ScoreTable_Limited2(string filepath, char list[21 * 21], int limit_l, int limit_h) {
	FILE *blosum62;
	char buf[256];
	int line = 0;
	vector<vector<int>> scorelist;
	fopen_s(&blosum62, filepath.c_str(), "rb");
	scorelist = vector<vector<int>>(21, vector<int>(21, '\0'));
	cout << "---------------------------------------------------------------\n\n";
	cout << filepath << "\n";
	cout << "---------------------------------------------------------------\n\n";
	while (fgets(buf, 256, blosum62) != NULL) {
		scorelist[line] = split_comma_i(buf);
		for (int i = 0; i < strlen(AAlist); ++i) {
			if (i >= 20) {
				list[line * 21 + i] = (char)limit_l;
			}
			else if (scorelist[line][i] < limit_l)
				list[line * 21 + i] = (char)limit_l;
			else if (scorelist[line][i] > limit_h)
				list[line * 21 + i] = (char)limit_h;
			else
				list[line * 21 + i] = (char)scorelist[line][i];

		}

		line++;
		if (line >= 20)
			break;
	}
	for (int i = 0; i <21; ++i) {
		list[21 * 20 + i] = limit_l;
	}
	for (int i = 0; i < 21; ++i) {
		for (int j = 0; j < 21; ++j) {
			cout << +list[21 * i + j] << ", ";
		}
		cout << endl;
	}
	


	fclose(blosum62);
	cout << "---------------------------------------------------------------\n\n";
}
void Import_Codon_Freq_P(string filepath, double *P, string codon) {
	FILE *stream;
	char buf[256];
	int line = 0;
	int readline = 2;
	if (codon == "NNK") {
		printf_s("Codon : NNK\n");
		readline = 4;
	}
	else {
		printf_s("Codon : NNN\n");
	}
	fopen_s(&stream, filepath.c_str(), "r");
	if (&stream != 0) {
		while (fgets(buf, 256, stream) != NULL) {
			if (line == readline) {
				vector<string> data = split_comma(buf);
				for (int i = 1; i < 22; ++i) {
					P[i - 1] = stod(data[i]);
				}
			}
			line++;
		}
		fclose(stream);
	}
	else {
		cout << "File open error !!: " << filepath << endl;
	}
}
void Import_Codon_Freq_Num(string filepath, int *P, string codon) {
	FILE *stream;
	char buf[256];
	int line = 0;
	int readline = 1;
	if (codon == "NNK") {
		printf_s("Codon : NNK\n");
		readline = 3;
	}
	else {
		printf_s("Codon : NNN\n");
	}
	fopen_s(&stream, filepath.c_str(), "r");
	if (&stream != 0) {
		while (fgets(buf, 256, stream) != NULL) {
			if (line == readline) {
				vector<string> data = split_comma(buf);
				for (int i = 1; i < 22; ++i) {
					P[i - 1] = stoi(data[i]);
				}
			}
			line++;
		}
		fclose(stream);
	}
	else {
		cout << "File open error !!: " << filepath << endl;
	}
}
inline string Revers_DNA(const string &seq) {
	string seq_r(seq.size(), '*');
	for (int i = 0; i < seq.size(); ++i)
		switch (seq[i]) {
		case 'A':
			seq_r[seq.size() - i] = 'T';
			break;
		case 'G':
			seq_r[seq.size() - i] = 'C';
			break;
		case 'C':
			seq_r[seq.size() - i] = 'G';
			break;
		case 'T':
			seq_r[seq.size() - i] = 'A';
			break;
		default:
			continue;
		}
	return seq_r;
}
inline char *Get_Complement_Sequence(const char *seq) {
	size_t size = strlen(seq);
	// cout << "DNA size:" << size << endl;;
	// cout << seq << endl;
	char *seq_r = new char[size];
	for (int i = 1; i < size; ++i)
		switch (seq[size - i]) {
		case 'A':
			seq_r[i - 1] = 'T';
			break;
		case 'G':
			seq_r[i - 1] = 'C';
			break;
		case 'C':
			seq_r[i - 1] = 'G';
			break;
		case 'T':
			seq_r[i - 1] = 'A';
			break;
		default:
			seq_r[i - 1] = '*';
			continue;
		}
	seq_r[size - 1] = '\0';
	// cout << seq_r << endl;
	return seq_r;
}

inline char *AA2Chr(string &src) {
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
inline vector<char *> StoC(const vector<string> &src) {
	vector<char *> vec;

	for (auto itr = 0; itr < src.size(); ++itr) {
		vec.push_back((char *)src[itr].c_str());
		cout << src[itr] << endl;
		cout << (char *)src[itr].c_str() << endl;
	}
	return vec;
}
vector<string> Get_File_path_in_dir(const string &dir_name,
	const string &extension) noexcept(false) {
	HANDLE hFind;
	WIN32_FIND_DATA win32fd; // defined at Windwos.h
	std::vector<std::string> file_names;

	//�g���q�̐ݒ�
	std::string search_name = dir_name + "\\*." + extension;

	hFind = FindFirstFile(search_name.c_str(), &win32fd);
	try {
		if (hFind == INVALID_HANDLE_VALUE) {
			throw std::runtime_error("file not found");
		}

		do {
			if (win32fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
			}
			else {
				file_names.push_back(win32fd.cFileName);
				// printf("%s\n", file_names.back().c_str());
			}
		} while (FindNextFile(hFind, &win32fd));

		FindClose(hFind);
	}
	catch (exception &e) {
		return file_names;
	}

	return file_names;
}
vector<string> Get_Dir_path_in_dir(string dir) {
	vector<string> dirlist;
	HANDLE hFind;
	WIN32_FIND_DATA win32fd;
	string filename = dir + "\\*";
	hFind = FindFirstFile(filename.c_str(), &win32fd);
	try {
		if (hFind == INVALID_HANDLE_VALUE) {
			throw std::runtime_error("file not found");
			return dirlist;
		}
		do {
			if (win32fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
				string buf = win32fd.cFileName;
				if (buf != "." && buf != "..") {
					if (buf != "�E" && buf != "�E�E") {
						//	printf("%s (DIR)\n", win32fd.cFileName);
						dirlist.push_back(win32fd.cFileName);
					}
				}
			}
			else {
				// printf("%s\n", win32fd.cFileName);
			}
		} while (FindNextFile(hFind, &win32fd));

		FindClose(hFind);
	}
	catch (exception &e) {
		return dirlist;
	}
	return dirlist;
}
inline string ExtractPathWithoutExt(const string &fn) {
	// http://www.slis.tsukuba.ac.jp/~fujisawa.makoto.fu/cgi-bin/wiki/index.php?string%A4%CB%A4%E8%A4%EB%CA%B8%BB%FA%CE%F3%BD%E8%CD%FD#accc1053
	string::size_type pos;
	if ((pos = fn.find_last_of(".")) == string::npos) {
		return fn;
	}

	return fn.substr(0, pos);
}
inline char *Importfile(string FP) {
	FILE *fpr;
	size_t size = GetFileSize(FP.c_str());
	// cout << FP << ": " << size / 1000000 << " Mbyte" << endl;
	char *buf = new (char[size + 1]);
	fopen_s(&fpr, FP.c_str(), "r");
	size_t rc = fread(buf, sizeof(char), size, fpr);
	if (rc == 0) {
		printf("File import error\n");
	}
	else {
		buf[rc] = '\0';
	}
	buf[size] = '\0';
	fclose(fpr);
	return buf;
}

inline vector<char *> Importfile_old(string FP) {
	FILE *fpr;
	vector<char *> vec;
	// mtx.lock();
	fopen_s(&fpr, FP.c_str(), "r");
	char str[256];
	if (fpr == NULL) {
		printf("%s file not open!\n", FP.c_str());
	}
	else {
		printf("%s file opened!\n", FP.c_str());
		while (fgets(str, 256, fpr) != NULL) {
			char *sp;
			sp = strchr(str, '\0');
			if (sp == NULL) {
				sp = strchr(str, '\n');
				if (sp == NULL) {
					sp = strchr(str, EOF);
				}
			}
			if (sp != NULL) {
				int len = (int)(sp - str);
				char *read = new char[len];
				strncpy_s(read, len, str, len - 1);
				vec.push_back(read);
			}

			// printf("%s", str);
		}
		fclose(fpr);
		mtx.unlock();
		return vec;
	}
	fclose(fpr);
	// mtx.unlock();
}

inline vector<Fasta> Import_Fasta(string dir, vector<string> list) {
	vector<Fasta> data;
	for (int i = 0; i < list.size(); ++i) {
		FILE *fp;
		char buf[2048] = { '\0' };
		string FP = dir + "\\" + list[i];
		fopen_s(&fp, FP.c_str(), "r");
		if (fp == NULL) {
			return data;
		}
		string seq_t;
		string name_t;
		int firstw = 0;

		while (fgets(buf, 2048, fp) != NULL) {
			if (buf[0] == '>') {
				// cout << buf << endl;
				if (firstw != 0) {
					Fasta *temp = new Fasta;
					temp->name = name_t;
					temp->seq = seq_t;
					data.push_back(*temp);
					delete temp;
					seq_t = "";
					name_t = buf;
					name_t.pop_back();
				}
				else {
					name_t = buf;
					name_t.pop_back();
					firstw = 1;
				}
			}
			else {
				seq_t = seq_t + buf;
				seq_t.pop_back();
				// cout << seq_t << endl;
			}
		}
		Fasta *temp = new (Fasta);
		temp->name = name_t;
		temp->seq = seq_t;
		data.push_back(*temp);
		delete temp;
		fclose(fp);
	}
	return data;
}
inline void Export_Fasta(string dir, vector<Fasta> &list) {
	FILE *fpw;
	fopen_s(&fpw, dir.c_str(), "w");
	for (auto n = 0; n < list.size(); ++n) {

		if (fpw == NULL) {
			printf("%s file not open!\n");
		}
		fprintf_s(fpw, "%s\n%s\n", (list[n].name).c_str(), (list[n].seq).c_str());
	}
	fclose(fpw);
}
inline void Export_char(const char *src, string dir, string filename) {
	FILE *fpw;
	Directory_check(dir);
	fopen_s(&fpw, (dir + "\\" + filename).c_str(), "a");
	if (fpw == NULL) {
		printf("%s file not open!\n", dir.c_str());
	}
	for (auto n = 0; n < strlen(src); ++n) {
		putc(src[n], fpw);
		putc(',', fpw);
	}
	fclose(fpw);
}
bool removeDirectory(string fileName) {
	bool retVal = true;
	string nextFileName;

	WIN32_FIND_DATA foundFile;

	HANDLE hFile = FindFirstFile((fileName + "\\*.*").c_str(), &foundFile);

	if (hFile != INVALID_HANDLE_VALUE) {
		do {
			// If a found file is . or .. then skip
			if (strcmp(foundFile.cFileName, ".") != 0 &&
				strcmp(foundFile.cFileName, "..") != 0) {
				// The path should be absolute path
				nextFileName = fileName + "\\" + foundFile.cFileName;

				// If the file is a directory
				if ((foundFile.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0) {
					removeDirectory(nextFileName.c_str());
					RemoveDirectory(nextFileName.c_str());
				}
				// If the file is a file
				else {
					DeleteFile(nextFileName.c_str());
				}
			}
		} while (FindNextFile(hFile, &foundFile) != 0);
	}

	FindClose(hFile);

	// Delete starting point itseft
	if (RemoveDirectory(fileName.c_str()) == 0)
		retVal = false;

	return retVal;
}
