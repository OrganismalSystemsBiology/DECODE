#include <DECODE.h>

#define FSEQ "GTTAACTTTAAGAAGG"
#define RSEQ "CGCTGCCGCTGCCGCA" //IonProton : Rseq
#define NOBC "CCTAATACGACTCA"
#define RNOBC "TGAGTCGTATTAGG"

struct seqdata
{
	string seq;
	size_t read;
};
bool sort_read(const seqdata &a, const seqdata &b)
{
	return a.read > b.read;
}
int CPU_Cores = 1;
//int Thread_max = 10;
int Thread_max = 1;
int IonProton_flag = 0;
size_t read_lines = 1000 * 1000;
size_t task = 0;
size_t fsize = 0;
fpos_t gfpos;
size_t Seq_Combine(int X, vector<string> filelist, vector<string> flist, string FPw, string FPr)
{
	size_t rread = 0;
	unordered_map<string, int> Result_map;
	char *data = Importfile(FPr + "/" + flist[X]);
	char *ps = data;
	char *pe = data;
	while (1)
	{
		pe = strchr(ps, ',');
		if (pe == NULL)
		{
			break;
		}
		size_t size = (uintptr_t)pe - (uintptr_t)ps;
		char *seq = new char[size + 1];
		for (int i = 0; i < size; ++i)
		{
			seq[i] = ps[i];
		}
		seq[size] = '\0';
		ps = strchr(pe, '\n');
		if (ps == NULL)
		{
			break;
		}
		pe++;
		size = (uintptr_t)ps - (uintptr_t)pe;
		char *read = new char[size + 1];
		for (int i = 0; i < size; ++i)
		{
			read[i] = pe[i];
		}
		read[size] = '\0';
		//cout << seq << " : " << read << endl;
		if (Result_map.count(seq) == 0)
		{
			Result_map[seq] = stoi(read);
		}
		else
		{
			Result_map[seq] += stoi(read);
		}
		rread += stoi(read);
		delete[] read;
		delete[] seq;
		ps++;
	}

	if (Result_map.size() > 0)
	{
		vector<Seq_read> results;
		for (auto itr = Result_map.begin(); itr != Result_map.end(); ++itr)
		{
			Seq_read *temp = new (Seq_read);
			temp->seq = itr->first;
			temp->read = itr->second;
			results.push_back(*temp);
			delete temp;
		}
		if (results.size() > 0)
		{
			sort(results.begin(), results.end());
			string to_file_name = FPw + '/' + flist[X];
			FILE *fpw;
			fopen_s(&fpw, to_file_name.c_str(), "w");
			setvbuf(fpw, NULL, _IOFBF, 1024 * 1024);
			for (size_t ite = 0; ite < results.size(); ++ite)
			{
				fprintf(fpw, "%s,%d\n", (results[ite].seq).c_str(), results[ite].read);
			}
			fclose(fpw);
			results.erase(results.begin());
			results.shrink_to_fit();
			assert(results.capacity() == results.size());
			vector<Seq_read>().swap(results);
		}
	}
	delete[] data;
	Result_map.clear();
	return rread;
}

inline int Seq_convert_Pep_char(const char *ps, const char *pe, char *pep, unordered_map<string, char> &codon_map, int FR)
{
	char codon[4] = {};
	size_t size = (uintptr_t)pe - (uintptr_t)ps;
	char *DNA = new char[size + 1];
	for (int i = 0; i < size; ++i)
	{
		DNA[i] = ps[i];
	}
	DNA[size] = '\0';
	//cout << "seq; " << DNA << endl;
	char *p;
	int len_beforATG = strlen(startseq);
	if (FR == 0)
	{
		p = strstr(DNA, startseq);
		if (p != NULL)
		{
			int length = strlen(p);
			if (length > (CDS_LENGTH * 3) + len_beforATG)
			{
				for (int i = 0; i < CDS_LENGTH; ++i)
				{
					try
					{
						codon[0] = p[3 * i + len_beforATG];
						codon[1] = p[3 * i + len_beforATG + 1];
						codon[2] = p[3 * i + len_beforATG + 2];
						codon[3] = '\0';
						if (codon[0] != '\0' && codon[1] != '\0' && codon[2] != '\0')
						{
							pep[i] = codon_map[codon];
						}
					}
					catch (exception &e)
					{
						pep[0] = '\0';
						delete[] DNA;
						return 0;
					}
				}
				pep[CDS_LENGTH] = '\0';
				delete[] DNA;
				return 1;
			}
		}
	}
	else if (FR == 1)
	{
		char *seqr = Get_Complement_Sequence(DNA);
		p = strstr(seqr, startseq);
		if (p != NULL)
		{
			int length = strlen(p);
			if (length > (CDS_LENGTH * 3) + len_beforATG)
			{
				for (int i = 0; i < CDS_LENGTH; ++i)
				{
					try
					{
						codon[0] = p[3 * i + len_beforATG];
						codon[1] = p[3 * i + len_beforATG + 1];
						codon[2] = p[3 * i + len_beforATG + 2];
						codon[3] = '\0';
						if (codon[0] != '\0' && codon[1] != '\0' && codon[2] != '\0')
						{
							pep[i] = codon_map[codon];
						}
					}
					catch (exception &e)
					{
						pep[0] = {};
						delete[] DNA;
						delete[] seqr;
						return 0;
					}
				}
			}
			pep[CDS_LENGTH] = '\0';
			delete[] DNA;
			delete[] seqr;
			return 1;
		}
		delete[] seqr;
	}
	delete[] DNA;
	return 0;
}
unordered_map<string, string> Import_Barcode(const string &FP)
{
	unordered_map<string, string> barcode_map;

	char buf[256];
	char *name;
	char *seq;
	char *next_token = NULL;
	FILE *fp;
	fopen_s(&fp, FP.c_str(), "r");
	while (fgets(buf, 256, fp) != NULL)
	{
		cout << buf;
		name = strtok_s(buf, ",", &next_token);
		seq = strtok_s(NULL, ",", &next_token);
		string s = seq;
		s.pop_back(); //remove \r
		s.pop_back(); //remove \n
		//size_t strLen = s.length();
		//string sb = s.substr(strLen - 6, 6);
		//barcode_map[sb] = name;
		barcode_map[s] = name;
	}

	return barcode_map;
}

inline string Beacode_check(char *ps, char *pe, unordered_map<string, string> &barcode_F,
	unordered_map<string, string> &barcode_R, int &FR)
{

	size_t size = (uintptr_t)pe - (uintptr_t)ps;
	char *DNA = new char[size + 1];
	for (int i = 0; i < size + 1; ++i)
	{
		DNA[i] = '\0';
	}
	char *p;
	char *c;
	for (int i = 0; i < size; ++i)
	{
		DNA[i] = ps[i];
	}
	//cout << "O: " << DNA << endl;
	p = strstr(DNA, NOBC);
	if (p != 0)
	{
		delete[] DNA;
		FR = 9;
		return "0";
	}
	p = strstr(DNA, RNOBC);
	if (p != 0)
	{
		delete[] DNA;
		FR = 9;
		return "0";
	}
	p = strstr(DNA, FSEQ);
	if (p != 0)
	{
		FR = 0;
	}
	else
	{
		p = strstr(DNA, RSEQ);
		{
			if (p != 0)
				FR = 1;
		}
	}
	for (auto ite = barcode_F.begin(); ite != barcode_F.end(); ++ite)
	{
		//cout << DNA << ":" << ite->first << endl;
		//string s = DNA;
		//if (s.find(ite->first) != string::npos){
		if (strstr(DNA, ite->first.c_str()) != NULL)
		{
			delete[] DNA;
			return ite->second;
		}
	}
	for (auto ite = barcode_R.begin(); ite != barcode_R.end(); ++ite)
	{
		//cout << DNA << ":" << ite->first << endl;
		//string s = DNA;
		//if (s.find(ite->first) != string::npos){
		if (strstr(DNA, ite->first.c_str()) != NULL)
		{
			delete[] DNA;
			return ite->second;
		}
	}
	delete[] DNA;
	return "0";
}

inline char *ImportFastq(string FP)
{
	cout << "ImportFastq" << endl;
	if (gfpos >= (fsize - 1))
	{
		char *b = new char[1];
		b[0] = '\0';
		return b;
	}

	size_t length = 4 * read_lines;
	fpos_t startpos = gfpos;
	FILE *fpr;
	fopen_s(&fpr, FP.c_str(), "r");
	//cout << "gfpos >= fsize - 2" << gfpos << " / " << fsize << endl;
	_fseeki64(fpr, startpos, SEEK_SET);
	char bufs[256 * 10];
	for (int i = 0; i < 4 * length; ++i)
	{
		if (fgets(bufs, 255 * 10, fpr) == NULL)
			break;
	}
	fgetpos(fpr, &gfpos);
	char *buf = new char[gfpos - startpos + 1];
	_fseeki64(fpr, startpos, SEEK_SET);
	size_t rc = fread(buf, sizeof(char), gfpos - startpos, fpr);
	buf[gfpos - startpos] = '\0';
	fclose(fpr);
	return buf;
}

void Export_map(unordered_map<string, unordered_map<string, int>> &Result, string fp)
{
	cout << "Export: " << fp << endl;
	for (auto itr = Result.begin(); itr != Result.end(); ++itr)
	{
		string fp2 = fp + "\\" + itr->first + ".csv";
		Directory_check(fp);
		FILE *fpw;
		fopen_s(&fpw, fp2.c_str(), "a");
		setvbuf(fpw, NULL, _IOFBF, 1024 * 1024);
		for (auto ite = itr->second.begin(); ite != itr->second.end(); ++ite)
		{
			fprintf(fpw, "%s,%d\n", (ite->first).c_str(), ite->second);
		}
		fclose(fpw);
	}
}
void Export_map2(vector<seqdata> &Result, string fp)
{
	//Directory_check(fp);
	cout << fp << endl;
	FILE *fpw;
	fopen_s(&fpw, fp.c_str(), "w+");
	for (auto itr = 0; itr < Result.size(); ++itr)
	{
		//cout << itr << "/" << Result.size() << "   ";
		//printf("%s,%ld\n", Result[itr].seq.c_str(), Result[itr].read);
		fprintf(fpw, "%s,%ld\n", Result[itr].seq.c_str(), Result[itr].read);
	}
	fclose(fpw);
}
inline void Seq2Pep(int X, string file, string &FP, string FPw, unordered_map<string, char> codon_map,
	unordered_map<string, string> barcode_F, unordered_map<string, string> barcode_R)
{
	while (1)
	{
		unordered_map<string, unordered_map<string, int>> Result_map;
		int i_ret = 0;
		char pep[CDS_LENGTH + 1] = {};
		mtx.lock();
		if (gfpos >= fsize - 2)
		{
			cout <<"gfpos >= fsize - 2" <<gfpos << " / " << fsize << endl;
			Result_map.clear();
			mtx.unlock();
			break;
		}
		char *data = ImportFastq(FP + "\\" + file);
		mtx.unlock();
		cout << "Import: " << ((float)gfpos) / 1000000000 << " / " << ((float)fsize) / 1000000000 << " Gbyte" << endl;
		char *ps = data;
		char *pe = data;
		ps = strchr(data, '\n');
		++ps;
		while (1)
		{
			if (ps == NULL) {
				break;
			}
			pe = strstr(ps, "\n+");
			if (pe == NULL) {
				break;
			}
			i_ret = 0;
			int FR = 0;
			string barcodename = Beacode_check(ps, pe, barcode_F, barcode_R, FR);
			//cout << barcodename << endl;
			if (barcodename != "0")
			{
				i_ret = Seq_convert_Pep_char(ps, pe, pep, codon_map, FR);
				if (i_ret == 1)
				{
					if (Result_map.count(barcodename) == 0)
					{
						Result_map[barcodename][pep] = 1;
					}
					else
					{
						if (Result_map[barcodename].count(pep) == 0)
						{
							Result_map[barcodename][pep] = 1;
						}
						else
						{
							Result_map[barcodename][pep] += 1;
						}
					}
				}
			}
			ps = strstr(pe, "\n@");
			if (ps == NULL) {
				cout << "ps == NULL" << endl;
				break;
			}
			++ps;
			ps = strchr(ps, '\n');
			if (ps == NULL) {
				cout << "ps == NULL" << endl;
				break;
			}
			++ps;
		}
		delete[] data;
		mtx2.lock();
		if (Result_map.size() > 0)
		{
			Export_map(Result_map, FPw);
		}
		mtx2.unlock();
		Result_map.clear();
		if (gfpos >= fsize - 2)
		{
			break;
		}
	}
	cout << "Seq2Pep end" << endl;
}

//for IonProton
inline void Seq2Pep2(string FP, vector<string> flist, string FPw, unordered_map<string, char> codon_map)
{
	int t = 0;
	unordered_map<string, int> Result_map;
	while (1)
	{
		mtx.lock();
		t = task;
		task++;
		mtx.unlock();
		if (t >= flist.size())
			break;
		char *data = Importfile(FP + "\\" + flist[t]);
		int i_ret = 0;
		Result_map.clear();
		char pep[CDS_LENGTH + 1] = {};
		char *ps = data;
		char *pe = data;
		ps = strchr(data, '\n');
		++ps;
		while (1)
		{
			if (ps == NULL)
				break;
			pe = strstr(ps, "\n+");
			if (pe == NULL)
				break;
			i_ret = 0;

			i_ret = Seq_convert_Pep_char(ps, pe, pep, codon_map, 1);
			if (i_ret == 1)
			{
				//cout << pep << endl;
				if (Result_map.count(pep) == 0)
				{
					Result_map[pep] = 1;
				}
				else
				{
					Result_map[pep] += 1;
				}
			}

			ps = strstr(pe, "\n@");
			if (ps == NULL)
				break;
			++ps;
			ps = strchr(ps, '\n');
			if (ps == NULL)
				break;
			++ps;
		}

		vector<seqdata> tmp;
		for (auto ite = Result_map.begin(); ite != Result_map.end(); ++ite)
		{
			seqdata buf;
			buf.seq = ite->first;
			buf.read = ite->second;
			tmp.push_back(buf);
		}
		sort(tmp.begin(), tmp.end(), sort_read);
		string fp_ext = ExtractPathWithoutExt(flist[t]);
		Export_map2(tmp, FPw + "\\" + fp_ext + ".csv");
		delete[] data;
		Result_map.clear();
	}
}
int main(int argc, char *argv[])
{
	cout << "argc:" << argc << endl;
	string FP = argv[1];
	string FPwd = argv[2];
	string Codon_FP = argv[3];
	cout << "Input: " << FP << endl;
	cout << "Export dir: " << FPwd << endl;
	//cout << "Codon file : " << Codon_FP << endl;
	Directory_check(FPwd);
	string barcode_dir = "";
	if (argc > 4)
	{
		barcode_dir = argv[4];
		cout << "Barcode: " << barcode_dir << endl;
		IonProton_flag = 0;
	}
	else
	{
		IonProton_flag = 1;
	}

	CPU_Cores = thread::hardware_concurrency() - 1;
	cout << " CPU: " << CPU_Cores << " cores\n\n";
	unordered_map<string, string> barcode_F;
	unordered_map<string, string> barcode_R;
	vector<string> barcodes;
	if (barcode_dir != "")
	{
		barcodes = Get_File_path_in_dir(barcode_dir, "csv");
		cout << barcodes[0] << endl;
		barcode_F = Import_Barcode(barcode_dir + "\\" + barcodes[0]);
		cout << barcodes[1] << endl;
		barcode_R = Import_Barcode(barcode_dir + "\\" + barcodes[1]);
		//for (auto i = barcode_F.begin(); i != barcode_F.end(); ++i)
		//	{
			//	cout << i->first << endl;
				//cout << i->first << ":" << i->second << endl
		//	}
	}
	unordered_map<string, char> codon_map = Import_Codon(Codon_FP.c_str());
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	cout << "Fastq to Peptide\n";
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	vector<string> filelist = Get_File_path_in_dir(FP, "fastq");
	sort(filelist.begin(), filelist.end());
	COUT_VEC_String(filelist);
	task = 0;
	if (IonProton_flag == 0)
	{
		for (int j = 0; j < filelist.size(); ++j)
		{
			gfpos = 0;
			fsize = GetFileSize((FP + "\\" + filelist[j]).c_str());
			cout << "File size: " << fsize << endl;
			string fp_ext = ExtractPathWithoutExt(filelist[j]);
			string FPw = FPwd + "\\" + fp_ext + "\\Peptide";
			string FPwt = FPwd + "\\" + fp_ext + "\\Peptide_temp";
			Directory_check(FPw);
			Directory_check(FPwt);
			gfpos = 0;
			//Thread_max = 1;
			if (CPU_Cores > Thread_max)
			{
				cout << " Start " << Thread_max << " threads\n\n";
				vector<thread> threads(Thread_max);
				for (int i = 0; i < Thread_max; ++i)
				{
					threads[i] = thread([&] { Seq2Pep(i, filelist[j], FP, FPwt, codon_map, barcode_F, barcode_R); });
				}
				for (auto &t : threads)
				{
					t.join();
				}
			}
			else if (CPU_Cores <= Thread_max)
			{
				cout << " Start " << CPU_Cores << " threads<<endl\n\n";
				vector<thread> threads(CPU_Cores);
				for (int i = 0; i < CPU_Cores; ++i)
				{
					threads[i] = thread([&] { Seq2Pep(i, filelist[j], FP, FPwt, codon_map, barcode_F, barcode_R); });
				}
				for (auto &t : threads)
				{
					t.join();
				}
			}
			vector<string> flist = Get_File_path_in_dir(FPwt, "csv");
			for (int n = 0; n < flist.size(); ++n)
			{
				int read = Seq_Combine(n, filelist, flist, FPw, FPwt);

				cout << flist[n] << "/" << read << endl;
				FILE *fpw2;
				fopen_s(&fpw2, (FPw + '\\' + "read_list.txt").c_str(), "a");
				string fname = flist[n];
				fprintf(fpw2, "%s,%d\n", fname.c_str(), read);
				fclose(fpw2);
			}
			for (int n = 0; n < flist.size(); ++n)
			{
				DelFile(FPwt + "\\" + flist[n]);
			}
			//DelDir(FPwt);
			task++;
		}
	}
	else
	{
		vector<thread> threads(CPU_Cores);
		for (int i = 0; i < CPU_Cores; ++i)
		{
			threads[i] = thread([&] { Seq2Pep2(FP, filelist, FPwd, codon_map); });
		}
		for (auto &t : threads)
		{
			t.join();
		}
	}
	//Step1();

	cout << "************************************" << endl
		<< endl;
	cout << "Fastq to Peptide Finish!" << endl;
	cout << "************************************" << endl
		<< endl;

	return 0;
}