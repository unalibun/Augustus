// standard C/C++ includes
#include <fstream>
#include <vector>
#include <iostream>
#include <map>

using namespace std;

#define FS_LEN          3
#define STOPCODON_LEN   3

#define SINGLE          0
#define TERMINAL        1
#define INITIAL         2
#define INTERNAL1       3
#define INTERNAL2       4

struct GenePrediction
{
    int start, stop;
    vector<pair<int,int>> intronArr, exonArr, CDSArr;
    vector<int> candidateFS, viableFS;
};

#include <sstream>
template <typename T>
std::string ToString(const T& val)
{
    ostringstream oss;
    oss << val;
    return oss.str();
}

void InsertArtificialFs(string& seq, string& slippery, int etype, GenePrediction& pred, string id)
{
    slippery = "";
    int len = slippery.length();
    bool reject;
    
    // single exons only so far 
    if(etype == SINGLE)
    {
        if(pred.intronArr.size() == 0 && pred.exonArr.size() == 1)
        {
            int i=pred.stop-1-3-3, j;
        
            while(i>pred.start+3)
            {
                // slow code : for producing tests only
                
                reject = false;
            
                if((seq[i] == 'a' || seq[i] == 'A') && (seq[i+1] == 'a' || seq[i+1] == 'A'))
                {
                    // make sure not within any intron
                    for(j=0;j<pred.intronArr.size() && !reject;++j)
                    {
                        if(i>=pred.intronArr[j].first-3 && i<=pred.intronArr[j].second)
                            reject = true;
                    }

                    if(len>0)
                    {
                        for(j=1;j<=len && !reject;++j)
                        {
                            if(i-j>=pred.start+3)
                            {
                                if(seq[i-j] != slippery[len-j])
                                    reject = true;
                            }
                            else
                                reject = true;
                        }                
                    }            
                }
                else
                    reject = true;
                    
                if(!reject)
                    pred.viableFS.push_back(i);

                i-=3;        
            
            }

        }
    }
    else if(etype == TERMINAL)
    {
        // terminal exons
        if(pred.exonArr.size()>1 && pred.intronArr.size() + 1 == pred.exonArr.size() )
        {
            int i=pred.stop-1-3-3, j;
        
            while(i>pred.intronArr.back().second+3)
            {
                // slow code : for testing only
                reject = false;
            
                if((seq[i] == 'a' || seq[i] == 'A') && (seq[i+1] == 'a' || seq[i+1] == 'A'))
                {
                    // make sure not within any intron
                    for(j=0;j<pred.intronArr.size() && !reject;++j)
                    {
                        if(i>=pred.intronArr[j].first-3 && i<=pred.intronArr[j].second)
                            reject = true;
                    }

                    if(len>0)
                    {
                        for(j=1;j<=len && !reject;++j)
                        {
                            if(i-j>=pred.start+3)
                            {
                                if(seq[i-j] != slippery[len-j])
                                    reject = true;
                            }
                            else
                                reject = true;
                        }                
                    }            
                }
                else
                    reject = true;
                    
                if(!reject)
                    pred.viableFS.push_back(i);
    
                i-=3;        
            }            
        }
    }
    else if(etype == INITIAL)
    {
        // initial exons
        if(pred.exonArr.size()>1 && pred.intronArr.size() + 1 == pred.exonArr.size() )
        {
            int i=pred.start-1+3+3, j;
            
            while(i<pred.intronArr.front().first-3)
            {
                // slow code : for testing only
                reject = false;
            
                if((seq[i] == 'a' || seq[i] == 'A') && (seq[i+1] == 'a' || seq[i+1] == 'A'))
                {
                    // make sure not within any intron
                    for(j=0;j<pred.intronArr.size() && !reject;++j)
                    {
                        if(i>=pred.intronArr[j].first-3 && i<=pred.intronArr[j].second)
                            reject = true;
                    }

                    if(len>0)
                    {
                        for(j=1;j<=len && !reject;++j)
                        {
                            if(i-j>=pred.start+3)
                            {
                                if(seq[i-j] != slippery[len-j])
                                    reject = true;
                            }
                            else
                                reject = true;
                        }                
                    }            
                }
                else
                    reject = true;
                    
                if(!reject)
                    pred.viableFS.push_back(i);
    
                i+=3;        
            }            
        }
    }
    else if(etype == INTERNAL1)
    {
        // internal exons
        if(pred.exonArr.size()>2 && pred.intronArr.size() + 1 == pred.exonArr.size() )
        {
            int i = 3 + pred.exonArr[1].first - 1 + (3-((pred.exonArr[0].second - pred.start+1) % 3)), j;
            
            while(i<pred.intronArr[1].first-3)
            {
                // slow code : for testing only
                reject = false;
            
                if((seq[i] == 'a' || seq[i] == 'A') && (seq[i+1] == 'a' || seq[i+1] == 'A'))
                {
                    // make sure not within any intron
                    for(j=0;j<pred.intronArr.size() && !reject;++j)
                    {
                        if(i>=pred.intronArr[j].first-3 && i<=pred.intronArr[j].second)
                            reject = true;
                    }

                    if(len>0)
                    {
                        for(j=1;j<=len && !reject;++j)
                        {
                            if(i-j>=pred.start+3)
                            {
                                if(seq[i-j] != slippery[len-j])
                                    reject = true;
                            }
                            else
                                reject = true;
                        }                
                    }            
                }
                else
                    reject = true;
                    
                if(!reject)
                    pred.viableFS.push_back(i);
    
                i+=3;        
            }            
        }
    }
    else if(etype == INTERNAL2)
    {
        // internal exons
        if(pred.exonArr.size()>3 && pred.intronArr.size() + 1 == pred.exonArr.size())
        {
            int i = 3 + 
            pred.exonArr[2].first - 1 + 
            (3-
            (pred.exonArr[0].second - pred.start+1 + pred.exonArr[1].second - pred.exonArr[1].first+1) % 3), j;

            while(i<pred.intronArr[2].first-3)
            {
                // slow code : for testing only
                reject = false;
            
                if((seq[i] == 'a' || seq[i] == 'A') && (seq[i+1] == 'a' || seq[i+1] == 'A'))
                {
                    // make sure not within any intron
                    for(j=0;j<pred.intronArr.size() && !reject;++j)
                    {
                        if(i>=pred.intronArr[j].first-3 && i<=pred.intronArr[j].second)
                            reject = true;
                    }

                    if(len>0)
                    {
                        for(j=1;j<=len && !reject;++j)
                        {
                            if(i-j>=pred.start+3)
                            {
                                if(seq[i-j] != slippery[len-j])
                                    reject = true;
                            }
                            else
                                reject = true;
                        }                
                    }            
                }
                else
                    reject = true;
                    
                if(!reject)
                    pred.viableFS.push_back(i);
    
                i+=3;        
            }            
        }
    }
 }

void ApplyFs(string filenameFasta, string slippery, int etype, map<string, vector<GenePrediction>>& predMap)
{
    string tmp, key, seq, filenameOut;
    bool stop;

    ifstream infile;
    infile.exceptions(ifstream::eofbit | istream::failbit | ifstream::badbit);
    
    try
    {
        infile.open(filenameFasta);
        infile >> key;
        key.erase(0, 1);

        while(infile.eof() == false)
        {
            stop = false;

            do
            {
                infile >> tmp;

                if(tmp.find(">") == string::npos)
                    seq += tmp;
                else
                    stop = true;
            }
            while(!infile.eof() && !stop);

            if(predMap.count(key) > 0)
                for(int i=0;i<predMap[key].size();++i)
                    InsertArtificialFs(seq, slippery, etype, predMap[key][i], key);

            seq.clear();

            if(!infile.eof())
            {
                key = tmp;
                key.erase(0, 1);
            }
        }

        infile.close();
    }
    catch(ifstream::failure e)
	{
        if(infile.eof()) // process the last unprocessed sequence
		{
            for(int i=0;i<predMap[key].size();++i)
                InsertArtificialFs(seq, slippery, etype, predMap[key][i], key);

            seq.clear();
        }
        
        if(infile.is_open())
			infile.close();
	}
}

enum Strand {STRAND_UNKNOWN, plusstrand, minusstrand, bothstrands};

Strand StrandChar(char chr)
{
    switch(chr)
    {
        case '+':
            return plusstrand;
        case '-':
            return minusstrand;
        default:
            return STRAND_UNKNOWN;
    }
}

void ReadPrediction(string filename, Strand strand, map<string, vector<GenePrediction>>& predMap)
{
    vector<GenePrediction> predArr;
    pair<map<string, vector<GenePrediction>>::iterator,bool> ret;

    ifstream infile;
    infile.exceptions(ifstream::eofbit | istream::failbit | ifstream::badbit);
    string tmp, key;	
    GenePrediction pred;
    bool stop;
    int start, end;
	
    try
    {
        infile.open(filename);
        
        while(!infile.eof())
        {
            infile >> tmp;
            
            if(tmp == "start")
            {

                infile >> tmp;
                if(tmp == "gene")
                {

                    infile >> tmp; // id gXXX expected
                    infile >> key; // id of sequence expected
                    
                    pred.start = pred.stop = -1;
                    pred.intronArr.clear();
                    pred.exonArr.clear();
                    pred.CDSArr.clear();

                    stop = false;
                    
                    while(tmp != "end" && !stop)
                    {
                        if(tmp == "intron")
                        {
                            // we do suppose introns to have been already sorted
                            infile >> start;
                            infile >> end;
                            pred.intronArr.push_back(make_pair(start, end));

                            infile >> tmp;  // skip score 
                            infile >> tmp;  // read strand

                            //if(StrandChar(tmp[0]) != strand)
                            if(tmp.find("+") == string::npos)
                                stop = true;    // overkill
                             
                        }
                        else if(tmp == "exon")
                        {
                            // we do suppose introns to have already sorted
                            infile >> start;
                            infile >> end;
                            pred.exonArr.push_back(make_pair(start, end));

                            infile >> tmp;  // skip score 
                            infile >> tmp;  // read strand

                            //if(StrandChar(tmp[0]) != strand)
                            if(tmp.find("+") == string::npos)
                                stop = true;    // overkill
                        }
                        else if(tmp == "CDS")
                        {
                            // we do suppose introns to have already sorted
                            infile >> start;
                            infile >> end;
                            pred.CDSArr.push_back(make_pair(start, end));

                            infile >> tmp;  // skip score 
                            infile >> tmp;  // read strand

                            //if(StrandChar(tmp[0]) != strand)
                            if(tmp.find("+") == string::npos)
                                stop = true;    // overkill
                          
                        }
                        else if(tmp == "exon")
                        {
                            // we do suppose introns to have already sorted
                            infile >> start;
                            infile >> end;
                            pred.exonArr.push_back(make_pair(start, end));

                            infile >> tmp;  // skip score 
                            infile >> tmp;  // read strand
                            
                            //if(StrandChar(tmp[0]) != strand)
                            if(tmp.find("+") == string::npos)
                                stop = true;    // overkill
                        }
                        else if(tmp == "start_codon")
                        {
                            infile >> pred.start;
                            infile >> tmp;  // skip end

                            infile >> tmp;  // skip score 
                            infile >> tmp;  // read strand

                            //if(StrandChar(tmp[0]) != strand)
                            if(tmp.find("+") == string::npos)
                                stop = true;    // overkill
                           
                        }
                        else if(tmp == "stop_codon")
                        {
                            infile >> pred.stop;
                            infile >> tmp;  // skip end
                            
                            infile >> tmp;  // skip score 
                            infile >> tmp;  // read strand

                            //if(StrandChar(tmp[0]) != strand)
                            if(tmp.find("+") == string::npos)
                                stop = true;    // overkill
                          
                        }
                        
                        infile >> tmp;                        
                    }       

                    if(!stop && pred.start>-1 && pred.stop>-1)
                    {
                        if(predMap.count(key) == 0)
                            predMap.insert(make_pair(key, predArr));
                        predMap[key].push_back(pred);

                    }
                    else
                    {
                        while(tmp != "end" )
                        {
                            infile >> tmp;
                        }                    
                    }
                }                
            }                
        }
        
        infile.close();
    }
    catch(ifstream::failure e)
	{
		if(infile.is_open())
			infile.close();
	}
}

void WritePrediction(int etype, map<string, vector<GenePrediction>>& predMap)
{
    ofstream outfile;
    outfile.exceptions(ofstream::eofbit | ostream::failbit | ofstream::badbit);
    string filename;

    try
    {
        map<string, vector<GenePrediction>>::iterator mapIt;        
        
        for(mapIt=predMap.begin();mapIt!=predMap.end();++mapIt)
        {
            filename = mapIt->first + "_plus";
            for(int i=0;i<mapIt->second.size();++i)
            {
                filename += "_g";
                filename += ToString(i);
                filename += "_i";
                filename += ToString(mapIt->second[i].intronArr.size());
            }
            filename += ".info";
            
            outfile.open(filename);

            for(int j, i=0;i<mapIt->second.size();++i)
            {
                if(!mapIt->second[i].viableFS.empty() && mapIt->second[i].intronArr.size() == mapIt->second[i].exonArr.size()-1)
                {
                    outfile << "GENE " << i << endl;

                    outfile << "startAt " << mapIt->second[i].start << " stopAt " << mapIt->second[i].stop << endl;
                    for(j=0;j<mapIt->second[i].intronArr.size();++j)
                        outfile << "intron["<< j << "]: " << mapIt->second[i].intronArr[j].first << ", " << mapIt->second[i].intronArr[j].second << endl;
                    for(j=0;j<mapIt->second[i].exonArr.size();++j)
                        outfile << "exon["<< j << "]: " << mapIt->second[i].exonArr[j].first << ", " << mapIt->second[i].exonArr[j].second << endl;
            
                    outfile << "viable fs: " << endl;

                    for(j=0;j<mapIt->second[i].viableFS.size();++j)
                    {
                        outfile << mapIt->second[i].viableFS[j] /*+j for multiple fs*/  << endl;
                    }
                
                    outfile << endl;
                }
            }
            outfile.close();
        }
    }
    catch(ifstream::failure e)
	{
		if(outfile.is_open())
			outfile.close();
	}
}

// pass by value and return seq : we do not care about inefficiency
string FS(string& seq, int fsAt)
{
    string res = seq;
    if(fsAt>-1 && fsAt<res.size())
        res.insert(res.begin() + fsAt, 't');
    
    return res;
}

void ReadFasta(string filename, string& seq)
{
    ifstream infile;
    infile.exceptions(ifstream::eofbit | istream::failbit | ifstream::badbit);
    string tmp;

    try
    {
        infile.open(filename);

        infile >> tmp;// skip key

        while(infile.eof() == false)
        {
            infile >> tmp;
            seq += tmp;
        }
        
        infile.close();
    }
    catch(ifstream::failure e)
	{
        if(infile.is_open())
			infile.close();
	}
}

void WriteFsFasta(map<string, vector<GenePrediction>>& predMap)
{
    string filename, seq;

    map<string, vector<GenePrediction>>::iterator mapIt;        
    
    for(mapIt=predMap.begin();mapIt!=predMap.end();++mapIt)
    {
        filename = "./splitFasta/" + mapIt->first + ".fasta";

        ReadFasta(filename, seq);

        if(!seq.empty())
        {
            for(int fsAt, geneAt=0;geneAt<mapIt->second.size();++geneAt)
            {
                for(fsAt=0;fsAt<mapIt->second[geneAt].viableFS.size();++fsAt)
                {
                    cout << ">" << mapIt->first << 
                    "_" << mapIt->second[geneAt].viableFS[fsAt] << endl;
                    cout << FS(seq, mapIt->second[geneAt].viableFS[fsAt]) << endl;


                }
            }

            seq.clear();
        }
    }
}

void WriteFs(map<string, vector<GenePrediction>>& predMap)
{
    string filename = "fs.hints";

    map<string, vector<GenePrediction>>::iterator mapIt;        
    
    ofstream outfile;
    outfile.exceptions(ofstream::eofbit | ostream::failbit | ofstream::badbit);

    try
    {
        outfile.open(filename);

        for(mapIt=predMap.begin();mapIt!=predMap.end();++mapIt)
        {
            for(int fsAt, geneAt=0;geneAt<mapIt->second.size();++geneAt)
            {
                for(fsAt=0;fsAt<mapIt->second[geneAt].viableFS.size();++fsAt)
                {
                    outfile << mapIt->second[geneAt].viableFS[fsAt] << endl;
                }
            }
        }

        outfile.close();
    }
    catch(ifstream::failure e)
	{
		if(outfile.is_open())
			outfile.close();
	}
}

void PrintUsage()
{
    cout << "usage:" << endl << "prepareTesters --strand=STRAND fastafile gfffile" << endl;
    cout << "fastafile contains the original sequences" << endl;
    cout << "gfffile contains Augustus predictions over original sequences" << endl;
    cout << "parameters:" << endl;
    cout << "--strand\tSTRAND can be + for plusstrand only, - for minusstrand only or . for bothstrand" << endl;
    cout << "--exon\tEXON can be single, initial, terminal, internal1 or internal2" << endl;
}

/*
 * main
 */
int main( int argc, char* argv[] )
{
    if(argc != 5)
        PrintUsage();
    else
    {
        string tmp = argv[1];
        
        if(tmp.length() != 10)
        {
            PrintUsage();
            return 0;
        }
        
        tmp.erase(9);
        
        if(tmp!="--strand=")
        {
            PrintUsage();
            return 0;
        }
        
        Strand strand = StrandChar(argv[1][9]);
        
        if(strand == STRAND_UNKNOWN)
        {
            PrintUsage();
            return 0;
        }

        int etype;
        tmp = argv[2];
        
        if(tmp == "--exon=single")
            etype = SINGLE;
        else if(tmp == "--exon=terminal")
            etype = TERMINAL;
        else if(tmp == "--exon=initial")
            etype = INITIAL;
        else if(tmp == "--exon=internal1")
            etype = INTERNAL1;
        else if(tmp == "--exon=internal2")
            etype = INTERNAL2;
        else
        {
            PrintUsage();
            return 0;
        }
        
        map<string, vector<GenePrediction>> predMap;
        ReadPrediction(argv[4], strand, predMap);
        ApplyFs(argv[3], "", etype, predMap);
        // WritePrediction(etype, predMap);
        WriteFsFasta(predMap);
        // WriteFs(predMap);
    }

    return 0;
}
