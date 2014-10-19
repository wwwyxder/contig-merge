#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <cstdlib>
#include <unistd.h>
#include <utility>
#include <cstdlib>
#include <sys/wait.h>
#include <errno.h>
#include <vector>
#include <algorithm>
#define __DEBUG_MODE_ON__
using namespace std;
struct MapTriad{
    int firstContigID;
    int secondContigID;
    int mappedPos;
};
struct MapTriadLess
{
    bool operator()(const MapTriad& lhs, const MapTriad& rhs) {
        return lhs.firstContigID < rhs.firstContigID;
    }
};
static map<string,string> config;
static vector<int> contigStartPoses;
static vector<int> contigLengthes;
static vector<bool> contigRemoved;
static vector<string> contigSeq;
static vector<MapTriad> mapTriads;
const static string tandemFile("localTandem.fa");
const static string readFile("localReads.fa");

static void readConfig()
{
    ifstream file("config.merge");
    string line;
    string key;
    string value;
    while(getline(file, line)) {
        stringstream sstrm(line);
        sstrm>>key>>value;
        config[key] = value;
    }
    file.close();
}
/*
 * getIdAndLength() :
 * parameter readHead is the head of a read, must be in format like ">103:1562", 
 * where 103 is the read ID(0 based), 1562 is the length of the read. 
 * return a std::pair <ID,length>
 */
static pair<int,int> getIdLength(string readHead)
{
    for(size_t i=0;i<readHead.size();++i) {
        if(readHead[i]=='>'||readHead[i]==':')
            readHead[i] = ' ';
    }
    stringstream sstrm(readHead);
    int id;
    int length;
    sstrm>>id>>length;
    return make_pair<int,int>(id,length);
}
static void tandem_cut()
{
    readConfig();
    if(config["inputFile"]=="") {
        cerr<<"input file no specified in config.merge"<<endl;
        exit(1);
    }
    int readLength=50;
    if(config["cutLength"]!="") {
        sscanf(config["cutLength"].c_str(),"%d",&readLength);
    }
    ifstream infile(config["inputFile"].c_str());
    ofstream tandem(tandemFile.c_str());
    ofstream cutreads(readFile.c_str());
    tandem<<">local\n";
    string line;
    int len = 50;
    int currentPos = 1;
    int IDchecker = 0;
    while(getline(infile,line)) {
        if(line.size()<=0) continue;
        if(line[0]=='>') { //this line is read head
            pair<int,int> position = getIdLength(line);
            if(IDchecker!=position.first) {
                cerr<<"read "<<line<<" has an error ID"<<endl;
                exit(1);
            }
            ++IDchecker;
            cutreads<<">"<<position.first<<":"<<position.second<<":"<<currentPos<<"\n";
            contigStartPoses.push_back(currentPos);
            contigLengthes.push_back(position.second);
            currentPos += position.second;
        } else { //this line is read sequence
            cutreads<<line.substr(0,readLength)<<"\n";
            size_t start = 0;
            while(start+len<=line.size()) {
                tandem<<line.substr(start,len)<<"\n";
                start += len;
                len = 50;
            }
            if(start<line.size()) {
                tandem<<line.substr(start);
                len = len - line.substr(start).size();
            }
            contigSeq.push_back(line);
            contigRemoved.push_back(0);
        }
    }
    infile.close();
    tandem.close();
    cutreads.close();
}
/*
 * find globalPos in [beg, end)
 */
static pair<int,int> _localBinaryFind(int beg, int end, int globalPos)
{
    if(end-beg==1) { //only one element in this interval
#ifdef __DEBUG_MODE_ON__
        if(globalPos-contigStartPoses[beg]>=contigLengthes[beg]) {
            cerr<<"error in binary find : found pos out of range. globalPos="<<globalPos<<endl;
            exit(1);
        }
#endif
        return make_pair<int,int>(beg, globalPos-contigStartPoses[beg]);
    }
    int middle = (beg+end)/2;
    if(contigStartPoses[middle]<=globalPos) {
        return _localBinaryFind(middle, end, globalPos);
    } else {
        return _localBinaryFind(beg, middle, globalPos);
    }
}
inline static pair<int,int> getMapIDPos(int globalPos)
{
    return _localBinaryFind(0,contigStartPoses.size(),globalPos);
}
/*
 * checkOverlap() :
 * if mt is a real overlap return true
 */
static bool embody(const MapTriad&);
static bool checkOverlap(const MapTriad& mt)
{
    for(size_t i=mt.mappedPos;
            i<contigSeq[mt.firstContigID].size() && i-mt.mappedPos<contigSeq[mt.secondContigID].size();
            ++i) {
        if(contigSeq[mt.firstContigID][i] != 
                contigSeq[mt.secondContigID][i-mt.mappedPos])
            return false;
    }
    return true;
}
static void make_candidate()
{
    pid_t pid = fork();
    if(pid<0) {
        perror("fork() error");
        exit(1);
    } else if(pid==0) { //child
        system("./makeCandidatePre.sh");
        exit(0);
    }
    //parent
    while(1) {
        int exit_child = waitpid(pid, NULL, 0);
        if(exit_child<0&&EINTR==errno) continue;
        break;
    }
    ifstream midfile("local.mid");
    ofstream outfile("local.candidate.pre");
    outfile<<"@READ"<<"\t"<<"READLEN"<<"\t"<<"MAPPED-READ"<<"\t"<<"MAPPED-POS-TO-READ"<<endl;
    string line;
    int readID, readLen, readGlobalPos, mappedGlobalPos;
    pair<int,int> p;
    MapTriad mt;
    while(getline(midfile, line)) {
        stringstream sstrm(line);
        sstrm>>readID>>readLen>>readGlobalPos>>mappedGlobalPos;
        p = getMapIDPos(mappedGlobalPos);
        if(p.second==0
            &&(contigLengthes[readID]>contigLengthes[p.first]||(contigLengthes[readID]==contigLengthes[p.first]&&readID<p.first))) {
            continue;
        }
        mt.firstContigID = p.first;
        mt.secondContigID = readID;
        mt.mappedPos = p.second;
        if(!checkOverlap(mt))
            continue;
        outfile<<readID<<"\t"<<readLen<<"\t"<<p.first<<"\t"<<p.second<<"\n";
        mapTriads.push_back(mt);
    }
    stable_sort(mapTriads.begin(), mapTriads.end(), MapTriadLess());
    midfile.close();
    outfile.close();
}
/*
 * embody() :
 * check if mt.secondContig is a substring of mt.firstContig at mt.mappedPos or not
 * return value : true - it is the substring; flase - not the substring;
 */
static bool embody(const MapTriad& mt)
{
    if(contigSeq[mt.firstContigID].size()-mt.mappedPos < contigSeq[mt.secondContigID].size()) {
        return false;
    }
    if( contigSeq[mt.firstContigID].substr(mt.mappedPos, contigSeq[mt.secondContigID].size())
            == contigSeq[mt.secondContigID] ) {
        return true;
    }
    return false;
}
static void contig_remove()
{
    for(vector<MapTriad>::size_type i=0;i<mapTriads.size();++i) {
        if(contigRemoved[ mapTriads[i].secondContigID ])
            continue;
//        if(contigRemoved[ mapTriads[i].firstContigID ]) {
//            contigRemoved[ mapTriads[i].secondContigID ] = true;
//            continue;
//        }
        if(embody(mapTriads[i])) {
            contigRemoved[ mapTriads[i].secondContigID ] = true;
        }
    }
}
static void write_overlap()
{
    map<int,int> newID;
    int cnter=0;
#ifdef __DEBUG_MODE_ON__
    cerr<<"contigRemoved.size():"<<contigRemoved.size()<<endl;
#endif
    for(vector<bool>::size_type i=0;i<contigRemoved.size();++i) {
        if(!contigRemoved[i]) {
            newID[i] = cnter;
            ++cnter;
        }
    }
    if(config["outputFile"]=="") {
        cerr<<"output file no specified in config.merge"<<endl;
        exit(1);
    }
    ofstream outFastaFile(config["outputFile"].c_str());
    for(vector<string>::size_type i=0;i<contigSeq.size();++i) {
        if(contigRemoved[i]) continue;
        outFastaFile<<">"<<newID[i]<<":"<<contigSeq[i].size()<<"\n";
        outFastaFile<<contigSeq[i]<<"\n";
    }
    outFastaFile.close();

    int minimumOL = 50;
    if(config["minimumOverlapLength"]!="") {
        stringstream sstrm(config["minimumOverlapLength"]);
        sstrm>>minimumOL;
    }
    ofstream overlapFile("local.overlap");
    for(vector<MapTriad>::size_type i=0;i<mapTriads.size();++i) {
        if(contigRemoved[ mapTriads[i].firstContigID ] 
                || contigRemoved[ mapTriads[i].secondContigID ])
            continue;
        int overlapLength 
            = contigSeq[ mapTriads[i].firstContigID ].size() - mapTriads[i].mappedPos;
        if(overlapLength<minimumOL)
            continue;
        overlapFile<<newID[ mapTriads[i].firstContigID ]
            <<"\t"<<newID[ mapTriads[i].secondContigID ]
            <<"\t"<<mapTriads[i].mappedPos
            <<"\t"<<overlapLength<<"\n";
    }
}
int main()
{
    for(map<string,string>::iterator iter=config.begin();iter!=config.end();++iter) {
        cout<<iter->first<<":"<<iter->second<<endl;
    }
    tandem_cut();
    make_candidate();
    contig_remove();
    write_overlap();
    return 0;
}

