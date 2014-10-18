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
#define __DEBUG_MODE_ON__
using namespace std;
static map<string,string> config;
static vector<int> readStartPoses;
static vector<int> readLengthes;
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
            readStartPoses.push_back(currentPos);
            readLengthes.push_back(position.second);
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
        if(globalPos-readStartPoses[beg]>=readLengthes[beg]) {
            cerr<<"error in binary find : found pos out of range. globalPos="<<globalPos<<endl;
            exit(1);
        }
#endif
        return make_pair<int,int>(beg, globalPos-readStartPoses[beg]);
    }
    int middle = (beg+end)/2;
    if(readStartPoses[middle]<=globalPos) {
        return _localBinaryFind(middle, end, globalPos);
    } else {
        return _localBinaryFind(beg, middle, globalPos);
    }
}
inline static pair<int,int> getMapIDPos(int globalPos)
{
    return _localBinaryFind(0,readStartPoses.size(),globalPos);
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
    cerr<<"make condidate-pre success"<<endl;
    ifstream midfile("local.mid");
    ofstream outfile("local.candidate");
    outfile<<"@READ"<<"\t"<<"READLEN"<<"\t"<<"MAPPED-READ"<<"\t"<<"MAPPED-POS-TO-READ"<<endl;
    string line;
    int readID, readLen, readGlobalPos, mappedGlobalPos;
    pair<int,int> p;
    while(getline(midfile, line)) {
        stringstream sstrm(line);
        sstrm>>readID>>readLen>>readGlobalPos>>mappedGlobalPos;
        p = getMapIDPos(mappedGlobalPos);
        if(p.second==0
            &&(readLengthes[readID]>readLengthes[p.first]||(readLengthes[readID]==readLengthes[p.first]&&readID<p.first))) {
            continue;
        }
        outfile<<readID<<"\t"<<readLen<<"\t"<<p.first<<"\t"<<p.second<<"\n";
    }
    midfile.close();
    outfile.close();
}
int main()
{
    for(map<string,string>::iterator iter=config.begin();iter!=config.end();++iter) {
        cout<<iter->first<<":"<<iter->second<<endl;
    }
    tandem_cut();
    make_candidate();
    return 0;
}

