#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <cstdlib>
using namespace std;
static map<string,string> config;
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
    while(getline(infile,line)) {
        if(line.size()<=0) continue;
        if(line[0]=='>') {
            cutreads<<line<<"\n";
        } else {
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
int main()
{
    for(map<string,string>::iterator iter=config.begin();iter!=config.end();++iter) {
        cout<<iter->first<<":"<<iter->second<<endl;
    }
    tandem_cut();
    return 0;
}

