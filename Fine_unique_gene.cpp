#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<map>
#include<utility>
#include<vector>
#include<iterator>

using namespace std;

int main()
{
    ifstream ac("./H3K56ac_gene");
    if(!ac){
        cout<<"H3K56ac cannot be opened!"<<endl;
        exit(1);
    }
    ifstream HMG("./H3K56HMG_gene");
    if(!HMG){
        cout<<"H3K56HMG cannot be opened!"<<endl;
        exit(1);
    }
    string element;
    map<string, int> H3K56ac;
    vector<string> H3K56ac_unique, H3K56HMG_unique, same;
    while(ac>>element) H3K56ac[element];
    while(HMG>>element){
        if(H3K56ac.count(element)) ++H3K56ac[element];
        else H3K56HMG_unique.push_back(element);
    }
    map<string, int>::const_iterator map_it = H3K56ac.begin();
    while(map_it != H3K56ac.end()){
        if((map_it->second) == 0)
            H3K56ac_unique.push_back(map_it->first);
        else
            same.push_back(map_it->first);
        ++map_it;
    }
    ac.close();
    HMG.close();

    ofstream out;
    vector<string>::const_iterator H3K56ac_it = H3K56ac_unique.begin(), H3K56HMG_it = H3K56HMG_unique.begin(), same_it = same.begin();
    out.open("./H3K56ac_unique");
    while(H3K56ac_it != H3K56ac_unique.end()){out<<*H3K56ac_it<<endl;++H3K56ac_it;}
    out.close();
    out.open("./H3K56HMG_unique");
    while(H3K56HMG_it != H3K56HMG_unique.end()){out<<*H3K56HMG_it<<endl;++H3K56HMG_it;}
    out.close();
    out.open("./same");
    while(same_it != same.end()){out<<*same_it<<endl;++same_it;}
    out.close();

	return 0;
}
