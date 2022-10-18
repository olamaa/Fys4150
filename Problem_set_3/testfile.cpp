#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include<iterator>
#include<string>
#include<fstream>
#include<iomanip>




int main(){

std::ofstream ofile;
std::string filename;
for (int i = 0;i<5;i++){
        filename = std::to_string(i) + ".txt";
        ofile.open(filename);
        ofile.close();
        ofile.open(filename, std::ofstream::app);
    for (int y =500;y<505;y++){
        
        ofile << y << std::endl;
        
    }
    ofile.close();
}

}