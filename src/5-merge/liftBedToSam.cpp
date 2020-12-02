//
//  main.cpp
//  AirLift
//
//  Created by Can Fırtına on 11.11.20.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
enum SamFields{qName, flag, rName, pos, mapQ, cigar, rNext, pNext, tLen, seq, qual, optional};

//https://stackoverflow.com/a/236803
template <typename Out>
void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

//https://stackoverflow.com/a/236803
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

void processBedReadName(std::string bedField, std::string& readName, int& pair){
    if(!bedField.substr(bedField.length()-2, 2).compare("\\2") ||
       !bedField.substr(bedField.length()-2, 2).compare("/2")){
        pair = 2;
        readName = bedField.substr(0, bedField.length()-2);
    }else if(!bedField.substr(bedField.length()-2, 2).compare("\\1") || !bedField.substr(bedField.length()-2, 2).compare("/1")){
        pair = 1;
        readName = bedField.substr(0, bedField.length()-2);
    }else{
        pair = 0;
        readName = bedField;
    }
}

void processBedLine(std::string bedLine, std::vector<std::string>& fields, int bedFieldCmpIndex, std::string& readName,
                    int& pair){
    
    std::istringstream bedLineStream(bedLine);
    std::string dataBed;
    int startPosBed;
    bedLineStream >> dataBed >> startPosBed;
    
    fields.push_back(dataBed);
    {std::ostringstream outputBedStream;
    outputBedStream << ++startPosBed;
    fields.push_back(outputBedStream.str());}
    while(!bedLineStream.eof()){
        bedLineStream >> dataBed;
        fields.push_back(dataBed);
    }
    
    processBedReadName(fields[bedFieldCmpIndex], readName, pair);
}

//argv[1] sam file
//argv[2] bed file
//argv[3] read name field in BED file -- 1-based
//argv[4] replace field in SAM file -- 1-based
//argv[5] replace field in BED file (this will replace the field in the SAM file) -- 1-based
int main(int argc, const char* argv[]) {
    
    std::string dash = "-";
    std::ifstream samFileStream;
    try{
        samFileStream.open(argv[1], ios::in);
    }catch(std::ios_base::failure e){
        std::cerr << "Could not open " << argv[1] << std::endl;
        return 1;
    }
    
    std::ifstream bedFileStream;
    bool isBedStreaming = true;
    if(dash.compare(argv[2])){
        isBedStreaming = false;
        try{
            bedFileStream.open(argv[2], ios::in);
        }catch(std::ios_base::failure e){
            std::cerr << "Could not open " << argv[2] << std::endl;
            return 1;
        }
    }
    
    int joinFieldBed = atoi(argv[3])-1;
    std::vector<std::string> replaceFieldSamStr = split(argv[4], ',');
    std::vector<std::string> replaceFieldBedStr = split(argv[5], ',');
    std::vector<int> replaceFieldSam;
    std::vector<int> replaceFieldBed;
    for(int i = 0; i < replaceFieldSamStr.size(); ++i){
        replaceFieldSam.push_back(atoi(replaceFieldSamStr[i].c_str())-1);
        replaceFieldBed.push_back(atoi(replaceFieldBedStr[i].c_str())-1);
    }
    
    if(joinFieldBed == 1 || joinFieldBed == 2){
        std::cerr << "Specified BED field is not supported to include read names" << std::endl;
        return 1;
    }
    
    std::string samLine;
    std::string bedLine;
    
    if(!std::getline((isBedStreaming)?cin:bedFileStream, bedLine)) return 1;
    std::vector<std::string> fieldsBed;
    std::string readNameBed;
    int pair = 0; //0 if not paired, otherwise 1 or 2 pair numbers
    processBedLine(bedLine, fieldsBed, joinFieldBed, readNameBed, pair);
    
    bool moreBedFields = true;
    //We assume regular SAM Format: https://en.wikipedia.org/wiki/SAM_(file_format)
    while(std::getline(samFileStream, samLine)){
        if(moreBedFields){
            std::istringstream samLineStream(samLine);
            std::ostringstream outputSamStream;
            std::vector<std::string> fieldsSam;
            
            std::string dataSam;
            while(!samLineStream.eof()){
                samLineStream >> dataSam;
                fieldsSam.push_back(dataSam);
            }
            
            int isSecond = atoi(fieldsSam[1].c_str()) & 128;
            
            //!match, time to replace...
            if(!fieldsSam[SamFields::qName].compare(readNameBed) &&
               (pair != 2 || isSecond)){
                for(int i = 0; i < replaceFieldSam.size(); ++i)
                    fieldsSam[replaceFieldSam[i]] = fieldsBed[replaceFieldBed[i]];
                    
                if(!std::getline((isBedStreaming)?cin:bedFileStream, bedLine)) moreBedFields = false;
                fieldsBed.clear();
                processBedLine(bedLine, fieldsBed, joinFieldBed, readNameBed, pair);
            }
            
            for (std::vector<std::string>::const_iterator i = fieldsSam.begin(); i != fieldsSam.end(); ++i)
                std::cout << *i << '\t';
        }else std::cout << samLine;
        
        std::cout << std::endl;
    }
    
    samFileStream.close();
    bedFileStream.close();
    return 0;
}

