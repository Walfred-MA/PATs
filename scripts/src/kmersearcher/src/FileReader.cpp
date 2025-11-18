//
//  FileReader.cpp
//  Ctyper2
//
//  Created by walfred on 6/24/25.
//

#include "FileReader.hpp"


void FileReader::Load()
{
    if (strlen(filepath)<2) return;
    
    fafile.open(filepath, ios_base::in);
    
    if (!fafile) {
        
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
        
        return;
    }
    
}

bool FileReader::nextLine(std::string &StrLine)
{
    return (bool)( getline(fafile,StrLine));
}

void FileReader::Close()
{
    
    fafile.close();
    
}

void FileReader::Reset()
{
    fafile.seekg(0, std::ios::beg);
}

