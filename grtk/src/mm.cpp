#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include "mm.hpp"

#include <iostream>

using namespace std;

MMFile::MMFile(string path) {
    struct stat st;
    stat(path.c_str(), &st);
    this->fileSize = st.st_size;

    handle = open(path.c_str(), O_RDONLY);
    this->data = mmap(0, fileSize, PROT_READ, 
            MAP_SHARED, handle, 0);
    if (this->data == MAP_FAILED) {
        throw 5;
    }
}

MMFile::~MMFile() {
    munmap(this->data, this->fileSize);
    close(handle);
}

void*
MMFile::operator[](size_t offset) {
    return ((char*) data) + offset;
}
