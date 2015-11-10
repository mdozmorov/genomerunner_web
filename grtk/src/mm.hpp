#include <string>

class MMFile {
    int handle;
    off_t fileSize;
    void *data;
public:
    MMFile(std::string path);
    ~MMFile();
    void* operator[](size_t offset);
};
