#include <vector>

/* Genomic region stuff */

typedef uint32_t pos_t;

/* A genomic region which has had the chrom string converted
 * to the proper BBI internal index */
struct Region {
    uint32_t chrom;
    pos_t start, end;
    bool overlaps(const Region& o);
};

struct BED {
    const char* chrom;
    uint32_t start, end;
    char* name;
    float score;
    char strand;
    const char* rest;

    bool overlaps(const BED& o);
};

/* 
 * BBI binary data format structs 
 */

enum MagicNumber {
    BigWigMagicNumber = 0x888FFC26,
    BigBEDMagicNumber = 0x8789F2EB,
    BTreeMagicNumber = 0x78CA8C91,
    RTreeMagicNumber = 0x2468ACE0
};

struct BBIHeader {
    uint32_t magic;
    uint16_t version, zoomLevels;
    uint64_t chromosomeTreeOffset, fullDataOffset, fullIndexOffset;
    uint16_t fieldCount, definedFieldCount;
    uint64_t autoSqlOffset, totalSummaryOffset;
    uint32_t uncompressBufSize;
} __attribute__((packed));

/* Shared between BTree and RTree */

struct TreeNode {
    bool isLeaf, reserved;
    uint16_t count;
} __attribute__((packed));

/* BTree */

struct BHeader {
    uint32_t magic, blockSize, keySize, valSize;
    uint64_t itemCount, reserved;
} __attribute__((packed));

struct Chromosome {
    std::string key;
    uint32_t id, length;
};

/* RTree */

struct RHeader {
    uint32_t magic, blockSize;
    uint64_t itemCount;
    uint32_t startChromIx, startBase, endChromIx, endBase;
    uint64_t endFileOffset;
    uint32_t itemsPerSlot, reserved;
} __attribute__((packed));

struct RNonLeaf {
    uint32_t startChromIx, startBase, endChromIx, endBase;
    uint64_t dataOffset;
} __attribute__((packed));

struct RLeaf {
    uint32_t startChromIx, startBase, endChromIx, endBase;
    uint64_t dataOffset, dataSize;
} __attribute__((packed));


struct BinaryWigHeader {
    uint32_t chromId, chromStart, chromEnd, itemStep, itemSpan;
    uint8_t type, reserved;
    uint16_t itemCount;
} __attribute__((packed));

struct BEDGraphDatum {
    uint32_t start, end;
    float value;
} __attribute__((packed));


/* 
 * BBI file API classes 
 */

struct Summary {
    uint32_t length, covered;
    float sum, mean0, mean;
};

class BBIFile {
protected:
    MMFile* handle;
    BBIHeader* header;
    std::vector<Chromosome> genome;
    char* uncBuf; // decompression buffer

    BBIFile(const char* path);
    virtual ~BBIFile();

    void readChromosomes(uint32_t keySize, char* pNode);
    void searchIndex(char* pNode, const Region& query, 
            std::vector<RLeaf>& hits);
    std::vector<RLeaf> searchIndex(const Region& query);
    std::vector<BED> search(const Region& query);

    virtual std::vector<BED> searchSection(const Region&, int) = 0;
public:
    std::vector<BED> search(const BED& query);
};

class BigWigFile : public BBIFile {
protected:
    std::vector<BED> searchSection(const Region&, int);
public:
    BigWigFile(const char* path);
    /* Compute coverage information about a particular region */
    Summary summary(const BED& bed);
};

class BigBEDFile : public BBIFile {
protected:
    std::vector<BED> searchSection(const Region& query, int avail);
public: 
    BigBEDFile(const char* path);
};
