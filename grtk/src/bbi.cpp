#include <stdint.h> 
#include <cstring>
#include <cstdio>
#include <cassert>

#include <zlib.h>

#include "mm.hpp"
#include "bbi.hpp"

using namespace std;

/* Read a struct from a pointer and increment the pointer */
template <typename T>
T* readPtr(char** ptr) {
    T* pT = (T*) *ptr;
    *ptr += sizeof(T);
    return pT;
}

/* Region methods */

bool 
Region::overlaps(const Region& o)
{
    // FIXME: figure out gt/lt or gte/lte
    return
        chrom == o.chrom &&
        end >= o.start &&
        start <= o.end;
}

bool 
BED::overlaps(const BED& o)
{
    return
        chrom == o.chrom &&
        end >= o.start &&
        start <= o.end;
}


/* BBI File implementation */

void 
BBIFile::readChromosomes(uint32_t keySize, char* pNode) {
    // Recursively read chromosome data from BTree
    TreeNode* node = readPtr<TreeNode>(&pNode);

    size_t elemSize = keySize + 8;
    char key[keySize];

    for (int i=0; i<node->count; i++) {
        memcpy(key, pNode, keySize);
        pNode += keySize;
        if (node->isLeaf) {
            Chromosome ch = {
                string(key, keySize),
                *(uint32_t*) pNode,
                *(uint32_t*) (pNode + 4)
            };
            genome.push_back(ch);
        } else {
            uint64_t childOffset = *(uint64_t*) pNode;
            char* pChild = (char*) (*handle)[childOffset];
            readChromosomes(keySize, pChild);
        }
        pNode += 8;
    }
}

// Search the R-Tree index for data nodes overlapping a query
void 
BBIFile::searchIndex(char* pNode, const Region& query,
        vector<RLeaf>& hits) {
    TreeNode* node = readPtr<TreeNode>(&pNode);
    
    // FIXME: Linear search
    for (int i=0; i<node->count; i++) {
        if (node->isLeaf) {
            RLeaf* e = readPtr<RLeaf>(&pNode);
            /* Replicate index entries that span
             * multiple chromosomes so that each entry
             * is on only one chromosome. */
            for (uint32_t chr=e->startChromIx; chr<=e->endChromIx; chr++) {
                pos_t start = chr == e->startChromIx ?
                    e->startBase : 0;
                pos_t end = chr == e->endChromIx ? e->endBase : 
                    genome[chr].length;
                Region rNode = {chr, start, end};
                if (rNode.overlaps(query)) {
                    hits.push_back(*e);
                    // Each query can only overlap one chromosome
                    break;
                }
            }
        } else {
            RNonLeaf* e = readPtr<RNonLeaf>(&pNode);
            searchIndex((char*) (*handle)[e->dataOffset],
                    query, hits);
        }
    }
}

vector<RLeaf> 
BBIFile::searchIndex(const Region& query) {
    void* rtreeOffset = (*handle)[header->fullIndexOffset];
    RHeader* rHeader = (RHeader*) rtreeOffset;
    assert(rHeader->magic == RTreeMagicNumber);

    vector<RLeaf> hits;
    searchIndex(((char*) rtreeOffset) + sizeof(RHeader), query, hits);
    return hits;
}

BBIFile::BBIFile(const char* path) {
    handle = new MMFile(string(path));
    header = (BBIHeader*) (*handle)[0];
    assert(header->magic==BigWigMagicNumber 
        || header->magic==BigBEDMagicNumber);
    uncBuf = new char[header->uncompressBufSize];

    // Read the chromosome names, IDs, and sizes from BTree
    void* btreeOffset = (*handle)[header->chromosomeTreeOffset];
    BHeader* bHeader = (BHeader*) btreeOffset;
    assert(bHeader->magic==BTreeMagicNumber);
    readChromosomes(bHeader->keySize, 
            (char*) btreeOffset + sizeof(BHeader));
}

BBIFile::~BBIFile() {
    delete handle;
    delete [] uncBuf;
}

vector<BED> 
BBIFile::search(const Region& query) {
    vector<BED> rs;
    vector<RLeaf> hits = searchIndex(query);

    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    for (RLeaf hit : hits) {
        int err;

        err = inflateInit2(&strm, 15 + 32);
        // if (err != Z_OK)
        strm.avail_in = hit.dataSize;
        strm.next_in = (Bytef*) (*handle)[hit.dataOffset];
        strm.avail_out = header->uncompressBufSize;
        strm.next_out = (Bytef*) uncBuf;
        inflate(&strm, Z_NO_FLUSH);

        for ( BED bed : searchSection(query, strm.avail_out)) {
            rs.push_back(bed);
        }
        inflateEnd(&strm);
    }
    return rs;
}

// FIXME: could use a map instead of linear search
// FIXME: restrict to strand
vector<BED> 
BBIFile::search(const BED& query) {
    uint32_t chrom = UINT_MAX;
    for (Chromosome ch : genome) {
        if (!strcmp(query.chrom,ch.key.c_str()))
            chrom = ch.id;
    }
    assert(chrom != UINT_MAX);
    Region q = {chrom, query.start, query.end};
    return search(q);
}

/* BigWig */

// FIXME: determine appropriate WIG type (fixed, bedgraph, etc) and convert.
// Currently assumes bedgraph
// This also assumes that the chromosome matches
vector<BED> 
BigWigFile::searchSection(const Region& query, int avail) {
    vector<BED> rs;
    char* pNode = uncBuf;
    BinaryWigHeader* header = readPtr<BinaryWigHeader>(&pNode);
    for (int i=0; i<header->itemCount; i++) {
        BEDGraphDatum* datum = readPtr<BEDGraphDatum>(&pNode);
        if (datum->start < query.end && query.start < datum->end) {
            rs.push_back(BED{genome[0].key.c_str(),datum->start,
                    datum->end,
                    NULL, datum->value, 0, 0});
        }
    }
    return rs;
}

BigWigFile::BigWigFile(const char* path) : BBIFile(path) {};

// FIXME: handle exons
Summary
BigWigFile::summary(const BED& q) {
    Summary s;
    s.sum = 0;
    s.mean0 = 0;
    s.mean = 0;
    s.covered = 0;
    s.length = q.end - q.start;
    vector<BED> rs = search(q);
    for (BED r : rs) {
        uint32_t covered = min(q.end, r.end) - max(q.start, r.start);
        s.covered += covered;
        s.sum += covered * r.score;
    }
    s.mean0 = s.sum ? 1.0 * s.sum / s.length : 0;
    s.mean = s.sum ? 1.0 * s.sum / s.covered : 0;
    return s;
}

/* BigBED */

vector<BED> 
BigBEDFile::searchSection(const Region& query, int avail) {
    vector<BED> rs;
    char* pNode = uncBuf;
    while (pNode < (uncBuf + avail)) {
        uint32_t chromId = *readPtr<uint32_t>(&pNode);
        uint32_t start = *readPtr<uint32_t>(&pNode);
        uint32_t end = *readPtr<uint32_t>(&pNode);
        char* rest = pNode;
        pNode = strchr(pNode, '\0') + 1;

        Region region = {chromId, start, end};
        if (region.overlaps(query)) {
            const char* chrom = genome[chromId].key.c_str();
            BED bed = {chrom, start, end, rest};
            rs.push_back(bed);
        }
    }
    return rs;
}

BigBEDFile::BigBEDFile(const char* path) : BBIFile(path) {};
