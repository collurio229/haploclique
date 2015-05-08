#ifndef VERTEXSET_H_
#define VERTEXSET_H_

#include <exception>

using namespace std;

class VertexSet {
private:
    unsigned int max_size_;
    unsigned int capacity_;
    int* bitset_;

    inline unsigned int index_(unsigned int id) { return id / (sizeof(int) * 8); }
    inline unsigned int offset_(unsigned int id) { return id % (sizeof(int) * 8); }
public:
    VertexSet(unsigned int size);
    VertexSet(const VertexSet& set);
    VertexSet(VertexSet&& set);
    virtual ~VertexSet();

    VertexSet& operator=(const VertexSet& set);
    VertexSet& operator=(VertexSet&& set);

    bool isempty() { return (size() == 0); }
    bool ismember(unsigned int id) { return bitset_[index_(id)] & (1 << offset_(id) ); }
    void set(unsigned int id) { bitset_[index_(id)] |= 1 << offset_(id); }
    void unset(unsigned int id) { bitset_[index_(id)] &= ~ (1 << offset_(id)); }
    void flip(unsigned int id) { bitset_[index_(id)] ^= 1 << offset_(id); }
    unsigned int size() {
        unsigned int size = 0;
        for(unsigned int i = 0; i < max_size_; i++) {
            if (this->ismember(i)) size++;
        }
        return size;
    }

    VertexSet& unite(const VertexSet& set);
    VertexSet& intersect(const VertexSet& set);
    VertexSet& complement(const VertexSet& set);

    bool operator[](unsigned int id) { return this->ismember(id); }
    VertexSet& operator|=(const VertexSet& set) { return this->unite(set); } 
    VertexSet& operator&=(const VertexSet& set) { return this->intersect(set); }
    VertexSet& operator-=(const VertexSet& set) { return this->complement(set); }
};

class UnmatchingSetCapacity : public exception {
    virtual const char* what() const throw() {
        return "The capacities of the sets didn't match";
    }
};

#endif // VERTEXSET_H_
