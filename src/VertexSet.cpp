#include "VertexSet.h"

VertexSet::VertexSet(unsigned int size) {
    max_size_ = size;
    capacity_ = size / (sizeof(int) * 8) + 1;    

    bitset_ = new int[capacity_]();    
}

VertexSet::VertexSet(const VertexSet& set)
: max_size_(set.max_size_), capacity_(set.capacity_) {

    bitset_ = new int[capacity_];
    for(unsigned int i = 0; i < capacity_; i++) {
        bitset_[i] = set.bitset_[i];
    }
}

VertexSet::VertexSet(VertexSet&& set) 
: max_size_(set.max_size_), capacity_(set.capacity_), bitset_(set.bitset_) {
    set.bitset_ = nullptr;
}

VertexSet::~VertexSet() {
    delete[] bitset_;
}

VertexSet& VertexSet::operator=(const VertexSet& set) {
    max_size_ = set.max_size_;
    
    if (capacity_ != set.capacity_) {
        capacity_ = set.capacity_;
        delete[] bitset_;
        bitset_ = new int[capacity_];
    }
    
    for (unsigned int i = 0; i < capacity_; i++) {
        bitset_[i] = set.bitset_[i];
    }

    return *this;
}

VertexSet& VertexSet::operator=(VertexSet&& set) {
    max_size_ = set.max_size_;
    
    capacity_ = set.capacity_;

    delete[] bitset_;
    bitset_ = set.bitset_;
    set.bitset_ = nullptr;

    return *this;
}

VertexSet& VertexSet::unite(const VertexSet& set) {
    if (capacity_ != set.capacity_) {
        throw UnmatchingSetCapacity();
    }
    for (unsigned int i = 0; i < capacity_; i++) {
        bitset_[i] |= set.bitset_[i];
    }

    return *this;
}

VertexSet& VertexSet::intersect(const VertexSet& set) {
    if (capacity_ != set.capacity_) {
        throw UnmatchingSetCapacity();
    }
    for (unsigned int i = 0; i < capacity_; i++) {
        bitset_[i] &= set.bitset_[i];
    }

    return *this;
}

VertexSet& VertexSet::complement(const VertexSet& set) {
    if (capacity_ != set.capacity_) {
        throw UnmatchingSetCapacity();
    }
    for (unsigned int i = 0; i < capacity_; i++) {
        bitset_[i] &= ~set.bitset_[i];
    }

    return *this;
}
