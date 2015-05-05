/* Copyright 2012-2014 Tobias Marschall and Armin Töpfer
 *
 * This file is part of HaploClique.
 *
 * HaploClique is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HaploClique is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HaploClique.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CLEVER_H_
#define CLEVER_H_

#include <set>
#include <list>
#include <boost/unordered_map.hpp>
 #include <boost/dynamic_bitset.hpp>

#include "Clique.h"
#include "CliqueFinder.h"

/** Implementation of the Maximal Clique Enumeration algorithm of CLEVER */
class CLEVER : public CliqueFinder {
private:
    typedef std::pair<unsigned int,size_t> length_and_index_t;
    std::set<length_and_index_t> alignments_by_length;
    alignment_id_t next_id;
    bool no_sort;
    void reorganize_storage();

    typedef struct {
        bool operator()(const Clique* c0,const Clique* c1) const {
            // __uint128_t bitInteger = 0;
            // int size = c0->getAlignmentSet().size();
            // for (int j = 0; j < size; j++) {
            //      if (c0->getAlignmentSet()[j]) {
            //         bitInteger |= (1 << j);
            //      }
            //  }
            //std::cerr << "SORT " << c0->leftmostSegmentStart() << "\t" << c1->leftmostSegmentStart();
            // if (c0->leftmostSegmentStart() < c1->leftmostSegmentStart()) {
                // std::cerr << "\t" << "!" << std::endl;
                // return 1;
            // }
            // int size = c0->getAlignmentSet().size();
            // if (size != c1->getAlignmentSet().size()) {
            //     throw "Bit vector are of different size.";
            // }
            // for (int i = 0; i < size; i++) {
            //     if (c0->getAlignmentSet()[i] != c1->getAlignmentSet()[i]) {
            //         if (c0->getAlignmentSet()[i] == 0) {
            //             return 0;
            //         } else {
            //             return 1;
            //         }
            //     }
            // }
            // return 0;
            return c0->getAlignmentSet() < c1->getAlignmentSet();
        }
    } clique_comp_t;

    typedef struct {
        bool operator()(const Clique* c0,const Clique* c1) const {
            return c0->getAlignmentSet() == c1->getAlignmentSet();
        }
    } clique_equal_t;

public:
    CLEVER(const EdgeCalculator& edge_calculator, CliqueCollector& clique_collector, const ReadGroups* read_groups, bool no_sort);
    virtual ~CLEVER();
    void addAlignment(std::auto_ptr<AlignmentRecord> ap);
    void finish();
};

#endif /* CLEVER_H_ */
