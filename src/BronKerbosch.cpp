#include "BronKerbosch.h"
#include <boost/unordered_set.hpp>
#include <boost/dynamic_bitset.hpp>

// using namespace boost;
using namespace std;

BronKerbosch::BronKerbosch(const EdgeCalculator& edge_calculator, CliqueCollector& clique_collector, const ReadGroups* read_groups, bool no_sort)
: CliqueFinder(edge_calculator, clique_collector, read_groups, no_sort), alignments_() {
    order_ = nullptr;
    vertices_ = nullptr;
    degree_map_ = new degree_map_t();
    actives_ = new list<adjacency_list_t*>();
    vertices_as_lists_ = new vector<adjacency_list_t*>();
}

BronKerbosch::~BronKerbosch() {
    if (cliques != nullptr) {
        finish();
    }
    for (auto&& alignment : alignments_) {
        delete alignment;
    }
    delete order_;
    delete vertices_;
}

void BronKerbosch::finish() {
    degeneracy_order();

    alignment_set_t R(alignment_count);
    alignment_set_t P(alignment_count);
    P.flip();
    alignment_set_t X(alignment_count);

    for (auto&& i : *order_) {

        bronkerbosch(R.set(i), P & (*vertices_)[i], X & (*vertices_)[i]);

        R.reset(i);
        P.reset(i);
        X.set(i);

    }

    for (auto clique_it = cliques->begin();clique_it!=cliques->end(); ++clique_it) {
    	Clique* clique = *clique_it;
        clique_collector.add(unique_ptr<Clique>(clique));
   }
	    delete cliques;
	    cliques = nullptr;
}

void BronKerbosch::degeneracy_order() {

    // Move all still active vertices into the degree map
    for (auto it = actives_->begin(); it != actives_->end();) {
        list<adjacency_list_t*>::size_type ind = (*it)->second.size();

        if(degree_map_->count(ind) == 0) {           
            degree_map_->emplace(ind, list<adjacency_list_t*>(1, *it) );
        } else {
            degree_map_->at(ind).push_back(*it);
        }
        
        it = actives_->erase(it);
    }

    // Delete isolated vertices
    degree_map_->erase(0);

    order_ = new list<size_t>();
    vertices_ = new vector<alignment_set_t>(alignment_count, alignment_set_t(alignment_count));

    while (not degree_map_->empty()) {
        auto it = degree_map_->begin();

        if (it->second.empty()) {
            degree_map_->erase(it);
            continue;
        }
        
        adjacency_list_t* list_vertex = it->second.front();
        it->second.pop_front();

        // Update degree map
        for (auto&& index : list_vertex->second) {
            adjacency_list_t* companion = (*vertices_as_lists_)[index];
            list<adjacency_list_t*>::size_type s = companion->second.size();

            (*vertices_)[list_vertex->first].set(companion->first);
            (*vertices_)[companion->first].set(list_vertex->first);

            companion->second.remove(list_vertex->first);

            degree_map_->at(s).remove(companion);

            if (degree_map_->count(s-1) == 0) {
                degree_map_->emplace(s-1, list<adjacency_list_t*>(1, companion) );                
            } else {
                degree_map_->at(s-1).push_front(companion);
            }
        }
        order_->push_back(list_vertex->first);

        delete list_vertex;
    }

    delete degree_map_;
    delete actives_;
    delete vertices_as_lists_;
}

alignment_set_t::size_type BronKerbosch::find_pivot(const alignment_set_t& P, const alignment_set_t& X) {
    alignment_set_t unit = P | X;
    alignment_set_t::size_type max_size = 0;
    alignment_set_t::size_type index = unit.find_first();

    for(auto i = unit.find_next(index); i != alignment_set_t::npos; i = unit.find_next(i)) {
        auto new_size = (P & (*vertices_)[i]).count();        

        if ( new_size > max_size ) {
            max_size = new_size;
            index = i;
        }
    }

    assert(index != alignment_set_t::npos);

    return index;
}

void BronKerbosch::bronkerbosch(alignment_set_t R, alignment_set_t P, alignment_set_t X) {
    if (P.none() and X.none()) {
        unique_ptr<alignment_set_t> ptr(new alignment_set_t(R) );
        cliques->push_back(new Clique(*this,  ptr) );
        return;
    }

    alignment_set_t::size_type pivot = find_pivot(P, X);
    alignment_set_t comp = P - (*vertices_)[pivot];

    for (auto i = comp.find_first(); i != alignment_set_t::npos; i = comp.find_next(i)) {
        bronkerbosch(R.set(i), P & (*vertices_)[i], X & (*vertices_)[i]);

        R.reset(i);
        P.reset(i);
        X.set(i);
    }
}

void BronKerbosch::addAlignment(std::unique_ptr<AlignmentRecord>& alignment_autoptr) {
	assert(alignment_autoptr.get() != 0);
	assert(cliques!=0);
	alignment_id_t id = next_id++;
	AlignmentRecord* alignment = alignment_autoptr.release();
	alignment->setID(id);
	coverage_monitor.addAlignment(*alignment);

	size_t index = alignment_count++;
	alignments_.push_back(alignment);

    adjacency_list_t* vertex = new adjacency_list_t(index, list<size_t>());

    vertices_as_lists_->push_back(vertex);

    for (auto it = actives_->begin(); it != actives_->end();) {
        AlignmentRecord* alignment2 = alignments_[(*it)->first];

        if (alignment->getIntervalStart() > alignment2->getIntervalEnd()) {
            list<adjacency_list_t*>::size_type ind = (*it)->second.size();

            if(degree_map_->count(ind) == 0) {           
                degree_map_->emplace(ind, list<adjacency_list_t*>(1, *it) );
            } else {
                degree_map_->at(ind).push_back(*it);
            }
            
            it = actives_->erase(it);
            continue;
        }

        if(edge_calculator.edgeBetween(*alignment, *alignment2) ) {

            if(second_edge_calculator != nullptr) {
                if(not second_edge_calculator->edgeBetween(*alignment, *alignment2) ) {
                    continue;
                }
            }

            // Draw edge between current alignment and old alignment in vert
            get<1>(*vertex).push_back((*it)->first);
            get<1>(**it).push_back(vertex->first);

   			if (edge_writer != 0) {
				edge_writer->addEdge(*alignment, *alignment2);
			}
        }
        it++;
    }

    actives_->push_back(vertex);

}
