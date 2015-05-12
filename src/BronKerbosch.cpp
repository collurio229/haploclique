#include "BronKerbosch.h"

using namespace boost;
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

    bronkerbosch(alignment_set_t(alignment_count), alignment_set_t(alignment_count).flip(), alignment_set_t(alignment_count));

    if (no_sort == 0) {
		sort(cliques.begin(), cliques.end(), clique_comp_t());
		// cliques.erase(std::unique(cliques.begin(), cliques.end(),clique_equal_t()), cliques.end());
	}

    for (auto clique_it = cliques->begin();clique_it!=cliques->end(); ++clique_it) {
    	Clique* clique = *clique_it;
        clique_collector.add(auto_ptr<Clique>(clique));
   }
	    delete cliques;
	    cliques = nullptr;
}

void BronKerbosch::degeneracy_order() {
    order_ = new list<size_t>();
    vertices_ = new vector<alignment_set_t>(alignment_count, alignment_set_t(alignment_count));

    while (not degree_map_->empty()) {
        auto it = degree_map_->begin();
        
        adjacency_list_t* list_vertex = it->second.front();
        it->second.pop_front();

        for (auto&& index : list_vertex->second) {
            size_type s = (*vertices_as_lists_)[index]->second.size();
            adjacency_list_t* companion = (*vertices_as_lists_)[index];

            vertices_[list_vertex->first].set(companion->first);
            vertices_[companion->first].set(list_vertex->first);

            companion->second.remove(list_vertex->first);

            degree_map_->at(s-1).push_back(companion);
            for (auto iterator = degree_map_->at(s).begin(); iterator != degree_map_->at(s).end();) {
                if( (*iterator)->first == companion->first) {
                    degree_map_->at(s).erase(iterator);
                    break;
                }
                iterator++;
            }
        }
        order_->push_back(list_vertex->first);

        delete list_vertex;
    }

    delete degree_map_;
    delete actives_;
    delete vertices_as_lists_;
}

size_type BronKerbosch::pivot(alignment_set_t& P, alignment_set_t& X) {
    unit = P | X;
    size_type max_size = 0;
    size_type index = alignment_set_t::npos;

    for(size_type i = unit.find_first(); i != alignment_set_t::npos; i = unit.find_next(i)) {
        size_type new_size = (P & (*vertices_)[i]).count();        

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
        cliques->push_back(new Clique(*this, R));
    }

    size_type pivot = pivot(P, X);
    comp = P - (*vertices_)[pivot];

    for (size_type i = comp.find_first(); i != alignment_set_t::npos; i = comp.find_next(i)) {
        bronkerbosch(R.set(i), P & (*vertices_)[i], X & (*vertices_)[i]);

        P.reset(i);
        X.set(i);
    }
}

void BronKerbosch::addAlignment(std::auto_ptr<AlignmentRecord> ap) {
	assert(alignment_autoptr.get() != 0);
	assert(cliques!=0);
	alignment_id_t id = next_id++;
	AlignmentRecord* alignment = alignment_autoptr.release();
	alignment->setID(id);
	coverage_monitor.addAlignment(*alignment);

	size_t index = alignment_count++;
	alignments.push_back(alignment);

    adjacency_list_t* vertex = new adjacency_list_t(index, list<size_t>());

    vertices_as_lists_->push_back(vertex);

    for (auto it = actives_->begin(); it != actives->end();) {
        alignment2* = alignments[(*it)->first];

        if (alignment->getIntervalStart() > alignment2->getIntervalEnd()) {
            size_type ind = (*it)->second.size();

            if(degree_map_->count(ind) == 0) {           
                degree_map_->emplace(ind, list(1, *it) );
            } else {
                degree_map->at(ind).push_back(*it);
            }
            
            it = actives_->erase(it);
            continue;
        }

        if(edge_calculator.edgeBetween(alignment, alignment2) ) {

            if(second_edge_calculator != 0) {
                if(not second_edge_calculator.edgeBetween(alignment, alignment2) ) {
                    continue;
                }
            }

            // Draw edge between current alignment and old alignment in vert
            get<1>(*vertex).push_back((*it)->first);
            get<1>(**it).push_back(vertex->first);
        }
        it++;
    }

}
