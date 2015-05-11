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

    min_key_ = 0;
    max_key_ = 0;
}

BronKerbosch::~BronKerbosch() {
    if (cliques != nullptr) {
        finish();
    }
}

void BronKerbosch::finish() {
    degeneracy_order();
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
        order_->push_back(list_vertex->second);
    }

    delete degree_map_;
    delete actives_;
    delete vertices_as_lists_;
}

void BronKerbosch::bronkerbosch(alignment_set_t R, alignment_set_t P, alignment_set_t X) {
    if (P.none() and X.none()) {
        cliques->push_back(new Clique(*this, R));
    }

    alignment_set_t pivot = choosePivot(P | X);

    for (v in P
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
