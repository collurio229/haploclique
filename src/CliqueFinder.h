#ifndef CLIQUEFINDER_H_
#define CLIQUEFINDER_H_

#include "CliqueCollector.h"
#include "AlignmentRecord.h"
#include "EdgeCalculator.h"
#include "CoverageMonitor.h"
#include "EdgeWriter.h"
#include "ReadGroups.h"


using namespace std;
using namespace boost;

class CliqueFinder {
protected:
    const EdgeCalculator & edge_calculator;
    CliqueCollector & clique_collector;
    CoverageMonitor coverage_monitor;
    bool no_sort;

    typedef std::list<Clique*> clique_list_t;
    clique_list_t *cliques;    

    size_t alignment_count;
    EdgeWriter *edge_writer;
    const EdgeCalculator *second_edge_calculator;
    alignment_id_t next_id;

    typedef struct {
        bool operator()(const Clique* c0,const Clique* c1) const {
            return c0->getAlignmentSet() < c1->getAlignmentSet();
        }
    } clique_comp_t;

    typedef struct {
        bool operator()(const Clique* c0,const Clique* c1) const {
            return c0->getAlignmentSet() == c1->getAlignmentSet();
        }
    } clique_equal_t;

public:
    CliqueFinder(const EdgeCalculator& edge_calculator, CliqueCollector& clique_collector, const ReadGroups* read_groups, bool no_sort) : edge_calculator(edge_calculator), clique_collector(clique_collector), coverage_monitor(read_groups), no_sort(no_sort) {
        cliques = new clique_list_t();
    	alignment_count = 0;
    	edge_writer = 0;
    	second_edge_calculator = 0;
        next_id = 0;
    }
    virtual ~CliqueFinder() {
        
    }

    virtual void addAlignment(std::unique_ptr<AlignmentRecord>& ap) = 0;

    virtual void finish() {
        if (edge_writer != 0) {
    		edge_writer->finish();
	    }
	    clique_list_t::iterator clique_it = cliques->begin();
	    for (;clique_it!=cliques->end(); ++clique_it) {
	    	Clique* clique = *clique_it;
	    	clique_collector.add(unique_ptr<Clique>(clique));
	    }
	    delete cliques;
	    cliques = 0;
    }

    virtual CoverageMonitor & getCoverageMonitor() { return coverage_monitor; }

    virtual const AlignmentRecord & getAlignmentByIndex(size_t index) const = 0;

    virtual void setEdgeWriter(EdgeWriter& edge_writer) { this->edge_writer = &edge_writer; }

    /** Add a second edge calculator: edges will only be drawn when both edge calculators agree
	 *  that the edge is present.
	 */
    virtual void setSecondEdgeCalculator(const EdgeCalculator* ec) { this->second_edge_calculator = ec; }
};

#endif
