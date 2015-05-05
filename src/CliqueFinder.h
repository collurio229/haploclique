#ifndef CLIQUEFINDER_H_
#define CLIQUEFINDER_H_

#include "CliqueCollector.h"
#include "AlignmentRecord.h"
#include "EdgeCalculator.h"
#include "CoverageMonitor.h"
#include "EdgeWriter.h"
#include "ReadGroups.h"

class CliqueFinder {
protected:
    const EdgeCalculator & edge_calculator;
    CliqueCollector & clique_collector;
    CoverageMonitor coverage_monitor;
    typedef std::list<Clique*> clique_list_t;
    clique_list_t *cliques;    
    size_t capacity;
    AlignmentRecord **alignments;
    size_t alignment_count;
    EdgeWriter *edge_writer;
    const EdgeCalculator *second_edge_calculator;
public:
    CliqueFinder(const EdgeCalculator& edge_calculator, CliqueCollector& clique_collector, const ReadGroups* read_groups) : edge_calculator(edge_calculator), clique_collector(clique_collector), coverage_monitor(read_groups) {
        cliques = new clique_list_t();
	    capacity = alignment_set_t::bits_per_block;
	    alignments = new AlignmentRecord*[capacity];
    	alignment_count = 0;
    	edge_writer = 0;
    	second_edge_calculator = 0;

    }
    virtual ~CliqueFinder() {
        for (size_t i=0; i<alignment_count; ++i) {
	    	delete alignments[i];
	    }
	    delete [] alignments;
    }

    virtual void addAlignment(std::auto_ptr<AlignmentRecord> ap) = 0;
    virtual void finish() = 0;
    virtual CoverageMonitor & getCoverageMonitor() { return coverage_monitor; }

    virtual const AlignmentRecord & getAlignmentByIndex(size_t index) const {
        assert(index<alignment_count);
    	return *(alignments[index]);
    }

    virtual void setEdgeWriter(EdgeWriter& edge_writer) { this->edge_writer = &edge_writer; }

    /** Add a second edge calculator: edges will only be drawn when both edge calculators agree
	 *  that the edge is present.
	 */
    virtual void setSecondEdgeCalculator(const EdgeCalculator* ec) { this->second_edge_calculator = ec; }
};

#endif
