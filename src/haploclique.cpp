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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <limits>
#include <cassert>
#include <ctime>
#include <algorithm>
#include <deque>

#include "docopt/docopt.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <api/BamReader.h>
#include <api/BamAlignment.h>

#include "AlignmentRecord.h"
#include "QuasispeciesEdgeCalculator.h"
#include "CliqueFinder.h"
#include "CLEVER.h"
#include "BronKerbosch.h"
#include "CliqueCollector.h"
#include "CoverageMonitor.h"
#include "HistogramBasedDistribution.h"
#include "AnyDistributionEdgeCalculator.h"
#include "ReadSetSignificanceTester.h"
#include "ReadSetZTester.h"
#include "ReadSetGenericTester.h"
#include "ReadGroups.h"
#include "ReadGroupAwareEdgeCalculator.h"
#include "ReadSetGroupWiseZTester.h"
#include "GaussianEdgeCalculator.h"

using namespace std;
using namespace boost;

static const char USAGE[] =
R"(haploclique predicts haplotypes from NGS reads.

Usage:
  haploclique bronkerbosch [options] <bamfile> [<output>]
  haploclique [clever] [options] <bamfile> [<output>]

  clever        use the original clever clique finder
  bronkerbosch  use the Bron-Kerbosch based clique finder

options:
  -e <file> --write_edges <file>              Write graph edges to file
  -q <num> --edge_quasi_cutoff_cliques <num>  edge calculator option
                                              [default: 0.99]
  -k <num> --edge_quasi_cutoff_mixed <num>    edge calculator option
                                              [default: 0.97]
  -g <num> --edge_quasi_cutoff_single <num>   edge calculator option
                                              [default: 0.95]
  -Q <num> --random_overlap_quality <num>     edge calculator option
                                              [default: 0.9]
  -m --frame_shift_merge                      Reads will be clustered with 
                                              single nucleotide insertions or
                                              deletions. Use for PacBio data.
  -o <num> --min_overlap_cliques <num>        edge calculator option
                                              [default: 0.9]
  -j <num> --min_overlap_single <num>         edge calculator option
                                              [default: 0.6]
  -A <file> --allel_frequencies <file>
  -I <file> --call_indels <file>              variant calling is not supported
                                              yet.
  -M <file> --mean_and_sd_filename <file>     Required for option -I
  -p <num> --indel_edge_sig_level <num>       [default: 0.2]
  -i <num> --iterations <num>                 Number of iterations.
                                              haploclique will stop if the 
                                              superreads converge.
                                              [default: -1]
  -f <num> --filter <num>                     Filter out reads with low
                                              frequency at the end.
                                              [default: 0.0]
  -n --no_singletons                          Filter out single read cliques
                                              after first iteration.
  -s <num> --significance <num>               Filter out reads with low
                                              frequency after every iteration.
                                              [default: 0.01]
)";

void usage() {
    cerr << USAGE;
    exit(1);
}

bool read_mean_and_sd(const string& filename, double* mean, double* sd) {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
    boost::char_separator<char> whitespace_separator(" \t");
    ifstream in(filename.c_str());
    string line;
    if (in.fail() || (!getline(in,line))) {
        return false;
    }
    tokenizer_t tokenizer(line,whitespace_separator);
    vector<string> tokens(tokenizer.begin(), tokenizer.end());
    if (tokens.size() != 2) {
        return false;
    }
    try {
        *mean = boost::lexical_cast<double>(tokens[0]);
        *sd = boost::lexical_cast<double>(tokens[1]);
    } catch(boost::bad_lexical_cast &){
        return false;
    }
    return true;
}

deque<AlignmentRecord*>* readBamFile(string filename, vector<string>& readNames) {
    typedef std::unordered_map<std::string, AlignmentRecord*> name_map_t;
    name_map_t names_to_reads;
    deque<AlignmentRecord*>* reads = new deque<AlignmentRecord*>;
    BamTools::BamReader bamreader;
    if (not bamreader.Open(filename)) {
        cerr << bamreader.GetFilename() << endl;
        throw std::runtime_error("Couldn't open Bamfile");
    }
    BamTools::BamAlignment alignment;

    while (bamreader.GetNextAlignment(alignment)) {
        if(names_to_reads.count(alignment.Name) > 0) {
            names_to_reads[alignment.Name]->pairWith(alignment);
            reads->push_back(names_to_reads[alignment.Name]);
            names_to_reads.erase(alignment.Name);
        } else {
            names_to_reads[alignment.Name] = new AlignmentRecord(alignment, readNames.size(), &readNames);
            readNames.push_back(alignment.Name);
        }
    }

    // Push all single-end reads remaining in names_to_reads into the reads vector
    for (const auto& i : names_to_reads) {
        reads->push_back(i.second);
    }

    bamreader.Close();

    auto comp = [](AlignmentRecord* r1, AlignmentRecord* r2) { return r1->getIntervalStart() < r2->getIntervalStart(); };

    sort(reads->begin(), reads->end(), comp);

    return reads;
}

int main(int argc, char* argv[]) {

    map<std::string, docopt::value> args
        = docopt::docopt(USAGE,
                         { argv + 1, argv + argc },
                         false);

    // PARAMETERS
    string bamfile = args["<bamfile>"].asString();
    string outfile = "quasispecies.fasta";
    if (args["<output>"]) outfile = args["<output>"].asString();
    string edge_filename = "";
    if (args["--write_edges"]) edge_filename = args["--write_edges"].asString();
    double edge_quasi_cutoff_cliques = stod(args["--edge_quasi_cutoff_cliques"].asString());
    double edge_quasi_cutoff_single = stod(args["--edge_quasi_cutoff_single"].asString());
    double edge_quasi_cutoff_mixed = stod(args["--edge_quasi_cutoff_mixed"].asString());
    double Q = stod(args["--random_overlap_quality"].asString());
    double overlap_cliques = stod(args["--min_overlap_cliques"].asString());
    double overlap_single = stod(args["--min_overlap_single"].asString());
    bool frameshift_merge = args["--frame_shift_merge"].asBool();
    string allel_frequencies_path = "";
    if (args["--allel_frequencies"]) allel_frequencies_path = args["--allel_frequencies"].asString();
    string mean_and_sd_filename = "";
    if (args["--mean_and_sd_filename"]) mean_and_sd_filename = args["--mean_and_sd_filename"].asString();
    string indel_edge_cutoff;
    double indel_edge_sig_level = stod(args["--indel_edge_sig_level"].asString());
    string indel_output_file = "";
    if (args["--call_indels"]) indel_output_file = args["--call_indels"].asString();
    int iterations = stoi(args["--iterations"].asString());
    double filter = stod(args["--filter"].asString());
    double significance = stod(args["--significance"].asString());
    bool filter_singletons = args["--no_singletons"].asBool();
    // END PARAMETERS

    bool call_indels = indel_output_file.size() > 0;
    if (call_indels && (mean_and_sd_filename.size() == 0)) {
        cerr << "Error: when using option -I, option -M must also be given." << endl;
        return 1;
    }

    //read allel frequency distributions
    std::map<int, double> simpson_map;
    //cerr << "PARSE PRIOR";
    cerr.flush();
    if (allel_frequencies_path.size() > 0) {
        ifstream ia(allel_frequencies_path.c_str());
        string ia_line;
        while (getline(ia, ia_line)) {
            std::vector<std::string> words;
            trim_right(ia_line);
            boost::split(words, ia_line, boost::is_any_of("\t"), boost::token_compress_on);

            std::vector<std::string> insertion_words;
            boost::split(insertion_words, words[0], boost::is_any_of("\\."), boost::token_compress_on);
            if (insertion_words.size() > 1) {
            } else {
                simpson_map[atoi(words[0].c_str())] = pow(atof(words[1].c_str()),2)+pow(atof(words[2].c_str()),2)+pow(atof(words[3].c_str()),2)+pow(atof(words[4].c_str()),2)+pow(atof(words[5].c_str()),2);
                //cerr << simpson_map[atoi(words[0].c_str())] << endl;
            }
        }
        ia.close();
    }
    //cerr << "PARSE PRIOR: done" << endl;


    clock_t clock_start = clock();
    EdgeCalculator* edge_calculator = nullptr;
    EdgeCalculator* indel_edge_calculator = nullptr;
    ReadGroups* read_groups = nullptr;
    unique_ptr<vector<mean_and_stddev_t> > readgroup_params(nullptr);
    edge_calculator = new QuasispeciesEdgeCalculator(Q, edge_quasi_cutoff_cliques, overlap_cliques, frameshift_merge, simpson_map, edge_quasi_cutoff_single, overlap_single, edge_quasi_cutoff_mixed);
    if (call_indels) {
        double insert_mean = -1.0;
        double insert_stddev = -1.0;
        if (!read_mean_and_sd(mean_and_sd_filename, &insert_mean, &insert_stddev)) {
            cerr << "Error reading \"" << mean_and_sd_filename << "\"." << endl;
            return 1;
        }
        cerr << "Null distribution: mean " << insert_mean << ", sd " <<  insert_stddev << endl;
        indel_edge_calculator = new GaussianEdgeCalculator(indel_edge_sig_level,insert_mean,insert_stddev);
    }
    std::ofstream* indel_os = nullptr;
    if (call_indels) {
        indel_os = new ofstream(indel_output_file.c_str());
    }
    CliqueCollector collector;
    CliqueFinder* clique_finder;
    if (args["bronkerbosch"].asBool()) {
        clique_finder = new BronKerbosch(*edge_calculator, collector, read_groups);
    } else {
        clique_finder = new CLEVER(*edge_calculator, collector, read_groups);
    }
    if (indel_edge_calculator != 0) {
        clique_finder->setSecondEdgeCalculator(indel_edge_calculator);
    }
    EdgeWriter* edge_writer = 0;
    ofstream* edge_ofstream = 0;
    if (edge_filename.size() > 0) {
        edge_ofstream = new ofstream(edge_filename.c_str());
        edge_writer = new EdgeWriter(*edge_ofstream);
        clique_finder->setEdgeWriter(*edge_writer);
    }
    ofstream* reads_ofstream = 0;

    vector<string> originalReadNames;
    deque<AlignmentRecord*>* reads = readBamFile(bamfile, originalReadNames);

    // Main loop
    int ct = 0;
    auto filter_fn = [&](unique_ptr<AlignmentRecord>& read) {
        return (ct == 1 and filter_singletons and read->getReadCount() <= 1) or (ct > 0 and read->getProbability() < significance);
    };


    while (ct != iterations) {
        clique_finder->initialize();

        while(not reads->empty()) {
            assert(reads->front() != nullptr);

            unique_ptr<AlignmentRecord> al_ptr(reads->front());

            reads->pop_front();

            if (filter_fn(al_ptr)) continue;

            clique_finder->addAlignment(al_ptr);
        }

        delete reads;

        clique_finder->finish();
        reads = collector.finish();

        setProbabilities(*reads);

        if (clique_finder->hasConverged()) break;
        cout << ct++ << ": " << reads->size() << endl;
    }

    // Filter superreads according to read probability
    if (filter > 0.0) {
        auto filter_fn = [&](AlignmentRecord* al) { return al->getProbability() < filter;};
        auto new_end_it = std::remove_if(reads->begin(), reads->end(), filter_fn);
        for (auto it = new_end_it; it != reads->end(); it++) {
            delete *it;
            *it = nullptr;
        }
        reads->erase(new_end_it, reads->end());
    }
    ofstream os(outfile, std::ofstream::out);
    setProbabilities(*reads);
    printReads(os, *reads);

    cout << "final: " << reads->size() << endl;

    for (auto&& r : *reads) {
        delete r;
    }

    delete reads;

    if (indel_os != nullptr) {
        indel_os->close();
        delete indel_os;
    }
    if (edge_calculator != nullptr) delete edge_calculator;
    if (edge_writer != nullptr) {
        delete edge_writer;
        delete edge_ofstream;
    }
    if (clique_finder != nullptr) delete clique_finder;
    if (reads_ofstream != nullptr) {
        delete reads_ofstream;
    }
    cout.precision(3);
    cout << std::fixed;
    double cpu_time = (double) (clock() - clock_start) / CLOCKS_PER_SEC;
    cout << "time: " <<  cpu_time << endl;
    return 0;
}
