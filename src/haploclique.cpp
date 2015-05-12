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
#include <limits>
#include <cassert>
#include <ctime>

#include "docopt/docopt.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <bamtools/api/BamReader.h>

#include "AlignmentRecord.h"
#include "QuasispeciesEdgeCalculator.h"
#include "CliqueFinder.h"
#include "CLEVER.h"
#include "BronKerbosch.h"
#include "CliqueWriter.h"
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

Usage: haploclique (--clever | --bronkerbosch) [options]

Reads alignments / alignment pairs from stdin and computes all cliques.
Format for single-end reads:
<read-name> <read-nr> <read-group> <phred-sum> <chromosome> <start> <end> <strand> <cigar> <seq> <qualities> <aln-prob>
Format for paired-end reads:
<read-name> <pair-nr> <read-group> <phred-sum1> <chromosome1> <start1> <end1> <strand1> <cigar1> <seq1> <qualities1> <phred-sum2> <chromosome2> <start2> <end2> <strand2> <cigar2> <seq2> <qualities2> <aln-pair-prob> <aln-pair-prob-inslength>

NOTE: Lines are assumed to be ordered by field 6 (start1).

--clever            Use the original CLEVER algorithm for enumerating maximal
                    cliques.
--bronkerbosch      Use the Bron-Kerbosch algorithm for enumerationg maximal
                    cliques.

<distribution-file> file with assumed internal segment length distribution.
                    in default mode, this is a file containing one line with
                    mean and standard deviation of the normal distribution
                    to be used. Such a file can be generated using option -m
                    of insert-length-histogram. Note that the
                    "internal segment" does NOT include the read (ends), i.e.
                    a fragment disjointly consists of two reads (read ends)
                    and an internal segment.

Options:
  -v --verbose                Be verbose: output additional statistics for
                              each variation.
  -w --min_aln_weight NUM     Minimum weight of alignment pairs to be
                              considered. [default: 0.0016]
  -l --max_insert_length NUM  Maximum insert length of alignments to be
                              considered (0=unlimited). [default: 50000]
  -c --max_coverage NUM       Maximum allowed coverage. If exceeded,
                              violating reads are ignored. The number of such
                              ignored reads is printed to stderr
                              (0=unlimited). [default: 200]
  -e --write_edges FILE       Write edges to file of given name.
  -f --fdr NUM                False discovery rate (FDR). [default: 0.1]
  -a --all                    Output all cliques instead of only the
                              significant ones. Cliques are not sorted and
                              last column (FDR) is not computed.
  -r --output_reads FILE      Output reads belonging to at least one
                              significant clique to the given filename (along
                              with their assignment to significant cliques).
  -C --output_coverage FILE   Output the coverage with considered insert
                              segments along the chromosome (one line per
                              position) to the given filename.
  -q --edge_quasi_cutoff_cliques NUM  End compatibility probability cutoff
                                      between error-corrected reads for
                                      quasispecies reconstruction.
                                      [default: 0.99]
  -k --edge_quasi_cutoff_mixed NUM    End compatibility probability cutoff
                                      between raw<->error-corrected reads for
                                      quasispecies reconstruction.
                                      [default: 0.97]
  -g --edge_quasi_cutoff_single NUM   End compatibility probability cutoff
                                      between raw<->raw reads for
                                      quasispecies reconstruction.
                                      [default:0.95]
  -Q --random_overlap_quality NUM     Probability that two random reads are
                                      equal at the same position.
                                      [default: 0.9]
  -m --frame_shift_merge              Reads will be clustered if one has
                                      single nucleotide deletions and
                                      insertions. Use for PacBio data sets.
  -o --min_overlap_cliques NUM        Minimum relative overlap between
                                      error-corrected reads. [default: 0.9]
  -j --min_overlap_single NUM         Minimum relative overlap between
                                      raw<->raw and raw<->error-corrected
                                      reads. [default: 0.6]
  -s --super_read_min_coverage NUM    Minimum coverage for super-read
                                      assembly. [default: 2]
  -A --allel_frequencies FILE         ???
  -I --call_indels FILE               Call indels from cliques. In this mode,
                                      the classical "CLEVER" edge criterion
                                      is used in addition to the new one.
                                      Filename to write indels to must be
                                      given as parameter.
  -M --mean_and_sd_filename FILE      Name of file with mean and standard
                                      deviation of insert size distribution
                                      (only required if option -I is used).
  -p --indel_edge_sig_level NUM       Significance level for "indel" edges
                                      criterion, see option -I (the lower the
                                      level, the more edges will be present).
                                      [default: 0.2]
  -t --time_limit NUM         Time limit for computation. If
                              exceeded, non processed reads will be
                              written to skipped.
  -N --no_sort                Do not sort new clique w.r.t. their bitsets.
  -S --suffix                 Suffix for clique names. Used for
                              parallelization.
  -L --minimal_super_read_length NUM
                              Minimal super read length. [default: 0]  
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

int main(int argc, char* argv[]) {

    map<std::string, docopt::value> args
        = docopt::docopt(USAGE,
                         { argv + 1, argv + argc },
                         true);
    // PARAMETERS
    double min_aln_weight = stod(args["--min_aln_weight"].asString());
    double max_insert_length = stod(args["--max_insert_length"].asString());
    long max_coverage = args["--max_coverage"].asLong();
    string edge_filename = "";
    if (args["--write_edges"]) edge_filename = args["--write_edges"].asString();
    double fdr = stod(args["--fdr"].asString());
    string reads_output_filename = "";
    if (args["--output_reads"]) reads_output_filename = args["--output_reads"].asString();
    string coverage_output_filename = "";
    if (args["--output_coverage"]) coverage_output_filename = args["--output_coverage"].asString();
    bool verbose = args["--verbose"].asBool();
    double edge_quasi_cutoff_cliques = stod(args["--edge_quasi_cutoff_cliques"].asString());
    double edge_quasi_cutoff_single = stod(args["--edge_quasi_cutoff_single"].asString());
    double edge_quasi_cutoff_mixed = stod(args["--edge_quasi_cutoff_mixed"].asString());
    double Q = stod(args["--random_overlap_quality"].asString());
    double overlap_cliques = stod(args["--min_overlap_cliques"].asString());
    double overlap_single = stod(args["--min_overlap_single"].asString());
    bool frameshift_merge = args["--frame_shift_merge"].asBool();
    long super_read_min_coverage = args["--super_read_min_coverage"].asLong();
    string allel_frequencies_path = "";
    if (args["--allel_frequencies"]) allel_frequencies_path = args["--allel_frequencies"].asString();
    string mean_and_sd_filename = "";
    if (args["--mean_and_sd_filename"]) mean_and_sd_filename = args["--mean_and_sd_filename"].asString();
    string indel_edge_cutoff;
    double indel_edge_sig_level = stod(args["--indel_edge_sig_level"].asString());
    string indel_output_file = "";
    if (args["--call_indels"]) indel_output_file = args["--call_indels"].asString();
    long time_limit = args["--time_limit"].asLong();
    bool no_sort = args["--no_sort"].asBool();
    string suffix = "";
    if (args["--suffix"]) suffix = args["--suffix"].asString();
    long minimal_superread_length = args["--minimal_superread_length"].asLong();
    // END PARAMETERS

    if (isatty(fileno(stdin))) {
        usage();
    }
    bool output_all = args["--all"].asBool();
    if (output_all && (reads_output_filename.size() > 0)) {
        cerr << "Error: options -a and -r are mutually exclusive." << endl;
        return 1;
    }
    bool call_indels = indel_output_file.size() > 0;
    if (call_indels && (mean_and_sd_filename.size() == 0)) {
        cerr << "Error: when using option -I, option -M must also be given." << endl;
        return 1;
    }

    std::map<string,string> clique_to_reads_map;
    ifstream tsv_stream("data_clique_to_reads.tsv");
    string tsv_stream_line;
    while (getline(tsv_stream, tsv_stream_line)) {
        std::vector<std::string> words;
        trim_right(tsv_stream_line);
        boost::split(words, tsv_stream_line, boost::is_any_of("\t"), boost::token_compress_on);
        clique_to_reads_map[words[0]] = words[1];
    }
    remove("data_clique_to_reads.tsv");

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
    EdgeCalculator* edge_calculator = 0;
    EdgeCalculator* indel_edge_calculator = 0;
    ReadSetSignificanceTester* significance_tester = 0;
    VariationCaller* variation_caller = 0;
    ReadGroups* read_groups = 0;
    auto_ptr<vector<mean_and_stddev_t> > readgroup_params(0);
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
        significance_tester = new ReadSetZTester(insert_mean, insert_stddev);
        variation_caller = new VariationCaller(insert_mean, *significance_tester);
    }
    std::ofstream* indel_os = 0;
    if (call_indels) {
        indel_os = new ofstream(indel_output_file.c_str());
    }
    CliqueWriter clique_writer(cout, variation_caller, indel_os, read_groups, false, output_all, fdr, verbose, super_read_min_coverage, frameshift_merge, suffix, minimal_superread_length);
    CliqueFinder* clique_finder = new CLEVER(*edge_calculator, clique_writer, read_groups, no_sort);
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
    if (reads_output_filename.size() > 0) {
        reads_ofstream = new ofstream(reads_output_filename.c_str());
        clique_writer.enableReadListOutput(*reads_ofstream);
    }

    size_t last_pos = 0;
    int n = 0;
    string line;
    // size_t skipped_by_weight = 0;
    // size_t skipped_by_length = 0;
    size_t skipped_by_coverage = 0;
    size_t valid_alignments = 0;
    size_t total_alignments = 0;
    cerr << "STATUS";
    cerr.flush();
    while (getline(cin, line)) {
        n += 1;
        total_alignments += 1;
        try {
            AlignmentRecord ap(line, clique_to_reads_map, read_groups);
            if (ap.getIntervalStart() < last_pos) {
                cerr << "Error: Input is not ordered by position (field 6)! Offending line: " << n << endl;
                return 1;
            }
            if (ap.isPairedEnd()) {
                if (ap.getChrom1().compare(ap.getChrom2()) != 0) continue;
                if (ap.getStrand1().compare(ap.getStrand2()) == 0) continue;
            }
            valid_alignments += 1;
            last_pos = ap.getIntervalStart();
            auto_ptr<AlignmentRecord> alignment_autoptr(new AlignmentRecord(ap));
            if (ap.isPairedEnd() && (max_insert_length > 0)) {
                if (alignment_autoptr->getInsertLength() > max_insert_length) {
                //skipped_by_length += 1;
                //continue;
                }
            }
            /*if (alignment_autoptr->getWeight() < min_aln_weight) {
            // cout << "Skipping alignment (weight): "  << alignment_autoptr->getName() << " weight: " << alignment_autoptr->getWeight() << endl;
                skipped << line;
                skipped_by_weight += 1;
                continue;
            }*/
            bool time = ((double) (clock() - clock_start) / CLOCKS_PER_SEC / 60) > time_limit;
            if (max_coverage > 0 || time) {
                if (clique_finder->getCoverageMonitor().probeAlignment(*alignment_autoptr) > (size_t) max_coverage || time) {
                // cout << "Skipping alignment (coverage): "  << alignment_autoptr->getName()  << endl;
                    skipped_by_coverage += 1;
                    continue;
                }
            }
            clique_finder->addAlignment(alignment_autoptr);
        } catch (std::runtime_error&) {
            cerr << "Error parsing input, offending line: " << n << endl;
            return 1;
        }
        cerr << "\rSTATUS: " << total_alignments;
        cerr.flush();
    }
    clique_finder->finish();
    clique_writer.finish();
    cerr << endl;

    if (indel_os != 0) {
        indel_os->close();
        delete indel_os;
    }
    if (edge_calculator != 0) delete edge_calculator;
    if (variation_caller != 0) delete variation_caller;
    if (significance_tester != 0) delete significance_tester;
    if (edge_writer != 0) {
        delete edge_writer;
        delete edge_ofstream;
    }
    if (clique_finder != 0) delete clique_finder;
    if (reads_ofstream != 0) {
        delete reads_ofstream;
    }
    double cpu_time = (double) (clock() - clock_start) / CLOCKS_PER_SEC;
    cerr << "Cliques/Uniques/CPU time:\t" << clique_writer.getPairedCount() << "/" << clique_writer.getSingleCount() << "/" << round(cpu_time) << endl;
    return 0;
}
