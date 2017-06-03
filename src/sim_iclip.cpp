// ======================================================================
// PureCLIP: capturing target-specific protein-RNA interaction footprints
// ======================================================================
// Copyright (C) 2017  Sabrina Krakau, Max Planck Institute for Molecular 
// Genetics
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Sabrina Krakau <krakau@molgen.mpg.de>
// ==========================================================================

#define SIM_ICLIP_PROFILE

#define SIM_ICLIP_ENABLE_PARALLELISM 1

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <iostream>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

#include <seqan/arg_parse.h>
#include <seqan/parallel.h>

#include "sim_iclip.h"

/*
~/Development/seqan-2.2/build/release/bin/sim_iclip -bam STAR/rep1/Aligned.R1.mapq3.bam -bai STAR/rep1/Aligned.R1.mapq3.bam.bai -ref /project/lincRNA_seq/bimorphic/genome/hg19_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -bs PUM2_FIMO/fimo_matches/fimo.so.intersectTranscript.d100.final.bed -bbs /project/lincRNA_seq/iCLIP/backCLIP/BackgroundTraining_19datasets.so.intersectTranscript.d100.final.bed -out sim/sim.SIM_SET4.tmp.openMP_test.bam -pd -ubs -baf 0.1 > sim/sim_iclip.SIM_SET4.openMP_test.log

Mod:
stream/iostream_bgzf.
numThreads = 1 

*/


using namespace seqan;

struct AppOptions
{
    CharString bamFileName;
    CharString baiFileName;
    CharString refFileName;
    CharString bsFileName;      // input binding sites
    CharString bbsFileName;      // input background binding sites
    CharString outFileName;
    CharString fldFileName;     // fragment length distr. 

    bool pullDown;
    unsigned noCrosslinkSites;  // for one binding site
    unsigned useRnCrosslinkSites;
    unsigned maxFragmentSize;
    unsigned minFragmentSize;
    unsigned fragmentSize;
    unsigned mincDNAlength;
    unsigned maxcDNAlength;
    unsigned maxReadLength;
    unsigned simReadLength;
    double bindingAffinity;
    double truncationRate;
    double otherTruncationRate;
    bool noRTStopAtBs;
    unsigned maxBSWidth;
    bool useBgBS;
    unsigned noBgCrosslinkSites;
    double bindingAffFactor;
    double bgTruncationRate;

    bool ignoreTargetSites;
    unsigned numThreads;
    CharString tempPath;
    int verbosity;

    AppOptions() :
        pullDown(true),             // if false, assume given pull down already
        noCrosslinkSites(4),  
        useRnCrosslinkSites(false),
        maxFragmentSize(11000),    // increases runtime dramatically -> ATTENTIONE: tradeoff  ("< 10% of introns are more than 11,000 bp in length", Sakharkar MK et al.)
        minFragmentSize(50),
        fragmentSize(75),
        mincDNAlength(30),
        maxcDNAlength(140),
        maxReadLength(50),
        simReadLength(50),
        bindingAffinity(1.0),
        truncationRate(0.7),
        otherTruncationRate(0.1),   // early RT stop downstream of crosslink sites due to other reasons (other remaining footprints) 
        noRTStopAtBs(true),        // no RT stop within binding site due to other reasons (than crosslink sites)
        maxBSWidth(8),
        useBgBS(false),
        noBgCrosslinkSites(10), 
        bindingAffFactor(1.0/10.0),
        bgTruncationRate(0.6),
        ignoreTargetSites(false),
        numThreads(1),
        verbosity(1)
    {}
};


ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("sim_iclip");
    // Set short description, version, and date.
    setShortDescription(parser, "Simulate iCLIP/eCLIP-seq data ");
    setVersion(parser, "0.1");
    setDate(parser, "Mai 2017");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] <-bam \\fIBAM FILE\\fP> <-bai \\fIBAI FILE\\fP> <-ref \\fIGENOME FILE\\fP>  <-bs \\fIBINDING REGIONS FILE\\fP> <-out \\fIOUTPUT BAM FILE\\fP> ");
    addDescription(parser, "Simulate iCLIP/eCLIP-seq data based on aligned RNA-seq data using a given list of binding sites.");

    // Files
    addOption(parser, ArgParseOption("b", "bam", "Input bam file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "bam", ".bam");
    setRequired(parser, "bam", true);

    addOption(parser, ArgParseOption("i", "bai", "Input bai file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "bai", ".bai");
    setRequired(parser, "bai", true);

    addOption(parser, ArgParseOption("r", "ref", "Input reference file.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "ref", ".fa .fasta");
    setRequired(parser, "ref", true);

    addOption(parser, ArgParseOption("bs", "bs", "Input file containing binding sites.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "bs", ".bed");
    setRequired(parser, "bs", true);

    addOption(parser, ArgParseOption("o", "out", "Output file containing reads.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "out", ".bam");
    setRequired(parser, "out", true);

    addSection(parser, "Options for simulating target signal");

    addOption(parser, ArgParseOption("fld", "fld", "File containing empirical fragment length distr.", ArgParseArgument::OUTPUT_FILE));
    setValidValues(parser, "fld", ".txt");
    addOption(parser, ArgParseOption("tr", "tr", "Truncation rate for target sites. Default: 0.7.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "tr", "0.0");
    setMaxValue(parser, "tr", "1.0");
    addOption(parser, ArgParseOption("cs", "cs", "Number of crosslink sites to simulate within each binding region. Default: 4.", ArgParseArgument::INTEGER));
    setMinValue(parser, "cs", "1");
    setMaxValue(parser, "cs", "6");
    addOption(parser, ArgParseOption("urc", "urc", "Use random number of crosslink sites for each binding site in range [1..cs]."));
    addOption(parser, ArgParseOption("fs", "fs", "Use fragment size (ignored if empirical distribution file is used).", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("ba", "ba", "Binding affinity for pull-down. Default: 1.0.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "ba", "0.005");
    setMaxValue(parser, "ba", "1.0");


    // Background noise
    addSection(parser, "Options for simulating background noise");

    addOption(parser, ArgParseOption("bbs", "bbs", "Input file containing background binding sites.", ArgParseArgument::INPUT_FILE));
    setValidValues(parser, "bbs", ".bed");
    addOption(parser, ArgParseOption("bcs", "bcs", "Number of crosslink sites to simulate within each background binding region, in range [cs..bcs].", ArgParseArgument::INTEGER));
    setMinValue(parser, "bcs", "1");
    setMaxValue(parser, "bcs", "50");
    addOption(parser, ArgParseOption("baf", "baf", "Factor to multiply background binding region score with to obtain binding affinity for pull-down.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "baf", "0.005");
    setMaxValue(parser, "baf", "0.2");
    addOption(parser, ArgParseOption("btr", "btr", "Truncation rate for background binding regions.", ArgParseArgument::DOUBLE));
    setMinValue(parser, "btr", "0.0");
    setMaxValue(parser, "btr", "1.0");
    addOption(parser, ArgParseOption("its", "its", "Ignore target sites (produce only background binding reads)."));


    addSection(parser, "General user options");

    addOption(parser, ArgParseOption("nt", "nt", "Thread number.", ArgParseArgument::INTEGER));
    addOption(parser, ArgParseOption("tmp", "tmp", "Path to directory to store intermediate files. ", ArgParseArgument::STRING));
    addOption(parser, ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, ArgParseOption("vv", "very-verbose", "Enable very verbose output."));
    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBsim_iclip\\fP \\fB-bam rna-seq.bam\\fP \\fB-bai rna-seq.bai\\fP \\fB-ref ref.fasta\\fP \\fB-bs binding_regions.bed\\fP \\fB-out sim.bam\\fP \\fB-nt 10\\fP ",
                "");
    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;


    getOptionValue(options.bamFileName, parser, "bam");
    getOptionValue(options.baiFileName, parser, "bai");
    getOptionValue(options.refFileName, parser, "ref");
    getOptionValue(options.bsFileName, parser, "bs");
    getOptionValue(options.bbsFileName, parser, "bbs");
    if (!empty(options.bbsFileName))
        options.useBgBS = true;
    getOptionValue(options.outFileName, parser, "out");
    getOptionValue(options.fldFileName, parser, "fld");

    getOptionValue(options.bindingAffinity, parser, "ba");
    getOptionValue(options.bindingAffFactor, parser, "baf");

    getOptionValue(options.truncationRate, parser, "tr");
    getOptionValue(options.noCrosslinkSites, parser, "cs");
    if (isSet(parser, "urc"))
        options.useRnCrosslinkSites = true;

    getOptionValue(options.fragmentSize, parser, "fs");

    getOptionValue(options.noBgCrosslinkSites, parser, "bcs");
    getOptionValue(options.bgTruncationRate, parser, "btr");

    if (isSet(parser, "its"))
        options.ignoreTargetSites = true;
    getOptionValue(options.numThreads, parser, "nt");
    getOptionValue(options.tempPath, parser, "tmp");
    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    return ArgumentParser::PARSE_OK;
}



int main(int argc, char const ** argv)
{
    // Parse the command line.
    ArgumentParser parser;
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    std::cout << "Simulate iCLIP/eCLIP-seq reads \n"
              << "===============\n\n";

    // Print the command line arguments back to the user.
    /*if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY\t" << options.verbosity << "\n\n";
    }*/


#ifdef _OPENMP
    omp_set_num_threads(options.numThreads);  
    std::cout << " Set number of openMP threads to " <<  options.numThreads << std::endl;
    //std::cout << " Get max. threads " <<  omp_get_max_threads()  << std::endl;
#endif

#ifdef SIM_ICLIP_PROFILE
    double timeStamp = sysTime();
#endif

    doIt(options);

#ifdef SIM_ICLIP_PROFILE
    Times::instance().time_all = sysTime() - timeStamp;
    std::cout << "  Time needed for all: " << Times::instance().time_all/60.0 << "min" << std::endl;
#endif

    return 0;
}






