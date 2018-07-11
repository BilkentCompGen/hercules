/** @file CommandLineOptions.h
 * @brief Parses command line options.
 *
 * It is used to parse command line options for both preprocessing and correction step.
 * Copyright © 2017 Can Firtina. All rights reserved.
 *
 * @author Can Firtina
 * @bug No bug currently
 */

#ifndef COMMAND_LINE_OPTIONS_H_
#define COMMAND_LINE_OPTIONS_H_

#include <iostream>
#include <vector>
#include <seqan/arg_parse.h>

#define PREPROCESS_FLAG 1
#define CORRECT_FLAG 2

/** @brief Struct for holding command line options and their values.
 *
 *
 *  @param CommandLineOptions All variables are defined below
 *  @see parseCommandOptions()
 *  @see parseInitialCommand()
 */
struct CommandLineOptions
{
    CommandLineOptions():
    filterSize(100), maxDeletion(10), maxInsertion(3), matchTransition(0.70), insertionTransition(0.25), shouldQuite(false),
    matchEmission(0.97), maxThread(1), mapQ(0), deletionTransitionFactor(2.5), shouldOutputCoverage(false), shouldCompress(true),
    maxCoverage(1), nonNCount(40), shouldApplyBloomFilter(false){}
    
    seqan::CharString longInputFile;
    seqan::CharString alignmentFile;
    std::vector<seqan::CharString> shortReads;
    seqan::CharString output;
    unsigned mapQ;
    unsigned filterSize;
    unsigned maxDeletion;
    unsigned maxInsertion;
    unsigned maxThread;
    unsigned nonNCount;
    unsigned maxCoverage;
    bool shouldCompress;
    bool shouldOutputCoverage;
    bool shouldApplyBloomFilter;
    bool shouldQuite;
    double matchTransition;
    double insertionTransition;
    double matchEmission;
    double deletionTransitionFactor;
    int phase;
};

/** @brief Parse values in order to run either preprocessing step or correction step.
 *
 *  @param options Stores parsed values
 *  @param argc Number of arguments specified while running Hercules
 *  @param argv Argument values array
 *  @param phase Preprocessing or Correction phase?
 *  @see parseInitialCommand()
 *  @return seqan::ArgumentParser::PARSE_OK if everything went well
 */
seqan::ArgumentParser::ParseResult
parseCommandOptions(CommandLineOptions& options, int argc, char const **argv, int phase){
    
    using namespace std;
    seqan::ArgumentParser parser("Hercules: A Profile HMM-based hybrid error correction algorithm for long reads");
    
    setVersion(parser, "0.1");
    setDate(parser, "November 2017");
    
    if(phase == PREPROCESS_FLAG || phase == CORRECT_FLAG){
        addOption(parser, seqan::ArgParseOption("li", "longInputFile", "fast{a,q} file which contains original long "
                                                "reads", seqan::ArgParseArgument::INPUT_FILE, "FILE"));
        setRequired(parser, "li");
    }
    
    if(phase == PREPROCESS_FLAG){
        addOption(parser, seqan::ArgParseOption("si", "shortRead", "Short reads file to align to the long reads. You "
                                                "may define as many short reads file as you wish with multiple -si "
                                                "options.", seqan::ArgParseArgument::INPUT_FILE, "FILE", true));
        setRequired(parser, "si");
        
        addOption(parser, seqan::ArgParseOption("o", "outputDir", "Preprocessing directory where the resulting files "
                                                "will be written. This directory **MUST** exist beforehand.",
                                                seqan::ArgParseArgument::OUTPUT_FILE, "FILE"));
        setRequired(parser, "o");
        
        addOption(parser, seqan::ArgParseOption("nonN", "nonN", "**Compressed** short read should have at least nonN "
                                                "many non-N characters not to be filtered out for the alignment phase.",
                                                seqan::ArgParseArgument::INPUT_FILE, "INT"));
        setDefaultValue(getOption(parser, "nonN"), options.nonNCount);
        
        addOption(parser, seqan::ArgParseOption("b", "bloomFilter", "Apply bloom filter to remove duplicates."));
        addOption(parser, seqan::ArgParseOption("nc", "noCompression", "Do not compress short reads. Reported short reads "
                                                "will only be filtered out according to nonN value and bloom filter if it is "
                                                "set. You need to set the same option in the correction phase as well."));
    }
    
    if(phase == CORRECT_FLAG){
        addOption(parser, seqan::ArgParseOption("ai", "alignmentFile", "{s,b}am file which contains alignments of short"
                                                " reads to long reads.", seqan::ArgParseArgument::INPUT_FILE, "FILE"));
        setRequired(parser, "ai");
        
        addOption(parser, seqan::ArgParseOption("si", "shortRead", "**Uncompressed** short read file created during the"
                                                " preprocessing step", seqan::ArgParseArgument::INPUT_FILE, "FILE"));
        setRequired(parser, "si");
        
        addOption(parser, seqan::ArgParseOption("o", "outputFile", "Output file to write the resulting reads",
                                                seqan::ArgParseArgument::OUTPUT_FILE, "FILE"));
        setRequired(parser, "o");
        
        addOption(parser, seqan::ArgParseOption("c", "outputCoverage", "If specified, Hercules creates another file "
                                                "within the same folder of corrected reads, which reports how much of "
                                                "a read is covered by short reads"));
        
        addOption(parser, seqan::ArgParseOption("q", "mapQ", "Minimum mapping quality for a long-short reads alignment "
                                                "to use in correction. Note that if multiple alignment specified, "
                                                "aligners report these mapping quality as 0",
                                                seqan::ArgParseArgument::INTEGER, "INT"));
        setDefaultValue(getOption(parser, "q"), options.mapQ);
        seqan::setMinValue(parser, "q", "0");
        seqan::setMaxValue(parser, "q", "255");
        
        addOption(parser, seqan::ArgParseOption("mc", "maxCoverage", "Maximum short read coverage per position of a long read. "
                                                "If provided, short reads are removed based on edit distance value. Otherwise, "
                                                " short reads are removed randomly. Setting 0 will use all short reads.",
                                                seqan::ArgParseArgument::INTEGER, "INT"));
        setDefaultValue(getOption(parser, "mc"), options.maxCoverage);
        seqan::setMinValue(parser, "mc", "0");
        
        addOption(parser, seqan::ArgParseOption("mf", "filterSize", "Filter size that allows calculation of at most mf "
                                                "many most probable transitions in each time step. This parameter is "
                                                "directly proportional to running time.",
                                                seqan::ArgParseArgument::INTEGER, "INT"));
        setDefaultValue(getOption(parser, "mf"), options.filterSize);
        seqan::setMinValue(parser, "mf", "1");
        
        addOption(parser, seqan::ArgParseOption("mi", "maxInsertion", "Maximum number of insertions in a row. This "
                                                "parameter is directly proportional to the running time.",
                                                seqan::ArgParseArgument::INTEGER, "INT"));
        setDefaultValue(getOption(parser, "mi"), options.maxInsertion);
        seqan::setMinValue(parser, "mi", "0");
        
        addOption(parser, seqan::ArgParseOption("md", "maxDeletion", "Maximum number of deletions in a row. This "
                                                "parameter is directly proportional to the running time.",
                                                seqan::ArgParseArgument::INTEGER, "INT"));
        setDefaultValue(getOption(parser, "md"), options.maxDeletion);
        seqan::setMinValue(parser, "md", "0");
        
        addOption(parser, seqan::ArgParseOption("trm", "matchTransition", "Initial transition probability to a match "
                                                "state. See --insertionTransition as well.",
                                                seqan::ArgParseArgument::DOUBLE, "FLOAT"));
        setDefaultValue(getOption(parser, "trm"), options.matchTransition);
        seqan::setMinValue(parser, "trm", "0");
        seqan::setMaxValue(parser, "trm", "1");
        
        addOption(parser, seqan::ArgParseOption("tri", "insertionTransition", "Initial transition probability to a "
                                                "insertion state. Note that: deletion transition probability = "
                                                "1 - (matchTransition + insertionTransition)",
                                                seqan::ArgParseArgument::DOUBLE, "FLOAT"));
        setDefaultValue(getOption(parser, "tri"), options.insertionTransition);
        seqan::setMinValue(parser, "tri", "0");
        seqan::setMaxValue(parser, "tri", "1");
        
        addOption(parser, seqan::ArgParseOption("df", "deletionTransitionFactor", "Factor of the polynomial "
                                                "distribution to calculate each deletion transition. Higher value "
                                                "favors less deletions.", seqan::ArgParseArgument::DOUBLE, "FLOAT"));
        setDefaultValue(getOption(parser, "df"), options.deletionTransitionFactor);
        seqan::setMinValue(parser, "df", "0");
        
        addOption(parser, seqan::ArgParseOption("emm", "matchEmission", "Initial emission probability of a match to a "
                                                "reference. Note that: mismatch emission probability = "
                                                "(1-matchEmission)/3", seqan::ArgParseArgument::DOUBLE, "FLOAT"));
        setDefaultValue(getOption(parser, "emm"), options.matchEmission);
        seqan::setMinValue(parser, "emm", "0");
        seqan::setMaxValue(parser, "emm", "1");
        
        addOption(parser, seqan::ArgParseOption("t", "thread", "Number of threads to use",
                                                seqan::ArgParseArgument::INTEGER, "INT"));
        setDefaultValue(getOption(parser, "t"), options.maxThread);
        seqan::setMinValue(parser, "t", "1");
        
        addOption(parser, seqan::ArgParseOption("nc", "noCompression", "Set this option if short and long reads are not compressed"
                                                "in the preprocessing step."));
        
        addOption(parser, seqan::ArgParseOption("nv", "noVerbose", "Hercules runs quitely with no informative output"));
    }
    
    const char** argNextOptions = argv+1;
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc-1, argNextOptions);
    if (res == seqan::ArgumentParser::PARSE_OK){
        
        options.phase = phase;
        
        if(phase == PREPROCESS_FLAG || phase == CORRECT_FLAG){
            getOptionValue(options.longInputFile, parser, "li");
        }
        
        if(phase == PREPROCESS_FLAG){
            int count = getOptionValueCount(parser, "si");
            const std::vector<std::string> shortreadvalues = getOptionValues(parser, "si");
            for(int i = 0; i < count; ++i) options.shortReads.push_back(shortreadvalues.at(i));
            getOptionValue(options.output, parser, "o");
            getOptionValue(options.nonNCount, parser, "nonN");
            options.shouldApplyBloomFilter = isSet(parser, "b");
            options.shouldCompress = !isSet(parser, "nc");
        }else if(phase == CORRECT_FLAG){
            getOptionValue(options.alignmentFile, parser, "ai");
            seqan::CharString shortFileStr;
            getOptionValue(shortFileStr, parser, "si");
            options.shortReads.push_back(shortFileStr);
            getOptionValue(options.output, parser, "o");
            options.shouldOutputCoverage = isSet(parser, "c");
            getOptionValue(options.mapQ, parser, "q");
            getOptionValue(options.maxCoverage, parser, "mc");
            getOptionValue(options.filterSize, parser, "mf");
            getOptionValue(options.maxInsertion, parser, "mi");
            getOptionValue(options.maxDeletion, parser, "md");
            getOptionValue(options.matchTransition, parser, "trm");
            getOptionValue(options.insertionTransition, parser, "tri");
            getOptionValue(options.deletionTransitionFactor, parser, "df");
            getOptionValue(options.matchEmission, parser, "emm");
            getOptionValue(options.maxThread, parser, "t");
            options.shouldCompress = !isSet(parser, "nc");
            options.shouldQuite = isSet(parser, "nv");
            
            if(options.matchTransition + options.insertionTransition > 1){
                std::cerr << "ERROR: (matchTransition + insertionTransition) cannot be larger than 1 whereas the sum "
                << "is now: " << options.matchTransition + options.insertionTransition << std::endl;
                return seqan::ArgumentParser::PARSE_ERROR;
            }
        }
    }
    
    return res;
}

/** @brief Determines if the user intends to run a preprocessing step or a correction step
 *
 *  @param options Stores parsed values
 *  @param argc Number of arguments specified while running Hercules
 *  @param argv Argument values array
 *  @see parseCommandOptions()
 *  @return seqan::ArgumentParser::PARSE_OK if everything went well
 */
seqan::ArgumentParser::ParseResult
parseInitialCommand(CommandLineOptions& options, int argc, char const** argv){
    
    using namespace std;
    seqan::ArgumentParser parser("Hercules: A Profile HMM-based hybrid error correction algorithm for long reads");
    
    setVersion(parser, "0.1");
    setDate(parser, "November 2017");

    addOption(parser, seqan::ArgParseOption("1", "preprocess", "Compresses the required reads and creates an proper "
                                            "fasta files and an index file for the long read. Created reads should be "
                                            "provided to a aligner"));
    addOption(parser, seqan::ArgParseOption("2", "correct", "Corrects the long reads using the alignment and index file"
                                            " created in the preprocess step"));
    
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, (argc > 2)?2:argc, argv);
    if (res == seqan::ArgumentParser::PARSE_OK){
        int phase = 0;
        if(isSet(parser, "preprocess"))
            phase |= PREPROCESS_FLAG;
        if(isSet(parser, "correct"))
            phase |= CORRECT_FLAG;
        
        if(phase == 0){
            std::cerr << "ERROR: You should specify a command" << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
        }else if(phase != PREPROCESS_FLAG && phase != CORRECT_FLAG){
            std::cerr << "ERROR: You cannot specify multiple commands or an unspecified command" << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
        }
        return parseCommandOptions(options, argc, argv, phase);
    }
    return res;
}

#endif
