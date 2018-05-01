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

#ifndef APPS_SIM_ICLIP_SIM_ICLIP_H_
#define APPS_SIM_ICLIP_SIM_ICLIP_H_


#include <iostream>
#include <fstream>
#include <seqan/random.h>
#include <seqan/bed_io.h>

using namespace seqan;


class Times
{
public:
    double time_all;

    static Times & instance()
    {
        static Times times;
        return times;
    }

private:
    Times() :
        time_all(0)
    {}
};


struct BindingSite
{
    unsigned beginPos;
    unsigned endPos;

    char strand;        // '+' / '-'

    double bindingAffinity;

    String<unsigned> crosslinkSites;    // within binding site
    String<double> crosslinkProbs;      // sum  up to one?

    BindingSite() :
        beginPos(0),
        endPos(0),
        strand('.'),
        bindingAffinity(0.0)
    {}
};

struct Stats
{
    unsigned bsCount;
    unsigned uncoveredBsCount;
    unsigned readCount;
    unsigned targetTruncCount;
    unsigned offTargetTruncCount;
    unsigned randomReadCount;

    Stats() :
        bsCount(0),
        uncoveredBsCount(0),
        readCount(0),
        targetTruncCount(0),
        offTargetTruncCount(0),
        randomReadCount(0)       
    {}
};


template <typename TSites, typename TStore, typename TOptions>
bool loadBindingSites(TSites &bindingSites, TStore &store, CharString fileName, bool const& useScore, TOptions &options)
{
    resize(bindingSites, length(store.contigNameStore));

    BedFileIn bedIn(toCString(fileName));
    BedRecord<Bed6> bedRecord;      

    unsigned contigId = 0;
    while (!atEnd(bedIn))
    {
        try
        {
            readRecord(bedRecord, bedIn);
        }
        catch (ParseError const & e)
        {
            std::cerr << "ERROR: input BED record is badly formatted. " << e.what() << "\n";
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy input BED record. " << e.what() << "\n";
        } 

        if (bedRecord.ref != store.contigNameStore[contigId])   // assuming same contig order ...
        {
            while (contigId < length(store.contigNameStore) && bedRecord.ref != store.contigNameStore[contigId]) ++contigId;
        }

        if (contigId < length(store.contigNameStore)) 
        {
            BindingSite b;
            b.beginPos = bedRecord.beginPos;
            b.endPos = bedRecord.endPos;
            // limit bs width (e.g. to target motif width) 
            if ((b.endPos - b.beginPos) > options.maxBSWidth) 
            {
                unsigned d = (b.endPos - b.beginPos) - options.maxBSWidth;
                b.beginPos = b.beginPos + std::floor((double)d/2.0);
                b.endPos = b.endPos - std::ceil((double)d/2.0);
            }
            b.strand = bedRecord.strand;
           
            if (useScore)
            {
                std::stringstream ss(toCString(bedRecord.score)); 
                int i;
                ss >> i;
                b.bindingAffinity = i * options.bindingAffFactor;         
                b.bindingAffinity = std::min(b.bindingAffinity, 1.0);
                b.bindingAffinity = std::max(b.bindingAffinity, 0.0);
            }
            else
            {
                b.bindingAffinity = options.bindingAffinity;
            }
            appendValue(bindingSites[contigId], b, Generous());                
        }
    } 
    return 0;
}


bool isEl(unsigned i, String<unsigned> list)
{
    for (unsigned j = 0; j < i; ++j)
    {
        if (list[i] == list[j]) 
        {
            return true;
        }
    }
    return false;
}

template<typename TRng, typename TOptions>
bool simulatePotentialCrosslinkSites(BindingSite &bindingSite, unsigned c, TRng &rng, TOptions &/*options*/)
{
    resize(bindingSite.crosslinkSites, c);
    resize(bindingSite.crosslinkProbs, c);

    // pick random sites and uniform probs for the beginning 
    for (unsigned i = 0; i < c; i++)
    {
        std::uniform_int_distribution<unsigned> dist(bindingSite.beginPos, bindingSite.endPos - 1);
        bindingSite.crosslinkSites[i] = dist(rng);
        // make sure to have c different positions
        while (isEl(i, bindingSite.crosslinkSites))
            bindingSite.crosslinkSites[i] = dist(rng);

        bindingSite.crosslinkProbs[i] = 1.0/((double)c);                
    }

    // TODO sequence bias ???
    // given base frequency at crosslink sites
    // draw for each crosslink site base type
    // then draw position out of all x-bases

    // assign individual prob. for each fragment to be crosslinked here
    
    // NOTE: base type influences crosslinking prob. / no. of crosslinks!
    // just choosing crosslink position dependent on base type not enough here!
    // draw potential crosslink sites with probability based 

    return 0;
}


template<typename TOptions>
bool adjustBamRecord(BamAlignmentRecord &newRecord, unsigned readLen, unsigned fragLen, TOptions &options)
{
    if (fragLen < options.mincDNAlength || fragLen > options.maxcDNAlength)
        return false;

    // make sure sequence cannot be used for mapping accidentally later
    clear(newRecord.seq);
    resize(newRecord.seq, readLen, 'A', Exact());

    clear(newRecord.qual);
    resize(newRecord.qual, readLen, '<', Exact());        // phred 60

    // set to only matches
    clear(newRecord.cigar);
    CigarElement<> cigarEl;
    cigarEl.operation = 'M';
    cigarEl.count = readLen;
    appendValue(newRecord.cigar, cigarEl);

    return true;
}


template<typename TRng, typename TOptions>
bool apply_iCLIPmodifications(Stats &stats, BamAlignmentRecord &newRecord, BamAlignmentRecord &record, BindingSite &bindingSite, unsigned fragLen, TRng &rng, bool isTarget, TOptions &options)
{
    newRecord = record;
    unsigned fragEndPos;    // 3'end
    if (!hasFlagRC(record))
        fragEndPos = record.beginPos + fragLen;
    else
        fragEndPos = record.beginPos + getAlignmentLengthInRef(record) - fragLen;

    // early RT stop caused by other footprints (uniform distribution)
    std::uniform_real_distribution<double> distEvent(0, 1);
    double x = distEvent(rng);
    if (x < options.otherTruncationRate)
    {
        // draw position within fragment:
        // if requested, exclude binding site and make sure, remaining fragment has length >= read length !

        // get bs length to subtract from sampling length (if bs overlaps with area to sample from)
        unsigned bsBuffer = 0;
        if (options.noRTStopAtBs && !isTarget)
        {
            if (!hasFlagRC(record))
            {
                if (bindingSite.endPos <= fragEndPos)
                    bsBuffer = bindingSite.endPos - bindingSite.beginPos;
                else if (bindingSite.beginPos <= fragEndPos) 
                    bsBuffer = record.beginPos + fragLen - bindingSite.beginPos;  
            }
            else
            {
                if (bindingSite.beginPos >= fragEndPos)  
                    bsBuffer = bindingSite.endPos - bindingSite.beginPos;
                else if (bindingSite.endPos >= fragEndPos)  
                    bsBuffer = bindingSite.endPos - fragEndPos;           
            }
        }      
        // draw position using length of area to sample from
        std::uniform_int_distribution<int> distPos(0, fragLen - bsBuffer); // NOTE use full length (extend over current length if necessary)
        unsigned x1 = distPos(rng);
        // project position to fragment
        if (!hasFlagRC(record))
        {
            if (x1 < (bindingSite.beginPos - record.beginPos))
                newRecord.beginPos = record.beginPos + x1;
            else
                newRecord.beginPos = record.beginPos + bsBuffer + x1;

            unsigned newFragLen = fragEndPos - newRecord.beginPos;
            unsigned readLen = std::min(options.simReadLength, newFragLen);
            // shorten sequences, qualities, adjust cigar, sequence invalid "AAAAA..."
            if(adjustBamRecord(newRecord, readLen, newFragLen, options))
            {
                SEQAN_OMP_PRAGMA(atomic)
                ++stats.offTargetTruncCount;
                return true; 
            }
            else
            {
                return false;
            }
        }      
        else
        {
            unsigned newEndPos; // read end pos
            if (x1 < (record.beginPos + getAlignmentLengthInRef(record) - bindingSite.endPos))  
                newEndPos = record.beginPos + getAlignmentLengthInRef(record) - x1;  
            else
                newEndPos = record.beginPos + getAlignmentLengthInRef(record) - bsBuffer - x1;   

            if (newEndPos - fragEndPos > options.simReadLength)             // still within fragment
                newRecord.beginPos = newEndPos - options.simReadLength;
            else                                                            // at end of fragment, shorter read
                newRecord.beginPos = fragEndPos;

            unsigned newFragLen = newEndPos - fragEndPos;
            unsigned readLen = newEndPos-newRecord.beginPos;
            // shorten sequences, qualities, adjust cigar, sequence invalid "AAAAA..."
            if (adjustBamRecord(newRecord, readLen, newFragLen, options))
            {
                SEQAN_OMP_PRAGMA(atomic)
                ++stats.offTargetTruncCount;
                return true; 
            }
            else
            {
                return false;
            }
        }      
    }

    // draw crosslink site(s) 
    unsigned crosslinkSite = bindingSite.crosslinkSites[0];
    double cumProb = 0;
    for (unsigned i = 0; i < length(bindingSite.crosslinkProbs); ++i) 
        cumProb += bindingSite.crosslinkProbs[i];
    std::uniform_real_distribution<double> distCS(0, cumProb);
    double x2 = distCS(rng); 

    cumProb = 0;
    for (unsigned i = 0; i < length(bindingSite.crosslinkProbs); ++i)
    {
        cumProb += bindingSite.crosslinkProbs[i];
        if (x2 < cumProb)        
        {
            crosslinkSite = bindingSite.crosslinkSites[i];
            break;
        }
    }

    // truncate with given prob.
    double tr;
    if (isTarget)
        tr = options.truncationRate;
    else
        tr = options.bgTruncationRate;

    if (x < options.otherTruncationRate + tr) 
    {
        if (!hasFlagRC(record) && ((signed)crosslinkSite > newRecord.beginPos)) // if no RT protection by protein, make sure read only shifted if cs downstream of start
        {
            newRecord.beginPos = crosslinkSite + 1;  

            unsigned newFragLen = fragEndPos - newRecord.beginPos;
            unsigned readLen = std::min(options.simReadLength, newFragLen);
            // shorten sequences, qualities, adjust cigar, sequence invalid "AAAAA..."
            if (adjustBamRecord(newRecord, readLen, newFragLen, options))    
            {
                SEQAN_OMP_PRAGMA(atomic)
                ++stats.targetTruncCount;
                return true; 
            }
            else
            {
                return false;
            }
        }
        else if (hasFlagRC(record) && (crosslinkSite < newRecord.beginPos + getAlignmentLengthInRef(record)))
        {
            unsigned newEndPos = crosslinkSite; 

            if (newEndPos - fragEndPos > options.simReadLength)             // still within fragment
                newRecord.beginPos = newEndPos - options.simReadLength;
            else                                                            // at end of fragment, shorter read
                newRecord.beginPos = fragEndPos;

            unsigned newFragLen = newEndPos - fragEndPos;
            unsigned readLen = newEndPos-newRecord.beginPos;
            // shorten sequences, qualities, adjust cigar, sequence invalid "AAAAA..."
            if (adjustBamRecord(newRecord, readLen, newFragLen, options))
            {
                SEQAN_OMP_PRAGMA(atomic)
                ++stats.targetTruncCount;
                return true; 
            }
            else
            {
                return false;
            }
        }
        // if read does not cover crosslink site (but selected for truncation), keep untruncated !?
    }
    
    // adjust length if fragment shorter than read
    unsigned readLen = std::min(options.simReadLength, fragLen);
    if (hasFlagRC(record))          // adjust 3' end of read if necessary
    {
        newRecord.beginPos = record.beginPos + getAlignmentLengthInRef(record) - readLen;
    }
    // shorten sequences, qualities, adjust cigar, sequence invalid "AAAAA..."
    return adjustBamRecord(newRecord, readLen, fragLen, options);                   // Note: size selection: ratio skewed! more loikely to be discarded than truncated cDNAs! 
    //return true;                                                                  // adjust? ensure final reads have 0.1 untruncated? TODO
}


template <typename TBedOut, typename TOptions>
void writeSim(TBedOut &outBedFile1, TBedOut &outBedFile2, CharString &contigName, BindingSite &bindingSite, unsigned readCount, TOptions &/*options*/)
{
    if (readCount > 0)
    {
        BedRecord<Bed6> bedRecord;
        bedRecord.ref = contigName;
 
        // output read count around binding sites 
        bedRecord.beginPos = bindingSite.beginPos;
        bedRecord.endPos = bindingSite.endPos;
        
        std::stringstream ss;
        ss << (int)readCount;
        bedRecord.name = ss.str();
        ss.str("");  
        ss.clear();  
        
        ss << (double)bindingSite.bindingAffinity;
        bedRecord.score = ss.str();
        ss.str("");  
        ss.clear();  
        bedRecord.strand = bindingSite.strand;
        writeRecord(outBedFile1, bedRecord);

        // output individual simulated crosslinkSites
        for (unsigned i = 0; i < length(bindingSite.crosslinkSites); ++i)
        {
            bedRecord.beginPos = bindingSite.crosslinkSites[i];    
            bedRecord.endPos = bindingSite.crosslinkSites[i] + 1;

            std::stringstream ss;
            ss << (int)readCount;
            bedRecord.name = ss.str();
            ss.str("");  
            ss.clear();  

            ss << (double)bindingSite.crosslinkProbs[i];
            bedRecord.score = ss.str();
            ss.str("");  
            ss.clear();  
            bedRecord.strand = bindingSite.strand;
            writeRecord(outBedFile2, bedRecord);
        }
    }
}
 

template <typename TRng, typename TOptions>
unsigned draw_fragLength(String<double> &fragLengthDistr, TRng &rng, TOptions &/*options*/)
{
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double x = dist(rng); 
    double cum = 0.0;

    for (unsigned i = 0; i < length(fragLengthDistr); ++i)
    {
        cum += fragLengthDistr[i];
        if (x <= cum)
            return i;
    }   
    return length(fragLengthDistr);     // should not happen 
}


template <typename TBamOut, typename TBedOut, typename TBamIn, typename TBai, typename TRng, typename TOptions>
bool process_bamRegion(Stats &stats, TBamOut &outBamFile, TBedOut &outBedFile1, TBedOut &outBedFile2, TBamIn &inFile, TBai &baiIndex, int const& rID, BindingSite &bindingSite, String<double> &fragLengthDistr, TRng &rng, bool isTarget, TOptions &options)
{
    int jump_beginPos;
    int jump_endPos; 
    if (bindingSite.strand == '+')
    {
        jump_beginPos = bindingSite.endPos - options.maxFragmentSize - 1;
        jump_endPos = bindingSite.beginPos + options.maxReadLength + 1;   // jumpToRegions() regarding left most position or whole reads?
    }
    else
    {
        jump_beginPos = bindingSite.endPos - options.maxReadLength - 1;                           
        jump_endPos = bindingSite.beginPos + options.maxFragmentSize + 1;   // ?
    }
    //beginPos -= 1;
    //endPos -= 1;

    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    if (!jumpToRegion(inFile, hasAlignments, rID, jump_beginPos, jump_endPos, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << jump_beginPos << ":" << jump_endPos << "\n";
        return false;
    }
    if (!hasAlignments)
    {
        if (options.verbosity > 1) std::cout << "WARNING: no alignments here " << jump_beginPos << ":" << jump_endPos << "\n";
        return false;  
    }


    // Seek linearly to the selected position.
    BamAlignmentRecord bamRecord;
    BamAlignmentRecord newBamRecord;
    CharString contigName;
    unsigned readCount = 0;
    while (!atEnd(inFile))
    {
        readRecord(bamRecord, inFile);

        // If we are on the next reference 
        if (bamRecord.rID == -1 || bamRecord.rID > rID || !hasFlagAllProper(bamRecord))
            break;

        unsigned l = options.fragmentSize;
        if (options.fldFileName != "")
            l = draw_fragLength(fragLengthDistr, rng, options);

        if (!hasFlagRC(bamRecord))          // Forward
        {
            if (bindingSite.strand == '-')
                continue;

            // If read starts behind start of the binding site (and thus also all following) RT protection
            if (isTarget && (bamRecord.beginPos > (signed)bindingSite.beginPos))
                break;

            // do not simulate RT protection by proteins for background binding
            if (!isTarget && (bamRecord.beginPos >= (signed)bindingSite.endPos))
                break;

            // If mate ends before end of binding site
            if (bamRecord.beginPos + (signed)l < (signed)bindingSite.endPos)
                continue;
        }
        else                                // Reverse  
        {
            if (bindingSite.strand == '+')
                continue;

            // If read starts that far away, that we know its mate also starts behind start of the binding site (and thus also all following) 
            if (bamRecord.beginPos > (signed)(bindingSite.beginPos + options.maxFragmentSize - options.maxReadLength))
                break;

            // If read ends before end of binding site  ... RT protection
            if (isTarget && (bamRecord.beginPos + (signed)getAlignmentLengthInRef(bamRecord) < (signed)bindingSite.endPos))
                continue;

            // do not simulate RT protection by proteins for background binding
            if (!isTarget && (bamRecord.beginPos + (signed)getAlignmentLengthInRef(bamRecord) <= (signed)bindingSite.beginPos))
                continue;

            // If mate starts behind start of binding site 
            if (bamRecord.beginPos + (signed)getAlignmentLengthInRef(bamRecord) - (signed)l > (signed)bindingSite.beginPos)     
                continue;
        }

        /*if (!hasFlagFirst(bamRecord)) // currently assuming reads already selected
        {
            std::cerr << "WARNING: not the first segment! R2 entries should be removed from BAM before. " << "\n";
            continue;
        }*/
        if (std::abs(bamRecord.tLen) < options.minFragmentSize || std::abs(bamRecord.tLen) > options.maxFragmentSize)
            continue;
        
        /////////////////////////////
        // apply iCLIP modifications
        /////////////////////////////
        // TODO shorten to maxReadLength

        // downSample corresponding to binding affinity
        std::uniform_real_distribution<double> distBA(0.0, 1.0);
        double x = distBA(rng); 
        if (x <= bindingSite.bindingAffinity)
        {
            //std::cout << "binding affinity : " << bindingSite.bindingAffinity  << "\n";
            if (apply_iCLIPmodifications(stats, newBamRecord, bamRecord, bindingSite, l, rng, isTarget, options))
            {
                writeRecord(outBamFile, newBamRecord);

                SEQAN_OMP_PRAGMA(atomic)
                ++readCount;
                
                contigName = getContigName(bamRecord, inFile);
            }
        }
    }
    writeSim(outBedFile1, outBedFile2, contigName, bindingSite, readCount, options);
    if (readCount > 0) 
    {
        SEQAN_OMP_PRAGMA(atomic)
        stats.readCount += readCount;
        SEQAN_OMP_PRAGMA(atomic)
        ++stats.bsCount;
    }

    return true;
}


template <typename TBamOut, typename TBamIn, typename TBai, typename TRng, typename TOptions>
bool process_bamRegion_subsample(Stats &stats, TBamOut &outBamFile, TBamIn &inFile, TBai &baiIndex, int const& rID, String<double> &fragLengthDistr, TRng &rng, TOptions &options)
{
    int jump_beginPos = 0;
    int jump_endPos = 500000000;    // TODO 

    // Jump the BGZF stream to this position.
    bool hasAlignments = false;
    if (!jumpToRegion(inFile, hasAlignments, rID, jump_beginPos, jump_endPos, baiIndex))
    {
        std::cerr << "ERROR: Could not jump to " << jump_beginPos << ":" << jump_endPos << "\n";
        return false;
    }
    if (!hasAlignments)
    {
        if (options.verbosity > 1) std::cout << "WARNING: no alignments here " << jump_beginPos << ":" << jump_endPos << "\n";
        return false;  
    }

    // Seek linearly to the selected position.
    BamAlignmentRecord bamRecord;
    BamAlignmentRecord newBamRecord;
    CharString contigName;
    unsigned readCount = 0;
    while (!atEnd(inFile))
    {
        readRecord(bamRecord, inFile);

        // If we are on the next reference 
        if (bamRecord.rID == -1 || bamRecord.rID > rID || !hasFlagAllProper(bamRecord))
            break;

        unsigned fragLen = options.fragmentSize;
        if (options.fldFileName != "")
            fragLen = draw_fragLength(fragLengthDistr, rng, options);
 
        /*if (!hasFlagFirst(bamRecord)) // currently assuming reads already selected
        {
            std::cerr << "ERROR: not the first segment! R2 entries should be removed from BAM before. " << "\n";
            continue;
        }*/     

        // downSample corresponding to subsampling rate
        std::uniform_real_distribution<double> distBA(0.0, 1.0);
        double x = distBA(rng); 
        if (x <= options.subsampleRate)
        {                
            newBamRecord = bamRecord;

            unsigned readLen = std::min(options.simReadLength, fragLen);
            // shorten sequences, qualities, adjust cigar, sequence invalid "AAAAA..."
            if (!adjustBamRecord(newBamRecord, readLen, fragLen, options))
                return false;  
  
            writeRecord(outBamFile, newBamRecord);

            SEQAN_OMP_PRAGMA(atomic)
            ++readCount;
            
            contigName = getContigName(bamRecord, inFile);
        }
    }
    if (readCount > 0) 
    {
        SEQAN_OMP_PRAGMA(atomic)
        stats.randomReadCount += readCount;
    }

    return true;
}



template <typename TOptions>
bool parse_fragLengthDistr(String<double> &fragLengthDistr, TOptions &options) {
    
    std::ifstream ifs(toCString(options.fldFileName));
    if (!ifs.good())
        std::cerr << "ERROR: Could not open fragment length distribution file!\n";
    std::string lineBuffer;

    String<unsigned> lengths;
    String<double> probs;
    while (std::getline(ifs, lineBuffer))
    {
        //std::cout << lineBuffer << std::endl;
        std::string value1;
        std::string value2;
        std::stringstream ss(lineBuffer);
        if (!std::getline(ss, value1, '\t'))
            std::cerr << "ERROR: could not read first value\n";
        if (!std::getline(ss, value2, '\t'))
            std::cerr << "ERROR: could not read second value\n";

        double d1 = std::strtod(value1.c_str(), NULL);
        double d2 = std::strtod(value2.c_str(), NULL);

        appendValue(lengths, d1);
        appendValue(probs, d2);
    }
    resize(fragLengthDistr, (lengths[length(lengths)-1]+1), 0.0, Exact());
    for (unsigned i = 0; i < length(lengths); ++i)
        fragLengthDistr[lengths[i]] = probs[i];
    
    return 0;
}


inline
const char * myTempFileName(std::string const suffix = "", std::string const tempPath = "")
{
    if (mkdir(tempPath.c_str(), 0777) == -1)
    {
        if(errno == EEXIST ) {
            //std::cout << "Directory already exists " << std::endl;
        } else {
            std::cerr << tempPath << std::endl;
            throw std::runtime_error("ERROR: Could not create directory");
        }
    }

    static char fileNameBuffer[1000];
    strcpy(fileNameBuffer, (tempPath + "SIM_ICLIP.XXXXXX" + suffix).c_str());

    int _tmp = mkstemps(fileNameBuffer, suffix.size());
    if (_tmp == -1)
        throw std::runtime_error("Could not create temp file");

    close(_tmp);    
    return fileNameBuffer;
}


inline bool exists_test(const CharString& fileName) {
    std::ifstream f(toCString(fileName));
    return f.good();
}


// TODO use for each crosslink sites probability of a fragment to be truncated here... (influence each other)
template <typename TOptions>
bool doIt(TOptions &options)
{
    // load reference
    std::cout << "Load reference ... "  << "\n";
    FragmentStore<> store;
    loadContigs(store, toCString(options.refFileName));

    // load given binding sites (assume same order as in BAM and reference file)
    std::cout << "Load target binding sites ... "  << "\n";
    String<String<BindingSite> > bindingSites;  
    loadBindingSites(bindingSites, store, options.bsFileName, false, options);

    String<String<BindingSite> > bgBindingSites; 
    if (options.useBgBS)
    {
        std::cout << "Load background binding sites ... " << options.bbsFileName << " " << true << "\n";
        loadBindingSites(bgBindingSites, store, options.bbsFileName, true, options);
    }

    String<double> fragLengthDistr;
    if (options.fldFileName != "")
        parse_fragLengthDistr(fragLengthDistr, options);

    ///////////////////////////
    // target binding sites
    ///////////////////////////
    Stats targetStats;

    // for each binding site  
    // NOTE: ignore additive binding here.... one fragment can cause reads at multiple binding sites
    const int SEED = 0;
    typedef std::mt19937 TRng;
    String<CharString> tempFileNames;
    resize(tempFileNames, length(store.contigStore));

    bool abort = false;
    if (!options.ignoreTargetSites)
    {
        std::cout << "  Simulate target signal ..." << std::endl;
#if SIM_ICLIP_ENABLE_PARALLELISM
        if (options.verbosity > 1) std::cout << "     in parallel ..." << std::endl;
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreads)) 
#endif  
        for (unsigned contigId = 0; contigId < length(bindingSites); ++contigId)
        {
            // Open BamFileIn for reading.
            if (options.verbosity > 1) std::cout << "Open Bam and Bai file ... "  << "\n";
            BamFileIn inFile;
            if (!open(inFile, toCString(options.bamFileName)))
            {
                std::cerr << "ERROR: Could not open " << options.bamFileName << " for reading.\n";
                abort = true;
            }
            BamHeader header;
            readHeader(header, inFile);
            // Read BAI index.
            BamIndex<Bai> baiIndex;
            if (!open(baiIndex, toCString(options.baiFileName)))
            {
                std::cerr << "ERROR: Could not read BAI index file " << options.baiFileName << "\n";
                abort = true;
            }

            // Temp. fileName
            CharString tempFileName;
            if (empty(options.tempPath))
            {
                SEQAN_OMP_PRAGMA(critical)
                append(tempFileName, SEQAN_TEMP_FILENAME());
                std::stringstream ss;
                ss << contigId;
                append(tempFileName, ss.str());
                append(tempFileName, ".bam");
            }
            else
            {
                SEQAN_OMP_PRAGMA(critical)
                tempFileName = myTempFileName(".bam", toCString(options.tempPath));
            }
            tempFileNames[contigId] = tempFileName;
            if (options.verbosity > 1) std::cout << "temp file Name: " << tempFileName << std::endl;

            // Bam output
            BamFileOut outBamFile(context(inFile), toCString(tempFileName));
            writeHeader(outBamFile, header);
            // Bed log output
            CharString tempBedName1 = prefix(tempFileName, length(tempFileName) - 4);
            append(tempBedName1, ".bs.bed");
            BedFileOut outBedFile1(toCString(tempBedName1));
            CharString tempBedName2 = prefix(tempFileName, length(tempFileName) - 4);
            append(tempBedName2, ".cs.bed");
            BedFileOut outBedFile2(toCString(tempBedName2));

            TRng rng(SEED);
            // Translate from contig name to rID.
            int rID = 0;
            if (!getIdByName(rID, contigNamesCache(context(inFile)), store.contigNameStore[contigId]))
            {
                std::cerr << "WARNING: Reference sequence named " << store.contigNameStore[contigId] << " not known.\n";
                //abort = true; 
            }
            else
            {
                for (unsigned i = 0; i < length(bindingSites[contigId]); ++i)
                {
                    // number of crosslink sites 
                    unsigned c = options.noCrosslinkSites;
                    if (options.useRnCrosslinkSites)
                    {
                        std::uniform_int_distribution<unsigned> dist(1, options.noCrosslinkSites);
                        c = dist(rng);
                        c = std::min((bindingSites[contigId][i].endPos-bindingSites[contigId][i].beginPos), c); // make sure fitting within binding region
                    }
                    // simulate crosslink sites
                    simulatePotentialCrosslinkSites(bindingSites[contigId][i], c, rng, options);
                    if(!process_bamRegion(targetStats, outBamFile, outBedFile1, outBedFile2, inFile, baiIndex, rID, bindingSites[contigId][i], fragLengthDistr, rng, true, options))
                    {
                        SEQAN_OMP_PRAGMA(atomic)
                            ++targetStats.uncoveredBsCount;
                    }
                }
            }
        }
    }
    if (abort) return 1;
    std::cout << "****************************** "  << "\n";
    std::cout << "Simulated " << targetStats.readCount << " reads for " << targetStats.bsCount << " target binding sites." << "\n";
    std::cout << "Simulated " << targetStats.targetTruncCount << " reads truncated at target crosslink sites. \n";
    std::cout << "Simulated " << targetStats.offTargetTruncCount << " reads truncated at off-target sites. \n";
    std::cout << "No. of sites uncovered by raw RNA-seq data:  " << targetStats.uncoveredBsCount << " ." << "\n";
    std::cout << "****************************** "  << "\n";

    ///////////////////////////
    // background binding sites
    ///////////////////////////
    Stats backgroundStats;
    String<CharString> bgTempFileNames;
    resize(bgTempFileNames, length(store.contigStore));
    if (options.useBgBS)    // !empty(options.bbsFileName))
    {
        std::cout << "  Simulate background signal ..." << std::endl;
        // for each binding site 
#if SIM_ICLIP_ENABLE_PARALLELISM
        if (options.verbosity > 1) std::cout << "     in parallel ..." << std::endl;
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 1) num_threads(options.numThreads)) 
#endif  
        for (unsigned contigId = 0; contigId < length(bgBindingSites); ++contigId)
        {
            // Open BamFileIn for reading.
            if (options.verbosity > 1) std::cout << "Open Bam and Bai file ... "  << "\n";
            BamFileIn inFile;
            if (!open(inFile, toCString(options.bamFileName)))
            {
                std::cerr << "ERROR: Could not open " << options.bamFileName << " for reading.\n";
                abort = true;
            }
            BamHeader header;
            readHeader(header, inFile);
            // Read BAI index.
            BamIndex<Bai> baiIndex;
            if (!open(baiIndex, toCString(options.baiFileName)))
            {
                std::cerr << "ERROR: Could not read BAI index file " << options.baiFileName << "\n";
                abort = true;    
            }

            // Temp. fileName
            CharString bgTempFileName;
            if (empty(options.tempPath))
            {
                SEQAN_OMP_PRAGMA(critical)
                append(bgTempFileName, SEQAN_TEMP_FILENAME());
                std::stringstream ss;
                ss << contigId;
                append(bgTempFileName, ss.str());
                append(bgTempFileName, ".bam");
            }
            else
            {
                SEQAN_OMP_PRAGMA(critical)
                bgTempFileName = myTempFileName(".bam", toCString(options.tempPath));   
            }
            bgTempFileNames[contigId] = bgTempFileName;
            if (options.verbosity > 1) std::cout << "temp file Name: " << bgTempFileName << std::endl;

            // Bam output
            BamFileOut outBamFile(context(inFile), toCString(bgTempFileName));
            writeHeader(outBamFile, header);
            // Bed log output
            CharString tempBedName3 = prefix(bgTempFileName, length(bgTempFileName) - 4);
            append(tempBedName3, ".bbs.bed");
            BedFileOut outBedFile3(toCString(tempBedName3));
            CharString tempBedName4 = prefix(bgTempFileName, length(bgTempFileName) - 4);
            append(tempBedName4, ".bcs.bed");
            BedFileOut outBedFile4(toCString(tempBedName4));

            TRng rng(SEED);
            // Translate from contig name to rID.
            int rID = 0;
            if (!getIdByName(rID, contigNamesCache(context(inFile)), store.contigNameStore[contigId]))
            {
                std::cerr << "WARNING: Reference sequence named " << store.contigNameStore[contigId] << " not known.\n";
                //abort = true;
            }
            else
            {
                for (unsigned i = 0; i < length(bgBindingSites[contigId]); ++i)
                {
                    // number of crosslink sites
                    // for the moment: make sure not less crosslink sites than for target sites (the fewer crosslink sites, the stronger signal at individual position)
                    unsigned c = options.noBgCrosslinkSites;
                    if (options.useRnCrosslinkSites)
                    {
                        std::uniform_int_distribution<unsigned> dist(options.noCrosslinkSites, options.noBgCrosslinkSites); 
                        c = dist(rng);
                        c = std::min((bgBindingSites[contigId][i].endPos-bgBindingSites[contigId][i].beginPos), c); // make sure fitting within binding region
                    }
                    // simulate crosslink sites
                    simulatePotentialCrosslinkSites(bgBindingSites[contigId][i], c, rng, options);
                    if(!process_bamRegion(backgroundStats, outBamFile, outBedFile3, outBedFile4, inFile, baiIndex, rID, bgBindingSites[contigId][i], fragLengthDistr, rng, false, options))
                    {
                        SEQAN_OMP_PRAGMA(atomic)
                            ++backgroundStats.uncoveredBsCount;
                    }
                }
            }
            // add random reads from RNA-seq data
            //process_bamRegion_subsample(backgroundStats, outBamFile, inFile, baiIndex, rID, fragLengthDistr, rng, options);
        }
        if (abort) return 1;

        std::cout << " Get max. threads " <<  omp_get_max_threads()  << std::endl;
        std::cout << "Simulated " << backgroundStats.readCount << " reads for " << backgroundStats.bsCount << " background binding sites." << "\n";
        std::cout << "Simulated " << backgroundStats.targetTruncCount << " reads truncated at target crosslink sites. \n";
        std::cout << "Simulated " << backgroundStats.offTargetTruncCount << " reads truncated at off-target sites. \n";
        std::cout << "Subsampled " << backgroundStats.randomReadCount << "random reads from RNA-seq data.\n";
        std::cout << "No. of sites uncovered by raw RNA-seq data:  " << backgroundStats.uncoveredBsCount << " ." << "\n";
        std::cout << "****************************** "  << "\n";
    }

    ////////////////////////////////////////////////////////////////
    // Append content of temp files to final output in contig order
    // Bam input just to get header
    BamFileIn inFile;
    if (!open(inFile, toCString(options.bamFileName)))
    {
        std::cerr << "ERROR: Could not open " << options.bamFileName << " for reading.\n";
        return 1;
    }
    BamHeader header;
    readHeader(header, inFile);
    // Bam output
    BamFileOut outBamFile(context(inFile), toCString(options.outFileName));
    writeHeader(outBamFile, header);
    // Target Bed log output
    CharString outBedName1 = prefix(options.outFileName, length(options.outFileName) - 4);
    append(outBedName1, ".bs.bed");
    BedFileOut outBedFile1(toCString(outBedName1));
    CharString outBedName2 = prefix(options.outFileName, length(options.outFileName) - 4);
    append(outBedName2, ".cs.bed");
    BedFileOut outBedFile2(toCString(outBedName2));
    // Background Bed log output
    CharString outBedName3 = prefix(options.outFileName, length(options.outFileName) - 4);
    append(outBedName3, ".bbs.bed");
    BedFileOut outBedFile3(toCString(outBedName3));
    CharString outBedName4 = prefix(options.outFileName, length(options.outFileName) - 4);
    append(outBedName4, ".bcs.bed");
    BedFileOut outBedFile4(toCString(outBedName4));

    for (unsigned contigId = 0; contigId < length(store.contigNameStore); ++contigId)
    {
        BamAlignmentRecord bamRecord;
        // Target BAM
        if (!options.ignoreTargetSites)
        {
            BamFileIn inFile;
            if (!open(inFile, toCString(tempFileNames[contigId])))
            {
                std::cerr << "ERROR: Could not open " << tempFileNames[contigId] << " for reading.\n";
                return 1;
            }
            readHeader(header, inFile);
            while (!atEnd(inFile))
            {
                readRecord(bamRecord, inFile);
                writeRecord(outBamFile, bamRecord);
            }
            BedRecord<seqan::Bed6> bedRecord;
            // Target BED
            BedFileIn bedFileIn1;
            CharString tempBedName1 = prefix(tempFileNames[contigId], length(tempFileNames[contigId]) - 4);
            append(tempBedName1, ".bs.bed");
            if (!open(bedFileIn1, toCString(tempBedName1)))
            {
                std::cerr << "ERROR: Could not open temporary bed file: " << tempBedName1 << "\n";
                return 1;
            }
            while (!atEnd(bedFileIn1))
            {
                readRecord(bedRecord, bedFileIn1);
                writeRecord(outBedFile1, bedRecord);
            }
            BedFileIn bedFileIn2;
            CharString tempBedName2 = prefix(tempFileNames[contigId], length(tempFileNames[contigId]) - 4);
            append(tempBedName2, ".cs.bed");
            if (!open(bedFileIn2, toCString(tempBedName2)))
            {
                std::cerr << "ERROR: Could not open temporary bed file: " << tempBedName2 << "\n";
                return 1;
            }
            while (!atEnd(bedFileIn2))
            {
                readRecord(bedRecord, bedFileIn2);
                writeRecord(outBedFile2, bedRecord);
            }
            // Remove temp files and check
            std::remove(toCString(tempFileNames[contigId]));
            if (exists_test(tempFileNames[contigId]))
            {
                std::cerr << "ERROR: Could open temporary bed file which should be deleted: " << tempFileNames[contigId]  << std::endl;
                return 1;
            }
            std::remove(toCString(tempBedName1));
            if (exists_test(tempBedName1))
            {
                std::cerr << "ERROR: Could open temporary bed file which should be deleted: " << tempBedName1  << std::endl;
                return 1;
            }
            std::remove(toCString(tempBedName2));
            if (exists_test(tempBedName2))
            {
                std::cerr << "ERROR: Could open temporary bed file which should be deleted: " << tempBedName2  << std::endl;
                return 1;
            }
        }

        // Background BAM
        if (options.useBgBS)
        {
            BamFileIn bgInFile;
            if (!open(bgInFile, toCString(bgTempFileNames[contigId])))
            {
                std::cerr << "ERROR: Could not open " << bgTempFileNames[contigId] << " for reading.\n";
                return 1;
            }
            readHeader(header, bgInFile);
            while (!atEnd(bgInFile))
            {
                readRecord(bamRecord, bgInFile);
                writeRecord(outBamFile, bamRecord);
            }
            BedRecord<seqan::Bed6> bedRecord;
            // Background BED
            BedFileIn bedFileIn3;
            CharString tempBedName3 = prefix(bgTempFileNames[contigId], length(bgTempFileNames[contigId]) - 4);
            append(tempBedName3, ".bbs.bed");
            if (!open(bedFileIn3, toCString(tempBedName3)))
            {
                std::cerr << "ERROR: Could not open temporary bed file: " << tempBedName3 << "\n";
                return 1;
            }
            while (!atEnd(bedFileIn3))
            {
                readRecord(bedRecord, bedFileIn3);
                writeRecord(outBedFile3, bedRecord);
            }
            BedFileIn bedFileIn4;
            CharString tempBedName4 = prefix(bgTempFileNames[contigId], length(bgTempFileNames[contigId]) - 4);
            append(tempBedName4, ".bcs.bed");
            if (!open(bedFileIn4, toCString(tempBedName4)))
            {
                std::cerr << "ERROR: Could not open temporary bed file: " << tempBedName4 << "\n";
                return 1;
            }
            while (!atEnd(bedFileIn4))
            {
                readRecord(bedRecord, bedFileIn4);
                writeRecord(outBedFile4, bedRecord);
            }
            // Remove temp files and check
            std::remove(toCString(bgTempFileNames[contigId]));
            if (exists_test(bgTempFileNames[contigId]))
            {
                std::cerr << "ERROR: Could open temporary bed file which should be deleted: " << bgTempFileNames[contigId]  << std::endl;
                return 1;
            }
            std::remove(toCString(tempBedName3));
            if (exists_test(tempBedName3))
            {
                std::cerr << "ERROR: Could open temporary bed file which should be deleted: " << tempBedName3  << std::endl;
                return 1;
            }
            std::remove(toCString(tempBedName4));
            if (exists_test(tempBedName4))
            {
                std::cerr << "ERROR: Could open temporary bed file which should be deleted: " << tempBedName4  << std::endl;
                return 1;
            }
        }
    }
    std::cout << "Simulation of iCLIP reads finished successfully. " << std::endl;

    return 0;
}


#endif
