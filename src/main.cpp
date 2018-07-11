/** @file main.cpp
 * @brief One sentence brief
 *
 * More details
 * In multiple lines
 * Copyright © 2016 Can Fırtına. All rights reserved.
 *
 * @author Can Firtina
 * @bug No bug hopefully
 */

#include <string>
#include <math.h>
#include <set>
#include <algorithm>
#include <stdio.h>
#include <vector>
#include <map>
#include <thread>
#include "CommandLineOptions.h"
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <mutex>
#include <condition_variable>
#include "bloom_filter.hpp"

#if defined(_WIN16) | defined(_WIN32) | defined(_WIN64)
#define SEPERATOR "\\"
#else
#define SEPERATOR "/"
#endif

#define MATCH_OFFSET(char, offset, insize) ((char+offset)*(insize+1))
#define INSERTION_OFFSET(char, offset, innumber, insize) ((char+offset)*(insize+1) + innumber) //innumber 1 based
#define GRAPH_SIZE(size, insize) ((size+1)*(insize+1) + 1)
#define END_STATE(size, insize) (GRAPH_SIZE(size, insize)-1)

using namespace seqan;

std::mutex indexMutex;
std::mutex outputFileMutex;

//#define DEBUG_HMM_DETAILED
//#define DEBUG_HMM_RESULT
//#define DEBUG_SHORT_READ_COMPLETION

enum Nucleotide {A, T, G, C, totalNuc};

/** @brief Parameters to be used in profile hidden Markov model graph
 *
 */
struct HMMParameters{
    
public:
    
    HMMParameters(unsigned filterSize, unsigned maxCoverage, unsigned maxDeletion, unsigned maxInsertion, double matchTransition,
                  double insertionTransition, double deletionTransitionFactor, double matchEmission): filterSize(filterSize),
    maxCoverage(maxCoverage), maxDeletion(maxDeletion), maxInsertion(maxInsertion), matchTransition(matchTransition),
    insertionTransition(insertionTransition), deletionTransitionFactor(deletionTransitionFactor), matchEmission(matchEmission){
        deletionTransition = 1.000 - (matchTransition + insertionTransition);
        mismatchEmission = (double)(1 - matchEmission)/3.00;
        insertionEmission = (double)1/3.00; //total nucleotide = 4; Emission prob for each except one
    }
    
    HMMParameters(const HMMParameters& cpy):
    filterSize(cpy.filterSize), maxCoverage(cpy.maxCoverage), maxDeletion(cpy.maxDeletion), maxInsertion(cpy.maxInsertion),
    matchTransition(cpy.matchTransition), insertionTransition(cpy.insertionTransition),
    deletionTransition(cpy.deletionTransition), deletionTransitionFactor(cpy.deletionTransitionFactor),
    matchEmission(cpy.matchEmission), mismatchEmission(cpy.mismatchEmission),
    insertionEmission(cpy.insertionEmission){}
    
    unsigned filterSize;
    unsigned maxDeletion;
    unsigned maxInsertion;
    unsigned maxCoverage;
    double matchTransition;
    double insertionTransition;
    double deletionTransition;
    double deletionTransitionFactor;
    double matchEmission;
    double mismatchEmission;
    double insertionEmission;
};

/** @brief Represents a state in the profile hidden Markov model graph
 *
 */
struct SequencingNode{
public:
    
    SequencingNode(){}
    SequencingNode(unsigned int index, unsigned int charIndex, unsigned int insize, char nuc, char nextNuc):
    index(index), charIndex(charIndex), nuc(nuc), isMatch((index%(insize+1) == 0)?true:false),
    isLastInsertion((insize > 0 && index%(insize+1) == insize)?true:false), nextNuc(nextNuc){}
    
    /** @brief Emission probability calculation
     *
     *  @param character Basepair to calculate the emission probability
     *  @param parameters Calculation is based on the specified parameters
     *  @return Emission probability of the given character (basepair) [0-1]
     */
    double getEmissionProb(const char& character, const HMMParameters& parameters){
        
        if(isMatch){//match state: either match or mismatch emission probability
            if(character == nuc || character == 'N')
                return parameters.matchEmission;
            return parameters.mismatchEmission;
        }else if(character == 'N' || character == nextNuc){
            //insertion state: no emission for 'N' or for the next character in long read
            return 0.00;
        }
        
        return parameters.insertionEmission; //insertion state
    }
    
    /** @brief Transition probability calculation from this state to the specified `toIndex` state
     *
     *  @param toIndex Specifies which state to transit from this state
     *  @param parameters Calculation is based on the specified parameters
     *  @return Transition probability from this state to the state indicated with `toIndex`
     */
    double transitionProbFromThisNode(const unsigned int& toIndex, HMMParameters& parameters){
        
        if(MATCH_OFFSET(charIndex, 1, parameters.maxInsertion) == toIndex){
            if(isLastInsertion){
                //match transition prob for last insertion state
                return parameters.matchTransition + parameters.insertionTransition;
            }
            return parameters.matchTransition; //match transition prob
        }
        
        if(index+1 == toIndex) return parameters.insertionTransition; //insertion transition
        
        //deletion transition calculations: normalized polynomial distribution
        int count = 0;
        int start = 1;
        double transitionProb = 0;
        for(int curDel = parameters.maxDeletion+1; curDel > 1; --curDel){
            if(MATCH_OFFSET(charIndex, curDel, parameters.maxInsertion) == toIndex)
                transitionProb = parameters.deletionTransition*start;
            count+=start;
            start*=parameters.deletionTransitionFactor;
        }
        
        return (count==0)?0:transitionProb/(double)count;
    }
    
    bool isMatchState() const { return isMatch; }
    bool isLastInsertionState() const { return isLastInsertion; }
    unsigned int getIndex() const { return index;}
    char getNucleotide() const { return nuc;}
    unsigned int getCharIndex() const { return charIndex; }
    
private:
    
    char nextNuc; //basepair in the next position
    bool isMatch; //is a match state?
    bool isLastInsertion;
    unsigned int index; //what is the index of this state in the graph
    unsigned int charIndex; //which character it corresponds to in the sequencing read
    char nuc; //basepair in the charIndex
};

struct ShortRead{
public:
    ShortRead(std::string shortRead, int pos, std::string cigar = ""):pos(pos), editDistance(-1){
        if(cigar.size() > 0){
            std::istringstream cigarStream(cigar);
            int curIndex = 0;
            int howmany;
            char type;
            while(cigarStream >> howmany >> type){
                if(type == 'S'){
                    shortRead.erase(curIndex, howmany);
                }else if(type == 'M' || type == 'I') curIndex += howmany;
            }
        }
        
        this->shortRead = shortRead;
    }
    
    std::string shortRead;
    int pos;
    int editDistance;
};

/*
 * Can be used as directed edge between from (outgoing edge source), and curState (incoming edge source).
 * operator< implemented so that std::set can evaulate this data structure
 */
struct TransitionInfoNode{
    
    int from, toState;
    TransitionInfoNode(int from, int toState):from(from), toState(toState){}
    
    bool operator==(const TransitionInfoNode& rhs) const{
        return from == rhs.from && toState == rhs.toState;
    }
    
    bool operator<(const TransitionInfoNode& rhs) const{
        return ((from < rhs.from) || (from == rhs.from && (toState < rhs.toState)));
    }
};

std::string compress(std::string& uncompressed, int& nonNCount){
    
    nonNCount = 0;
    std::string compressed = "";
    for(int curChar = 0; curChar < uncompressed.size(); ++curChar){
        compressed.append(&uncompressed[curChar], 1);
        if(uncompressed[curChar] != 'N') nonNCount++;
        while(curChar + 1 < uncompressed.size() && uncompressed[curChar] == uncompressed[curChar+1]) ++curChar;
    }
    
    return compressed;
}

/*
 * Starting from startValues, until endValues, finds the indices that has the greatest values in values array and puts
 * in maxValues array. If there cannot be more than maxValuesSize different indices, rest of maxValues left as -1
 * Return parameter: maxValues array
 * Return value:
 */
template <typename T>
void findMaxValues(const T* values, bool* selectedIndices, const int startValues, const int endValues,
                     const int maxValuesSize){
    
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int> >,
                        std::greater<std::pair<double, int> > > maxValues;
    
    for(int curState = 0; curState < endValues - startValues; ++curState) {
        if(maxValues.size() < maxValuesSize)
            maxValues.push(std::pair<double, int>(values[curState], curState));
        else if(maxValues.top().first < values[curState]){
            maxValues.pop();
            maxValues.push(std::pair<double, int>(values[curState], curState));
        }
    }

    while(!maxValues.empty()){
        selectedIndices[maxValues.top().second] = true;
        maxValues.pop();
    }
}

/*
 * Using the information provided with node and numberOfDeletions, inserts the transitions that should be made afterward
 * If you are going to change the transition structure, change it from there. These are imaginary edges in the graph.
 */
void insertNewForwardTransitions(std::vector<TransitionInfoNode>* transitionSet, const SequencingNode& node,
                                 const int numberOfDeletions, const int maxInsertion){
    //next insertion
    if(!node.isLastInsertionState())
        transitionSet->push_back(TransitionInfoNode(node.getIndex(), node.getIndex()+1));
    
    //match and deletions
    for(int curOffset = 1; curOffset <= numberOfDeletions+1; ++curOffset)
        transitionSet->push_back(TransitionInfoNode(node.getIndex(), MATCH_OFFSET(node.getCharIndex(), curOffset,
                                                                                  maxInsertion)));
}

/*
 * @graph: graph structure for HMM. Built based on the read to correct. There are multiple states representing the same
 * character in the input read. Provides transition probabilities to each other states from a state and emission
 * probabilities given a nucleotide for each state.
 * @forwardMatrix: Return value of the function. Forward-likelihood values are computed and written on this matrix.
 * First dimension is time, second dimension is a state in @graph
 * @shortRead: Read that aligns to a specific position of the read to be corrected (so aligns to @graph as well)
 * @startPosition: Represent for which character in the long read the forward likelihood starts to be computed. It's a
 * 1-based value. So if the forward likelihood will start with the first character, this value is 1
 * @maxDistanceOnLongRead: Which character (or which state in @graph) is the last character to calculate its
 * forward-likelihood value. This is again 1 based, and should be greater than @startPosition
 * @longReadLength: length of the read to be corrected
 * @shortReadLength: length of the read that aligns to the longer read to correct
 */
int fillForwardMatrix(SequencingNode* graph, HMMParameters parameters, double* calculatedTransitionProbs,
                      double** forwardMatrix, const char* shortRead, int const startPosition, int maxDistanceOnLongRead,
                      const int longReadLength, int const shortReadLength){
    
    if(startPosition < 0 || startPosition >= maxDistanceOnLongRead) return -1;
    if(GRAPH_SIZE(longReadLength, parameters.maxInsertion) < maxDistanceOnLongRead)
        maxDistanceOnLongRead = GRAPH_SIZE(longReadLength, parameters.maxInsertion);
    
    //which transitions requested from the previous time
    std::vector<TransitionInfoNode>* curTrSet = new std::vector<TransitionInfoNode>;
    //which transitions to be made for the next time
    std::vector<TransitionInfoNode>* nextTrSet = new std::vector<TransitionInfoNode>;
    std::vector<TransitionInfoNode>* tmpTrSet;
    
    bool* allowedParentStates = new bool[maxDistanceOnLongRead-startPosition+1];
    std::vector<bool> hasStateBeenProcessedBefore;
    hasStateBeenProcessedBefore.resize(maxDistanceOnLongRead-startPosition+1);
    std::fill_n(hasStateBeenProcessedBefore.begin(), hasStateBeenProcessedBefore.size(), false);
    std::fill_n(allowedParentStates, maxDistanceOnLongRead-startPosition+1, false);
    
    //1-initialization (t = 1)
    int curTime = 0; //represents the current time (1...T)
    insertNewForwardTransitions(curTrSet, graph[startPosition], parameters.maxDeletion, parameters.maxInsertion);
    for(int curTransition = 0; curTransition < curTrSet->size(); ++curTransition){
        int matchoff = MATCH_OFFSET(graph[curTrSet->at(curTransition).from].getCharIndex(), 0, parameters.maxInsertion);
        //0->insertion, 1-> match, 2,3...->deletions
        int transitionIndex = (curTrSet->at(curTransition).toState - matchoff)/(parameters.maxInsertion+1);
        forwardMatrix[0][curTrSet->at(curTransition).toState - startPosition] += calculatedTransitionProbs[transitionIndex]*
        graph[curTrSet->at(curTransition).toState].getEmissionProb(shortRead[curTime], parameters);
        
        insertNewForwardTransitions(nextTrSet, graph[curTrSet->at(curTransition).toState], parameters.maxDeletion,
                                    parameters.maxInsertion);
    }
    curTrSet->clear();
    
    //find the most likely states that should be allowed to make the next transitions
    findMaxValues(forwardMatrix[0], allowedParentStates, startPosition, maxDistanceOnLongRead, parameters.filterSize);
    
    tmpTrSet = curTrSet;
    curTrSet = nextTrSet;
    nextTrSet = tmpTrSet;
    
    //2-recursion (1 < t <= T)
    while (curTime < shortReadLength-1 && !curTrSet->empty()){
        curTime++;
        for(int curTransition = 0; curTransition < curTrSet->size(); ++curTransition){
            const TransitionInfoNode& frontTr = curTrSet->at(curTransition);
            if(allowedParentStates[frontTr.from - startPosition] && frontTr.toState < maxDistanceOnLongRead) {
                int matchoff = MATCH_OFFSET(graph[frontTr.from].getCharIndex(), 0, parameters.maxInsertion);
                //0->insertion, 1-> match, 2,3...->deletions
                int transitionIndex = (frontTr.toState - matchoff)/(parameters.maxInsertion+1);
                forwardMatrix[curTime][frontTr.toState-startPosition] +=
                    forwardMatrix[curTime-1][frontTr.from-startPosition]*
                    calculatedTransitionProbs[transitionIndex]*
                    graph[frontTr.toState].getEmissionProb(shortRead[curTime],parameters);
                
                if(!hasStateBeenProcessedBefore[frontTr.toState - startPosition]){
                    insertNewForwardTransitions(nextTrSet, graph[frontTr.toState], parameters.maxDeletion,
                                                parameters.maxInsertion);
                    hasStateBeenProcessedBefore[frontTr.toState - startPosition] = true;
                }
            }
        }
        curTrSet->clear();
        
        std::fill_n(hasStateBeenProcessedBefore.begin(), hasStateBeenProcessedBefore.size(), false);
        std::fill_n(allowedParentStates, maxDistanceOnLongRead-startPosition+1, false);
        
        //find the most likely states that should be allowed to make the next transitions
        findMaxValues(forwardMatrix[curTime], allowedParentStates, startPosition, maxDistanceOnLongRead,
                      parameters.filterSize);
        tmpTrSet = curTrSet;
        curTrSet = nextTrSet;
        nextTrSet = tmpTrSet;
    }
    
    int max = -1;
    for(int curStateForMax = 0; curStateForMax < maxDistanceOnLongRead-startPosition+1; ++curStateForMax){
        if(allowedParentStates[curStateForMax] &&
           (max == -1 || forwardMatrix[curTime][curStateForMax] > forwardMatrix[curTime][max])){
            max = curStateForMax;
        }
    }
    
    delete curTrSet; delete nextTrSet; delete[] allowedParentStates;
    
    return max + startPosition;
}

/*
 * Using the information provided with node and numberOfDeletions, inserts the transitions that should be made after
 * this node.
 * If you are going to change the transition structure, change it from there. These are imaginary edges in the graph.
 */
void insertNewBackwardTransitions(std::vector<TransitionInfoNode>* transitionSet, const SequencingNode& node,
                                  const int numberOfDeletions, const int maxInsertion){
    
    if(node.getCharIndex() == 0) return;
    
    if(node.isMatchState()){
        for(int curBeforeOffset = 0; curBeforeOffset <= numberOfDeletions; ++curBeforeOffset){
            int offset = -1*curBeforeOffset - 1;
            transitionSet->push_back(TransitionInfoNode(node.getIndex(), MATCH_OFFSET(node.getCharIndex(), offset,
                                                                                      maxInsertion))); //deletion
            for(int curInsertion = 1; curInsertion <= maxInsertion; ++curInsertion)
                transitionSet->push_back(TransitionInfoNode(node.getIndex(),
                                                            INSERTION_OFFSET(node.getCharIndex(),offset, curInsertion,
                                                                             maxInsertion)));//match
        }
    }else{
        transitionSet->push_back(TransitionInfoNode(node.getIndex(), node.getIndex()-1));//match
    }
}

/*
 * @graph: graph structure for HMM. Built based on the read to correct. There are multiple states representing the same
 * character in the input read. Provides transition probabilities to each other states from a state and emission
 * probabilities given a nucleotide for each state.
 * @backwardMatrix: Return value of the function. Backward-likelihood values are computed and written on this matrix.
 * First dimension is time, second dimension is a state in @graph
 * @shortRead: Read that aligns to a specific position of the read to be corrected (so aligns to @graph as well)
 * @startPosition: Represents where the backward likelihood computation will start. It is the character number on the
 * long read. 1-based value. So if it starts from the end of the read, the value is @longReadLength
 * @maxDistanceOnLongRead: Which character is the last character to compute its likelihood. It should be lower than
 * @startPosition because it will start from @startPosition and go backwards until this value. 1-based
 * @longReadLength: length of the read to be corrected
 * @shortReadLength: length of the read that aligns to the longer read to correct
 */
bool fillBackwardMatrix(SequencingNode* graph, HMMParameters parameters, double* calculatedTransitionProbs,
                        double** backwardMatrix, const char* shortRead, int startPosition, int maxDistanceOnLongRead,
                        const int longReadLength, const int shortReadLength){
    
    if(maxDistanceOnLongRead < 0 || maxDistanceOnLongRead >= startPosition) return false;
    if(GRAPH_SIZE(longReadLength, parameters.maxInsertion) < startPosition)
        startPosition = GRAPH_SIZE(longReadLength, parameters.maxInsertion);
    
    //which transitions requested from the previous time
    std::vector<TransitionInfoNode>* curTrSet = new std::vector<TransitionInfoNode>;
    //which transitions to be made for the next time
    std::vector<TransitionInfoNode>* nextTrSet = new std::vector<TransitionInfoNode>;
    std::vector<TransitionInfoNode>* tmpTrSet;
    
    bool* allowedParentStates = new bool[startPosition - maxDistanceOnLongRead + 1];
    std::vector<bool> hasStateBeenProcessedBefore;
    hasStateBeenProcessedBefore.resize(startPosition - maxDistanceOnLongRead + 1);
    std::fill_n(hasStateBeenProcessedBefore.begin(), hasStateBeenProcessedBefore.size(), false);
    std::fill_n(allowedParentStates, startPosition - maxDistanceOnLongRead + 1, false);
    
    //1-initialization
    int curTime = shortReadLength-1; //@curTime value is 0-based. So 0th index is the first character [T....1]
    insertNewBackwardTransitions(curTrSet, graph[startPosition], parameters.maxDeletion, parameters.maxInsertion);
    for(int curTransition = 0; curTransition < curTrSet->size(); ++curTransition){
        if(curTrSet->at(curTransition).toState > maxDistanceOnLongRead){
            int matchoff = MATCH_OFFSET(graph[curTrSet->at(curTransition).toState].getCharIndex(), 0,
                                        parameters.maxInsertion);
            //0->insertion, 1-> match, 2,3...->deletions
            int transitionIndex = (curTrSet->at(curTransition).from - matchoff)/(parameters.maxInsertion+1);
            backwardMatrix[curTime][curTrSet->at(curTransition).toState - maxDistanceOnLongRead] +=
            calculatedTransitionProbs[transitionIndex];
            
            insertNewBackwardTransitions(nextTrSet, graph[curTrSet->at(curTransition).toState], parameters.maxDeletion,
                                         parameters.maxInsertion);
        }
    }
    curTrSet->clear();
    
    findMaxValues(backwardMatrix[curTime], allowedParentStates, maxDistanceOnLongRead, startPosition,
                  parameters.filterSize);
    tmpTrSet = curTrSet;
    curTrSet = nextTrSet;
    nextTrSet = tmpTrSet;
    
    //2-recursion
    while (curTime > 0 && !curTrSet->empty()){
        curTime--;
        for(int curTransition = 0; curTransition < curTrSet->size(); ++curTransition){
            const TransitionInfoNode& frontTr = curTrSet->at(curTransition);
            if(allowedParentStates[frontTr.from - maxDistanceOnLongRead] && frontTr.toState > maxDistanceOnLongRead) {
                int matchoff = MATCH_OFFSET(graph[frontTr.toState].getCharIndex(), 0, parameters.maxInsertion);
                //0->insertion, 1-> match, 2,3...->deletions
                int transitionIndex = (frontTr.from - matchoff)/(parameters.maxInsertion+1);
                backwardMatrix[curTime][frontTr.toState - maxDistanceOnLongRead] +=
                    backwardMatrix[curTime+1][frontTr.from - maxDistanceOnLongRead]*
                    calculatedTransitionProbs[transitionIndex]*
                    graph[frontTr.from].getEmissionProb(shortRead[curTime+1], parameters);
                
                if(!hasStateBeenProcessedBefore[frontTr.toState - maxDistanceOnLongRead]){
                    insertNewBackwardTransitions(nextTrSet, graph[frontTr.toState], parameters.maxDeletion,
                                                 parameters.maxInsertion);
                    hasStateBeenProcessedBefore[frontTr.toState - maxDistanceOnLongRead] = true;
                }
            }
        }
        curTrSet->clear();
        
        std::fill_n(hasStateBeenProcessedBefore.begin(), hasStateBeenProcessedBefore.size(), false);
        
        std::fill_n(allowedParentStates, startPosition - maxDistanceOnLongRead + 1, false);
        
        findMaxValues(backwardMatrix[curTime], allowedParentStates, maxDistanceOnLongRead, startPosition,
                      parameters.filterSize);
        tmpTrSet = curTrSet;
        curTrSet = nextTrSet;
        nextTrSet = tmpTrSet;
    }
    
    delete curTrSet; delete nextTrSet; delete[] allowedParentStates;
    
    return true;
}

int backtraceWithViterbi(SequencingNode* graph, HMMParameters parameters, double** transitionProbs,
                         std::pair<double, char>* emissionProbs, int numberOfStates, int longReadLength,
                         std::vector<std::map<int, int> >& backtrace){
    
    std::vector<bool> allowedParentStates;
    allowedParentStates.resize(numberOfStates);
    std::fill_n(allowedParentStates.begin(), allowedParentStates.size(), false);
    
    std::set<int>* curTrSet = new std::set<int>();
    std::set<int>* nextTrSet = new std::set<int>();
    std::set<int>* tmpTrSet;
    
    std::map<int, double>* prevViterbi = new std::map<int, double>();
    std::map<int, double>* curViterbi = new std::map<int, double>();
    std::map<int, double>* tmpViterbi;
    std::map<int, double>::iterator itVit;
    
    //initialization step t = 1 (curTime = 0)
    std::vector<TransitionInfoNode> curStateTransitions;
    insertNewForwardTransitions(&curStateTransitions, graph[0], parameters.maxDeletion, parameters.maxInsertion);
    backtrace.push_back(std::map<int, int >());
    
    int curTime = 0;
    prevViterbi->clear();
    while(!curStateTransitions.empty()){
        
        int next = curStateTransitions.begin()->toState;
        int transitionIndex = next/(parameters.maxInsertion+1); //0->insertion, 1-> match, 2,3...->deletions
        
        if(transitionProbs[0][transitionIndex] > 0.001){
            (*prevViterbi)[next] = log10(transitionProbs[0][transitionIndex]) + log10(emissionProbs[next].first);
            curTrSet->insert(next);
            allowedParentStates[next] = true;
            backtrace[curTime][next] = curStateTransitions.begin()->from; //0
        }
        curStateTransitions.erase(curStateTransitions.begin());
    }
    curStateTransitions.clear();
    curTime++;
    
    bool finished = false;
    int maxFinalStateTime = 0;
    double maxFinalStateValue = INT_MIN;
    while(!finished && curTime < (parameters.maxInsertion+1)*longReadLength){
        
        //new storage for the current time
        curViterbi->clear();
        backtrace.push_back(std::map<int, int>());
        
        for(std::set<int>::iterator itSet = (*curTrSet).begin(); itSet != curTrSet->end(); ++itSet){
            //calculate the viterbi values from this state to the states that it has transitions
            int fromState = (*itSet);
            //if this state is one of the major ones from the previous time
            if(allowedParentStates[fromState]){
                
                insertNewForwardTransitions(&curStateTransitions, graph[fromState], parameters.maxDeletion,
                                            parameters.maxInsertion);
                for(int toStateIndex = 0; toStateIndex < curStateTransitions.size(); ++toStateIndex){
                    
                    int toState = curStateTransitions[toStateIndex].toState;
                    int matchoff = MATCH_OFFSET(graph[fromState].getCharIndex(), 0, parameters.maxInsertion);
                    //0->insertion, 1-> match, 2,3...->deletions
                    int transitionIndex = (toState - matchoff)/(parameters.maxInsertion+1);
                    
                    if(toState < numberOfStates && transitionProbs[fromState][transitionIndex] > 0.001){
                        
                        double newViterbi = (*prevViterbi)[fromState] +
                            log10(transitionProbs[fromState][transitionIndex]) + log10(emissionProbs[toState].first);
                        
                        itVit = curViterbi->find(toState);
                        if(itVit == curViterbi->end() || (*curViterbi)[toState] < newViterbi){
                            //regular viterbi here
                            (*curViterbi)[toState] = newViterbi;
                            backtrace[curTime][toState] = fromState;
                            (*nextTrSet).insert(toState);
                        }
                    }
                }
                
                curStateTransitions.clear();
            }
        }
        
        std::fill_n(allowedParentStates.begin(), allowedParentStates.size(), false);
        (*curTrSet).clear();
        
        //find the most possible states that will contribute to the viterbi of the next time
        std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int> >,
                            std::greater<std::pair<double, int> > > maxValues;
        for(std::set<int>::iterator itSet = (*nextTrSet).begin(); itSet != (*nextTrSet).end(); ++itSet) {
            if(maxValues.size() < parameters.filterSize)
                maxValues.push(std::pair<double, int>((*curViterbi)[*itSet], *itSet));
            else if(maxValues.top().first < (*curViterbi)[*itSet]){
                maxValues.pop();
                maxValues.push(std::pair<double, int>((*curViterbi)[*itSet], *itSet));
            }
        }
        
        while(!maxValues.empty()){
            allowedParentStates.operator[](maxValues.top().second) = true;
            maxValues.pop();
        }
        
        if(allowedParentStates[numberOfStates-1] &&
           (maxFinalStateTime == 0 || (*curViterbi)[numberOfStates-1] > maxFinalStateValue)){
            maxFinalStateTime = curTime;
            maxFinalStateValue = (*curViterbi)[numberOfStates-1];
        }
        
        //stop condition: we've reached final state and max valued final state is "way" left behind the current time
        if(maxFinalStateTime > 0 && curTime > maxFinalStateTime + 100) finished = true;
        else curTime++;
        
        //swaps
        tmpTrSet = curTrSet; curTrSet = nextTrSet; nextTrSet = tmpTrSet;
        tmpViterbi = curViterbi; curViterbi = prevViterbi; prevViterbi = tmpViterbi;
    }
    
    if(maxFinalStateTime == 0)
        std::cout << "Problem with viterbi, a long read will not be corrected, original one to be reportded" << std::endl;
    
    curTrSet->clear();
    nextTrSet->clear();
    allowedParentStates.clear();
    delete curTrSet; delete nextTrSet;
    delete curViterbi; delete prevViterbi;
    
    return maxFinalStateTime; //this is the time where backtrace will start
}

/*
 * @longRead: read to be corrected. Graph is built based on this read. For each character of this read, there are
 * multiple states (match, insertion) representing what to do with this character or after.
 * @shortReads: reads that align to the long read. These reads are used to compute the forward and backward likelihood
 * values. So these are considered as observations on HMM.
 * @startPositions: This is 1-based value. Represents for which character of the long read a short reads starts aligning
 * @maxTransition: Indicates how many transitions can be made for each f/b likelihood computations.
 * @longReadLength: length of the long read
 * @shortReadLength: length of the short reads.
 * @noOfShortReads: size of @shortReads
 * @cigarString: not supported yet but this is a hint for the correction. Cigar string for the current short-long read
 * alignment
 */
void correctLongRead(HMMParameters parameters, const std::string& longRead, std::vector<ShortRead>& alignedShortReads,
                     std::string& correctedRead){

    //There are x2 states for each character in the long read. Additionally, there are start state, insertion state
    //before the first match state, and end state. So total size is (@longReadLength*2 + 3)
    int hmmGraphSize = (int)GRAPH_SIZE(longRead.size(), parameters.maxInsertion);
    SequencingNode* sequencingGraph = new SequencingNode[hmmGraphSize];
    
    //constructing hidden states according to the long read to be corrected
    sequencingGraph[0] = SequencingNode(0, 0, parameters.maxInsertion, '\0', longRead[0]);
    for(int curIn = 1; curIn <= parameters.maxInsertion; ++curIn)
        sequencingGraph[curIn] = SequencingNode(curIn, 0, parameters.maxInsertion, '\0', longRead[0]);
    for(int curCharacter = 1; curCharacter < (int)longRead.size(); ++curCharacter){
        sequencingGraph[MATCH_OFFSET(curCharacter, 0, parameters.maxInsertion)] =
            SequencingNode(MATCH_OFFSET(curCharacter, 0, parameters.maxInsertion), curCharacter, parameters.maxInsertion,
                           longRead[curCharacter-1], longRead[curCharacter]);
        for(int curIn = 1; curIn <= parameters.maxInsertion; ++curIn)
            sequencingGraph[INSERTION_OFFSET(curCharacter, 0, curIn, parameters.maxInsertion)] =
                SequencingNode(INSERTION_OFFSET(curCharacter, 0, curIn, parameters.maxInsertion), curCharacter,
                               parameters.maxInsertion, longRead[curCharacter-1], longRead[curCharacter]);
    }

    sequencingGraph[MATCH_OFFSET((int)longRead.size(), 0, parameters.maxInsertion)] =
            SequencingNode(MATCH_OFFSET((int)longRead.size(), 0, parameters.maxInsertion), (int)longRead.size(), parameters.maxInsertion,
                           longRead[(int)longRead.size()-1], '\0');
        for(int curIn = 1; curIn <= parameters.maxInsertion; ++curIn)
            sequencingGraph[INSERTION_OFFSET((int)longRead.size(), 0, curIn, parameters.maxInsertion)] =
                SequencingNode(INSERTION_OFFSET((int)longRead.size(), 0, curIn, parameters.maxInsertion), (int)longRead.size(),
                               parameters.maxInsertion, longRead[(int)longRead.size()-1], '\0');

    sequencingGraph[END_STATE(longRead.size(), parameters.maxInsertion)] =
        SequencingNode((int)END_STATE(longRead.size(),parameters.maxInsertion), (int)longRead.size()+1, parameters.maxInsertion,
                       '\0', '\0');
    
    std::vector<bool> shouldUseShort;
    shouldUseShort.resize(alignedShortReads.size());
    std::fill(shouldUseShort.begin(), shouldUseShort.end(), true);
    if(parameters.maxCoverage > 0){
        std::vector<std::vector<int>> alignedCount;
        alignedCount.resize(longRead.size());
        for(int curAlign=0; curAlign < alignedShortReads.size(); ++curAlign)
            alignedCount[alignedShortReads[curAlign].pos-1].push_back(curAlign);
        
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int> >,
                            std::greater<std::pair<int, int> > > selectedShorts;
        for(int curPos=0; curPos < longRead.size(); ++curPos){
            if(alignedCount[curPos].size() > parameters.maxCoverage){
                for(int i=0; i < alignedCount[curPos].size(); ++i)
                    selectedShorts.push(std::pair<int, int>(alignedShortReads[alignedCount[curPos][i]].editDistance,
                                                            alignedCount[curPos][i]));
                
                int covCount = 0;
                while(!selectedShorts.empty()){
                    if(covCount < parameters.maxCoverage){ selectedShorts.pop(); covCount++;}
                    else {
                        shouldUseShort[selectedShorts.top().second] = false;
                        selectedShorts.pop();
                    }
                }
            }
        }
    }
    
    std::vector<TransitionInfoNode> curStateTransitions;
    std::vector<TransitionInfoNode> stateTransitions; //all possible transitions from a state to a state.
    insertNewForwardTransitions(&stateTransitions, sequencingGraph[0], parameters.maxDeletion, parameters.maxInsertion);
    //pre calculated initial transition values for any state in hmm
    double* calculatedTransitionProbs = new double[stateTransitions.size()];
    for(int curTr = 0; curTr < stateTransitions.size(); ++curTr)
        calculatedTransitionProbs[curTr] = sequencingGraph[0].transitionProbFromThisNode(stateTransitions[curTr].toState,
                                                                                         parameters);
    
    //transitionProbs[from state][to state transition]
    double** transitionProbs = new double*[hmmGraphSize];
    //emissionProbs[state][character index]
    double** emissionProbs = new double*[hmmGraphSize];
    double* curStateTransitionLikelihood = new double[stateTransitions.size()];
    double* curStateEmissionProbs = new double[totalNuc];
    //how many times a state has been processed
    int* stateProcessedCount = new int[hmmGraphSize]; std::fill_n(stateProcessedCount, hmmGraphSize, 0);
    int** transitionProcessedCount = new int*[hmmGraphSize];
    
    for(int curState = 0; curState < hmmGraphSize; ++curState){
        transitionProcessedCount[curState] = new int[stateTransitions.size()];
        std::fill_n(transitionProcessedCount[curState], stateTransitions.size(), 0);
        transitionProbs[curState] = new double[stateTransitions.size()];
        emissionProbs[curState] = new double[totalNuc];
        std::fill_n(transitionProbs[curState], stateTransitions.size(), 0.0);
        std::fill_n(emissionProbs[curState], totalNuc, 0.0);
    }
    
    for(int curShortRead = 0; curShortRead < (int)alignedShortReads.size(); ++curShortRead){
        if(shouldUseShort[curShortRead]){
            //to which character should correction extend at maximum
            int shortLength = (int)alignedShortReads[curShortRead].shortRead.size();
            int maxTransition = shortLength + shortLength/3 + 1;
            int maxDistanceOnLongRead = std::min(alignedShortReads[curShortRead].pos + maxTransition, (int)longRead.size());
            //states prior to this wont be processed. offset value is to ignore these states
            int offset = MATCH_OFFSET(alignedShortReads[curShortRead].pos, -1, parameters.maxInsertion);
            //maximum number of states to be processed
            int fbMatrixSize = MATCH_OFFSET(maxDistanceOnLongRead, 1, parameters.maxInsertion) - offset;
            
            int j; //j value in i,j transitions
            double** forwardMatrix = new double*[shortLength]; //forwardMatrix[time][state]
            for(int i = 0; i < shortLength; i++){
                forwardMatrix[i] = new double[fbMatrixSize];
                std::fill_n(forwardMatrix[i], fbMatrixSize, 0);
            }
            double** backwardMatrix = new double*[shortLength]; //backwardMatrix[time][state]
            for(int i = 0; i < shortLength; i++){
                backwardMatrix[i] = new double[fbMatrixSize];
                std::fill_n(backwardMatrix[i], fbMatrixSize, 0);
            }
            
    #if defined(DEBUG_SHORT_READ_COMPLETION)
            std::cout << "Short Read (" << curShortRead << ") f/b calculation starts from " <<
            alignedShortReads[curShortRead].pos << "th to " << maxDistanceOnLongRead << " positions:" << std::endl;
            clock_t curShortClock = clock();
    #endif
            
            int startForBackward = fillForwardMatrix(sequencingGraph,parameters,calculatedTransitionProbs, forwardMatrix,
                                                     alignedShortReads[curShortRead].shortRead.c_str(),offset,
                                                     MATCH_OFFSET(maxDistanceOnLongRead, 1, parameters.maxInsertion),
                                                     (int)longRead.size(),shortLength);
            
            if(startForBackward != -1 && sequencingGraph[startForBackward].getCharIndex() > alignedShortReads[curShortRead].pos){
                startForBackward = MATCH_OFFSET(sequencingGraph[startForBackward].getCharIndex(), 1, parameters.maxInsertion);
                
                if(fillBackwardMatrix(sequencingGraph, parameters, calculatedTransitionProbs, backwardMatrix,
                                      alignedShortReads[curShortRead].shortRead.c_str(), startForBackward, offset,
                                      (int)longRead.size(), shortLength)){
                    //updating probabilities wrt the f/b matrices computed just now
                    for(int curState = INSERTION_OFFSET(alignedShortReads[curShortRead].pos, -1, 1, parameters.maxInsertion);
                        curState < startForBackward; ++curState){
                        
                        if(sequencingGraph[curState].isLastInsertionState())
                            //for the last insertion state, the insertion probs change so that it wont have insertion trans.
                            calculatedTransitionProbs[1] = parameters.matchTransition + parameters.insertionTransition;
                        
                        if(curState-offset < fbMatrixSize){
                            int matchoff = MATCH_OFFSET(sequencingGraph[curState].getCharIndex(), 0, parameters.maxInsertion);
                            std::fill_n(curStateTransitionLikelihood, stateTransitions.size(), 0.0);
                            std::fill_n(curStateEmissionProbs, totalNuc, 0.0);
                            insertNewForwardTransitions(&curStateTransitions, sequencingGraph[curState],
                                                        parameters.maxDeletion, parameters.maxInsertion);
                            
                            for(int t = 0; t < shortLength; ++t){
                                //transition probabilities
                                if(t < shortLength-1){
                                    for(int curTr = 0; curTr < curStateTransitions.size(); ++curTr){
                                        if(curStateTransitions.at(curTr).toState - offset < fbMatrixSize){
                                            j = curStateTransitions[curTr].toState;
                                            //0->insertion, 1-> match, 2,3...->deletions
                                            int transitionIndex = (j - matchoff)/(parameters.maxInsertion+1);
                                            curStateTransitionLikelihood[transitionIndex] +=
                                            forwardMatrix[t][curState-offset]*calculatedTransitionProbs[transitionIndex]*
                                            sequencingGraph[j].getEmissionProb(alignedShortReads[curShortRead].shortRead[t+1],
                                                                               parameters)*backwardMatrix[t+1][j-offset];
                                            
                                        }
                                    }
                                }
                                
                                //emission probabilities
                                char emitChar = (alignedShortReads[curShortRead].shortRead[t] != 'N')?
                                alignedShortReads[curShortRead].shortRead[t]:
                                (sequencingGraph[curState].isMatchState())?sequencingGraph[curState].getNucleotide():'\0';
                                Nucleotide chosenNuc = (emitChar == 'A' || emitChar == 'a')?A:
                                (emitChar == 'T' || emitChar == 't')?T:
                                (emitChar == 'G' || emitChar == 'g')?G:
                                (emitChar == 'C' || emitChar == 'c')?C:totalNuc;
                                if(chosenNuc < totalNuc)
                                    curStateEmissionProbs[chosenNuc] += forwardMatrix[t][curState-offset]*backwardMatrix[t][curState-offset];
                            }
                            curStateTransitions.clear();
                            double totalEmissionProbs = curStateEmissionProbs[A] + curStateEmissionProbs[T] +
                                                        curStateEmissionProbs[G] + curStateEmissionProbs[C];
                            double totalTransitionLikelihoods = 0;
                            double processedTransitionProb = 0;
                            for(int i = (sequencingGraph[curState].isLastInsertionState())?1:0; i < stateTransitions.size(); ++i){
                                if(curStateTransitionLikelihood[i] > 0 || i == 0){
                                    totalTransitionLikelihoods += curStateTransitionLikelihood[i];
                                    processedTransitionProb += calculatedTransitionProbs[i];
                                }
                            }
                            
                            if(totalEmissionProbs != 0 && curState < startForBackward){
                                if(totalTransitionLikelihoods != 0){
                                    for(int i = (sequencingGraph[curState].isLastInsertionState())?1:0; i < stateTransitions.size(); ++i){
                                        if(curStateTransitionLikelihood[i] > 0 || i == 0){
                                            transitionProbs[curState][i] += (curStateTransitionLikelihood[i]/totalTransitionLikelihoods)*
                                                                            processedTransitionProb;
                                            transitionProcessedCount[curState][i]++;
                                        }
                                    }
                                }
                                
                                emissionProbs[curState][A] += curStateEmissionProbs[A]/totalEmissionProbs;
                                emissionProbs[curState][T] += curStateEmissionProbs[T]/totalEmissionProbs;
                                emissionProbs[curState][G] += curStateEmissionProbs[G]/totalEmissionProbs;
                                emissionProbs[curState][C] += curStateEmissionProbs[C]/totalEmissionProbs;
                                
                                stateProcessedCount[curState]++;
                                
                            }
                        }
                        
                        if(sequencingGraph[curState].isLastInsertionState())
                            //for the last insertion state, the insertion probs change so that it wont have insertion transition.
                            //putting it back to normal now
                            calculatedTransitionProbs[1] = parameters.matchTransition;
                    }
                }
            }
            
    #if defined(DEBUG_SHORT_READ_COMPLETION)
            printf("Time taken for using the current short read: %.2fs\n", (double)(clock() - curShortClock)/CLOCKS_PER_SEC);
    #endif
            
            for(int i = 0; i < shortLength; ++i) { delete[] backwardMatrix[i]; delete[] forwardMatrix[i];}
            delete[] backwardMatrix;
            delete[] forwardMatrix;
        }
    }
    
    for(int curState = 0; curState < hmmGraphSize; ++curState){
        if(sequencingGraph[curState].isLastInsertionState()){
            //for the last insertion state, the insertion probs change so that it wont have insertion transition
            calculatedTransitionProbs[1] = parameters.matchTransition + parameters.insertionTransition;
        }
        if(stateProcessedCount[curState] > 0){ //if this state ever processed then its probs may need to be updated
            for(int curTransition = (sequencingGraph[curState].isLastInsertionState())?1:0; curTransition < stateTransitions.size();
                ++curTransition){
                transitionProbs[curState][curTransition] = (transitionProcessedCount[curState][curTransition] > 0)?
                transitionProbs[curState][curTransition]/transitionProcessedCount[curState][curTransition]:
                calculatedTransitionProbs[curTransition];
            }
            
            emissionProbs[curState][A] /= stateProcessedCount[curState];
            emissionProbs[curState][T] /= stateProcessedCount[curState];
            emissionProbs[curState][G] /= stateProcessedCount[curState];
            emissionProbs[curState][C] /= stateProcessedCount[curState];
        }else{ //initial probs to be set unless this state has been processed
            for(int curTr = (sequencingGraph[curState].isLastInsertionState())?1:0; curTr < stateTransitions.size(); ++curTr)
                transitionProbs[curState][curTr] = calculatedTransitionProbs[curTr];
            
            emissionProbs[curState][A] = sequencingGraph[curState].getEmissionProb('A', parameters);
            emissionProbs[curState][T] = sequencingGraph[curState].getEmissionProb('T', parameters);
            emissionProbs[curState][G] = sequencingGraph[curState].getEmissionProb('G', parameters);
            emissionProbs[curState][C] = sequencingGraph[curState].getEmissionProb('C', parameters);
        }
        
        if(sequencingGraph[curState].isLastInsertionState()){
            //for the last insertion state, the insertion probs change so that it wont have insertion transition. putting it
            //back to normal now
            calculatedTransitionProbs[1] = parameters.matchTransition;
        }
    }
    
    delete[] stateProcessedCount;
    delete[] curStateTransitionLikelihood;
    delete[] curStateEmissionProbs;
    
    std::pair<double, char>* maxEmissionProbs = new std::pair<double, char>[hmmGraphSize];
    for(int curState = 0; curState < hmmGraphSize; ++curState){
        if(emissionProbs[curState][A] >= std::max(emissionProbs[curState][T],
                                                  std::max(emissionProbs[curState][G], emissionProbs[curState][C]))){
            maxEmissionProbs[curState] = std::make_pair(emissionProbs[curState][A], 'A');
        }else if(emissionProbs[curState][T] >= std::max(emissionProbs[curState][A],
                                                        std::max(emissionProbs[curState][G], emissionProbs[curState][C]))){
            maxEmissionProbs[curState] = std::make_pair(emissionProbs[curState][T], 'T');
        }else if(emissionProbs[curState][G] >= std::max(emissionProbs[curState][A],
                                                        std::max(emissionProbs[curState][T], emissionProbs[curState][C]))){
            maxEmissionProbs[curState] = std::make_pair(emissionProbs[curState][G], 'G');
        }else{
            maxEmissionProbs[curState] = std::make_pair(emissionProbs[curState][C], 'C');
        }
    }
    
    std::vector<std::map<int, int>> backtrace;
    int timeToBackTrace = backtraceWithViterbi(sequencingGraph, parameters, transitionProbs, maxEmissionProbs, hmmGraphSize,
                                               (int)longRead.size(), backtrace);
    
    for(int i = 0; i < hmmGraphSize; ++i){
        delete[] transitionProbs[i]; delete[] emissionProbs[i]; delete[] transitionProcessedCount[i];
    }
    delete[] transitionProbs; delete[] emissionProbs; delete[] transitionProcessedCount;
    
    correctedRead.clear();
    correctedRead.resize(timeToBackTrace);
    int lastIndexUsed = (int)END_STATE(longRead.size(), parameters.maxInsertion);
    for(int curChar = timeToBackTrace-1; curChar >= 0; --curChar){
        int curIndex = backtrace[curChar+1][lastIndexUsed];
        lastIndexUsed = curIndex;

        correctedRead[curChar] = maxEmissionProbs[curIndex].second;
    }
    
#if defined(DEBUG_HMM_RESULT) || defined(DEBUG_HMM_DETAILED)
    std::cout << correctedRead << std::endl;
#endif
    
    delete[] maxEmissionProbs;
    delete[] sequencingGraph;
    delete[] calculatedTransitionProbs;
}

void fillCompressedReadIndexArray(std::string uncompressed, std::vector<int>& index){
    
    index.clear();
    int curCompressedChar = 0;
    int curUncompressedChar = 0;
    index.push_back(0);
    while(curUncompressedChar+1 < uncompressed.size() && uncompressed[curUncompressedChar] == uncompressed[curUncompressedChar+1]){
        ++index.at(curCompressedChar);
        ++curUncompressedChar;
    }
    
    ++curUncompressedChar;
    ++curCompressedChar;
    while(curUncompressedChar < uncompressed.size()){
        index.push_back(index.at(curCompressedChar-1));
        while(curUncompressedChar+1 < uncompressed.size() &&
              uncompressed[curUncompressedChar] == uncompressed[curUncompressedChar+1]){
            ++index.at(curCompressedChar);
            ++curUncompressedChar;
        }
        
        ++curUncompressedChar;
        ++curCompressedChar;
    }
    
    for(int i = (int)index.size()-1; i > 0; --i) index[i] = index[i-1];
    index.at(0) = 0;
}

void correctReadThreadPool(HMMParameters parameters, FaiIndex& longIndex, FaiIndex& shortIndex, uint64_t size, uint& readCount,
                           int mapQ, BamFileIn& alignmentFileIn, BamAlignmentRecord& record, std::fstream& correctedReadStream,
                           std::fstream& coverageStream, bool shouldCalculateCoverage, bool shouldCompress, bool shouldQuite){
    
    int index = 0; //0-based. refers to the current read id
    CharString tagKey = "NM";
    //calculate coverage as you go
    while(index < size){
        bool shouldCorrect = false; //is there at least one alignment for the current read?
        int coverageCount = 0;
        int alignmentCount = 0;
        CharString curSeq; //original long read sequence
        CharString curSeqId; //original long read id
        std::string correctedSeq; //corrected sequence (if there is at least one alignment or just the original one)
        std::vector<bool> coverage; //which positions covered?
        std::vector<ShortRead> shortReads; //aligned short reads
        
        {//block for obtaining new index for the next long read
            std::lock_guard<std::mutex> lk(indexMutex);
            index = readCount;
            if(index < size){
                //read the current long read's original sequence
                readSequence(curSeq, longIndex, index);
                //read the current long read's original id
                curSeqId = sequenceName(longIndex, index);
                int seqId = -1;

                //id of the last "unproccessed" record in the alignment file. Id-names 0-based as well
                if(!atEnd(alignmentFileIn))
                	seqId = std::stoi(toCString(getContigName(record, alignmentFileIn)));
                
                //did we hit the current alignment or there is no alignment for the current long read?
                //since the alignment file is **sorted** we can assume there is no alignment if index+1 < seqId
                if(seqId < 0 || index < seqId) shouldCorrect = false;
                else if(index == seqId){ //there is an alignment. so read the aligned short reads.
                    coverage.resize(length(curSeq));
                    std::fill(coverage.begin(), coverage.end(), false);
                    int rId = record.rID;
                    std::string curSeqStd = toCString(curSeq);
                    std::vector<int> compressedPosIndex;
                    if(shouldCompress){
                        compressedPosIndex.resize(curSeqStd.size());
                        fillCompressedReadIndexArray(curSeqStd, compressedPosIndex);
                    }
                    while(rId == record.rID){
                        if(!hasFlagUnmapped(record) && record.mapQ >= mapQ){
                            std::stringstream seqStream;
                            String<Dna> alignedString;
                            int qId = std::stoi(toCString(record.qName));
                            readSequence(alignedString, shortIndex, qId);
                            if(hasFlagRC(record)) reverseComplement(alignedString);
                            seqStream << alignedString;
                            
                            if(seqStream.str().size() > 0){
                                std::stringstream cigarStream;
                                for(unsigned i  = 0; i < length(record.cigar); ++i){
                                    cigarStream << record.cigar[i].count;
                                    cigarStream << record.cigar[i].operation;
                                }
                                
                                int pos = shouldCompress?compressedPosIndex[record.beginPos] + record.beginPos + 1:
                                          record.beginPos + 1;
                                shortReads.push_back(ShortRead(seqStream.str(), pos, cigarStream.str()));
                                BamTagsDict tagsDict(record.tags);
                                unsigned id;
                                if(findTagKey(id, tagsDict, tagKey)){
                                    int32_t tagValue = 0;
                                    if (extractTagValue(tagValue, tagsDict, id))
                                        shortReads.back().editDistance = tagValue;
                                }
                                if(shouldCalculateCoverage){
                                    int offset = shortReads.back().pos-1;
                                    for(int j = 0; j < shortReads.back().shortRead.size() && offset + j < coverage.size(); ++j){
                                        if(!coverage.at(offset + j)) coverageCount++;
                                        coverage.at(offset + j) = true;
                                    }
                                    alignmentCount++;
                                }
                            }
                        }
                        
                        if(!atEnd(alignmentFileIn)) readRecord(record, alignmentFileIn);
                        else rId = -1;
                    }
                    shouldCorrect = true;
                }
                if(shouldCalculateCoverage){
                    coverageStream << toCString(curSeqId) << " " << (double)coverageCount/length(curSeq) << " " <<
                    alignmentCount << std::endl;
                }
                readCount++; //refers to the next long read if left any
                if(!shouldQuite && readCount%1000 == 0)
                    std::cout << ((double)readCount/size)*100 << "% " << std::endl;
            }
        }//end of obtaining index block
        
        if(index < size){
            if(shouldCorrect){
                std::string corr;
                correctLongRead(parameters, toCString(curSeq), shortReads, corr);
                correctedSeq = corr;
            }else correctedSeq = toCString(curSeq);
            
            //write the result
            {
                std::lock_guard<std::mutex> lk(outputFileMutex);
                SeqFileOut correctedReadsOut(correctedReadStream, Fasta());
                writeRecord(correctedReadsOut, curSeqId, correctedSeq);
            }
        }
        shouldCorrect = false;
    }//end of while
}

bool correctionPhase(HMMParameters& parameters, std::string longReadFile, std::string shortReadFile, std::string alignmentFile,
                     std::string outputFile, unsigned mapQ, int thread, bool shouldOutputCoverage, bool shouldCompress,
                     bool shouldQuite){
    
    BamFileIn alignmentFileIn;
    BamAlignmentRecord record;
    if (!open(alignmentFileIn, alignmentFile.c_str())){
        std::cerr << "ERROR: Could not open " << alignmentFile << std::endl;
        return false;
    }
    try{
        BamHeader header;
        readHeader(header, alignmentFileIn);
        if(!atEnd(alignmentFileIn)) readRecord(record, alignmentFileIn);
    }catch(Exception const & e){
        std::cerr << "ERROR: " << e.what() << std::endl;
        return false;
    }
    
    FaiIndex longIndex;
    if(!build(longIndex, longReadFile.c_str())){
        std::cerr << "ERROR: Could not build FAI index for file " << longReadFile << std::endl;
        return false;
    }
    uint64_t size = numSeqs(longIndex);
    uint readCount = 0;
    
    FaiIndex shortIndex;
    if(!build(shortIndex, shortReadFile.c_str())){
        std::cerr << "ERROR: Could not build FAI index for file " << shortReadFile << std::endl;
        return false;
    }
    
    std::fstream correctedReadStream;
    try{
        correctedReadStream.open(outputFile.c_str(), std::fstream::out | std::fstream::trunc);
    }catch(std::ios_base::failure e){
        std::cerr << "Could not open " << outputFile << std::endl;
        return false;
    }
    
    std::fstream coverageStream;
    if(shouldOutputCoverage){
        std::string coverageFile = outputFile + ".coverage";
        try{
            coverageStream.open(coverageFile.c_str(), std::fstream::out | std::fstream::trunc);
            coverageStream << "#read_id coverage_fraction number_of_aligned_reads" << std::endl;
        }catch(std::ios_base::failure e){
            std::cerr << "Could not open " << coverageFile << ". Won't report coverage." << std::endl;
            shouldOutputCoverage = false;
        }
    }
    
    if(!shouldQuite) std::cout << "Correction has begun..." << std::endl;
    std::vector<std::thread> threads;
    for(int i = 0; i < thread; ++i){
        threads.push_back(std::thread(correctReadThreadPool, parameters, std::ref(longIndex), std::ref(shortIndex), size,
                                      std::ref(readCount), mapQ, std::ref(alignmentFileIn), std::ref(record),
                                      std::ref(correctedReadStream), std::ref(coverageStream), shouldOutputCoverage,
                                      shouldCompress, shouldQuite));
    }
    for(int i = 0; i < thread; ++i) threads[i].join();
    threads.clear();
    correctedReadStream.close();
    if(!shouldQuite) std::cout << std::endl << "Results have been written under " << outputFile << std::endl;
    
    return true;
}

bool reFactorReads(std::string inputFile, std::string outputFile, bool shouldCompress, bool rename = false){
    
    SeqFileIn seqFileIn;
    SeqFileOut compressedFileOut;
    if(!open(seqFileIn, inputFile.c_str())){
        std::cerr << "ERROR: Could not open the file.: " << inputFile << std::endl;
        return false;
    }else if(!open(compressedFileOut, outputFile.c_str())){
        std::cerr << "ERROR: Could not open " << outputFile << " to save compressed reads" << std::endl;
        return false;
    }
    
    int count = 0; //0-based ids
    CharString id; //read id
    Dna5String seq; //read sequence
    std::string seqStr; //read sequence in std
    int nonN; //count of non-n characters in the compressed read
    try{
        while(!atEnd(seqFileIn)){
            readRecord(id, seq, seqFileIn);
            assign(seqStr, seq);
            if(rename){
                if(shouldCompress)
                    writeRecord(compressedFileOut, count++, compress(seqStr, nonN));
                else
                    writeRecord(compressedFileOut, count++, seqStr);
            }
            else{
                if(shouldCompress)
                    writeRecord(compressedFileOut, id, compress(seqStr, nonN));
                else
                    writeRecord(compressedFileOut, id, seqStr);
            }
            
        }
    }catch (Exception const & e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return false;
    }
    
    return true;
}

std::string getBasename(std::string fullPath){
    size_t baseStart = fullPath.rfind(SEPERATOR, fullPath.length());
    if( baseStart != std::string::npos)
        return fullPath.substr(baseStart+1, fullPath.length() - baseStart);
    return fullPath;
}

bool preprocessPhase(std::string longReadFile, std::vector<seqan::CharString> shortReads, std::string outputDirectory,
                     int minNonN, bool shouldCompress, bool applyBloomFilter){
    
    //compressing and writing compressed long reads
    std::string reFactoredLongRead(outputDirectory);
    reFactoredLongRead.append(SEPERATOR);
    if(shouldCompress) reFactoredLongRead.append("compressed_" + getBasename(longReadFile));
    else reFactoredLongRead.append("renamed_" + getBasename(longReadFile));
    if(shouldCompress && !reFactorReads(longReadFile, reFactoredLongRead, true, true)) return false;
    else if(!shouldCompress && !reFactorReads(longReadFile, reFactoredLongRead, false, true)) return false;
    
    //renamed original short reads will be placed here
    std::string reshapedShortRead(outputDirectory);
    reshapedShortRead.append(SEPERATOR);
    reshapedShortRead.append("short.fasta");
    //renamed compressed short reads will be placed here
    std::string reshapedCompressedShortRead(outputDirectory);
    reshapedCompressedShortRead.append(SEPERATOR);
    reshapedCompressedShortRead.append("compressed_short.fasta");
    
    SeqFileOut reshapedShortReadOut;
    SeqFileOut reshapedCompressedShortReadOut;
    if(!open(reshapedShortReadOut, reshapedShortRead.c_str())){
        std::cerr << "ERROR: Could not open " << reshapedShortRead << std::endl;
        return false;
    }else if(shouldCompress && !open(reshapedCompressedShortReadOut, reshapedCompressedShortRead.c_str())){
        std::cerr << "ERROR: Could not open " << reshapedCompressedShortRead << std::endl;
        return false;
    }
    
    int curShortReadId = 0; //0-based ids
    int nonN = 0; //count of non-N characters in compressed short read
    CharString id; //short read id
    Dna5String seq; //short read sequence
    std::string seqStr; //short read std sequence
    
    bloom_parameters bloomParams;
    if(applyBloomFilter) bloomParams.compute_optimal_parameters();
    bloom_filter filter(bloomParams);
    
    //for each specified short read file:
    //read a short read & compress it & write original and compressed if compressed has more than minNonN non 'N' characters
    for(int curShort = 0; curShort < shortReads.size(); ++curShort){
        std::string curShortFile = toCString(shortReads[curShort]);
        SeqFileIn seqFileIn;
        if(!open(seqFileIn, curShortFile.c_str())){
            std::cerr << "ERROR: Could not open the file.: " << curShortFile << std::endl;
            return false;
        }

        try{
            while(!atEnd(seqFileIn)){
                readRecord(id, seq, seqFileIn);
                assign(seqStr, seq);
                std::string compressed;
                if(shouldCompress) compressed = compress(seqStr, nonN);
                else{
                    nonN = 0;
                    for(int i = 0; i < seqStr.size(); ++i)
                        if(seqStr[i] != 'N') nonN++;
                }
                if(nonN >= minNonN && (!applyBloomFilter || !filter.contains(compressed))){
                    writeRecord(reshapedShortReadOut, curShortReadId, seqStr);
                    if(shouldCompress) writeRecord(reshapedCompressedShortReadOut, curShortReadId, compressed);
                    if(applyBloomFilter) filter.insert(compressed);
                    curShortReadId++;
                }
            }
        }catch (Exception const & e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
            return false;
        }
        close(seqFileIn);
    }
    
    std::cout << "Preprocessing is done. Now you have two external work before running correction phase:" << std::endl;
    if(shouldCompress){
        std::cout << "(1)Align compressed short reads to compressed long reads using any aligner. Note that multiple alignment"
        " option of the aligner should be enabled to benefit from compression." << std::endl;
        std::cout << "Compressed short reads file is at: " << reshapedCompressedShortRead << std::endl;
        std::cout << "Compressed long reads file is at: " << reFactoredLongRead << std::endl;
    }else{
        std::cout << "(1)Align short reads to long reads using any aligner." << std::endl;
        std::cout << "Short reads file is at: " << reshapedShortRead << std::endl;
        std::cout << "Long reads file is at: " << reFactoredLongRead << std::endl;
    }
    std::cout << "(2)Use \"samtools sort -m 8G -l 0 | samtools rmdup -S\" commands to remove the duplicates in your alignment file (e.g. see utils/afteralignment.sh or utils/runBowtieRmDup.sh)" << std::endl;
    std::cout << "(3)Run correction phase of Hercules using: " << std::endl;
    std::cout << "-The alignment file that you removed its duplicates in (2)" << std::endl;
    std::cout << "-Original long reads file at " << longReadFile << std::endl;
    std::cout << "-Short reads file at: " << reshapedShortRead << std::endl;
    
    close(reshapedShortReadOut);
    close(reshapedCompressedShortReadOut);
    
    return true;
}

int main(int argc, const char** argv) {
    
    CommandLineOptions options;
    seqan::ArgumentParser::ParseResult parseRes = parseInitialCommand(options, argc, argv);
    
    if(parseRes != seqan::ArgumentParser::PARSE_OK)
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;

    if(options.phase == PREPROCESS_FLAG){
        std::cout << "Preprocessing has begun..." << std::endl;
        if(!preprocessPhase(toCString(options.longInputFile), options.shortReads, toCString(options.output), options.nonNCount,
           options.shouldCompress, options.shouldApplyBloomFilter))
            std::cout << "Failure on preprocessing step. Check out the errors." << std::endl;
    }else if(options.phase == CORRECT_FLAG){
        HMMParameters parameters(options.filterSize, options.maxCoverage, options.maxDeletion, options.maxInsertion,
                                 options.matchTransition, options.insertionTransition, options.deletionTransitionFactor,
                                 options.matchEmission);
        
        if(!options.shouldQuite){
            std::string coverageAnswer = options.shouldOutputCoverage?"Yes":"No";
            std::cout << "Long Read: " << toCString(options.longInputFile) << std::endl << "Short read: " <<
            toCString(options.shortReads.front()) << std::endl << "Alignment File: " << toCString(options.alignmentFile) <<
            std::endl << "Output file: " << options.output << std::endl << "Will output coverage?: " << coverageAnswer <<
            std::endl << "Min mapping quality: " << options.mapQ << std::endl << "Max Coverage: " << parameters.maxCoverage <<
            std::endl << "Filter size: " << parameters.filterSize << std::endl << "Maximum insertion: " << parameters.maxInsertion
            << std::endl << "Maximum deletion: " << parameters.maxDeletion << std::endl << "Match transition probability: " <<
            parameters.matchTransition << std::endl << "Insertion transition probability: " << parameters.insertionTransition <<
            std::endl << "Deletion transition probability: " << parameters.deletionTransition << std::endl <<
            "Deletion transition factor: " << parameters.deletionTransitionFactor << std::endl << "Match emission probability: "
            << parameters.matchEmission << std::endl << "Mismatch emission probability: " << parameters.mismatchEmission <<
            std::endl << "Insertion emission probability: " << parameters.insertionEmission << std::endl << "Max thread: " <<
            options.maxThread << std::endl;
        }
        
        correctionPhase(parameters, toCString(options.longInputFile), toCString(options.shortReads.front()),
                        toCString(options.alignmentFile), toCString(options.output), options.mapQ, options.maxThread,
                        options.shouldOutputCoverage, options.shouldCompress, options.shouldQuite);
    }
    
    return 0;
}
