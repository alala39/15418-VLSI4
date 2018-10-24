#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <assert.h>
#include "mpi.h"
#include "wire_route.h"

#define root 0
#define tag 0

inline int max (int a, int b){return a >= b ? a : b;}

inline int min (int a, int b){return a <= b ? a : b;}

int *globalWires;

int *costs;

// Initialize problem
// Init cost array in root and broadcast to other processes
void init(int numRows, int numCols)
{
    //TODO Implement code here
    costs = calloc(numRows*numCols,sizeof(int));
    //int *costs = calloc(numRows*numCols,sizeof(int));
    //MPI_Bcast(costs,numRows*numCols,MPI_INT,root,MPI_COMM_WORLD);
    //Initializes the cost array. Need to redo this as it creates one for each process.
    //Just call this if the procID is 0 (root)
    //return costs;
}

// Initialize a given wire
static inline void initWire(int wireIndex, int x1, int y1, int x2, int y2,wire_t *wireArray)
{
    //TODO Implement code here
    //wireArray[wireIndex] = {x1,y1,y1,y2,undef,undef,undef,undef};
    wireArray[wireIndex].x1 = x1;
    wireArray[wireIndex].y1 = y1;
    wireArray[wireIndex].x2 = x2;
    wireArray[wireIndex].y2 = y2;
    wireArray[wireIndex].b1x = undef;
    wireArray[wireIndex].b1y = undef;
    wireArray[wireIndex].b2x = undef;
    wireArray[wireIndex].b2y = undef;
}

// Return number of rows
static inline int getNumRows()
{
    //TODO Implement code here
    return 0;
}

// Return number of cols
static inline int getNumCols()
{
    //TODO Implement code here
    return 0;
}

// Return delta
static inline int getDelta()
{
	//TODO Implement code here
    return 0;
}

// Return number of wires
static inline int getNumWires()
{
    //TODO Implement code here
    return 0;
}

// Get cost array entry
static inline int getCost(int row, int col,int numCols)
{
    //TODO Implement code here
    return costs[row * numCols + col];
}

// Get a wire placement. Returns number of points (should be 2-4 for 0-2 bends).
static inline int getWire(int wireIndex, int* x1, int* y1, int* x2, int* y2, int* x3, int* y3, int* x4, int* y4,wire_t *wireArray)
{
    //TODO Implement code here
    *x1 = wireArray[wireIndex].x1;
    *y1 = wireArray[wireIndex].y1;
    *x2 = wireArray[wireIndex].x2;
    *y2 = wireArray[wireIndex].y2;
    *x3 = wireArray[wireIndex].b1x;
    *y3 = wireArray[wireIndex].b1y;
    *x4 = wireArray[wireIndex].b2x;
    *y4 = wireArray[wireIndex].b2y;
    if(*x3==undef && *y3==undef)
    {
        return 2;
    }
    if(*x4==undef && *y4==undef)
    {
        return 3;
    }
    else
    {
        return 4;
    }
}

void removeEqualXWire(int x, int y1,int y2,int *costs,int numRows,int numCols,int by)
{
    int start = min(y1,y2);
    int end = max(y1,y2);
    while(start <= end)
    {
        if(costs[start * numCols + x] && start!=by)
        {
            costs[start * numCols + x]-=1;
        }
        start+=1;
    }
}

void removeEqualYWire(int x1,int x2,int y,int *costs,int numRows,int numCols,int bx)
{
    int start = min(x1,x2);
    int end = max(x1,x2);
    while(start <= end)
    {
        if(costs[y * numCols + start] && start!=bx)
        {
            costs[y * numCols + start]-=1;
        }
        start+=1;
    }
}

void placeXWire(int x,int y1,int y2,int *costs,int numRows,int numCols)
{
    int start = min(y1,y2);
    int end = max(y1,y2);
    while(start<=end)
    {
        costs[start * numCols + x]+=1;
        start+=1;
    }
}

void placeYWire(int x1,int x2,int y,int *costs,int numRows,int numCols)
{
    int start = min(x1,x2);
    int end = max(x1,x2);
    while(start <= end)
    {
        costs[y * numCols + start]+=1;
        start+=1;
    }
}

void placeWire(wire_t *wire,int *costs,int numRows,int numCols)
{
    int x1 = wire->x1;
    int y1 = wire->y1;
    int x2 = wire->x2;
    int y2 = wire->y2;
    int b1x = wire->b1x;
    int b1y = wire->b1y;
    int b2x = wire->b2x;
    int b2y = wire->b2y;

    if(b1x!=undef && b1y!=undef)
    {
        if(x1==b1x)
        {
            placeXWire(x1,y1,b1y,costs,numRows,numCols);
        }
        else
        {
            //printf("x1: %d y1: %d bx1: %d by1: %d\n",x1,y1,b1x,b1y);
            assert(y1==b1y);
            placeYWire(x1,b1x,y1,costs,numRows,numCols);
        }
        if(b2x!=undef && b2y!=undef)
        {
            if(b1x==b2x)
            {
                placeXWire(b1x,b1y,b2y,costs,numRows,numCols);
            }
            else
            {
                assert(b1y==b2y);
                placeYWire(b1x,b2x,b1y,costs,numRows,numCols);
            }
            if(b2x==x2)
            {
                placeXWire(x2,b2y,y2,costs,numRows,numCols);
            }
            else
            {
                assert(b2y==y2);
                placeYWire(b2x,x2,y2,costs,numRows,numCols);
            }
        }
        else
        {
            if(b1x==x2)
            {
                placeXWire(x2,b1y,y2,costs,numRows,numCols);
            }
            else
            {
                assert(b1y==y2);
                placeYWire(b1x,x2,y2,costs,numRows,numCols);
            }
        }
    }
    else
    {
        if(x1==x2)
        {
            placeXWire(x1,y1,y2,costs,numRows,numCols);
        }
        else
        {
            assert(y1==y2);
            placeYWire(x1,x2,y1,costs,numRows,numCols);
        }
    }
    if(b1x!=undef && b1y!=undef)
    {
        costs[b1y * numCols + b1x]-=1;
    }
    if(b2x!=undef && b2y!=undef)
    {
        costs[b2y * numCols + b2x]-=1;
    }
}

int costPerLine(int x1,int x2, int y1, int y2,int *costs,int numRows,int numCols)
{
    int result = 0;
    if(x1==x2)
    {
        int start = min(y1,y2);
        int end = max(y1,y2);
        while(start <= end)
        {
            if(costs[start * numCols + x1]>result)
            {
                result = costs[start * numCols + x1];
            }
            start+=1;
        }
    }
    else
    {
        assert(y1==y2);
        int start = min(x1,x2);
        int end = max(x1,x2);
        while(start <= end)
        {
            if(costs[y1 * numCols + start]>result)
            {
                result = costs[y1 * numCols + start];
            }
            start+=1;
        }
    }
    return result;
}

int costEqualX(int x,int y1,int y2,int *costs,int numRows,int numCols,int delta,wire_t *wire)
{
    int result = 0;
    int start = min(y1,y2);
    int end = max(y1,y2);
    result = costPerLine(x,x,y1,y2,costs,numRows,numCols);
    wire_t w = {x,y1,x,y2,undef,undef,undef,undef};
    *wire = w;
    if(y1==y2 || delta==0)
    {
        return result;
    }
    int curMin = -1;
    for(int i = 1;i<=delta/2;i++)
    {
        if(x+i>=numCols || x-i < 0)
        {
            continue;
        }
        int res1 = max(costPerLine(x,x+i,start,start,costs,numRows,numCols),max(costPerLine(x+i,x+i,start,end,costs,numRows,numCols),costPerLine(x+i,x,end,end,costs,numRows,numCols)));
        int res2 = max(costPerLine(x,x-i,start,start,costs,numRows,numCols),max(costPerLine(x-i,x-i,start,end,costs,numRows,numCols),costPerLine(x-i,x,end,end,costs,numRows,numCols)));
        result = min(result,min(res1,res2));
        if(y1==46 && x+i==46)
        {
            //printf("This is wrong\n");
        }
        else if(y1==46 && x-i==46)
        {
            //printf("This is also wrong\n");
        }
        if(result==res1 && (result<=curMin || curMin==-1))
        {
            wire_t w = {x,y1,x,y2,x+i,y1,x+i,y2};
            *wire = w;
            curMin = res1;
        }
        else if(result==res2 && (result<=curMin || curMin==-1))
        {
            wire_t w = {x,y1,x,y2,x-i,y1,x-i,y2};
            *wire = w;
            curMin = res2;
        }
        else
        {
            curMin = result;
        }
    }
    return curMin;
}

int costEqualY(int x1,int x2,int y,int *costs,int numRows,int numCols,int delta,wire_t *wire)
{
    int result = 0;
    int start = min(x1,x2);
    int end = max(x1,x2);
    result = costPerLine(x1,x2,y,y,costs,numRows,numCols);
    wire_t w = {x1,y,x2,y,undef,undef,undef,undef};
    *wire = w;
    if(x1==x2 || delta==0)
    {
        return result;
    }
    int curMin = -1;
    for(int i = 1;i<=delta/2;i++)
    {
        if(y+i>=numRows || y-i < 0)
        {
            continue;
        }
        int res1 = max(costPerLine(start,start,y,y+i,costs,numRows,numCols),max(costPerLine(start,end,y+i,y+i,costs,numRows,numCols),costPerLine(end,end,y+i,y,costs,numRows,numCols)));
        int res2 = max(costPerLine(start,start,y,y-i,costs,numRows,numCols),max(costPerLine(start,end,y-i,y-i,costs,numRows,numCols),costPerLine(end,end,y-i,y,costs,numRows,numCols)));
        result = min(result,min(res1,res2));
        if(y==46 && x1==46)
        {
            //printf("This is not the right place\n");
        }
        if(result==res1 && (result<=curMin || curMin==-1))
        {
            wire_t w = {x1,y,x2,y,x1,y+i,x2,y+i};
            *wire = w;
            curMin = res1;
        }
        else if(result==res2 && (result<=curMin || curMin==-1))
        {
            wire_t w = {x1,y,x2,y,x1,y-i,x2,y-i};
            *wire = w;
            curMin = res2;
        }
        else
        {
            curMin = result;
        }
    }
    return curMin;
}

int checkAllHorizontal(int x1,int x2,int y1,int y2,int *costs,int numRows,int numCols,int delta,wire_t *wire)
{
    int startX = min(x1,x2);
    int endX = max(x1,x2);
    int lowerBound = startX - (delta/2);
    int upperBound = endX + (delta/2);
    int startY;
    int endY;
    int curMin = -1;
    if(startX==x1)
    {
        startY = y1;
        endY = y2;
    }
    else
    {
        assert(startX==x2);
        startY = y2;
        endY = y1;
    }
    for(int i = lowerBound;i<=upperBound;i++)
    {
        if(i>=numCols || i < 0)
        {
            continue;
        }
        int rec = max(costPerLine(startX,i,startY,startY,costs,numRows,numCols),max(costPerLine(i,i,startY,endY,costs,numRows,numCols),costPerLine(i,endX,endY,endY,costs,numRows,numCols)));
        if(startY==46 && i==46)
        {
            //printf("This is not the right place please do not\n");
        }
        if(rec<=curMin || curMin==-1)
        {
            if(i!=startX && i!=endX)
            {
                curMin = rec;
                wire_t w = {startX,startY,endX,endY,i,startY,i,endY};
                *wire = w;
            }
            else if(i==startX)
            {
                curMin = rec;
                wire_t w = {startX,startY,endX,endY,i,endY,undef,undef};
                *wire = w;
            }
            else if(i==endX)
            {
                curMin = rec;
                wire_t w = {startX,startY,endX,endY,i,startY,undef,undef};
                *wire = w;
            }
        }
    }
    return curMin;
}

int checkAllVertical(int x1,int x2,int y1,int y2,int *costs, int numRows,int numCols,int delta,wire_t *wire)
{
    int startY = min(y1,y2);
    int endY = max(y1,y2);
    int lowerBound = startY - (delta/2);
    int upperBound = endY + (delta/2);
    int startX;
    int endX;
    int curMin = -1;
    if(startY==y1)
    {
        startX = x1;
        endX = x2;
    }
    else
    {
        assert(startY==y2);
        startX = x2;
        endX = x1;
    }
    for(int i = lowerBound;i<=upperBound;i++)
    {
        if(i>=numRows || i<0)
        {
            continue;
        }
        int rec = max(costPerLine(startX,startX,startY,i,costs,numRows,numCols),max(costPerLine(startX,endX,i,i,costs,numRows,numCols),costPerLine(endX,endX,i,endY,costs,numRows,numCols)));
        if(startY==46 && startX==46)
        {
            //printf("Please not here\n");
        }
        if(rec<=curMin || curMin==-1)
        {
            if(i!=startY && i!=endY)
            {
                curMin = rec;
                wire_t w = {startX,startY,endX,endY,startX,i,endX,i};
                *wire = w;
            }
            else if(i==startY)
            {
                curMin = rec;
                wire_t w = {startX,startY,endX,endY,endX,i,undef,undef};
                *wire = w;
            }
            else if(i==endY)
            {
                curMin = rec;
                wire_t w = {startX,endX,startY,endY,startX,i,undef,undef};
                *wire = w;
            }
        }
    }
    return curMin;
}

void mergeCost(int *aggCost,int *newCost,int numRows,int numCols)
{
    for(size_t i = 0;i<numRows;i++)
    {
        for(size_t j = 0;j<numCols;j++)
        {
            //printf("MERGE THOSE COSTS\n");
            aggCost[i * numCols + j]+=newCost[i * numCols + j];
        }
    }
    //This seqential and is super slow
}

//The function below will be called after each iteration. All non root nodes will send cost array to root node who will consolidate arrays into one.
void combineCostArrays(int procID,int *costs,int numRows,int numCols,int nproc)
{
    if(procID!=root)
    {
        MPI_Send(costs,numRows * numCols,MPI_INT,root,tag,MPI_COMM_WORLD);
    }
    else
    {
        MPI_Status status;
        int *nodeCost = calloc(numRows * numCols,sizeof(int));
        for(size_t node = 1;node<nproc;node++)
        {
            //TODO NEED TO COMBINE COST ARRAYS HERE
            MPI_Recv(nodeCost,numRows * numCols,MPI_INT,node,tag,MPI_COMM_WORLD,&status);
            mergeCost(costs,nodeCost,numRows,numCols);  //Merges the cost arrays of the node
            //and the root
            //Tells node that the root received the data
        }
    }
    MPI_Bcast(costs,numRows * numCols,MPI_INT,root,MPI_COMM_WORLD);
    //The line above broadcasts the udpated cost array to each node
}


//TODO NEED TO WRITE FUNCTION THAT COMBINES COST ARRAYS FOR EACH NODE AT THE ROOT NODE
//NODES WILL BROADCAST THEIR COST ARRAYS AFTER EACH ITERATION TO BEGIN WITH

// Perform computation, including reading/writing output files
void compute(int procID, int nproc, char* inputFilename, double prob, int numIterations)
{
    int numRows;
    int numCols;
    int delta;
    int wiresPerProc;
    int numWires;
    int curCostWireH;
    int curCostWireV;

    MPI_Datatype wireMPI;
    MPI_Type_contiguous(8,MPI_INT,&wireMPI);
    MPI_Type_commit(&wireMPI);

    wire_t *wireArray = NULL;
    if(procID==root)
    {
        wireArray = readInput(inputFilename,procID,nproc,&numRows,&numCols,&delta,&wiresPerProc,&numWires);
    }
    wire_t *wireSubArray = calloc(wiresPerProc,sizeof(wire_t));
    MPI_Scatter(wireArray,wiresPerProc,wireMPI,wireSubArray,wiresPerProc,wireMPI,root,MPI_COMM_WORLD);  //Root sends data to other nodes using scatter
    init(numRows,numCols);
    //TODO Implement code here
    //TODO Decide which processors should be reading/writing files

    for(int iter = 0;iter<numIterations;iter++)
    {
        double probVal = ((double) rand() / (RAND_MAX));
        if(probVal > 1-prob)
        {
            continue;
        }
        //Search time
        for(int i = 0;i<wiresPerProc;i++)
        {
            int x1 = wireSubArray[i].x1;
            int y1 = wireSubArray[i].y1;
            int x2 = wireSubArray[i].x2;
            int y2 = wireSubArray[i].y2;
            int b1x = wireSubArray[i].b1x;
            int b1y = wireSubArray[i].b1y;
            int b2x = wireSubArray[i].b2x;
            int b2y = wireSubArray[i].b2y;
            int curMin = 0;
            if(x1==x2 || y1==y2)  //Points in a straight line
            {
                if(x1==x2)
                {
                    removeEqualXWire(x1,y2,y2,costs,numRows,numCols,undef);
                    wire_t EqX = {x1,y1,x2,y2,undef,undef,undef,undef};
                    curCostWireH = costEqualX(x1,y1,y2,costs,numRows,numCols,delta,&EqX);
                    placeWire(&EqX,costs,numRows,numCols);
                }
                else
                {
                    assert(y1==y2);
                    removeEqualYWire(x1,x2,y1,costs,numRows,numCols,undef);
                    wire_t EqY = {x1,y1,x2,y2,undef,undef,undef,undef};
                    curCostWireV = costEqualY(x1,x2,y1,costs,numRows,numCols,delta,&EqY);
                    placeWire(&EqY,costs,numRows,numCols);
                }
            }
            else  //TODO DEAL WITH BENDS NOW
            {
                if(x1==b1x)
                {
                    removeEqualXWire(x1,y1,b1y,costs,numRows,numCols,undef);
                }
                else if(y1==b1y)
                {
                    removeEqualYWire(x1,b1x,y1,costs,numRows,numCols,undef);
                }
                if(b2x!=undef && b2y!=undef)
                {
                    if(b1x==b2x)
                    {
                        removeEqualXWire(b1x,b1y,b2y,costs,numRows,numCols,b1y);
                    }
                    else if(b1y==b2y)
                    {
                        removeEqualYWire(b1x,b2x,b1y,costs,numRows,numCols,b1x);
                    }
                    if(b2x==x2)
                    {
                        removeEqualXWire(b2x,b2y,y2,costs,numRows,numCols,b2y);
                    }
                    else if(b2y==y2)
                    {
                        removeEqualYWire(b2x,x2,y2,costs,numRows,numCols,b2x);
                    }
                }
                else if(b1x==x2)
                {
                    removeEqualXWire(x2,b1y,y2,costs,numRows,numCols,b1y);
                }
                else if(b1y==y2)
                {
                    removeEqualYWire(b1x,x2,y2,costs,numRows,numCols,b1x);
                }
                wire_t H = {x1,y1,x2,y2,undef,undef,undef,undef};
                wire_t V = {x1,y1,x2,y2,undef,undef,undef,undef};

                curCostWireH = checkAllHorizontal(x1,x2,y1,y2,costs,numRows,numCols,delta,&H);
                curCostWireV = checkAllVertical(x1,x2,y1,y2,costs,numRows,numCols,delta,&V);
                curMin = min(curCostWireH,curCostWireV);
                if(curMin==curCostWireH)
                {
                    wireSubArray[i] = H;
                    placeWire(&H,costs,numRows,numCols);
                }
                else
                {
                    assert(curMin==curCostWireV);
                    wireSubArray[i] = V;
                    placeWire(&V,costs,numRows,numCols);
                }
            }
        }
        combineCostArrays(procID,costs,numRows,numCols,nproc);
        //Communicating the info to the root so the root can merge arrays
    }
    //TODO NEED TO MERGE WIRE ARRAYS SO THAT THE ROOT CAN WRITE TO OUTPUT
    //COMBINEWIREARRAY
    MPI_Gather(wireSubArray,wiresPerProc,wireMPI,wireArray,wiresPerProc,wireMPI,root,MPI_COMM_WORLD);

    if(procID==root)  //Root has the cost and wire array so it writes out
    {
        writeCost(inputFilename, nproc,numRows,numCols);
        writeOutput(inputFilename, nproc,numRows,numCols,delta,wiresPerProc,wireArray);
    }
}

// Read input file
wire_t *readInput(char* inputFilename,int procID,int nProcs,int *numRows,int *numCols,int *delta,int *wiresPerProc,int *numWires)
{
    FILE* fp = fopen(inputFilename, "r");
    //int numRows;
    //int numCols;
    //int delta;
    int x1;
    int y1;
    int x2;
    int y2;
    if (fp == NULL) {
        fprintf(stderr, "Failed to read input file %s\n", inputFilename);
        exit(-1);
    }
    if (fscanf(fp, "%d %d", numCols, numRows) != 2) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if ((*numRows <= 0) || (*numCols <= 0)) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if (fscanf(fp, "%d", delta) != 1) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if (fscanf(fp, "%d", numWires) != 1) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if (numWires <= 0) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    //init(*numRows, *numCols, *delta, numWires,procID);
    
    //int *globalWires = calloc(numWires,sizeof(wire_t));
    //Global array that contains all wires, will contain mergered wires at root

    //The above code creates a wire struct type for MPI

    *wiresPerProc = (*numWires + nProcs - 1)/nProcs;
    wire_t *wireArray = calloc(*numWires,sizeof(wire_t));
    for(int i = 0;i<(*numWires);i++)
    {
        if(fscanf(fp,"%d %d %d %d",&x1,&y1,&x2,&y2)!=4)
        {
            fprintf(stderr,"Invalid input file format\n");
            exit(-1);
        }
        initWire(i,x1,y1,x2,y2,wireArray);
    }
    //wire_t *wireSubArray = calloc(*wiresPerProc,sizeof(wire_t));
    //MPI_Scatter(wireArray,*wiresPerProc,wireMPI,wireSubArray,*wiresPerProc,wireMPI,root,MPI_COMM_WORLD);
    fclose(fp);
    return wireArray;
}

// Write cost array file based on input filename
void writeCost(char* inputFilename, int nproc,int numRows,int numCols)
{
    //Only need root to do this
    char* dname = strdup(inputFilename);
    char* bname = strdup(inputFilename);
    char* costsFilename = malloc(strlen(inputFilename) + 100);
    FILE* fp;
    //int numRows = getNumRows();
    //int numCols = getNumCols();
    int row;
    int col;
    assert(costsFilename != NULL);
    sprintf(costsFilename, "%s/costs_%s_%d.txt", dirname(dname), basename(bname), nproc);
    fp = fopen(costsFilename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Failed to write costs file %s\n", costsFilename);
        exit(-1);
    }
    fprintf(fp, "%d %d\n", numCols, numRows);
    for (row = 0; row < numRows; row++) {
        for (col = 0; col < numCols; col++) {
            fprintf(fp, "%d ", getCost(row, col,numCols));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    free(costsFilename);
    free(bname);
    free(dname);
}

// Write wire output file based on input filename
void writeOutput(char* inputFilename, int nproc,int numRows,int numCols,int delta,int numWires,
        wire_t *wireArray)
{
    char* dname = strdup(inputFilename);
    char* bname = strdup(inputFilename);
    char* outputFilename = malloc(strlen(inputFilename) + 100);
    FILE* fp;
    //int numRows = getNumRows();
    //int numCols = getNumCols();
    //int delta = getDelta();
    //int numWires = getNumWires();
    int wireIndex;
    int numPoints;
    int x1;
    int y1;
    int x2;
    int y2;
    int x3;
    int y3;
    int x4;
    int y4;
    assert(outputFilename != NULL);
    sprintf(outputFilename, "%s/output_%s_%d.txt", dirname(dname), basename(bname), nproc);
    fp = fopen(outputFilename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Failed to write output file %s\n", outputFilename);
        exit(-1);
    }
    fprintf(fp, "%d %d\n", numCols, numRows);
    fprintf(fp, "%d\n", delta);
    fprintf(fp, "%d\n", numWires);
    for (wireIndex = 0; wireIndex < numWires; wireIndex++) {
        numPoints = getWire(wireIndex, &x1, &y1, &x2, &y2, &x3, &y3, &x4, &y4,wireArray);
        switch (numPoints) {
            case 2:
                fprintf(fp, "%d %d %d %d\n", x1, y1, x2, y2);
                break;

            case 3:
                fprintf(fp, "%d %d %d %d %d %d\n", x1, y1, x2, y2, x3, y3);
                break;

            case 4:
                fprintf(fp, "%d %d %d %d %d %d %d %d\n", x1, y1, x2, y2, x3, y3, x4, y4);
                break;

            default:
                assert(0); // invalid number of points
                break;
        }
    }
    fclose(fp);
    free(outputFilename);
    free(bname);
    free(dname);
}
