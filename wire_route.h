#ifndef _WIRE_ROUTE_H
#define _WIRE_ROUTE_H
#define undef -1

typedef struct wire_t{
    int x1;
    int y1;
    int x2;
    int y2;
    int b1x;
    int b1y;
    int b2x;
    int b2y;
}wire_t;

// Perform computation, including reading/writing output files
void compute(int procID, int nproc, char* inputFilename, double prob, int numIterations);

// Read input file
wire_t *readInput(char* inputFilename,int procID,int nProcs,int *numRows,int *numCols,int *delta,int *wiresPerProc);

// Write cost array file based on input filename
void writeCost(char* inputFilename, int nproc,int numRows,int numCols,int procID);

// Write wire output file based on input filename
void writeOutput(char* inputFilename, int nproc,int numRows,int numCols,int delta,int numWires,
        wire_t *wireArray);

#endif // _WIRE_ROUTE_H
