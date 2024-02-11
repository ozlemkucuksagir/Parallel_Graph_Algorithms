//Ozlem KUCUKSAGIR
//19050111021

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define N 701


// Structure representing the graph
typedef struct {
    int vertices[N][N];  // Adjacency matrix to represent edges between vertices
    int visited[N];      // Array to keep track of visited vertices during traversal
} Graph;

// Function to initialize the graph with random edges
void initGraph(Graph *graph) {
    srand(time(NULL));

    for (int i = 0; i < N; ++i) {
        graph->visited[i] = 0;  // Mark all vertices as not visited
        for (int j = i + 1; j < N; ++j) {
            // Generate random edges with increased probability (1/3 chance)
            if (rand() % 3 == 0) {
                graph->vertices[i][j] = 1;  // Set edge from i to j
                graph->vertices[j][i] = 1;  // Set edge from j to i (undirected graph)
            }
        }
    }
}

typedef struct {
    double processingTime;
} PerformanceMetrics;


PerformanceMetrics measureTimeWithMetrics(void (*func)(Graph*, int), Graph *graph, int startNode) {
    double start = omp_get_wtime();
    func(graph, startNode);          // Execute the specified traversal algorithm
    double end = omp_get_wtime();

    PerformanceMetrics metrics;
    metrics.processingTime = end - start;  // Calculate processing time
    return metrics;
}
// Function to calculate speedup
double calculateSpeedup(double sequentialTime, double parallelTime) {
    if (parallelTime == 0) {
        // Avoid division by zero
        return 0.0;
    }
    return sequentialTime / parallelTime;
}

// Function to calculate efficiency
double calculateEfficiency(double speedup, int numThreads) {
    if (numThreads == 0) {
        // Avoid division by zero
        return 0.0;
    }
    return speedup / numThreads;
}


// Sequential Depth-First Search (DFS) traversal
void dfsSequential(Graph *graph, int startNode) {
    int stack[N];
    int top = -1;

    stack[++top] = startNode;  // Push the starting node onto the stack

    while (top >= 0) {
        int node = stack[top--];  // Pop a node from the stack

        if (!graph->visited[node]) {
            graph->visited[node] = 1;  // Mark the node as visited

            // Explore neighbors and push unvisited neighbors onto the stack
            for (int i = N - 1; i >= 0; --i) {
                if (graph->vertices[node][i] && !graph->visited[i]) {
                    stack[++top] = i;  // Push unvisited neighbor onto the stack
                }
            }
        }
    }
}

// Parallel Depth-First Search (DFS) traversal using OpenMP
void dfsParallel(Graph *graph, int startNode) {
#pragma omp parallel
    {
        int stack[N];
        int top = -1;

#pragma omp single
        stack[++top] = startNode;  // Push the starting node onto the stack

        while (top >= 0) {
            int node;

#pragma omp critical
            {
                node = stack[top--];  // Pop a node from the stack
            }

            if (!graph->visited[node]) {
#pragma omp critical
                {
                    graph->visited[node] = 1;  // Mark the node as visited
                }

#pragma omp for nowait
                // Explore neighbors and push unvisited neighbors onto the stack
                for (int i = N - 1; i >= 0; --i) {
                    if (graph->vertices[node][i] && !graph->visited[i]) {
#pragma omp critical
                        {
                            stack[++top] = i;  // Push unvisited neighbor onto the stack
                        }
                    }
                }
            }
        }
    }
}

// Sequential Breadth-First Search (BFS) traversal
void bfsSequential(Graph *graph, int startNode) {
    int queue[N];
    int front = -1, rear = -1;

    queue[++rear] = startNode;  // Enqueue the starting node

    while (front != rear) {
        int node = queue[++front];  // Dequeue a node

        if (!graph->visited[node]) {
            graph->visited[node] = 1;  // Mark the node as visited

            // Enqueue unvisited neighbors
            for (int i = 0; i < N; ++i) {
                if (graph->vertices[node][i] && !graph->visited[i]) {
                    queue[++rear] = i;  // Enqueue unvisited neighbor
                }
            }
        }
    }
}

// Parallel Breadth-First Search (BFS) traversal using OpenMP
void bfsParallel(Graph *graph, int startNode) {
#pragma omp parallel
    {
        int queue[N];
        int front = -1, rear = -1;

#pragma omp single
        queue[++rear] = startNode;  // Enqueue the starting node

        while (front != rear) {
            int node;

#pragma omp critical
            {
                node = queue[++front];  // Dequeue a node
            }

            if (!graph->visited[node]) {
#pragma omp critical
                {
                    graph->visited[node] = 1;  // Mark the node as visited
                }

#pragma omp for nowait
                // Enqueue unvisited neighbors
                for (int i = 0; i < N; ++i) {
                    if (graph->vertices[node][i] && !graph->visited[i]) {
#pragma omp critical
                        {
                            queue[++rear] = i;  // Enqueue unvisited neighbor
                        }
                    }
                }
            }
        }
    }
}


// Main function
int main() {
    Graph graph;
    initGraph(&graph);  // Initialize the graph with random edges

    int startNode = 0;   // Starting node for traversal
    int numThreads;

#pragma omp parallel
    {
#pragma omp master
        {
            numThreads = omp_get_num_threads();  // Get the number of threads
        }
    }

    printf("Number of Threads: %d\n", numThreads);


    // Measure and print execution times for DFS and BFS along with speedup and efficiency metrics
    PerformanceMetrics sequentialDFSResult = measureTimeWithMetrics(dfsSequential, &graph, startNode);
    PerformanceMetrics parallelDFSResult = measureTimeWithMetrics(dfsParallel, &graph, startNode);

    // Reset visited array for BFS
    for (int i = 0; i < N; ++i) {
        graph.visited[i] = 0;
    }

    PerformanceMetrics sequentialBFSResult = measureTimeWithMetrics(bfsSequential, &graph, startNode);
    PerformanceMetrics parallelBFSResult = measureTimeWithMetrics(bfsParallel, &graph, startNode);

    // Print results
    printf("Sequential DFS Time: %f s\n", sequentialDFSResult.processingTime);
    printf("Parallel DFS Time: %f s\n", parallelDFSResult.processingTime);
    printf("Sequential BFS Time: %f s\n", sequentialBFSResult.processingTime);
    printf("Parallel BFS Time: %f s\n", parallelBFSResult.processingTime);

    // Print speedup and efficiency metrics
    printf("DFS Speedup: %f\n", calculateSpeedup(sequentialDFSResult.processingTime,parallelDFSResult.processingTime));

    printf("BFS Speedup: %f\n", calculateSpeedup(sequentialBFSResult.processingTime,parallelBFSResult.processingTime));

    printf("DFS Efficiency: %f\n", calculateEfficiency(calculateSpeedup(sequentialDFSResult.processingTime,parallelDFSResult.processingTime),numThreads));

    printf("BFS Efficiency: %f\n", calculateEfficiency(calculateSpeedup(sequentialBFSResult.processingTime,parallelBFSResult.processingTime),numThreads));

    return 0;
}
