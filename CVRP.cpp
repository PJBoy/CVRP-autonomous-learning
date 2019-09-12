// The capacitated vehicle routing problem //
/*
Truck begins/ends at depot.
Trucks all have same fixed capacity.
Customers have demand.
Minimise distance travelled by all trucks.
Any number of trucks.
Output cost and route to "best-solution.txt" (look at assignment description)

Thus:
    truck_i = v_0, U_j v_j, v_0
    for all i != j:
        (truck_i \ v_0) ^ (truck_j \ v_0) = {}
    minimise:
        sum_i distance(truck_i)
    subject to:
        sum_{i != 0} demand(v_i) <= truckCapacity
*/

/*
Things to consider:
Crossover rate                   | Medium
Crossover method                 | Partially matched crossover
Mutation rate                    | Low
Mutation method                  | Two swap
Fitness function                 | Sum of edges
Population growth                | 0
Percentage population to replace | 100%
Representation of chromosome     | Permutation vertex list with repeating depots
Method of parent selection       | Roulette wheel selection (for now)
Penalty                          | Weighted penalty for a chromosome violating the subject constraint
Repair                           | Repair routes whose trucks exceed capacity if sufficiently many trucks are violating
*/


#include <array>
#include <chrono>
#include <fstream>
#include <functional>
#include <forward_list>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>


const unsigned
    input_n_vertices       = 250,
    input_truckCapacity    = 500,
    n_vertices             = input_n_vertices,
    n_genes                = n_vertices - 1,
    truckCapacity          = input_truckCapacity,
    n_chromosomes          = 1 << 6,
    minParents             = n_chromosomes >> 3,
    windowSize             = 1 << 11;


typedef unsigned                                                         i_Vertex_t, i_evaluation_t, i_clusters_t, i_cluster_t;
typedef float                                                            distance_measure;
typedef i_Vertex_t                                                       allele_t;
typedef std::vector<allele_t>                                            cluster_t, chromosome_vertex_t;
typedef std::vector<cluster_t>                                           clusters_t;
typedef std::vector<i_cluster_t>                                         cluster_order_t;
typedef std::vector<cluster_order_t>                                     chromosome_t;
typedef std::array<std::array<distance_measure, n_vertices>, n_vertices> distances_t;
typedef std::array<chromosome_t, n_chromosomes>                          population_t;
typedef std::array<chromosome_vertex_t, n_chromosomes>                   population_vertex_t;


struct Evaluations;

struct Evaluation
{
    distance_measure value;
    Evaluations* parent;
    i_evaluation_t index;

    Evaluation();
    Evaluation(Evaluations* parent, i_evaluation_t index);
    Evaluation& operator=(distance_measure value_in);

#define operator(op)                               \
    bool operator op (distance_measure rhs) const;
    operator(<)
    operator(<=)
    operator(>)
    operator(>=)
    operator(==)
    operator(!=)
#undef operator
};

struct Evaluations
{
    // static double globalMinVariance, globalMaxVariance;

    std::array<Evaluation, n_chromosomes> evaluations;

    i_evaluation_t min_i, max_i;
    distance_measure min, max;
    // distance_measure sum, sumSquare;

    // const unsigned n{n_chromosomes};

    Evaluations()
    {
        reset();
        for (i_evaluation_t i(0); i < evaluations.size(); ++i)
            evaluations[i] = {this, i};
    }

    const Evaluation& operator[](i_evaluation_t i) const
    {
        return evaluations[i];
    }

    Evaluation& operator[](i_evaluation_t i)
    {
        return evaluations[i];
    }

    void reset()
    {
        min = std::numeric_limits<distance_measure>::max();
        max = 0;
        // sum = 0;
        // sumSquare = 0;
    }

    void update(i_evaluation_t index, distance_measure value)
    {
        if (value < min)
        {
            min = value;
            min_i = index;
        }
        if (value > max)
        {
            max = value;
            max_i = index;
        }

        // sum += value;
        // sumSquare += value * value;
    }

#if 0
    double mutationRate() const
    {
        // variance = sum^n_{i=1} v_i^2 / (n - 1) - n / (n - 1) (sum^n_{i=1} v_i / n)^2
        const double variance((sumSquare - sum*sum / n) / (n - 1));
        globalMinVariance = std::min(globalMinVariance, variance);
        globalMaxVariance = std::max(globalMaxVariance, variance);
        const double variance_scaled((variance - globalMinVariance) / (globalMaxVariance - globalMinVariance));
        const double variableMutationRange(maximumMutationRate - minimumMutationRate);
        
        return minimumMutationRate + variableMutationRange * (1.0 - std::min(1.0, variance_scaled * varianceScalar));
    }
#endif
};

// double Evaluations::globalMinVariance = std::numeric_limits<double>::max();
// double Evaluations::globalMaxVariance = 0;

// struct Evaluation
#if 1
Evaluation::Evaluation() = default;

Evaluation::Evaluation(Evaluations* parent, i_evaluation_t index)
    : parent(parent), index(index)
{}

Evaluation& Evaluation::operator=(distance_measure value_in)
{
    value = value_in;
    parent->update(index, value);
    return *this;
}

#define operator(op)                                      \
bool Evaluation::operator op (distance_measure rhs) const\
{                                                         \
    return value op rhs;                                  \
}
operator(<)
operator(<=)
operator(>)
operator(>=)
operator(==)
operator(!=)
#undef operator
#endif


struct Vertex
{
    i_Vertex_t i;
    int X, Y;
    unsigned demand;
};

typedef std::array<Vertex, input_n_vertices> input_vertices_t;

const input_vertices_t input_vertices
{{
    {1,   -33, 33,  0},
    {2,   -99, -97, 6},
    {3,   -59, 50,  72},
    {4,   0,   14,  93},
    {5,   -17, -66, 28},
    {6,   -69, -19, 5},
    {7,   31,  12,  43},
    {8,   5,   -41, 1},
    {9,   -12, 10,  36},
    {10,  -64, 70,  53},
    {11,  -12, 85,  63},
    {12,  -18, 64,  25},
    {13,  -77, -16, 50},
    {14,  -53, 88,  57},
    {15,  83,  -24, 1},
    {16,  24,  41,  66},
    {17,  17,  21,  37},
    {18,  42,  96,  51},
    {19,  -65, 0,   47},
    {20,  -47, -26, 88},
    {21,  85,  36,  75},
    {22,  -35, -54, 48},
    {23,  54,  -21, 40},
    {24,  64,  -17, 8},
    {25,  55,  89,  69},
    {26,  17,  -25, 93},
    {27,  -61, 66,  29},
    {28,  -61, 26,  5},
    {29,  17,  -72, 53},
    {30,  79,  38,  8},
    {31,  -62, -2,  24},
    {32,  -90, -68, 53},
    {33,  52,  66,  13},
    {34,  -54, -50, 47},
    {35,  8,   -84, 57},
    {36,  37,  -90, 9},
    {37,  -83, 49,  74},
    {38,  35,  -1,  83},
    {39,  7,   59,  96},
    {40,  12,  48,  42},
    {41,  57,  95,  80},
    {42,  92,  28,  22},
    {43,  -3,  97,  56},
    {44,  -7,  52,  43},
    {45,  42,  -15, 12},
    {46,  77,  -43, 73},
    {47,  59,  -49, 32},
    {48,  25,  91,  8},
    {49,  69,  -19, 79},
    {50,  -82, -14, 79},
    {51,  74,  -70, 4},
    {52,  69,  59,  14},
    {53,  29,  33,  17},
    {54,  -97, 9,   19},
    {55,  -58, 9,   44},
    {56,  28,  93,  5},
    {57,  7,   73,  37},
    {58,  -28, 73,  100},
    {59,  -76, 55,  62},
    {60,  41,  42,  90},
    {61,  92,  40,  57},
    {62,  -84, -29, 44},
    {63,  -12, 42,  37},
    {64,  51,  -45, 80},
    {65,  -37, 46,  60},
    {66,  -97, 35,  95},
    {67,  14,  89,  56},
    {68,  60,  58,  56},
    {69,  -63, -75, 9},
    {70,  -18, 34,  39},
    {71,  -46, -82, 15},
    {72,  -86, -79, 4},
    {73,  -43, -30, 58},
    {74,  -44, 7,   73},
    {75,  -3,  -20, 5},
    {76,  36,  41,  12},
    {77,  -30, -94, 3},
    {78,  79,  -62, 8},
    {79,  51,  70,  31},
    {80,  -61, -26, 48},
    {81,  6,   94,  3},
    {82,  -19, -62, 52},
    {83,  -20, 51,  99},
    {84,  -81, 37,  29},
    {85,  7,   31,  12},
    {86,  52,  12,  50},
    {87,  83,  -91, 98},
    {88,  -7,  -92, 4},
    {89,  82,  -74, 56},
    {90,  -70, 85,  24},
    {91,  -83, -30, 33},
    {92,  71,  -61, 45},
    {93,  85,  11,  98},
    {94,  66,  -48, 4},
    {95,  78,  -87, 36},
    {96,  9,   -79, 72},
    {97,  -36, 4,   26},
    {98,  66,  39,  71},
    {99,  92,  -17, 84},
    {100, -46, -79, 21},
    {101, -30, -63, 99},
    {102, -42, 63,  33},
    {103, 20,  42,  84},
    {104, 15,  98,  74},
    {105, 1,   -17, 93},
    {106, 64,  20,  25},
    {107, -96, 85,  39},
    {108, 93,  -29, 42},
    {109, -40, -84, 77},
    {110, 86,  35,  68},
    {111, 91,  36,  50},
    {112, 62,  -8,  42},
    {113, -24, 4,   71},
    {114, 11,  96,  85},
    {115, -53, 62,  78},
    {116, -28, -71, 64},
    {117, 7,   -4,  5},
    {118, 95,  -9,  93},
    {119, -3,  17,  18},
    {120, 53,  -90, 38},
    {121, 58,  -19, 29},
    {122, -83, 84,  81},
    {123, -1,  49,  4},
    {124, -4,  17,  23},
    {125, -82, -3,  11},
    {126, -43, 47,  86},
    {127, 6,   -6,  2},
    {128, 70,  99,  31},
    {129, 68,  -29, 54},
    {130, -94, -30, 87},
    {131, -94, -20, 17},
    {132, -21, 77,  81},
    {133, 64,  37,  72},
    {134, -70, -19, 10},
    {135, 88,  65,  50},
    {136, 2,   29,  25},
    {137, 33,  57,  71},
    {138, -70, 6,   85},
    {139, -38, -56, 51},
    {140, -80, -95, 29},
    {141, -5,  -39, 55},
    {142, 8,   -22, 45},
    {143, -61, -76, 100},
    {144, 76,  -22, 38},
    {145, 49,  -71, 11},
    {146, -30, -68, 82},
    {147, 1,   34,  50},
    {148, 77,  79,  39},
    {149, -58, 64,  6},
    {150, 82,  -97, 87},
    {151, -80, 55,  83},
    {152, 81,  -86, 22},
    {153, 39,  -49, 24},
    {154, -67, 72,  69},
    {155, -25, -89, 97},
    {156, -44, -95, 65},
    {157, 32,  -68, 97},
    {158, -17, 49,  79},
    {159, 93,  49,  79},
    {160, 99,  81,  46},
    {161, 10,  -49, 52},
    {162, 63,  -41, 39},
    {163, 38,  39,  94},
    {164, -28, 39,  97},
    {165, -2,  -47, 18},
    {166, 38,  8,   3},
    {167, -42, -6,  23},
    {168, -67, 88,  19},
    {169, 19,  93,  40},
    {170, 40,  27,  49},
    {171, -61, 56,  96},
    {172, 43,  33,  58},
    {173, -18, -39, 15},
    {174, -69, 19,  21},
    {175, 75,  -18, 56},
    {176, 31,  85,  67},
    {177, 25,  58,  10},
    {178, -16, 36,  36},
    {179, 91,  15,  84},
    {180, 60,  -39, 59},
    {181, 49,  -47, 85},
    {182, 42,  33,  60},
    {183, 16,  -81, 33},
    {184, -78, 53,  62},
    {185, 53,  -80, 70},
    {186, -46, -26, 79},
    {187, -25, -54, 98},
    {188, 69,  -46, 99},
    {189, 0,   -78, 18},
    {190, -84, 74,  55},
    {191, -16, 16,  75},
    {192, -63, -14, 94},
    {193, 51,  -77, 89},
    {194, -39, 61,  13},
    {195, 5,   97,  19},
    {196, -55, 39,  19},
    {197, 70,  -14, 90},
    {198, 0,   95,  35},
    {199, -45, 7,   76},
    {200, 38,  -24, 3},
    {201, 50,  -37, 11},
    {202, 59,  71,  98},
    {203, -73, -96, 92},
    {204, -29, 72,  1},
    {205, -47, 12,  2},
    {206, -88, -61, 63},
    {207, -88, 36,  57},
    {208, -46, -3,  50},
    {209, 26,  -37, 19},
    {210, -39, -67, 24},
    {211, 92,  27,  14},
    {212, -80, -31, 18},
    {213, 93,  -50, 77},
    {214, -20, -5,  28},
    {215, -22, 73,  72},
    {216, -4,  -7,  49},
    {217, 54,  -48, 58},
    {218, -70, 39,  84},
    {219, 54,  -82, 58},
    {220, 29,  41,  41},
    {221, -87, 51,  98},
    {222, -96, -36, 77},
    {223, 49,  8,   57},
    {224, -5,  54,  39},
    {225, -26, 43,  99},
    {226, -11, 60,  83},
    {227, 40,  61,  54},
    {228, 82,  35,  86},
    {229, -92, 12,  2},
    {230, -93, -86, 14},
    {231, -66, 63,  42},
    {232, -72, -87, 14},
    {233, -57, -84, 55},
    {234, 23,  52,  2},
    {235, -56, -62, 18},
    {236, -19, 59,  17},
    {237, 63,  -14, 22},
    {238, -13, 38,  28},
    {239, -19, 87,  3},
    {240, 44,  -84, 96},
    {241, 98,  -17, 53},
    {242, -16, 62,  15},
    {243, 3,   66,  36},
    {244, 26,  22,  98},
    {245, -38, -81, 78},
    {246, 70,  80,  92},
    {247, 17,  -35, 65},
    {248, 96,  -83, 64},
    {249, -77, 80,  43},
    {250, -14, 44,  50}
}};

const i_Vertex_t depot(input_vertices[0].i);
const int
    input_vertices_min(-100),
    input_vertices_max(100),
    input_vertices_range(input_vertices_max - input_vertices_min);


class SlidingWindow
{
    using T = distance_measure;

    std::array<T, windowSize> window;
    unsigned size{0};
    unsigned next{0};

public:
    void push(T value)
    {
        window[next] = value;
        next = (next + 1) % windowSize;
        size = std::min(size + 1, windowSize);
    }

    T mean() const
    {
        return
            size < windowSize
            ? std::numeric_limits<T>::max()
            : std::accumulate(std::begin(window), std::end(window), T()) / windowSize;
    }
};


template<typename T> 
std::pair<const T&, const T&> minmax(const T& a, const T& b)
{
    return
        b < a
        ? std::make_pair(b, a)
        : std::make_pair(a, b);
}


template<typename BidirIt>
void reverse(BidirIt first, BidirIt last)
{
    while (first != last && first != --last)
        std::iter_swap(first++, last);
}


template<typename RandomIt, typename URNG>
void shuffle(RandomIt first, RandomIt last, URNG&& g)
{
    typedef typename std::iterator_traits<RandomIt>::difference_type diff_t;
    typedef typename std::make_unsigned<diff_t>::type udiff_t;
    typedef typename std::uniform_int_distribution<udiff_t> distr_t;
    typedef typename distr_t::param_type param_t;
 
    distr_t D;
    for (diff_t i(last - first - 1); i; --i)
        std::swap(first[i], first[D(g, param_t(0, i))]);
}


template<typename InputIt, typename T>
InputIt find(InputIt first, InputIt last, const T& value)
{
    for (; first != last; ++first)
        if (*first == value)
            return first;
    return last;
}


template<typename InputIt, typename T>
InputIt find_circular(InputIt begin, InputIt end, InputIt first, const T& value)
{
    InputIt i(find(first, end, value));
    return i != end ? i : find(begin, first, value);
}


template<class OutputIt, class Size, class T>
OutputIt iota_n(OutputIt first, Size count, T value)
{
    for (Size i(count); i > 0; --i)
    {
        *first++ = value;
        ++value;
    }
    return first;
}


distances_t distances;

void initialiseDistances()
{
    // Calculate pairwise euclidean distances
    for (unsigned y(0); y < distances.size(); ++y)
        for (unsigned x(0); x < distances.size(); ++x)
        {
            const distance_measure
                d_x(distance_measure(input_vertices[y].X - input_vertices[x].X)),
                d_y(distance_measure(input_vertices[y].Y - input_vertices[x].Y));
            distances[y][x] = std::sqrt(d_x * d_x + d_y * d_y);
        }
}


void population_clusterToVertex(population_vertex_t& population_vertex, const population_t& population, const clusters_t& clusters)
{
    for (unsigned i(0); i < population.size(); ++i)
    {
        const chromosome_t& chromosome(population[i]);
        chromosome_vertex_t& chromosome_vertex(population_vertex[i]);

        chromosome_vertex = {};

        for (unsigned gene(0); gene < chromosome.size(); ++gene)
            for (i_cluster_t i : chromosome[gene])
                chromosome_vertex.push_back(clusters[gene][i]);
    }
}


void population_vertexToCluster(const population_vertex_t& population_vertex, population_t& population, clusters_t& clusters)
{
    // Use first chromosome to set up the clusters
    const chromosome_vertex_t& chromosome_vertex(population_vertex.front());

    // For genes until demand > capacity, make that a cluster,
    clusters = {};
    {
        cluster_t cluster;
        unsigned clusterDemand(0);
        allele_t from(depot);
        for (allele_t to : chromosome_vertex)
        {
            const unsigned demand(input_vertices[to - 1].demand);

            if (clusterDemand + demand > truckCapacity)
            {
                clusters.push_back(cluster);
                cluster = {};
                clusterDemand = 0;
            }

            cluster.push_back(to);
            clusterDemand += demand;
        }
        clusters.push_back(cluster);
    }

    chromosome_t chromosome;

    for (const cluster_t& cluster : clusters)
    {
        cluster_order_t clusterOrder;
        iota_n(std::back_inserter(clusterOrder), cluster.size(), 0);
        chromosome.push_back(clusterOrder);
    }

    population.fill(chromosome);
}


namespace NS_Cluster
{
    const double
        // minimumMutationRate = 0.950,
        // maximumMutationRate = 0.999,
        // varianceScalar      = 4,
        mutationRate_shuffle = 0.5,
        mutationRate_reverse = 0.99,
        mutationRate_swap = 0.99,
        mutationRate_exchangeCustomers = 0.0001;


    typedef unsigned                                                         i_Vertex_t, i_evaluation_t, i_clusters_t, i_cluster_t;
    typedef float                                                            distance_measure;
    typedef i_Vertex_t                                                       allele_t;
    typedef std::vector<allele_t>                                            cluster_t;
    typedef std::vector<cluster_t>                                           clusters_t;
    typedef std::vector<i_cluster_t>                                         cluster_order_t;
    typedef std::vector<cluster_order_t>                                     chromosome_t;
    typedef std::vector<chromosome_t>                                        parents_t;
    typedef std::array<std::array<distance_measure, n_vertices>, n_vertices> distances_t;
    typedef std::array<chromosome_t, n_chromosomes>                          population_t;
    typedef population_t                                                     offspring_t;


    const distance_measure
        margin_leastFit = -1 + 0.0001f, // Minimum -1. The greater, the more likely the least fit parents will be selected, the lesser, the more likely that fit parents will be selected
        satisfictionTheshold = 0.001f,
        penalty = 1,
        penalty_n_invalidTrucks = 1.01f,
        distanceThreshold = 15;


    void initialPopulation(population_t& population, clusters_t& clusters)
    {
        std::pair<i_Vertex_t, unsigned> verticesSpace[input_vertices_range][input_vertices_range] = {};
        for (const Vertex& vertex : input_vertices)
            verticesSpace[vertex.Y - input_vertices_min][vertex.X - input_vertices_min] = std::make_pair(vertex.i, vertex.demand);

        cluster_t cluster;
        unsigned clusterDemand(0);
        for (unsigned y(0); y < 5; ++y)
            for (unsigned x(y & 1 ? 5 : 0); (y & 1) ? x-->0 : x < 5; (y & 1) ? x : ++x)
                for (unsigned r(1); r <= input_vertices_range / 5; ++r)
                {
                    const unsigned yy(y * (input_vertices_range / 5) + r - 1);

                    for (unsigned xx(x * input_vertices_range / 5); xx < x * input_vertices_range / 5 + r; ++xx)
                    {
                        const unsigned demand(verticesSpace[yy][xx].second);
                        if (demand)
                        {
                            if (clusterDemand + demand > truckCapacity)
                            {
                                clusters.push_back(cluster);
                                cluster = {};
                                clusterDemand = 0;
                            }

                            cluster.push_back(verticesSpace[yy][xx].first);
                            clusterDemand += demand;
                        }
                    }

                    const unsigned xx(x * (input_vertices_range / 5) + r - 1);

                    for (unsigned yy(y * input_vertices_range / 5); yy + 1 < y * input_vertices_range / 5 + r; ++yy)
                    {
                        const unsigned demand(verticesSpace[yy][xx].second);
                        if (demand)
                        {
                            if (clusterDemand + demand > truckCapacity)
                            {
                                clusters.push_back(cluster);
                                cluster = {};
                                clusterDemand = 0;
                            }

                            cluster.push_back(verticesSpace[yy][xx].first);
                            clusterDemand += demand;
                        }
                    }
                }

        clusters.push_back(cluster);

        //std::default_random_engine UPRNG(std::random_device{}());
        chromosome_t chromosome;

        for (const cluster_t& cluster : clusters)
        {
            cluster_order_t clusterOrder;
            iota_n(std::back_inserter(clusterOrder), cluster.size(), 0);
            //shuffle(std::begin(clusterOrder), std::end(clusterOrder), UPRNG);
            chromosome.push_back(clusterOrder);
        }

        population.fill(chromosome);
    }


    distance_measure evaluateCluster(const cluster_order_t& cluster_order, const cluster_t& cluster)
    {
        distance_measure distance(0);
        allele_t from(depot);

        for (i_cluster_t i_to : cluster_order)
        {
            const allele_t to(cluster[i_to]);
            distance += distances[from - 1][to - 1];
            from = to;
        }

        distance += distances[from - 1][depot - 1];

        return distance;
    }


    distance_measure evaluateChromosome(const chromosome_t& chromosome, const clusters_t& clusters)
    {
        distance_measure distance(0);

        for (i_clusters_t gene(0); gene < chromosome.size(); ++gene)
            distance += evaluateCluster(chromosome[gene], clusters[gene]);

        return distance;
    }


    Evaluations& evaluatePopulation(const population_t& population, const clusters_t& clusters)
    {
        static Evaluations evaluations;
        evaluations.reset();

        const unsigned n(population.size());
#pragma omp parallel for
        for (int i = 0; unsigned(i) < n; ++i)
            evaluations[i] = evaluateChromosome(population[i], clusters);

        return evaluations;
    }


    parents_t selectParents(const population_t& population, const Evaluations& evaluations)
    {
        static std::default_random_engine UPRNG(std::random_device{}());

        const distance_measure
            min(evaluations.min),
            max(evaluations.max);

        std::uniform_real_distribution<distance_measure> URD
            (
                min,
                max + (max - min) * margin_leastFit
                );
        std::function<distance_measure()> random_real(std::bind(URD, UPRNG));

        parents_t selectedParents;
        for (;;)
        {
            for (unsigned i(0); i < population.size(); ++i)
                if (evaluations[i] <= random_real())
                    selectedParents.push_back(population[i]);
            if (selectedParents.size() >= minParents)
                return selectedParents;
        }
    }


    offspring_t& reproduce(const parents_t parents, clusters_t& clusters, const distance_measure rateOfChange)
    {
        // want |offspring| = |population|, but |parents| <= |population|
        // so lets choose pairs of parents at random until there are enough offspring

        static offspring_t offspring;

        std::default_random_engine UPRNG(std::random_device{}());

        std::uniform_int_distribution<unsigned> UID_parent(0, parents.size() - 1);
        std::function<unsigned()> random_parent(std::bind(UID_parent, UPRNG));

        std::uniform_int_distribution<unsigned> UID_offspring(0, offspring.size() - 1);
        std::function<unsigned()> random_offspring(std::bind(UID_offspring, UPRNG));

        std::bernoulli_distribution BD_fair;
        std::function<bool()> coin_flip(std::bind(BD_fair, UPRNG));

        std::bernoulli_distribution BD_mutation_shuffle(mutationRate_shuffle);
        std::function<bool()> coin_mutation_shuffle(std::bind(BD_mutation_shuffle, UPRNG));

        std::bernoulli_distribution BD_mutation_swap(mutationRate_swap);
        std::function<bool()> coin_mutation_swap(std::bind(BD_mutation_swap, UPRNG));

        std::bernoulli_distribution BD_mutation_reverse(mutationRate_reverse);
        std::function<bool()> coin_mutation_reverse(std::bind(BD_mutation_reverse, UPRNG));

        std::bernoulli_distribution BD_mutation_exchangeCustomers(mutationRate_exchangeCustomers / rateOfChange);
        std::function<bool()> coin_mutation_exchangeCustomers(std::bind(BD_mutation_exchangeCustomers, UPRNG));


        // (Possibly) exchange some customers between two clusters
        if (coin_mutation_exchangeCustomers())
        {
            std::uniform_int_distribution<unsigned> UID_cluster(0, clusters.size() - 1);
            std::function<unsigned()> random_cluster(std::bind(UID_cluster, UPRNG));

            cluster_t& cluster1(clusters[random_cluster()]);

            allele_t cluster1allele(cluster1[0]);
            i_clusters_t i_cluster(random_cluster());
            while (distances[clusters[i_cluster][0] - 1][cluster1allele - 1] > distanceThreshold * rateOfChange)
                i_cluster = random_cluster();
            cluster_t& cluster2(clusters[i_cluster]);

            std::uniform_int_distribution<unsigned> UID_customer(0, std::min(cluster1.size(), cluster2.size()) - 1);
            std::function<unsigned()> random_customer(std::bind(UID_customer, UPRNG));

            for (;;)
            {
                while (coin_flip())
                    std::iter_swap(std::begin(cluster1) + random_customer(), std::begin(cluster2) + random_customer());

                unsigned demand1(0);
                for (allele_t allele : cluster1)
                    demand1 += input_vertices[allele - 1].demand;

                if (demand1 > truckCapacity)
                    continue;

                unsigned demand2(0);
                for (allele_t allele : cluster2)
                    demand2 += input_vertices[allele - 1].demand;

                if (demand2 <= truckCapacity)
                    break;
            }
        }


        const unsigned n(offspring.size());
#pragma omp parallel for
        for (int i = 0; unsigned(i) < n; i += 2)
        {
            const chromosome_t
                &mum(parents[random_parent()]),
                &dad(parents[random_parent()]);

            chromosome_t
                offspring1 = mum,
                offspring2 = dad;


            // Swap approximately half of the clusters
            for (unsigned i(0); i < std::min(offspring1.size(), offspring2.size()); ++i)
                if (coin_flip())
                    std::iter_swap(std::begin(offspring1) + i, std::begin(offspring2) + i);

            std::uniform_int_distribution<unsigned> UID_cluster(0, std::min(offspring1.size(), offspring2.size()) - 1);
            std::function<unsigned()> random_cluster(std::bind(UID_cluster, UPRNG));

            // (Possibly) randomise a cluster order
            if (coin_mutation_shuffle())
            {
                const unsigned i(random_cluster());
                cluster_order_t& cluster_order(coin_flip() ? offspring1[i] : offspring2[i]);
                shuffle(std::begin(cluster_order), std::end(cluster_order), UPRNG);
            }

            // (Possibly) swap the order of two adjacent customers in a cluster
            if (coin_mutation_swap())
            {
                const unsigned i(random_cluster());
                cluster_order_t& cluster_order(coin_flip() ? offspring1[i] : offspring2[i]);
                std::uniform_int_distribution<unsigned> UID_gene(0, cluster_order.size() - 1);
                std::function<unsigned()> random_gene(std::bind(UID_gene, UPRNG));

                const unsigned
                    gene(random_gene()),
                    nextGene((gene + 1) % cluster_order.size());

                std::iter_swap(std::begin(cluster_order) + gene, std::begin(cluster_order) + nextGene);
            }

            // (Possibly) reverse the order of some contiguous interval in a cluster
            if (coin_mutation_reverse())
            {
                const unsigned i(random_cluster());
                cluster_order_t& cluster_order(coin_flip() ? offspring1[i] : offspring2[i]);
                std::uniform_int_distribution<unsigned> UID_gene(0, cluster_order.size() - 1);
                std::function<unsigned()> random_gene(std::bind(UID_gene, UPRNG));
                unsigned l, r;
                std::tie(l, r) = minmax(random_gene(), random_gene());
                reverse(std::begin(cluster_order) + l, std::begin(cluster_order) + r);
            }

            // Handle case of odd population size
            if (unsigned(i + 1) < offspring.size())
            {
                offspring[i] = offspring1;
                offspring[i + 1] = offspring2;
            }
            else
                offspring[i] = coin_flip() ? offspring1 : offspring2;
        }

#if 0
        std::bernoulli_distribution BD_mutation(mutationRate);
        std::function<bool()> coin_mutation(std::bind(BD_mutation, UPRNG));

        std::bernoulli_distribution BD_mutation_reverse(mutationRate_reverse);
        std::function<bool()> coin_mutation_reverse(std::bind(BD_mutation_reverse, UPRNG));

        while (coin_mutation())
        {
            chromosome_t& mutant(offspring[random_offspring()]);
            std::uniform_int_distribution<unsigned> UID_gene(0, mutant.size() - 1);
            std::function<unsigned()> random_gene(std::bind(UID_gene, UPRNG));
            unsigned l, r;
            std::tie(l, r) = minmax(random_gene(), random_gene());

            if (coin_mutation_reverse())
                std::iter_swap(std::begin(mutant) + l, std::begin(mutant) + r);
            else
                reverse(std::begin(mutant) + l, std::begin(mutant) + r);
        }
#endif

        return offspring;
    }


    void outputSolution(const chromosome_t& solution, const clusters_t& clusters, const distance_measure distance)
    {
        // Output everything
        std::cout <<
            "login pj12469 53904\n"
            "name Patrick Johnston\n"
            //         ---  0---|--- 10---|--- 20---|--- 30---|--- 40---|--- 50---|--- 60---|--- 70---|--- 80---|--- 90---|
            //         123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
            "algorithm Genetic algorithm\n"
            "cost " << distance << '\n';

        for (unsigned gene(0); gene < solution.size(); ++gene)
        {
            std::cout << "1->";

            for (i_cluster_t i : solution[gene])
                std::cout << clusters[gene][i] << "->";

            std::cout << "1\n";
        }
    }


    void run(population_t& population, clusters_t& clusters, distance_measure& min_global, chromosome_t& solution_global, clusters_t& solution_clusters_global)
    {
        // Satisfiction condition
        SlidingWindow differences;
        distance_measure min_prev(0);

        for (unsigned i(0); differences.mean() >= satisfictionTheshold; ++i)
        {
            // evaluate population
            Evaluations& evaluations(evaluatePopulation(population, clusters));

            const distance_measure min(evaluations.min);
            if (min < min_global)
            {
                min_global = min;
                solution_global = population[evaluations.min_i];
                solution_clusters_global = clusters;
            }
            differences.push(std::abs(min - min_prev));
            min_prev = min;

            if (i % (n_chromosomes >> 0) == 0)
                std::cout << min << ".." << evaluations.max << " | " << differences.mean() << '\n';

            // select parents
            parents_t parents(selectParents(population, evaluations));

            // reproduce using genetic operators
            population = reproduce(parents, clusters, differences.mean());
        }

        Evaluations& evaluations(evaluatePopulation(population, clusters));
        const distance_measure min(evaluations.min);
        if (min < min_global)
        {
            min_global = min;
            solution_global = population[evaluations.min_i];
            solution_clusters_global = clusters;
        }
    }
}


namespace NS_Vertex
{
    const double
        mutationRate = 0.9;


    typedef unsigned                                                         i_Vertex_t, i_population_t;
    typedef float                                                            distance_measure;
    typedef i_Vertex_t                                                       allele_t;
    typedef std::vector<allele_t>                                            chromosome_t;
    typedef std::array<chromosome_t, n_chromosomes>                          population_t;
    typedef std::array<distance_measure, n_chromosomes>                      evaluation_t;
    typedef std::array<std::array<distance_measure, n_vertices>, n_vertices> distances_t;
    typedef std::vector<chromosome_t>                                        parents_t;
    typedef population_t                                                     offspring_t;


    const distance_measure
        margin_leastFit = -1 + 0.00001f, // Minimum -1. The greater, the more likely the least fit parents will be selected, the lesser, the more likely that fit parents will be selected
        satisfictionTheshold = 0.001f,
        penalty = 1,
        penalty_n_invalidTrucks = 1.01f;


    population_t& initialPopulation()
    {
        // Add depot after every fifth gene

        // Populate initial ordered chromosome
        static chromosome_t id_list;
        static bool id_list_init(false);
        if (!id_list_init)
        {
            id_list_init = true;
            id_list.reserve(n_genes);
            iota_n(std::back_inserter(id_list), n_genes, depot + 1); // starting at two because the begin/end is fixed
        }

        // Make lots of permutations
        static population_t population;
        std::default_random_engine UPRNG(std::random_device{}());
        
        const unsigned n(population.size());
#pragma omp parallel for
        for (int i = 0; unsigned(i) < n; ++i)
        {
            population[i] = id_list;
            shuffle(std::begin(population[i]), std::end(population[i]), UPRNG);
        }

        return population;
    }


    distance_measure evaluateChromosome(const chromosome_t& chromosome)
    {
        distance_measure distance(0);
        {
            unsigned demandFulfilled(0);
            allele_t from(depot);

            // Pre:
            //     First demand considered is of first allele == to
            //     First distance considered is depot->[first allele] == from->to
            // In:
            //     Add from->(depot->)to to distance
            //     Add to.demand to demandFulfilled
            // Post:
            //     Last demand considered is of last allele == from
            //     Last distance considered is [prev. last allele]->[last allele] == _->from
            for (allele_t to : chromosome)
            {
                const distance_measure
                    distance_from_to(distances[from - 1][to - 1]),
                    distance_from_depot_to(distances[from - 1][depot - 1] + distances[depot - 1][to - 1]);

                if (demandFulfilled + input_vertices[to - 1].demand > truckCapacity)
                {
                    distance += distance_from_depot_to;
                    demandFulfilled = input_vertices[to - 1].demand;
                }
                else
                {
                    distance += distance_from_to;
                    demandFulfilled += input_vertices[to - 1].demand;
                }

                from = to;
            }

            distance += distances[from - 1][depot - 1];
        }

        return distance;
    }


    Evaluations& evaluatePopulation(population_t& population)
    {
        static Evaluations evaluations;
        evaluations.reset();

        const unsigned n(population.size());
#pragma omp parallel for
        for (int i = 0; unsigned(i) < n; ++i)
            evaluations[i] = evaluateChromosome(population[i]);

        return evaluations;
    }


    parents_t selectParents(const population_t& population, const Evaluations& evaluations)
    {
        static std::default_random_engine UPRNG(std::random_device{}());

        const distance_measure
            min(evaluations.min),
            max(evaluations.max);

        std::uniform_real_distribution<distance_measure> URD
        (
            min,
            max + (max - min) * margin_leastFit
        );
        std::function<distance_measure()> random_real(std::bind(URD, UPRNG));

        parents_t selectedParents;
        for (;;)
        {
            for (unsigned i(0); i < population.size(); ++i)
                if (evaluations[i] <= random_real())
                    selectedParents.push_back(population[i]);
            if (selectedParents.size() >= minParents)
                return selectedParents;
        }
    }


    offspring_t& reproduce(const parents_t parents)
    {
        // want |offspring| = |population|, but |parents| <= |population|
        // so lets choose pairs of parents at random until there are enough offspring

        static offspring_t offspring;

        std::default_random_engine UPRNG(std::random_device{}());

        std::uniform_int_distribution<unsigned> UID_parent(0, parents.size() - 1);
        std::function<unsigned()> random_parent(std::bind(UID_parent, UPRNG));

        std::uniform_int_distribution<unsigned> UID_offspring(0, offspring.size() - 1);
        std::function<unsigned()> random_offspring(std::bind(UID_offspring, UPRNG));

        std::bernoulli_distribution BD_fair;
        std::function<bool()> coin_flip(std::bind(BD_fair, UPRNG));

        const unsigned n(offspring.size());
#pragma omp parallel for
        for (int i = 0; unsigned(i) < n; i += 2)
        {
            const chromosome_t
                &mum(parents[random_parent()]),
                &dad(parents[random_parent()]);

            chromosome_t
                offspring1 = mum,
                offspring2 = dad;

            std::uniform_int_distribution<unsigned> UID_gene(0, std::min(mum.size(), dad.size()));
            std::function<unsigned()> random_gene(std::bind(UID_gene, UPRNG));

            unsigned l, r;
            std::tie(l, r) = minmax(random_gene(), random_gene());

            for (unsigned i(l); i < r; ++i)
            {
                i_Vertex_t
                    allele1(mum[i]),
                    allele2(dad[i]);

                std::iter_swap
                (
                    find_circular(std::begin(offspring1), std::end(offspring1), std::begin(offspring1) + random_gene(), allele1),
                    find_circular(std::begin(offspring1), std::end(offspring1), std::begin(offspring1) + random_gene(), allele2)
                );
                std::iter_swap
                (
                    find_circular(std::begin(offspring2), std::end(offspring2), std::begin(offspring2) + random_gene(), allele1),
                    find_circular(std::begin(offspring2), std::end(offspring2), std::begin(offspring2) + random_gene(), allele2)
                );
            }

            // Handle case of odd population size
            if (unsigned(i + 1) < offspring.size())
            {
                offspring[i] = offspring1;
                offspring[i + 1] = offspring2;
            }
            else
                offspring[i] = coin_flip() ? offspring1 : offspring2;
        }

        std::bernoulli_distribution BD_mutation(mutationRate);
        std::function<bool()> coin_mutation(std::bind(BD_mutation, UPRNG));

        while (coin_mutation())
        {
            chromosome_t& mutant(offspring[random_offspring()]);
            std::uniform_int_distribution<unsigned> random_gene(0, mutant.size() - 1);
            std::swap(mutant[random_gene(UPRNG)], mutant[random_gene(UPRNG)]);
        }

        return offspring;
    }


    auto verboseChromosome(const chromosome_t& chromosome_in) -> std::pair<chromosome_t, distance_measure>
    {
        chromosome_t chromosome_out;
        distance_measure distance(0);
        {
            unsigned demandFulfilled(0);
            allele_t from(depot);

            // Pre:
            //     First demand considered is of first allele == to
            //     First distance considered is depot->[first allele] == from->to
            //     First allele appended is of first allele = to
            // In:
            //     Add from->(depot->)to to distance
            //     Add to.demand to demandFulfilled
            //     Append (depot, )to to chromosome_out
            // Post:
            //     Last demand considered is of last allele == from
            //     Last distance considered is [prev. last allele]->[last allele] == _->from
            //     Last allele appended is last allele == from
            for (allele_t to : chromosome_in)
            {
                const distance_measure
                    distance_from_to(distances[from - 1][to - 1]),
                    distance_from_depot_to(distances[from - 1][depot - 1] + distances[depot - 1][to - 1]);

                if (demandFulfilled + input_vertices[to - 1].demand > truckCapacity)
                {
                    chromosome_out.insert(std::end(chromosome_out), {depot, to});
                    distance += distance_from_depot_to;
                    demandFulfilled = input_vertices[to - 1].demand;
                }
                else
                {
                    chromosome_out.push_back(to);
                    distance += distance_from_to;
                    demandFulfilled += input_vertices[to - 1].demand;
                }

                from = to;
            }

            distance += distances[from - 1][depot - 1];
        }

        return{chromosome_out, distance};
    }


    void outputSolution(const chromosome_t& solution, distance_measure distance)
    {
        // Output everything
        std::cout <<
            "login pj12469 53904\n"
            "name Patrick Johnston\n"
            //         ---  0---|--- 10---|--- 20---|--- 30---|--- 40---|--- 50---|--- 60---|--- 70---|--- 80---|--- 90---|
            //         123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|
            "algorithm Genetic algorithm\n"
            "cost " << distance << '\n';

        std::cout << "1->";
        for (allele_t allele : solution)
            if (allele == depot)
                std::cout << "1\n1->";
            else
                std::cout << allele << "->";
        std::cout << "1\n";
    }


    void run(population_t& population, distance_measure& min_global, chromosome_t& solution_global)
    {
        // Satisfiction condition
        SlidingWindow differences;
        distance_measure min_prev(0);

        for (unsigned i(0); differences.mean() >= satisfictionTheshold; ++i)
        {
            // evaluate population
            Evaluations& evaluations(evaluatePopulation(population));

            const distance_measure min(evaluations.min);
            if (min < min_global)
            {
                min_global = min;
                solution_global = population[evaluations.min_i];
            }
            differences.push(std::abs(min - min_prev));
            min_prev = min;

            if (i % (n_chromosomes >> 2) == 0)
                std::cout << min << ".." << evaluations.max << '\n';

            // select parents
            parents_t parents(selectParents(population, evaluations));

            // reproduce using genetic operators
            population = reproduce(parents);
        }

        // choose chromosome with best evaluation
        Evaluations& evaluations(evaluatePopulation(population));
        const distance_measure min(evaluations.min);
        if (min < min_global)
        {
            min_global = min;
            solution_global = population[evaluations.min_i];
        }
    }
}


int main()
{
    const auto doNotPrint(
        //std::ios::goodbit
        std::ios::badbit
    );

    std::cout.setstate(doNotPrint);
    using clock = std::chrono::high_resolution_clock;

    auto clock_start(clock::now());
    std::cout.sync_with_stdio(false);

    // Construct distance matrix
    initialiseDistances();

    // Construct initial chromosomes
    static population_vertex_t population_vertex;
    static population_t population;
    static clusters_t clusters;
    NS_Cluster::initialPopulation(population, clusters);
    distance_measure min_global(10000), min_global_vertex(10000);
    chromosome_t solution_global;
    chromosome_vertex_t solution_vertex_global;
    clusters_t solution_clusters_global;

    for (unsigned i(2); i; --i)
    {
        NS_Cluster::run(population, clusters, min_global, solution_global, solution_clusters_global);

        std::cout << "Switching to vertex list representation\n" << std::flush;
        population_clusterToVertex(population_vertex, population, clusters);
        NS_Vertex::run(population_vertex, min_global_vertex, solution_vertex_global);
        
        if (i == 1)
            break;

        std::cout << "Switching to cluster reference representation\n" << std::flush;
        population_vertexToCluster(population_vertex, population, clusters);
    }

    std::cout.clear();
    if (min_global <= min_global_vertex)
        NS_Cluster::outputSolution(solution_global, solution_clusters_global, min_global);
    else
        NS_Vertex::outputSolution(solution_vertex_global, min_global_vertex);
    std::cout.setstate(doNotPrint);

    auto clock_end(clock::now());
    std::cout << std::chrono::duration<float>(clock_end - clock_start).count() << " (s)\n";
}