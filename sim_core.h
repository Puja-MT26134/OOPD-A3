#pragma once

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <cmath>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <future>
#include <memory>
#include <filesystem>

// All the settings of the simulation live here
struct Params
{
    int N = 1000;             // grid size
    int days = 200;           // number of days to simulate
    uint64_t seed = 42;       // random seed
    double beta = 0.08;       // transmission probability
    int d = 7;                // days infectious
    double mu = 0.01;         // death probability
    int I0 = 10;              // initially infected
    double v = 0.0;           // vaccinated fraction
    double mask = 1.0;        // mask factor (1 = off)
    double lambda = 0.3;      // lockdown factor
    double thetaOn = 0.02;    // start lockdown above this
    double thetaOff = 0.005;  // end lockdown below this
    int threads = 1;          // number of threads (used in Part B)
    int chunks = 1;           // number of chunks (used in Part B)
    int snapshotEvery = 0;    // snapshot every K days (0 = off)
    bool useLockdown = false; // lockdown on or off
    bool fixedDays = false;   // true = never stop early
    std::string outFile = "out.csv";
    std::string partition = "strips";      // strips, blocks or interleaved
    std::string mode = "sequential";       // sequential, thread-per-day or pool
    std::string snapshotDir = "out/snaps"; // folder for the images
    int ioLatency = 0;                     // extra delay per write, in ms
        std::string snapshotMode = "blocking";   // blocking or async
    int maxPending = 4;                      // Q: most writes allowed at once
};

Params parseArgs(int argc, char *argv[])
{
    Params p; // starts with all the default values

    for (int i = 1; i < argc; i++)
    {
        std::string flag = argv[i];

        // every flag needs a value after it
        if (i + 1 >= argc)
        {
            std::cout << "Missing value for " << flag << "\n";
            exit(1);
        }
        std::string value = argv[i + 1];
        i++; // skip the value, we have already read it

        if (flag == "--N")
            p.N = std::stoi(value);
        else if (flag == "--days")
            p.days = std::stoi(value);
        else if (flag == "--seed")
            p.seed = std::stoull(value);
        else if (flag == "--beta")
            p.beta = std::stod(value);
        else if (flag == "--d")
            p.d = std::stoi(value);
        else if (flag == "--mu")
            p.mu = std::stod(value);
        else if (flag == "--I0")
            p.I0 = std::stoi(value);
        else if (flag == "--v")
            p.v = std::stod(value);
        else if (flag == "--mask")
            p.mask = std::stod(value);
        else if (flag == "--lambda")
            p.lambda = std::stod(value);
        else if (flag == "--theta-on")
            p.thetaOn = std::stod(value);
        else if (flag == "--theta-off")
            p.thetaOff = std::stod(value);
        else if (flag == "--threads")
            p.threads = std::stoi(value);
        else if (flag == "--chunks")
            p.chunks = std::stoi(value);
        else if (flag == "--snapshot-every")
            p.snapshotEvery = std::stoi(value);
        else if (flag == "--out")
            p.outFile = value;
        else if (flag == "--lockdown")
            p.useLockdown = (std::stoi(value) != 0);
        else if (flag == "--fixed-days")
            p.fixedDays = (std::stoi(value) != 0);
        else if (flag == "--partition")
            p.partition = value;
        else if (flag == "--mode")
            p.mode = value;
        else if (flag == "--snapshot-dir")
            p.snapshotDir = value;
        else if (flag == "--io-latency")
            p.ioLatency = std::stoi(value);
                    else if (flag == "--snapshot-mode") p.snapshotMode = value;
        else if (flag == "--max-pending") p.maxPending = std::stoi(value);
        else
        {
            std::cout << "Unknown flag: " << flag << "\n";
            exit(1);
        }
    }
    return p;
}

// The five possible states of a person
enum State
{
    S,
    I,
    R,
    V,
    D
};

// One person in the grid
struct Cell
{
    int state = S; // S, I, R, V or D
    int c = 0;     // days infected (only matters when state is I)
};

// Count how many of the 8 neighbours of cell (i, j) are infected
int countInfectedNeighbours(const std::vector<Cell> &grid, int N, int i, int j)
{
    int count = 0;
    for (int di = -1; di <= 1; di++)
    {
        for (int dj = -1; dj <= 1; dj++)
        {
            if (di == 0 && dj == 0)
                continue; // skip the cell itself

            int ni = i + di; // neighbour row
            int nj = j + dj; // neighbour column

            // cells outside the grid do not exist (no wrap-around)
            if (ni < 0 || ni >= N || nj < 0 || nj >= N)
                continue;

            if (grid[ni * N + nj].state == I)
                count++;
        }
    }
    return count;
}
// Mixes the bits of a number so the result looks random
uint64_t mix(uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Returns a random number between 0 and 1, decided only by the four inputs
double rand01(uint64_t seed, int day, int cell, int stream)
{
    uint64_t h = mix(seed);
    h = mix(h ^ (uint64_t)day);
    h = mix(h ^ (uint64_t)cell);
    h = mix(h ^ (uint64_t)stream);
    // keep the top 53 bits and turn them into a number in [0, 1)
    return (h >> 11) * (1.0 / 9007199254740992.0);
}

// Counts of each state in a grid
struct Counts
{
    long long s = 0, i = 0, r = 0, v = 0, d = 0;
    uint64_t checksum = 0;
};

// Count everybody in the whole grid
Counts countAll(const std::vector<Cell> &grid)
{
    Counts t;
    for (const Cell &cell : grid)
    {
        if (cell.state == S)
            t.s++;
        else if (cell.state == I)
            t.i++;
        else if (cell.state == R)
            t.r++;
        else if (cell.state == V)
            t.v++;
        else
            t.d++;
    }
    return t;
}

// Set up the starting grid: vaccinate first, then infect
void initGrid(std::vector<Cell> &grid, const Params &p)
{
    int total = p.N * p.N;

    // 1. Vaccinate a fraction v of the people
    int toVaccinate = (int)(p.v * total + 0.5);
    int done = 0;
    int attempt = 0;
    while (done < toVaccinate)
    {
        int cell = (int)(rand01(p.seed, -1, attempt, 2) * total);
        attempt++;
        if (grid[cell].state == S)
        { // skip if already vaccinated
            grid[cell].state = V;
            done++;
        }
    }

    // 2. Infect I0 susceptible people (never more than are available)
    int susceptible = total - toVaccinate;
    int toInfect = p.I0;
    if (toInfect > susceptible)
        toInfect = susceptible;

    done = 0;
    attempt = 0;
    while (done < toInfect)
    {
        int cell = (int)(rand01(p.seed, -1, attempt, 3) * total);
        attempt++;
        if (grid[cell].state == S)
        { // only infect susceptible people
            grid[cell].state = I;
            grid[cell].c = 0;
            done++;
        }
    }
}

// Add one person's state to a Counts box
void countCell(Counts &t, int id, const Cell &cell)
{
    if (cell.state == S)
        t.s++;
    else if (cell.state == I)
        t.i++;
    else if (cell.state == R)
        t.r++;
    else if (cell.state == V)
        t.v++;
    else
        t.d++;

    // checksum: position, state and counter all matter
    uint64_t key = ((uint64_t)id << 24) | ((uint64_t)cell.state << 20) | (uint64_t)(cell.c & 0xFFFFF);
    t.checksum += mix(key);
}

// One piece of the grid that a task works on
struct Chunk
{
    int rowStart, rowEnd, rowStep;
    int colStart, colEnd;
};

// Divide the grid into chunks, using the scheme in p.partition
std::vector<Chunk> makeChunks(const Params &p)
{
    std::vector<Chunk> chunks;
    int N = p.N;
    int C = p.chunks;
    if (C < 1)
        C = 1;

    if (p.partition == "strips")
    {
        for (int c = 0; c < C; c++)
        {
            Chunk ch;
            ch.rowStart = c * N / C;
            ch.rowEnd = (c + 1) * N / C;
            ch.rowStep = 1;
            ch.colStart = 0;
            ch.colEnd = N;
            chunks.push_back(ch);
        }
    }
    else if (p.partition == "interleaved")
    {
        for (int c = 0; c < C; c++)
        {
            Chunk ch;
            ch.rowStart = c;
            ch.rowEnd = N;
            ch.rowStep = C;
            ch.colStart = 0;
            ch.colEnd = N;
            chunks.push_back(ch);
        }
    }
    else if (p.partition == "blocks")
    {
        // g x g blocks, where g*g is at least C
        int g = (int)std::ceil(std::sqrt((double)C));
        for (int bi = 0; bi < g; bi++)
        {
            for (int bj = 0; bj < g; bj++)
            {
                Chunk ch;
                ch.rowStart = bi * N / g;
                ch.rowEnd = (bi + 1) * N / g;
                ch.rowStep = 1;
                ch.colStart = bj * N / g;
                ch.colEnd = (bj + 1) * N / g;
                chunks.push_back(ch);
            }
        }
    }
    else
    {
        std::cout << "Unknown partition: " << p.partition << "\n";
        exit(1);
    }
    return chunks;
}

// Update one chunk of the grid for one day.
// Reads from cur, writes only to next, returns counts of what it wrote.
Counts updateChunk(const std::vector<Cell> &cur, std::vector<Cell> &next,
                   const Params &p, int day, const double P[9],
                   const Chunk &ch)
{
    Counts t;
    int N = p.N;

    for (int i = ch.rowStart; i < ch.rowEnd; i += ch.rowStep)
    {
        for (int j = ch.colStart; j < ch.colEnd; j++)
        {
            int id = i * N + j;
            Cell newCell = cur[id]; // by default, nothing changes

            if (cur[id].state == S)
            {
                // susceptible: may get infected
                int k = countInfectedNeighbours(cur, N, i, j);
                if (rand01(p.seed, day, id, 0) < P[k])
                {
                    newCell.state = I;
                    newCell.c = 0;
                }
            }
            else if (cur[id].state == I)
            {
                // infected: one more day sick
                if (cur[id].c + 1 == p.d)
                {
                    // the illness ends today: die or recover
                    if (rand01(p.seed, day, id, 1) < p.mu)
                        newCell.state = D;
                    else
                        newCell.state = R;
                    newCell.c = 0;
                }
                else
                {
                    newCell.c = cur[id].c + 1;
                }
            }
            // R, V and D never change, so newCell is already correct

            next[id] = newCell;
            countCell(t, id, newCell);
        }
    }
    return t;
}

// Add the counts of b into a
void addCounts(Counts &a, const Counts &b)
{
    a.s += b.s;
    a.i += b.i;
    a.r += b.r;
    a.v += b.v;
    a.d += b.d;
    a.checksum += b.checksum;
}

// A pool of worker threads that run tasks which return Counts
class ThreadPool
{
public:
    // start numThreads workers (once)
    ThreadPool(int numThreads)
    {
        for (int i = 0; i < numThreads; i++)
        {
            workers.push_back(std::thread([this]()
                                          { workerLoop(); }));
        }
    }

    // tell the workers to stop, then wait for them
    ~ThreadPool()
    {
        {
            std::lock_guard<std::mutex> lock(mtx);
            stopping = true;
        }
        cv.notify_all();
        for (std::thread &w : workers)
        {
            w.join();
        }
    }

    // add a task to the waiting list, get a future for its result
    std::future<Counts> submit(std::function<Counts()> job)
    {
        auto task = std::make_shared<std::packaged_task<Counts()>>(job);
        std::future<Counts> result = task->get_future();
        {
            std::lock_guard<std::mutex> lock(mtx);
            tasks.push([task]()
                       { (*task)(); });
        }
        cv.notify_one(); // wake up one sleeping worker
        return result;
    }

private:
    // what every worker does: wait for a task, run it, repeat
    void workerLoop()
    {
        while (true)
        {
            std::function<void()> job;
            {
                std::unique_lock<std::mutex> lock(mtx);
                cv.wait(lock, [this]()
                        { return stopping || !tasks.empty(); });
                if (stopping && tasks.empty())
                    return;
                job = tasks.front();
                tasks.pop();
            }
            job(); // run the task outside the lock
        }
    }

    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex mtx;
    std::condition_variable cv;
    bool stopping = false;
};

// Simulate one whole day and swap the grids
Counts simulateDay(std::vector<Cell> &cur, std::vector<Cell> &next,
                   const Params &p, int day, double betaEff,
                   const std::vector<Chunk> &chunks, ThreadPool &pool)
{
    // P[k] = chance of infection with k infected neighbours
    double P[9];
    for (int k = 0; k <= 8; k++)
    {
        P[k] = 1.0 - std::pow(1.0 - betaEff, k);
    }

    Counts total;

    if (p.mode == "sequential")
    {
        for (const Chunk &ch : chunks)
        {
            Counts t = updateChunk(cur, next, p, day, P, ch);
            addCounts(total, t);
        }
    }
    else if (p.mode == "thread-per-day")
    {
        // each thread writes its result into its own slot
        std::vector<Counts> results(chunks.size());
        std::vector<std::thread> workers;

        for (size_t c = 0; c < chunks.size(); c++)
        {
            workers.push_back(std::thread([&, c]()
                                          { results[c] = updateChunk(cur, next, p, day, P, chunks[c]); }));
        }

        // wait for every thread to finish
        for (std::thread &w : workers)
        {
            w.join();
        }

        // the main thread adds up the counts
        for (const Counts &t : results)
        {
            addCounts(total, t);
        }
    }

    else if (p.mode == "pool")
    {
        std::vector<std::future<Counts>> futures;

        // submit one task per chunk
        for (size_t c = 0; c < chunks.size(); c++)
        {
            futures.push_back(pool.submit([&, c]()
                                          { return updateChunk(cur, next, p, day, P, chunks[c]); }));
        }

        // wait for each task and add up its counts
        for (std::future<Counts> &f : futures)
        {
            Counts t = f.get();
            addCounts(total, t);
        }
    }

    else
    {
        std::cout << "Unknown mode: " << p.mode << "\n";
        exit(1);
    }
    cur.swap(next); // tomorrow becomes today
    return total;
}

// Decide tomorrow's betaEff from today's totals.
// 'lockdown' remembers whether a lockdown is active right now.
double nextBeta(const Params &p, const Counts &t, bool &lockdown)
{
    double fraction = (double)t.i / ((double)p.N * p.N);

    if (p.useLockdown)
    {
        if (!lockdown && fraction > p.thetaOn)
            lockdown = true;
        else if (lockdown && fraction < p.thetaOff)
            lockdown = false;
    }

    double b = p.beta * p.mask;
    if (lockdown)
        b = b * p.lambda;
    return b;
}

// One line of the output file
struct DayRow
{
    int day;
    Counts t;
    double beta;
};

// Save the grid as a PPM image, one colour per state. Returns true if it worked.
bool writePPM(const std::vector<Cell> &grid, int N,
              const std::string &filename, int latencyMs)
{
    // imitate slow storage
    if (latencyMs > 0)
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(latencyMs));
    }

    std::ofstream f(filename, std::ios::binary);
    if (!f)
        return false;

    f << "P6\n"
      << N << " " << N << "\n255\n";

    std::vector<unsigned char> pixels(N * N * 3);
    for (int id = 0; id < N * N; id++)
    {
        unsigned char red = 0, green = 0, blue = 0;
        if (grid[id].state == S)
        {
            red = 60;
            green = 160;
            blue = 60;
        } // green
        else if (grid[id].state == I)
        {
            red = 220;
            green = 30;
            blue = 30;
        } // red
        else if (grid[id].state == R)
        {
            red = 50;
            green = 90;
            blue = 220;
        } // blue
        else if (grid[id].state == V)
        {
            red = 230;
            green = 210;
            blue = 40;
        } // yellow
        else
        {
            red = 20;
            green = 20;
            blue = 20;
        } // dark (dead)
        pixels[id * 3] = red;
        pixels[id * 3 + 1] = green;
        pixels[id * 3 + 2] = blue;
    }
    f.write((const char *)pixels.data(), pixels.size());
    return f.good();
}

// Wait for the oldest pending write, count it if it failed, remove it
void finishOne(std::deque<std::future<bool>>& pending, int& failedWrites) {
    bool ok = false;
    try {
        ok = pending.front().get();
    } catch (...) {
        ok = false;
    }
    if (!ok) failedWrites++;
    pending.pop_front();
}

// The short result of one finished simulation
struct Summary {
    double peakFraction = 0;     // highest infected fraction
    int peakDay = 0;             // day when the peak happened
    long long totalInfected = 0; // everybody who was ever infected
    long long totalDeaths = 0;
    int days = 0;                // number of days simulated
};

// Run one whole simulation and return only the summary (no CSV, no images)
Summary runSimulation(const Params &p)
{
    std::vector<Cell> cur(p.N * p.N);
    std::vector<Cell> next(p.N * p.N);
    initGrid(cur, p);

    std::vector<Chunk> chunks = makeChunks(p);
    ThreadPool pool(p.mode == "pool" ? p.threads : 0);

    Summary s;
    Counts t;
    bool lockdown = false;
    double betaEff = p.beta * p.mask;

    for (int day = 0; day < p.days; day++)
    {
        t = simulateDay(cur, next, p, day, betaEff, chunks, pool);
        s.days = day + 1;

        double fraction = (double)t.i / ((double)p.N * p.N);
        if (fraction > s.peakFraction)
        {
            s.peakFraction = fraction;
            s.peakDay = day;
        }

        if (!p.fixedDays && t.i == 0)
            break;
        betaEff = nextBeta(p, t, lockdown);
    }

    // people who were ever infected = still infected + recovered + dead
    s.totalInfected = t.i + t.r + t.d;
    s.totalDeaths = t.d;
    return s;
}