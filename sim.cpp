#include "sim_core.h"

int main(int argc, char *argv[])
{
    Params p = parseArgs(argc, argv);

    // Set up the two grids
    std::vector<Cell> cur(p.N * p.N);
    std::vector<Cell> next(p.N * p.N);
    initGrid(cur, p);

    std::vector<DayRow> rows;
    bool lockdown = false;
    double betaEff = p.beta * p.mask;

    std::vector<Chunk> chunks = makeChunks(p);

    if (p.snapshotEvery > 0)
    {
        std::filesystem::create_directories(p.snapshotDir);
    }
    int failedWrites = 0;

        std::deque<std::future<bool>> pending;   // writes running in the background
    int Q = p.maxPending;
    if (Q < 1) Q = 1;
    if (p.snapshotMode != "blocking" && p.snapshotMode != "async") {
        std::cout << "Unknown snapshot mode: " << p.snapshotMode << "\n";
        return 1;
    }

    // the pool is created once; it only has threads in "pool" mode
    ThreadPool pool(p.mode == "pool" ? p.threads : 0);

    // Only this loop is timed
    auto start = std::chrono::steady_clock::now();

    for (int day = 0; day < p.days; day++)
    {
        Counts t = simulateDay(cur, next, p, day, betaEff, chunks, pool);

        DayRow row;
        row.day = day;
        row.t = t;
        row.beta = betaEff;
        rows.push_back(row);

                // blocking snapshot: the simulation waits until the file is written
        /*if (p.snapshotEvery > 0 && day % p.snapshotEvery == 0) {
            std::string name = p.snapshotDir + "/snap_" + std::to_string(day) + ".ppm";
            if (!writePPM(cur, p.N, name, p.ioLatency)) failedWrites++;
        }*/

                if (p.snapshotEvery > 0 && day % p.snapshotEvery == 0) {
            std::string name = p.snapshotDir + "/snap_" + std::to_string(day) + ".ppm";

            if (p.snapshotMode == "blocking") {
                // the simulation waits until the file is written
                if (!writePPM(cur, p.N, name, p.ioLatency)) failedWrites++;
            }
            else {
                // 1. clear writes that have already finished
                while (!pending.empty() &&
                       pending.front().wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
                    finishOne(pending, failedWrites);
                }
                // 2. if Q writes are still running, wait for the oldest
                while (!pending.empty() && (int)pending.size() >= Q) {
                    finishOne(pending, failedWrites);
                }
                // 3. copy the grid and write the copy in the background
                auto copy = std::make_shared<std::vector<Cell>>(cur);
                int N = p.N;
                int latency = p.ioLatency;
                pending.push_back(std::async(std::launch::async, [copy, N, name, latency]() {
                    return writePPM(*copy, N, name, latency);
                }));
            }
        }

        // stop when nobody is infected (unless fixed-days is on)
        if (!p.fixedDays && t.i == 0)
            break;

        // today's totals decide tomorrow's betaEff
        betaEff = nextBeta(p, t, lockdown);
    }

        // wait for all writes that are still running (this is part of the timed run)
    while (!pending.empty()) {
        finishOne(pending, failedWrites);
    }

    auto end = std::chrono::steady_clock::now();
    double seconds = std::chrono::duration<double>(end - start).count();

    // Write the CSV file (not timed)
    std::ofstream out(p.outFile);
    if (!out)
    {
        std::cout << "Cannot open " << p.outFile << "\n";
        return 1;
    }
    out << std::setprecision(10);
    out << "day,S,I,R,V,D,beta_eff,checksum\n";
    for (const DayRow &row : rows)
    {
        out << row.day << "," << row.t.s << "," << row.t.i << ","
            << row.t.r << "," << row.t.v << "," << row.t.d << ","
            << row.beta << "," << row.t.checksum << "\n";
    }

    std::cout << "Simulated " << rows.size() << " days in "
              << seconds << " seconds";
              if(p.snapshotEvery > 0) std::cout << ", failed writes: " << failedWrites;
              std::cout << "\n";
    return 0;
}