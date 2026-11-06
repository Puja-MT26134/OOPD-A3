#include "sim_core.h"

// one run of the sweep
struct Job
{
    int id;
    double beta;
    double v;
    uint64_t seed;
};

// a run that has been started (a future cannot be copied, only moved)
struct Running
{
    Job job;
    std::future<Summary> result;
};

int main(int argc, char *argv[])
{
    // settings of the sweep (can be changed from the command line)
    int N = 100;
    int days = 300;
    int numSeeds = 3;
    int workers = 4;       // R: how many runs go at the same time
    int threadsPerRun = 1; // T: threads inside each run
    int maxRuns = 0;       // 0 = run everything
    double timeLimit = 0;  // seconds, 0 = no limit
    std::string outFile = "sweep_results.csv";

    for (int i = 1; i < argc; i += 2)
    {
        if (i + 1 >= argc)
        {
            std::cout << "Missing value for " << argv[i] << "\n";
            return 1;
        }
        std::string flag = argv[i];
        std::string value = argv[i + 1];

        if (flag == "--N") N = std::stoi(value);
        else if (flag == "--days") days = std::stoi(value);
        else if (flag == "--seeds") numSeeds = std::stoi(value);
        else if (flag == "--workers") workers = std::stoi(value);
        else if (flag == "--threads") threadsPerRun = std::stoi(value);
        else if (flag == "--max-runs") maxRuns = std::stoi(value);
        else if (flag == "--time-limit") timeLimit = std::stod(value);
        else if (flag == "--out") outFile = value;
        else
        {
            std::cout << "Unknown flag: " << flag << "\n";
            return 1;
        }
    }
    if (workers < 1) workers = 1;

    // every combination of beta, v and seed
    double betas[] = {0.04, 0.06, 0.08, 0.10, 0.12};
    double vaccs[] = {0.0, 0.2, 0.4, 0.6};

    std::vector<Job> jobs;
    for (double b : betas)
    {
        for (double v : vaccs)
        {
            for (int s = 1; s <= numSeeds; s++)
            {
                Job j;
                j.id = (int)jobs.size();
                j.beta = b;
                j.v = v;
                j.seed = s;
                jobs.push_back(j);
            }
        }
    }
    if (maxRuns > 0 && maxRuns < (int)jobs.size())
    {
        jobs.resize(maxRuns);
    }
    int total = (int)jobs.size();

    // results file
    std::ofstream out(outFile);
    if (!out)
    {
        std::cout << "Cannot open " << outFile << "\n";
        return 1;
    }
    out << std::setprecision(8);
    out << "run,beta,v,seed,peak_fraction,peak_day,total_infected,total_deaths,days\n";

    std::vector<Running> running;
    std::vector<bool> done(total, false);
    int nextJob = 0;
    int finished = 0;
    bool timeIsUp = false;

    auto startTime = std::chrono::steady_clock::now();

    while (true)
    {
        // has the time limit run out?
        if (timeLimit > 0 && !timeIsUp)
        {
            double elapsed = std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - startTime).count();
            if (elapsed >= timeLimit)
            {
                timeIsUp = true;
                std::cout << "Time limit reached: no new runs will be started\n";
            }
        }

        // start new runs while there is room (not after the time is up)
        while (!timeIsUp && nextJob < total && (int)running.size() < workers)
        {
            Params p;
            p.N = N;
            p.days = days;
            p.beta = jobs[nextJob].beta;
            p.v = jobs[nextJob].v;
            p.seed = jobs[nextJob].seed;
            if (threadsPerRun > 1)
            {
                p.mode = "pool";
                p.threads = threadsPerRun;
                p.chunks = threadsPerRun;
            }

            Running r;
            r.job = jobs[nextJob];
            r.result = std::async(std::launch::async, [p]()
                                  { return runSimulation(p); });
            running.push_back(std::move(r));
            nextJob++;
        }

        // nothing running and nothing more to start: we are finished
        if (running.empty()) break;

        // handle the runs that have finished, in the order they finish
        bool anyFinished = false;
        for (size_t k = 0; k < running.size();)
        {
            if (running[k].result.wait_for(std::chrono::seconds(0)) == std::future_status::ready)
            {
                Job j = running[k].job;
                try
                {
                    Summary s = running[k].result.get();
                    out << j.id << "," << j.beta << "," << j.v << "," << j.seed << ","
                        << s.peakFraction << "," << s.peakDay << ","
                        << s.totalInfected << "," << s.totalDeaths << ","
                        << s.days << std::endl;
                    done[j.id] = true;
                    finished++;
                    std::cout << finished << "/" << total << " runs done (run " << j.id
                              << ": beta=" << j.beta << " v=" << j.v
                              << " seed=" << j.seed << ")\n";
                }
                catch (...)
                {
                    std::cout << "Run " << j.id << " failed\n";
                }
                running.erase(running.begin() + k);
                anyFinished = true;
            }
            else
            {
                k++;
            }
        }

        // nothing finished yet: rest for a moment
        if (!anyFinished)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
    }

    double seconds = std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - startTime).count();
    std::cout << "Finished " << finished << " of " << total
              << " runs in " << seconds << " seconds\n";

    // report the runs that are missing
    if (finished < total)
    {
        std::cout << "Missing runs:\n";
        for (const Job &j : jobs)
        {
            if (!done[j.id])
            {
                std::cout << "  run " << j.id << ": beta=" << j.beta
                          << " v=" << j.v << " seed=" << j.seed << "\n";
            }
        }
    }
    return 0;
}