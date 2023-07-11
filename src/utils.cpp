#include <sstream>  // oss

#include "utils.h"
#include "ArbitraryPrecision/ap/ap.hpp"

ap_uint<PRECISON> arbitrary_precision_integer_pow(
    const int base, const int exponent)
{
    if (exponent == 0)
    {
        ap_uint<PRECISON> val = 1;
        return val;   
    }
    ap_uint<PRECISON> val = base;
    for (int ii=0; ii<exponent - 1; ii++)
    {
        val *= base;
    }
    return val;
}

long long ipow(long long base, long long exp)
{
    assert(base > 0);
    assert(exp > 0);
    long long result = 1;
    for (;;)
    {
        if (exp & 1){result *= base;}
        exp >>= 1;
        if (!exp){break;}
        base *= base;
    }
    return result;
}

namespace state
{
    void get_neighbors_(ap_uint<PRECISON> *neighbors, ap_uint<PRECISON> n,
        int bitLength)
    {   
        for (int b=0; b<bitLength; b++)
        {
            ap_uint<PRECISON> one = 1;
            neighbors[bitLength - 1 - b] = n ^ (one << b);
        }
    }

    ap_uint<PRECISON> flip_bit(const ap_uint<PRECISON> state, const int k, const int bitLength)
    {
        const int bit = bitLength - k;
        const ap_uint<PRECISON> one = 1;
        return (state ^ (one << (bit - 1)));
    }

    void arbitrary_precision_integer_from_int_array_(
        const int *config, const int N, ap_uint<PRECISON> &res)
    {
        res = 0;
        for (int ii=N-1; ii>=0; ii--)
        {
            const ap_uint<PRECISON> config_val = config[N - ii - 1];
            const ap_uint<PRECISON> power_val = arbitrary_precision_integer_pow(2, ii);
            res += config_val * power_val;
        }
    }

    void int_array_from_arbitrary_precision_integer_(
        int *config, const int N, const ap_uint<PRECISON> &integer)
    {
        ap_uint<PRECISON> _integer = integer;
        for (int ii=0; ii<N; ii++)
        {
            ap_uint<PRECISON> remainder = _integer % 2;
            _integer = _integer / 2;
            config[N - ii - 1] = int(remainder);
        }
    }

    std::string string_rep_from_arbitrary_precision_integer(const ap_uint<PRECISON> current_state, const unsigned int N)
    {
        int* binary_array = 0;
        binary_array = new int [N];
        int_array_from_arbitrary_precision_integer_(binary_array, N, current_state);

        std::string s = "";
        for (int ii=0; ii<N; ii++)
        {
            s += std::to_string(binary_array[ii]);
        }

        delete[] binary_array;

        return s;
    }
}


namespace parameters
{

    bool _key_exists(const json inp, const std::string key)
    {
        const std::string_view string_key = key;
        if (inp.contains(std::string(string_key))){return true;}
        return false;
    }

    void log_json(const json inp)
    {
        printf("------------------- config.json loaded -------------------\n");
        for (auto it=inp.begin(); it!=inp.end(); ++it)
        {
            std::cout << it.key() << " : " << it.value() << std::endl;
        }   
    }

    void log_parameters(const SimulationParameters p)
    {
        
        printf("--------------- simulation parameters set ----------------\n");
        printf("log10_N_timesteps        \t\t\t= %i\n", p.log10_N_timesteps);
        printf("N_timesteps              \t\t\t= %lli\n", p.N_timesteps);
        printf("N_spins                  \t\t\t= %i\n", p.N_spins);
        printf("beta                     \t\t\t= %.05f\n", p.beta);
        printf("beta_critical            \t\t\t= %.05f\n", p.beta_critical);
        printf("landscape                \t\t\t= %s\n", p.landscape.c_str());
        printf("dynamics                 \t\t\t= %s\n", p.dynamics.c_str());
        printf("memory                   \t\t\t= %lli\n", p.memory);
        printf("energetic threshold      \t\t\t= %.05f\n", p.energetic_threshold);
        printf("entropic attractor       \t\t\t= %.05f\n", p.entropic_attractor);
        printf("valid_entropic_attractor \t\t\t= %i\n", p.valid_entropic_attractor);
        printf("grid_size                \t\t\t= %i\n", p.grid_size);
        printf("dw                       \t\t\t= %.05f\n", p.dw);
        printf("n_tracers_per_MPI_rank   \t\t\t= %i\n", p.n_tracers_per_MPI_rank);
        if (p.use_manual_seed)
        {
            printf("manual seed              \t\t\t= %i\n", p.seed);
        }
        printf("PRECISON                 \t\t\t= %i\n", PRECISON);
        printf("----------------------------------------------------------\n");
    }

    SimulationParameters get_parameters(const json inp)
    {
        SimulationParameters p;

        // Run assertions
        const int N_required_parameters = 6;
        const std::string required_parameters[N_required_parameters] = {
            "log10_N_timesteps",
            "N_spins",
            "landscape",
            "beta",
            "dynamics",
            "n_tracers_per_MPI_rank"
        };
        bool error = false;
        for (int ii=0; ii<N_required_parameters; ii++)
        {
            if (!_key_exists(inp, required_parameters[ii]))
            {
                printf("Key %s not found", required_parameters[ii].c_str());
                error = true;
            }
        }
        if (error){throw std::runtime_error("At least one required config key not provided");}

        // Required parameters

        p.log10_N_timesteps = inp["log10_N_timesteps"];
        p.N_timesteps = ipow(10, int(p.log10_N_timesteps));

        // Assert N_timesteps
        if (p.N_timesteps < 0)
        {
            throw std::runtime_error("Somehow N_timesteps is negative");
        }

        p.N_spins = int(inp["N_spins"]);
        p.landscape = inp["landscape"];

        // Assert landscape
        if (
            p.landscape != "EREM"
            && p.landscape != "GREM"
        )
        {
            throw std::runtime_error("Invalid landscape");
        }

        // Set beta critical
        if (p.landscape == "EREM"){p.beta_critical = 1.0;}

        // This is ~sqrt(2 ln 2)
        else{p.beta_critical = 1.177410022515475;}

        // The provided beta is actually beta/betac
        p.beta = float(inp["beta"]) * p.beta_critical;

        // Dynamics
        p.dynamics = inp["dynamics"];
        if (p.dynamics != "standard" && p.dynamics != "gillespie")
        {
            throw std::runtime_error("Invalid dynamics");
        }

        // Memory status
        if (_key_exists(inp, "memory"))
        {
            p.memory = int(inp["memory"]);
            if (p.memory < -1 || p.memory == 0)
            {
                throw std::runtime_error("Invalid memory value");
            }
        }
        else
        {
            p.memory = pow(2, 28);
        }
        
        p.n_tracers_per_MPI_rank = int(inp["n_tracers_per_MPI_rank"]);

        // Get the energy barrier information
        double et, ea;
        if (p.landscape == "EREM")
        {
            et = -1.0 / p.beta_critical * log(p.N_spins);
            ea = 1.0 / (p.beta - p.beta_critical)
                * log((2.0 * p.beta_critical - p.beta) / p.beta_critical);

            if (p.beta >= 2.0 * p.beta_critical | p.beta <= p.beta_critical)
            {
                ea = 1e16;  // Set purposefully invalid value instead of nan or inf
                p.valid_entropic_attractor = false;
            }
        }
        else if (p.landscape == "GREM")
        {
            et = -sqrt(2.0 * p.N_spins * log(p.N_spins));
            ea = -p.N_spins * p.beta / 2.0;
        }
        else
        {
            throw std::runtime_error("Invalid landscape");
        }

        p.energetic_threshold = et;
        p.entropic_attractor = ea;

        // Not-required parameters

        if (_key_exists(inp, "seed"))
        {
            p.seed = inp["seed"];
            p.use_manual_seed = true;
        }

        p.grid_size = inp.value("grid_size", p.grid_size);
        p.dw = inp.value("dw", p.dw);

        p.seed = inp.value("seed", 0);
        p.use_manual_seed = inp.value("use_manual_seed", false);

        if (_key_exists(inp, "calculate_inherent_structure_observables"))
        {
            p.calculate_inherent_structure_observables = inp["calculate_inherent_structure_observables"];
        }

        return p;
    }

    FileNames get_filenames(const int ii)
    {
        std::string ii_str = std::to_string(ii);
        ii_str.insert(ii_str.begin(), 8 - ii_str.length(), '0');

        FileNames fnames;

        // Energy
        fnames.energy = "data/" + ii_str + "_energy.txt";
        fnames.energy_IS = "data/" + ii_str + "_energy_IS.txt";

        // Ridges
        fnames.ridge_E = "data/" + ii_str + "_ridge_E.txt";
        fnames.ridge_S = "data/" + ii_str + "_ridge_S.txt";

        // Misc
        fnames.cache_size = "data/" + ii_str + "_cache_size.txt";
        fnames.acceptance_rate = "data/" + ii_str + "_acceptance_rate.txt";
        fnames.walltime_per_waitingtime = "data/" + ii_str + "_walltime_per_waitingtime.txt";

        fnames.ii_str = ii_str;
        fnames.grids_directory = "grids";
        return fnames;
    }
}


void make_directories()
{
    system("mkdir data");
    system("mkdir grids");
}


namespace grids
{

    void make_energy_grid_logspace(const int log10_timesteps, const int n_gridpoints)
    {
        std::vector <long long> v;
        v.push_back(0.0);
        const double delta = ((double) log10_timesteps) / ((double) n_gridpoints);
        for (int ii=0; ii<n_gridpoints + 1; ii++)
        {
            int val = int(pow(10, ((double) ii * delta)));
            v.push_back(val);
        }

        v.erase(unique(v.begin(), v.end()), v.end());

        const std::string fname = "grids/energy.txt";
        FILE* outfile = fopen(fname.c_str(), "w");
        for (int ii=0; ii<v.size(); ii++)
        {
            fprintf(outfile, "%lli\n", v[ii]);
        }
        fclose(outfile);
    }

    void make_pi_grids(const int log10_timesteps, const double dw, const int n_gridpoints)
    {
        std::vector <long long> v1;
        std::vector <long long> v2;
        const int nMC = int(pow(10, log10_timesteps));
        const int tw_max = int(nMC / (dw + 1.0));
        const double delta = ((double) log10(tw_max)) / ((double) n_gridpoints);

        int _v1;
        for (int ii=0; ii<n_gridpoints + 1; ii++)
        {
            _v1 = int(pow(10, ((double) ii * delta)));
            v1.push_back(_v1);
        }
        v1.erase(unique(v1.begin(), v1.end()), v1.end());

        int _v2;
        for (int ii=0; ii<v1.size(); ii++)
        {
            _v2 = int(v1[ii] * (dw + 1.0));
            v2.push_back(_v2);
        }

        const std::string fname1 = "grids/pi1.txt";
        const std::string fname2 = "grids/pi2.txt";
        FILE* outfile1 = fopen(fname1.c_str(), "w");
        FILE* outfile2 = fopen(fname2.c_str(), "w");

        for (int ii=0; ii<v1.size(); ii++)
        {
            fprintf(outfile1, "%lli\n", v1[ii]);
            fprintf(outfile2, "%lli\n", v2[ii]);
        }

        fclose(outfile1);
        fclose(outfile2);
    }

    void load_long_long_grid_(std::vector<long long> &grid, const std::string loc)
    {
        std::ifstream myfile (loc);
        std::string line;
        if (myfile.is_open())
        {
            while (getline(myfile, line))
            {
                grid.push_back(stoll(line));
            }
            myfile.close();
        }
    }

}


namespace time_utils
{
    std::string get_datetime()
    {
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
        return oss.str();
    }

    double get_time_delta(const std::chrono::time_point<std::chrono::high_resolution_clock> start)
    {
        std::chrono::time_point<std::chrono::high_resolution_clock> stop = std::chrono::high_resolution_clock::now();
        const auto dur = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        return std::chrono::duration<double>(dur).count();
    }
}



