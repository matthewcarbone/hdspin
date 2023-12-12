#include <sstream>  // oss

#include "utils.h"
#include "ArbitraryPrecision/ap/ap.hpp"


namespace utils
{

double mean_vector(const std::vector<double> v)
{
    const double sum = std::accumulate(v.begin(), v.end(), 0.0);
    return sum / v.size();
}

double variance_vector(const std::vector<double> v)
{
    const double mean = mean_vector(v);
    const double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    return sq_sum / v.size() - mean * mean;
}

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


void get_neighbors_(ap_uint<PRECISON> *neighbors, ap_uint<PRECISON> n,
    const unsigned int bitLength)
{   
    for (int b=0; b<bitLength; b++)
    {
        ap_uint<PRECISON> one = 1;
        neighbors[bitLength - 1 - b] = n ^ (one << b);
    }
}

ap_uint<PRECISON> flip_bit(const ap_uint<PRECISON> state, const unsigned int k, const unsigned int bitLength)
{
    const int bit = bitLength - k;
    const ap_uint<PRECISON> one = 1;
    return (state ^ (one << (bit - 1)));
}

void arbitrary_precision_integer_from_int_array_(
    const unsigned int *config, const unsigned int N, ap_uint<PRECISON> &res)
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
    unsigned int *config, const unsigned int N, const ap_uint<PRECISON> &integer)
{
    ap_uint<PRECISON> _integer = integer;
    for (unsigned int ii=0; ii<N; ii++)
    {
        ap_uint<PRECISON> remainder = _integer % 2;
        _integer = _integer / 2;
        config[N - ii - 1] = int(remainder);
    }
}

std::string string_rep_from_arbitrary_precision_integer(const ap_uint<PRECISON> current_state, const unsigned int N)
{
    unsigned int* binary_array = 0;
    binary_array = new unsigned int [N];
    int_array_from_arbitrary_precision_integer_(binary_array, N, current_state);

    std::string s = "";
    for (int ii=0; ii<N; ii++)
    {
        s += std::to_string(binary_array[ii]);
    }

    delete[] binary_array;

    return s;
}

bool _key_exists(const json inp, const std::string key)
{
    const std::string_view string_key = key;
    if (inp.contains(std::string(string_key))){return true;}
    return false;
}

void print_json(const json jrep)
{
    std::cout << "hdspin - simulation parameters" << std::endl;
    std::cout << std::setw(4) << jrep << std::endl;
}


json simulation_parameters_to_json(const utils::SimulationParameters p)
{
    json j = {
        {"log10_N_timesteps", p.log10_N_timesteps},
        {"N_timesteps", p.N_timesteps},
        {"N_spins", p.N_spins},
        {"beta", p.beta},
        {"beta_critical", p.beta_critical},
        {"landscape", p.landscape},
        {"dynamics", p.dynamics},
        {"memory", p.memory},
        {"energetic_threshold", p.energetic_threshold},
        {"entropic_attractor", p.entropic_attractor},
        {"valid_entropic_attractor", p.valid_entropic_attractor},
        {"grid_size", p.grid_size},
        {"dw", p.dw},
        {"n_tracers", p.n_tracers},
        {"use_manual_seed", p.use_manual_seed},
        {"seed", p.seed},
        {"calculate_inherent_structure_observables", p.calculate_inherent_structure_observables},
        {"PRECISON", PRECISON}
    };

    return j;
}

void json_to_file(const json jrep, const std::string& filename)
{
    std::ofstream outputFile(filename);
    if (outputFile.is_open())
    {
        outputFile << std::setw(4) << jrep;
        outputFile.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}


utils::FileNames get_filenames(const unsigned int ii)
{
    std::string ii_str = std::to_string(ii);
    ii_str.insert(ii_str.begin(), 8 - ii_str.length(), '0');

    utils::FileNames fnames;

    // Energy
    fnames.energy = "data/" + ii_str + "_energy.txt";
    fnames.energy_IS = "data/" + ii_str + "_energy_IS.txt";
    fnames.max_energy = "data/" + ii_str + "_max_energy.txt";

    // Ridges
    fnames.ridge_E = "data/" + ii_str + "_ridge_E.txt";
    fnames.ridge_S = "data/" + ii_str + "_ridge_S.txt";

    // Psi
    fnames.psi_config = "data/" + ii_str + "_psi_config.txt";
    fnames.psi_basin_E = "data/" + ii_str + "_psi_basin_E.txt";
    fnames.psi_basin_S = "data/" + ii_str + "_psi_basin_S.txt";

    // Pi
    fnames.Pi_config = "data/" + ii_str + "_Pi_config.txt";
    fnames.Pi_basin_E = "data/" + ii_str + "_Pi_basin_E.txt";
    fnames.Pi_basin_S = "data/" + ii_str + "_Pi_basin_S.txt";

    // Misc
    fnames.cache_size = "data/" + ii_str + "_cache_size.txt";
    fnames.acceptance_rate = "data/" + ii_str + "_acceptance_rate.txt";
    fnames.walltime_per_waitingtime = "data/" + ii_str + "_walltime_per_waitingtime.txt";

    fnames.ii_str = ii_str;
    fnames.grids_directory = "grids";
    return fnames;
}



void make_directories()
{
    system("mkdir data");
    system("mkdir grids");
}


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
