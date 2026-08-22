/**
 * @file Main_3D.cpp
 * @brief CPU entry for 3D Euler MVP.
 */
#include "Euler3D.hpp"

#include <cstdlib>
#include <iostream>
#include <sys/stat.h>

static void ensure_parent_dirs(const std::string &path)
{
    auto pos = path.find_last_of('/');
    if (pos == std::string::npos)
        return;
    std::string dir = path.substr(0, pos);
    std::string cur;
    for (size_t i = 0; i < dir.size(); ++i)
    {
        cur.push_back(dir[i]);
        if (dir[i] == '/' || i + 1 == dir.size())
        {
            if (cur.empty() || cur == "/")
                continue;
            mkdir(cur.c_str(), 0755);
        }
    }
}

int main(int argc, char **argv)
{
    const std::string cfg_path =
        (argc > 1) ? argv[1] : "json_Files/Run_3D_Sod_LLF_smoke.json";
    try
    {
        euler3d::Config cfg = euler3d::load_config(cfg_path);
        cfg.use_cuda = false;
        ensure_parent_dirs(cfg.vtk_out);
        return euler3d::run_host(cfg);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Fatal: " << e.what() << "\n";
        return 1;
    }
}
