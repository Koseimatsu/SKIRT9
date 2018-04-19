/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "StaticPaths.hpp"
#include "FatalError.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include <mutex>
#include <unordered_map>

////////////////////////////////////////////////////////////////////

namespace
{
    // flag becomes true if the static paths have been initialized
    std::once_flag _initialized;

    // the resource paths:  <filename, complete_path>
    std::unordered_map<string,string> _resourcePaths;

    // relative paths to check for presence of built-in resources
    const char* _intpaths[] = { "../../../git/SKIRT/resources", "../../../../git/SKIRT/resources" };
    const int _Nintpaths = sizeof(_intpaths) / sizeof(const char*);

    // relative paths to check for presence of external resources
    const char* _extpaths[] = {"../../../resources", "../../../../resources" };
    const int _Nextpaths = sizeof(_extpaths) / sizeof(const char*);

    // recursively searches the given directory and its subdirectories for any files
    // and adds the corresponding canonical paths to the resource dictionary
    void findResourcePathsIn(string directory)
    {
        // search the current level
        for (const string& filename : System::filesInDirectory(directory))
        {
            if (!_resourcePaths.count(filename))
                _resourcePaths.emplace(filename, System::canonicalPath(StringUtils::joinPaths(directory, filename)));
        }

        // search the subdirectories
        for (const string& subdirname : System::dirsInDirectory(directory))
        {
            findResourcePathsIn(StringUtils::joinPaths(directory, subdirname));
        }
    }

    // populates the dictionary with resource paths, or throws an error if there is a problem
    void findResourcePaths()
    {
        // get the executable path (or the empty string in case of failure)
        string executableFilePath = System::executablePath();
        if (executableFilePath.empty()) throw FATALERROR("Could not determine path to executable");

        // get the location of the executable (i.e. the path to the containing directory)
        string executableDirPath = StringUtils::dirPath(executableFilePath);

        // iterate over the relative paths for built-in resources
        for (int i=0; i<_Nintpaths; i++)
        {
            string test = StringUtils::joinPaths(executableDirPath, _intpaths[i]);
            if (System::isDir(test)) findResourcePathsIn(test);
        }

        // at this point we should have found at least some internal resources
        if (_resourcePaths.empty())
            throw FATALERROR("Could not locate built-in resources relative to '" + executableDirPath + "'");

        // iterate over the relative paths for external resources
        for (int i=0; i<_Nextpaths; i++)
        {
            string test = StringUtils::joinPaths(executableDirPath, _extpaths[i]);
            if (System::isDir(test)) findResourcePathsIn(test);
        }
    }
}

////////////////////////////////////////////////////////////////////

string StaticPaths::resource(string name)
{
    // initialize the resource paths if needed
    std::call_once(_initialized, findResourcePaths);

    if (!_resourcePaths.count(name))
    {
        throw FATALERROR("Could not locate resource '" + name + "'"
                         "\nDownload additional resources from www.skirt.ugent.be");
    }
    return _resourcePaths.at(name);
}

////////////////////////////////////////////////////////////////////