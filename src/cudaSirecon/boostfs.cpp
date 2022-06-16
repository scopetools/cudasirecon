#include <boost/filesystem.hpp>
#include <regex>   // requires c++11
#include <iostream>
#include <string>
#include <algorithm> // sort

static boost::filesystem::path outputDir;

std::vector<std::string> gatherMatchingFiles(std::string target_path, std::string pattern)
{
  // Check if the pattern is specified as a full file name (i.e.; if '.tif' is in the name)
  size_t p1 = pattern.rfind(".tif");
  size_t p2 = pattern.rfind(".TIF");
  if (p1 != std::string::npos || p2 != std::string::npos) {
    // if no matches were found, rfind() returns 'string::npos'
    // do not append ".*tif" at the end of pattern
  }
  else
    pattern.append(".*\\.[tT][iI][fF]");
  pattern.insert(0, ".*");  // '.' is the wildcard in Perl regexp; '*' just means "repeat".

  const std::regex my_filter(pattern);

  std::vector< std::string > all_matching_files;

  boost::filesystem::directory_iterator end_itr; // Constructs the end iterator.
  for (boost::filesystem::directory_iterator i(target_path); i != end_itr; ++i) {

    // Skip if not a file
    if( !boost::filesystem::is_regular_file( i->status() ) ) continue;

    std::smatch what;

    auto fname = i->path().string(); //somehow this local variable is necessary with VS2017
    // With GCC 7.5 this variable is not needed; in-place "i->path().string()" is enough

    // Skip if no match
    if( !std::regex_match(fname/*i->path().string()*/, what, my_filter) ) continue;

    // File matches, store it
    all_matching_files.push_back(fname /*i->path().string()*/);
  }

  // sort file names so that earlier time points will be processed first:
  sort(all_matching_files.begin(), all_matching_files.end());


  // Create output subfolder "decon/" just under the data folder:
  outputDir = target_path;
  outputDir /= "GPUsirecon";

  if (! boost::filesystem::exists(outputDir) )
    boost::filesystem::create_directory(outputDir);

  return all_matching_files;
}


std::string makeOutputFilePath(std::string inputFileName, std::string insert)
{
  boost::filesystem::path inputpath(inputFileName);
  boost::filesystem::path outputpath(outputDir);

  std::string basename = inputpath.filename().string();
  int pos = basename.find_last_of(".tif");
  basename.insert(pos - 3, insert);

  outputpath /= basename;

  std::cout << "Output: " << outputpath.string() << '\n';
  return outputpath.string();
}
