#include <sys/stat.h>
#include <iostream> 
#include <string> 
#include <vector> 
#include <algorithm> 
#include <stdexcept> 

#include <cstdio> 
#include <cstdlib> 

#include <boost/regex.hpp>

int fexist(char*);

class ExecCommand { 
public: 
  ExecCommand( const std::string &cmd ) { 
    FILE *pfd = popen( cmd.c_str(), "r" ); 
    if ( pfd == 0 ) { 
      throw std::runtime_error( "Command or process could not be executed." ); 
    } 
    while (!feof(pfd)) { 
      char buf[ 1024 ] = {0}; 
      if (fgets(buf, sizeof(buf), pfd) > 0) { 
	m_results.push_back( std::string(buf) ); 
      } 
    } 
    pclose( pfd ); 
  } 

  const std::vector<std::string>& getResults() const { return m_results; } 
  int getTag() {
    int result = 0;
    std::string tag = "None";
    const char* pattern = "Sticky Tag:\\s+V(\\d+)_(\\d+)_(\\d+)\\s+\\(revision:";
    boost::regex re(pattern);
    boost::cmatch matches;

    for(unsigned int i=0; i<m_results.size(); i++) {
      if (boost::regex_search(m_results[i].c_str(), matches, re)) {
	result = atoi((matches[1]+matches[2]+matches[3]).c_str());
	break;
      }
    }
    
    return result;
  }
  
private:
  std::vector<std::string> m_results; 
}; 
