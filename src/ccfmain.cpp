/*
 * ccfmain.cpp
 *
 *  Created on: 15-oct-2018
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "readmatrix.h"
#include <lemon/arg_parser.h>
#include <fstream>

std::map<std::string, DoubleVector> parseClustering(const std::string& filename)
{
  std::map<std::string, DoubleVector> clusterMap;
  std::ifstream inFile(filename);
  std::string line;
  getline(inFile, line);
  while (inFile.good())
  {
    getline(inFile, line);
    if (line.empty()) continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of("\t"));
    DoubleVector d;
    for (int i = 2; i < s.size(); ++i)
    {
      d.push_back(boost::lexical_cast<double>(s[i]));
    }
    clusterMap[s[0]] = d;
  }
  
  return clusterMap;
}

int main(int argc, char** argv)
{
  std::string purityString;
  std::string clusteringFilename;
  lemon::ArgParser ap(argc, argv);
  ap.other("filename")
    .refOption("p", "Sample purities", purityString)
    .refOption("c", "Clustering filename", clusteringFilename);
  ap.parse();
  
  StringVector purityStringVector;
  boost::split(purityStringVector, purityString, boost::is_any_of(","));
  DoubleVector purityVector;
  if (!purityString.empty())
  {
    for (const std::string& s : purityStringVector)
    {
      double purity = boost::lexical_cast<double>(s);
      purityVector.push_back(purity);
    }
  }
    
  if (ap.files().empty())
  {
    std::cerr << "Error: missing input filename" << std::endl;
    return 1;
  }
  
  if (ap.files().empty())
  {
    std::cerr << "Error: missing input filename" << std::endl;
    return 1;
  }
  
  std::ifstream inR(ap.files()[0].c_str());
  if (!inR.good())
  {
    std::cerr << "Error: could not open '" << ap.files()[0] << "' for reading" << std::endl;
    return 1;
  }
  
  ReadMatrix R;
  try
  {
    inR >> R;
  }
  catch (std::runtime_error& e)
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  inR.close();
  
  const int m = R.getNrSamples();
  const int n = R.getNrCharacters();
  
  if (m != purityVector.size())
  {
    std::cerr << "Incorrect number of purity values" << std::endl;
    return 1;
  }
  
  std::vector<std::vector<IntDoublePairVector> > ccfMatrix(m, std::vector<IntDoublePairVector>(n));
  for (int i = 0; i < n; ++i)
  {
    for (int p = 0; p < m; ++p)
    {
      ccfMatrix[p][i] = R.computeSangerCCFs(R.getRef(p, i),
                                            R.getVar(p, i),
                                            purityVector[p],
                                            R.getCopyNumberStates(p, i));
    }
  }
  
  IntMatrix bestIdx(m, IntVector(n, -1));
  if (!clusteringFilename.empty())
  {
    std::map<std::string, DoubleVector> clusterMap = parseClustering(clusteringFilename);

    for (int i = 0; i < n; ++i)
    {
      std::string best_cluster;
      double best_val = std::numeric_limits<double>::max();

      for (const auto& kv : clusterMap)
      {
        double val = 0;

        for (int p = 0; p < m; ++p)
        {
          double min_val_p = std::numeric_limits<double>::max();
//          int best_idx_p = -1;
          for (int idx = 0; idx < ccfMatrix[p][i].size(); ++idx)
          {
            double val_p = pow(kv.second[p] - ccfMatrix[p][i][idx].second, 2.);
            if (val_p < min_val_p)
            {
              min_val_p = val_p;
//              best_idx_p = idx;
            }
          }
//          bestIdx[p][i] = best_idx_p;
          val += min_val_p;
        }
        
        val = sqrt(val);
        
        if (val < best_val)
        {
          best_val = val;
          best_cluster = kv.first;
          
          double val2 = 0;
          for (int p = 0; p < m; ++p)
          {
            double min_val_p = std::numeric_limits<double>::max();
            int best_idx_p = -1;
            for (int idx = 0; idx < ccfMatrix[p][i].size(); ++idx)
            {
              double val_p = pow(kv.second[p] - ccfMatrix[p][i][idx].second, 2.);
              if (val_p < min_val_p)
              {
                min_val_p = val_p;
                best_idx_p = idx;
              }
            }
            bestIdx[p][i] = best_idx_p;
            val2 += min_val_p;
          }
          val2 = sqrt(val2);
        }
      }

      std::cerr << R.indexToCharacter(i) << "\t" << best_cluster << std::endl;
    }
  }
  
  std::cout << "SNV";
  for (int p = 0; p < m; ++p)
  {
    std::cout << "\t" << R.indexToSample(p);
  }
  for (int p = 0; p < m; ++p)
  {
    std::cout << "\t" << "n_chr_" << R.indexToSample(p);
  }
  std::cout << std::endl;
  for (int i = 0; i < n; ++i)
  {
    std::cout << R.indexToCharacter(i);
    for (int p = 0; p < m; ++p)
    {
      std::cout << "\t" << ccfMatrix[p][i][bestIdx[p][i]].second;
    }
    for (int p = 0; p < m; ++p)
    {
      std::cout << "\t" << ccfMatrix[p][i][bestIdx[p][i]].first;
    }
    std::cout << std::endl;
  }
  
  return 0;
}
