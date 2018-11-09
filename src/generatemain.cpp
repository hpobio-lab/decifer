/*
 * generate.cpp
 *
 *  Created on: 19-oct-2017
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "generativemodel.h"
#include <lemon/arg_parser.h>
#include <fstream>

int main(int argc, char** argv)
{
  int seed = 0;
  int m;
  int n;
  int k;
  int coverage = 200;
  int maxXY = 3;
  int maxCN = 3;
  std::string outputDirectory = "./";
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("s", "Random number generator seed (default: 0)", seed, false)
    .refOption("C", "Target coverage", coverage, true)
    .refOption("maxXY", "Maximum number of maternal/paternal copies (default: 3)", maxXY, false)
    .refOption("maxCN", "Maximum number of copy number events (default: 3)", maxCN, false)
    .refOption("k", "Number of clusters", k, true)
    .refOption("m", "Number of samples", m, true)
    .refOption("n", "Number of SNVs", n, true)
    .refOption("o", "Output directory (default: '.')", outputDirectory);
  ap.parse();
  
  g_rng = std::mt19937(seed);
  
  GenerativeModel genModel(m, n, k, maxXY, maxCN, coverage);
  genModel.generate();
  
  char buf[1024];
  snprintf(buf, 1024, "%s/m%d_n%d_k%d_C%d_s%d",
           outputDirectory.c_str(), m, n, k, coverage, seed);
  
  std::ofstream outR((std::string(buf) + ".input.tsv").c_str());
  outR << genModel.R();
  outR.close();

  std::ofstream outSNV((std::string(buf) + ".SNV.tsv").c_str());
  genModel.writeMutationProperties(outSNV);
  outSNV.close();
  
  std::ofstream outS((std::string(buf) + ".S.true").c_str());
  outS << genModel.S();
  outS.close();
  
  std::ofstream outPhyloT((std::string(buf) + ".phyloT.true").c_str());
  outPhyloT << genModel.phyloT();
  outPhyloT.close();
  
  std::ofstream outT((std::string(buf) + ".T.dot").c_str());
  genModel.writeDOT(outT);
  outT.close();
  
  std::ofstream outPhyloTDOT((std::string(buf) + ".phyloT.dot").c_str());
  genModel.phyloT().writeDOT(genModel.S(), outPhyloTDOT);
  outPhyloTDOT.close();
  
  std::ofstream outPurities((std::string(buf) + ".purities").c_str());
  genModel.writeSamplePurities(outPurities);
  outPurities.close();
  
  std::ofstream outF((std::string(buf) + ".F.true").c_str());
  genModel.writeFandPi(outF);
  outF.close();
  
  std::ofstream outClustering((std::string(buf) + ".clustering.true").c_str());
  genModel.writeClustering(outClustering);
  outClustering.close();
  
//  for (int i = 0; i < n; ++i)
//  {
//    snprintf(buf, 1024, "%s/m%d_n%d_k%d_C%d_s%d.S%d.dot",
//             outputDirectory.c_str(), m, n, k, coverage, seed, i);
//
//    std::ofstream outS_i(buf);
//    genModel.S()[i].writeDOT(m, outS_i);
//  }
  
  return 0;
}
