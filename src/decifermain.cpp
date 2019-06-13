/*
 * decifermain.cpp
 *
 *  Created on: 20-oct-2018
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include <lemon/arg_parser.h>
#include <fstream>
#include <lemon/time_measure.h>
#include "precluster.h"

#ifdef CPLEX
  #include "emcplex.h"
  #include "softclusterilpcplex.h"
  #include "softclusterlpcplex.h"
  #include "softclusterlpcplexext.h"
  #include "incrementalsolver.h"
  typedef EMCplex EMAlg;
  typedef ClusterIlpCplex ClusterIlpAlg;
  typedef HardClusterIlpCplex HardClusterIlpAlg;
  typedef SoftClusterIlpCplex SoftClusterIlpAlg;
  typedef SoftClusterLpCplex SoftClusterLpAlg;

  #include "minclusterilpcplex.h"
  typedef MinClusterIlpCplex MinClusterIlpAlg;
#else
  #include "emgurobi.h"
  typedef EMGurobi EMAlg;
  typedef ClusterIlpGurobi ClusterIlpAlg;
  typedef HardClusterIlpGurobi HardClusterIlpAlg;

  #include "minclusterilpgurobi.h"
  typedef MinClusterIlpGurobi MinClusterIlpAlg;
#endif

typedef std::pair<double, BoolTensor> Solution;

Solution runEM(const ReadMatrix& R,
               const Solver::ClusterStatisticType statType,
               int k,
               int nrSegments,
               int seed,
               int nrThreads,
               int timeLimit,
               int maxIterations,
               int nrRestarts,
               const std::string& outputPrefix,
               int method,
               const BoolTensor& prevSolution,
               bool verbose,
               int downsampleSNVs,
               int localTimeLimit,
               int memoryLimit,
               bool forceTruncal,
               double betaBin)
{
  int new_seed = seed;
  if (method == 1)
  {
    new_seed = -1;
  }
  else if (method == 2)
  {
    new_seed = -2;
  }
  else if (method == 5)
  {
    new_seed = 0;
  }
  
  lemon::Timer timer;
  BoolTensor initY = prevSolution;
  
  char buf[1024];
  snprintf(buf, 1024, "%s_k%d", outputPrefix.c_str(), k);
  
  double obj = -std::numeric_limits<double>::max();
  try {
    int rr = 0;
    while (true)
    {
      if (localTimeLimit == -1 && rr == nrRestarts)
      {
        break;
      }
      else if (localTimeLimit != -1 && timer.realTime() > localTimeLimit)
      {
        break;
      }
      
      EMAlg em(R, k, nrSegments, statType, betaBin, forceTruncal);
      em.init();
      em.initHotStart(initY);
      g_rng.seed(seed + rr);
      if (em.solve(rr,
                   new_seed,
                   maxIterations,
                   nrThreads,
                   timeLimit,
                   verbose,
                   downsampleSNVs,
                   memoryLimit))
      {
        if (obj < em.getLogLikelihood())
        {
          obj = em.getLogLikelihood();
          std::ofstream outSNV((std::string(buf) + ".SNV.tsv").c_str());
          em.writeMutationProperties(outSNV);
          outSNV.close();
          
          std::ofstream outT((std::string(buf) + ".T.res").c_str());
          outT << em.getPosteriorStateTrees();
          outT.close();
          
          std::ofstream outD((std::string(buf) + ".D.res").c_str());
          outD << k << " #clusters" << std::endl;
          outD << R.getNrSamples() << " #samples" << std::endl;
          for (int j = 0; j < k; ++j)
          {
            for (int p = 0; p < R.getNrSamples(); ++p)
            {
              if (p > 0)
              {
                outD << " ";
              }
              outD << em.getD(p, j);
            }
            outD << std::endl;
          }
          for (int j = 0; j < k; ++j)
          {
            outD << em.getPi(j) << std::endl;
          }
          outD.close();
          
          initY = em.getInitY();
        }
        std::cerr << std::endl;
      }
      ++rr;
    }
  }
  catch (std::runtime_error& e)
  {
    std::cerr << e.what() << std::endl;
    obj = -std::numeric_limits<double>::max();
  }
#ifndef CPLEX
  catch (GRBException& e)
  {
    std::cerr << e.getMessage() << std::endl;
    obj = -std::numeric_limits<double>::max();
  }
#endif

  initY = prevSolution;
  return std::make_pair(obj, initY);
}

Solution runCluster(const ReadMatrix& R,
                    const Solver::ClusterStatisticType statType,
                    int k,
                    double alpha,
                    int nrThreads,
                    int timeLimit,
                    const std::string& outputPrefix,
                    const BoolTensor& prevSolution,
                    bool verbose,
                    int memoryLimit,
                    bool forceTruncal,
                    double betaBin)
{
  ClusterIlpAlg solver(R, k, alpha, statType, betaBin, forceTruncal);
  solver.init();
  solver.initHotStart(prevSolution);
  
  double logLikelihood = -std::numeric_limits<double>::max();
  if (solver.solve(nrThreads, timeLimit, verbose, memoryLimit))
  {
    char buf[1024];
    snprintf(buf, 1024, "%s_k%d", outputPrefix.c_str(), k);
    
    std::ofstream outSNV((std::string(buf) + ".SNV.tsv").c_str());
    solver.writeMutationProperties(outSNV);
    outSNV.close();
    
    std::ofstream outT((std::string(buf) + ".T.res").c_str());
    outT << solver.getStateTrees();
    outT.close();
    
    std::ofstream outD((std::string(buf) + ".D.res").c_str());
    outD << k << " #clusters" << std::endl;
    outD << R.getNrSamples() << " #samples" << std::endl;
    for (int j = 0; j < k; ++j)
    {
      for (int p = 0; p < R.getNrSamples(); ++p)
      {
        if (p > 0)
        {
          outD << " ";
        }
        outD << solver.getD(p, j);
      }
      outD << std::endl;
    }
    for (int j = 0; j < k; ++j)
    {
      outD << solver.getPi(j) << std::endl;
    }
    outD.close();
    
    logLikelihood = solver.getLogLikelihood();
  }
  
  return std::make_pair(logLikelihood, solver.getSolY());
}

Solution runSoftCluster(const ReadMatrix& R,
                        const Solver::ClusterStatisticType statType,
                        int k,
                        int nrSegments,
                        int nrThreads,
                        int timeLimit,
                        const std::string& outputPrefix,
                        const BoolTensor& prevSolution,
                        bool verbose,
                        int memoryLimit,
                        bool forceTruncal,
                        double betaBin,
                        const IntMatrix& preClustering)
{
  SoftClusterLpAlg solver(R, k, log(nrSegments)/log(2), statType, betaBin, forceTruncal, false);
  solver.init();
  solver.initHotStart(prevSolution);
  if (!preClustering.empty())
  {
    solver.initPreClusteringConstraints(preClustering);
  }
  
  double logLikelihood = -std::numeric_limits<double>::max();
  if (solver.solve(nrThreads, timeLimit, verbose, memoryLimit))
  {
    char buf[1024];
    snprintf(buf, 1024, "%s_k%d", outputPrefix.c_str(), k);
    
    std::ofstream outSNV((std::string(buf) + ".SNV.tsv").c_str());
    solver.writeMutationProperties(outSNV);
    outSNV.close();
    
    std::ofstream outT((std::string(buf) + ".T.res").c_str());
    outT << solver.getPosteriorStateTrees();
    outT.close();
    
    std::ofstream outD((std::string(buf) + ".D.res").c_str());
    outD << k << " #clusters" << std::endl;
    outD << R.getNrSamples() << " #samples" << std::endl;
    for (int j = 0; j < k; ++j)
    {
      for (int p = 0; p < R.getNrSamples(); ++p)
      {
        if (p > 0)
        {
          outD << " ";
        }
        outD << solver.getD(p, j);
      }
      outD << std::endl;
    }
    for (int j = 0; j < k; ++j)
    {
      outD << solver.getPi(j) << std::endl;
    }
    outD.close();
    
    logLikelihood = solver.getLogLikelihood();
  }
  
  return std::make_pair(logLikelihood, BoolTensor());
}

Solution runHardCluster(const ReadMatrix& R,
                        const Solver::ClusterStatisticType statType,
                        int k,
                        int nrSegments,
                        int nrThreads,
                        int timeLimit,
                        const std::string& outputPrefix,
                        const BoolTensor& prevSolution,
                        bool verbose,
                        int memoryLimit,
                        bool forceTruncal,
                        double betaBin)
{
  HardClusterIlpAlg solver(R, k, nrSegments, statType, betaBin, forceTruncal, false);
  solver.init();
  solver.initHotStart(prevSolution);
  
  double logLikelihood = -std::numeric_limits<double>::max();
  if (solver.solve(nrThreads, timeLimit, verbose, memoryLimit))
  {
    char buf[1024];
    snprintf(buf, 1024, "%s_k%d", outputPrefix.c_str(), k);
    
    std::ofstream outSNV((std::string(buf) + ".SNV.tsv").c_str());
    solver.writeMutationProperties(outSNV);
    outSNV.close();
    
    std::ofstream outT((std::string(buf) + ".T.res").c_str());
    outT << solver.getStateTrees();
    outT.close();
    
    std::ofstream outD((std::string(buf) + ".D.res").c_str());
    outD << k << " #clusters" << std::endl;
    outD << R.getNrSamples() << " #samples" << std::endl;
    for (int j = 0; j < k; ++j)
    {
      for (int p = 0; p < R.getNrSamples(); ++p)
      {
        if (p > 0)
        {
          outD << " ";
        }
        outD << solver.getD(p, j);
      }
      outD << std::endl;
    }
    for (int j = 0; j < k; ++j)
    {
      outD << solver.getPi(j) << std::endl;
    }
    outD.close();
    
    logLikelihood = solver.getLogLikelihood();
  }
  
  return std::make_pair(logLikelihood, solver.getSolY());
}

Solution runIncrementalSolver(const ReadMatrix& R,
                              const Solver::ClusterStatisticType statType,
                              int k,
                              int nrSegments,
                              int nrThreads,
                              int timeLimit,
                              const std::string& outputPrefix,
                              const BoolTensor& prevSolution,
                              bool verbose,
                              int memoryLimit,
                              bool forceTruncal,
                              double betaBin,
                              int nrRestarts,
                              int nrDownSampledSNVs)
{
  IncrementalSolver solver(R, k, nrSegments, statType, betaBin, forceTruncal);
  solver.init();
  
  double logLikelihood = -std::numeric_limits<double>::max();
  bool success = false;
//  if (nrRestarts == -1)
  {
    success = solver.solve(nrThreads, timeLimit, verbose, memoryLimit);
  }
//  else
//  {
//    success = solver.solve(nrThreads, timeLimit, verbose, memoryLimit, nrDownSampledSNVs, nrRestarts);
//  }
  if (success)
  {
    char buf[1024];
    snprintf(buf, 1024, "%s_k%d", outputPrefix.c_str(), k);
    
    std::ofstream outSNV((std::string(buf) + ".SNV.tsv").c_str());
    solver.writeMutationProperties(outSNV);
    outSNV.close();
    
    std::ofstream outT((std::string(buf) + ".T.res").c_str());
    outT << solver.getPosteriorStateTrees();
    outT.close();
    
    std::ofstream outD((std::string(buf) + ".D.res").c_str());
    outD << k << " #clusters" << std::endl;
    outD << R.getNrSamples() << " #samples" << std::endl;
    for (int j = 0; j < k; ++j)
    {
      for (int p = 0; p < R.getNrSamples(); ++p)
      {
        if (p > 0)
        {
          outD << " ";
        }
        outD << solver.getD(p, j);
      }
      outD << std::endl;
    }
    for (int j = 0; j < k; ++j)
    {
      outD << solver.getPi(j) << std::endl;
    }
    outD.close();
    
    logLikelihood = solver.getLogLikelihood();
  }
  
  return std::make_pair(logLikelihood, BoolTensor());
}

int getNrParameters(const ReadMatrix& R,
                    const Solver::ClusterStatisticType statType,
                    int k,
                    int method,
                    double betaBin)
{
  int res = R.getNrSamples() * (k+1);
  switch (method)
  {
    case 0:
    case 1:
    case 2:
    case 6:
      break;
    case 3:
    case 4:
      {
        Solver solver(R, k, 0, statType, betaBin, false);
        solver.init();
        
        for (int i = 0; i < R.getNrCharacters(); ++i)
        {
          res += k * solver.getScriptT(i).size() * k;
        }
      }
      break;
    default:
      std::cerr << "Error: invalid method" << std::endl;
      return -1;
  }
  return res;
}

Solution run(const ReadMatrix& R,
             const Solver::ClusterStatisticType statType,
             int k,
             int nrSegments,
             int seed,
             int nrThreads,
             int timeLimit,
             int maxIterations,
             int nrRestarts,
             const std::string& outputPrefix,
             int method,
             double alpha,
             const BoolTensor& prevSolution,
             bool verbose,
             int downsampleSNVs,
             int globalTimeLimit,
             int memoryLimit,
             bool forceTruncal,
             double betaBin,
             const IntMatrix& preClustering)
{
  switch (method)
  {
    case 0:
    case 1:
    case 2:
      return runEM(R, statType, k,
                   nrSegments, seed,
                   nrThreads, timeLimit,
                   maxIterations, nrRestarts,
                   outputPrefix, method,
                   prevSolution,
                   verbose, downsampleSNVs,
                   globalTimeLimit,
                   memoryLimit,
                   forceTruncal,
                   betaBin);
    case 3:
      return runCluster(R, statType, k, alpha,
                        nrThreads, timeLimit,
                        outputPrefix,
                        prevSolution,
                        verbose,
                        memoryLimit,
                        forceTruncal,
                        betaBin);
    case 4:
      return runHardCluster(R, statType, k,
                            nrSegments, nrThreads,
                            timeLimit, outputPrefix,
                            prevSolution,
                            verbose,
                            memoryLimit,
                            forceTruncal,
                            betaBin);
    case 5:
      return runEM(R, statType, k,
                   nrSegments, seed,
                   nrThreads, timeLimit,
                   maxIterations, nrRestarts,
                   outputPrefix, method,
                   prevSolution,
                   verbose, downsampleSNVs,
                   globalTimeLimit,
                   memoryLimit,
                   forceTruncal,
                   betaBin);
    case 6:
      return runSoftCluster(R, statType, k,
                            nrSegments, nrThreads,
                            timeLimit, outputPrefix,
                            prevSolution,
                            verbose,
                            memoryLimit,
                            forceTruncal,
                            betaBin,
                            preClustering);
    case 7:
      return runIncrementalSolver(R, statType, k,
                                  nrSegments, nrThreads,
                                  timeLimit, outputPrefix,
                                  prevSolution,
                                  verbose,
                                  memoryLimit,
                                  forceTruncal,
                                  betaBin,
                                  nrRestarts,
                                  downsampleSNVs);
    default:
      std::cerr << "Error: invalid method" << std::endl;
      return std::make_pair(-std::numeric_limits<double>::max(), BoolTensor());
  }
}

int main(int argc, char** argv)
{
  int seed = 0;
  std::string outputPrefix = "./output";
  int min_k = -1;
  int globalTimeLimit = -1;
  int k = -1;
  int nrSegments = 50;
  int nrRestarts = 10;
  int maxIterations = 50;
  int nrThreads = -1;
  int timeLimit = -1;
  int method = 2;
  double alpha = 0.25;
  bool bic = false;
  bool verbose = false;
  int downsampleSNVs = 25;
  int memoryLimit = -1;
  int clusterStatistic = 1;
  std::string stateTreeFilename;
  double precisionBetaBin = -1;
  bool forceTruncal;
  std::string preClusteringFilename = "-";

  lemon::ArgParser ap(argc, argv);
  ap.refOption("s", "Random number generator seed (default: 0)", seed)
    .refOption("a", "Confidence interval (CI) width used for hard clustering with CI (default: 0.25)", alpha)
    .refOption("bic", "Use Bayesian Information Criterion for model selection of #clusters", bic)
    .refOption("truncal", "Force the presence of a dominant truncal cluster in the solution", forceTruncal)
    .refOption("d", "Downsample SNVs for EM initialization (default: 25)", downsampleSNVs)
    .refOption("min_k", "Specify minimum number of clusters (only used in BIC mode, default: -1 => let ILP decide)", min_k)
    .refOption("betaBin", "Beta binomial precision parameter (default: -1, binomial model)", precisionBetaBin)
    .refOption("C", "Clustering statistic:\n" \
                    "     0 -- Cancer Cell Fraction (CCF)\n" \
                    "     1 -- Descendant Cell Fraction (DCF, default)", clusterStatistic)
    .refOption("m", "Clustering method:\n" \
                    "     0 -- Expectation-Maximization using copy-neutral SNV k-Means initialization\n" \
                    "     1 -- Expectation-Maximization using hard clustering with confidence intervals\n" \
                    "     2 -- Expectation-Maximization using hard clustering with probabilistic likelihood function (default)\n" \
                    "     3 -- Hard clustering with confidence intervals\n" \
                    "     4 -- Hard clustering with probabilistic likelihood function\n" \
                    "     5 -- Expectation-Maximization using hard clustering with distance-based likelihood function", method, false)
    .refOption("k", "Number of clusters", k, true)
    .refOption("i", "Maximum number of iterations during EM (default: 50)", maxIterations)
    .refOption("t", "Number of threads (default: -1, limited by CPU)", nrThreads)
    .refOption("tl", "Time limit in seconds (default: -1, unlimited)", timeLimit)
    .refOption("TL", "Global time limit in seconds (default: -1, unlimited)", globalTimeLimit)
    .refOption("N", "Number of segments (default: 50)", nrSegments, false)
    .refOption("r", "Number of restarts (default: 10)", nrRestarts, false)
    .refOption("o", "Output prefix (default: './output')", outputPrefix)
    .refOption("v", "Verbose (default: 0)", verbose, false)
    .refOption("ML", "Memory limit", memoryLimit)
    .refOption("S", "State tree file", stateTreeFilename, false)
    .refOption("P", "Pre clustering filename (default: no preclustering)", preClusteringFilename)
    .other("filename");
  ap.parse();
  
  if (!(0 <= clusterStatistic && clusterStatistic < 2))
  {
    std::cerr << "Error: invalid clustering statistic specified" << std::endl;
    return 1;
  }
  
  Solver::ClusterStatisticType statType = static_cast<Solver::ClusterStatisticType>(clusterStatistic);
  
  g_rng = std::mt19937(seed);
  
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
  
  bool writeStateTrees = false;
  if (!stateTreeFilename.empty())
  {
    std::ifstream inS(stateTreeFilename.c_str());
    if (!inS.good())
    {
      std::cerr << "Error: could not open '" << stateTreeFilename << "' for reading. Will generate state tree file." << std::endl;
      writeStateTrees = true;
    }
    else
    {
      StateGraph::readStateTrees(inS);
    }
  }
  
  std::cerr << "Arguments: ";
  for (int i = 1; i < argc; ++i)
  {
    std::cerr << " " << argv[i];
  }
  std::cerr << std::endl;
  std::cerr << "Input:      instance with n = " << R.getNrCharacters() << " SNVs and m = " << R.getNrSamples() << " samples" << std::endl;
  
  IntMatrix preClustering;
  if (preClusteringFilename.empty())
  {
    PreCluster pc(R);
    pc.run(nrSegments, statType, precisionBetaBin, nrThreads, timeLimit, verbose, memoryLimit);
    
    preClustering = pc.getPreClustering();
  }
  else if (preClusteringFilename != "-")
  {
    std::ifstream inC(preClusteringFilename.c_str());
    if (!inC.good())
    {
      std::cerr << "Error: failed to open '" << preClusteringFilename << "' for reading" << std::endl;
      return 1;
    }
    else
    {
      while (inC.good())
      {
        std::string line;
        getline(inC, line);
        
        StringVector s;
        boost::split(s, line, boost::is_any_of("\t "));
        
        IntVector cc;
        for (const std::string& str : s)
        {
          if (str != "")
          {
            cc.push_back(boost::lexical_cast<int>(str));
          }
        }
        
        if (!cc.empty())
        {
          preClustering.push_back(cc);
        }
      }
    }
  }
  
  std::cerr << "Algorithm:  ";
  switch (method)
  {
    case 0:
      std::cerr << "Expectation-Maximization using copy-neutral SNV k-Means initialization";
      break;
    case 1:
      std::cerr << "Expectation-Maximization using hard clustering with confidence intervals";
      break;
    case 2:
      std::cerr << "Expectation-Maximization using hard clustering with probabilistic likelihood function";
      break;
    case 3:
      std::cerr << "Hard clustering with confidence intervals";
      break;
    case 4:
      std::cerr << "Hard clustering with probabilistic likelihood function";
      break;
    case 5:
      std::cerr << "Expectation-Maximization using hard clustering with distance-based likelihood function";
      break;
    case 6:
      std::cerr << "Exact EM algorithm";
      break;
    case 7:
      std::cerr << "Incremental algorithm";
      break;
  }
  std::cerr << std::endl;
  
  if (bic)
  {
    int nrObservations = R.getNumberOfObservations();
    std::cerr << "Using BIC -- #Observations: " << nrObservations << std::endl;
    
    int minK = 0;
    if (min_k == -1)
    {
      MinClusterIlpAlg minCluster(R, k, statType, forceTruncal, precisionBetaBin);
      minCluster.init();
      minCluster.solve(nrThreads, timeLimit);
      minK = minCluster.getMinK();
    }
    else
    {
      minK = min_k;
    }
    
    BoolTensor initY;
    
    std::ofstream outBIC((outputPrefix + ".BIC.tsv").c_str());
    outBIC << "k\tloglikelihood\tobservations\tparameters\tBIC" << std::endl;
    std::cerr << "Minimum #clusters: " << minK << std::endl;
    int localTimeLimit = globalTimeLimit == -1 ? -1 : globalTimeLimit / (k - minK + 1);
    for (int kk = minK; kk <= k; ++kk)
    {
      Solution sol = run(R, statType, kk, nrSegments,
                         seed, nrThreads,
                         timeLimit, maxIterations,
                         nrRestarts, outputPrefix,
                         method, alpha, initY, verbose,
                         downsampleSNVs, localTimeLimit,
                         memoryLimit, forceTruncal, precisionBetaBin, preClustering);
      initY = sol.second;
      int nrParameters = getNrParameters(R, statType, kk, method, precisionBetaBin);
      double b = log(nrObservations) * nrParameters - 2 * sol.first;
      std::cerr << "k = " << kk << " -- Log likelihood " << sol.first << " -- BIC " << b << std::endl;
      outBIC << kk << "\t" << sol.first << "\t" << nrObservations << "\t" << nrParameters << "\t" << b << std::endl;
    }
    outBIC.close();
  }
  else
  {
    Solution sol = run(R, statType, k, nrSegments,
                       seed, nrThreads,
                       timeLimit, maxIterations,
                       nrRestarts, outputPrefix,
                       method, alpha, BoolTensor(), verbose,
                       downsampleSNVs, globalTimeLimit,
                       memoryLimit, forceTruncal,
                       precisionBetaBin, preClustering);
    std::cerr << "Log likelihood " << sol.first << std::endl;
  }
  
  if (writeStateTrees)
  {
    std::ofstream outS(stateTreeFilename);
    StateGraph::writeStateTrees(outS);
    outS.close();
  }

  return 0;
}
