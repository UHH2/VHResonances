cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
scramv1 b clean; scram build --ignore-arch; scramv1 b -j 20
cd -
