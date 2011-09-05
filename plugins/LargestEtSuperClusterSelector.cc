#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SortCollectionSelector.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/EtSCComparator.h"

 typedef ObjectSelector<
           SortCollectionSelector<
             reco::SuperClusterCollection, 
             GreaterBySCEt<reco::SuperCluster> 
           > 
         > LargestEtSuperClusterSelector;

DEFINE_FWK_MODULE( LargestEtSuperClusterSelector );
