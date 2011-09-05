#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SortCollectionSelector.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "CommonTools/Utils/interface/EtComparator.h"

 typedef ObjectSelector<
           SortCollectionSelector<
             reco::GsfElectronCollection, 
             GreaterByEt<reco::GsfElectron> 
           > 
         > LargestEtGsfElectronSelector;

DEFINE_FWK_MODULE( LargestEtGsfElectronSelector );
