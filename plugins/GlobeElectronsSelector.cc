/* \class GenParticleSelector
* 
* Configurable GenParticle Selector
*
* \author: Luca Lista, INFN
*
*/
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
 
typedef SingleObjectSelector<reco::GsfElectronCollection, StringCutObjectSelector<reco::GsfElectron> > GlobeElectronSelector;

DEFINE_FWK_MODULE(GlobeElectronSelector);
