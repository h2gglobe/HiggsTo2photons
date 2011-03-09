#!/bin/bash

python looper.py /vols/cms02/mjarvis/ntuples/h2gred/Born_10
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/Born_25
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/Born_250
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/Box_10
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/Box_25
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/Box_250
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/DYEEM1020
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/DYEEM20
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/QCD40
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/GJet20
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/VBF130
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/GGH130
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/WZTTH130
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/VBF120
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/GGH120
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/WZTTH120
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/VBF115
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/GGH115
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/WZTTH115
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/Run2010A
python looper.py /vols/cms02/mjarvis/ntuples/h2gred/Run2010B

hadd -f hist/Data.root hist/Run2010A_hist.root hist/Run2010B_hist.root
