import FWCore.ParameterSet.Config as cms

hggPhotonIDCuts = cms.PSet(

#  ----------------------------sob value=0.0002     end of iteration 8          Loose  ---------------------------
      cutsubleadisosumoet0 = cms.vdouble(     8.2,       4.1,       5.4,       2.6),
   cutsubleadisosumoetbad0 = cms.vdouble(      67,        69,        85,       7.2),
    cutsubleadtrkisooetom0 = cms.vdouble(     7.5,       4.5,       5.2,       2.5),
          cutsubleadsieie0 = cms.vdouble(  0.0112,    0.0102,     0.029,     0.028),
         cutsubleadhovere0 = cms.vdouble(    0.09,     0.089,     0.101,     0.073),
             cutsubleadr90 = cms.vdouble(    0.94,      0.31,      0.92,      0.29),
  cutsublead_drtotk_25_990 = cms.vdouble(    0.26,     0.029,    0.0062,    0.0055),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.940332        fake=0.0848119        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0004     end of iteration 6          Medium  ---------------------------
      cutsubleadisosumoet1 = cms.vdouble(     6.4,       3.2,       3.4,       2.2),
   cutsubleadisosumoetbad1 = cms.vdouble(      64,      10.8,        13,       3.5),
    cutsubleadtrkisooetom1 = cms.vdouble(     6.4,       3.4,       3.8,       2.1),
          cutsubleadsieie1 = cms.vdouble(  0.0109,      0.01,     0.029,     0.028),
         cutsubleadhovere1 = cms.vdouble(   0.089,     0.079,      0.09,     0.061),
             cutsubleadr91 = cms.vdouble(    0.94,      0.32,      0.94,      0.29),
  cutsublead_drtotk_25_991 = cms.vdouble(    0.98,     0.029,    0.0109,    0.0111),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.918779        fake=0.0596949        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0008     end of iteration 6          Tight  ---------------------------
      cutsubleadisosumoet2 = cms.vdouble(     4.7,       2.8,       2.5,      1.46),
   cutsubleadisosumoetbad2 = cms.vdouble(      62,       5.2,       7.3,       2.5),
    cutsubleadtrkisooetom2 = cms.vdouble(     4.7,       2.9,       3.8,      1.63),
          cutsubleadsieie2 = cms.vdouble(  0.0107,    0.0099,     0.028,     0.027),
         cutsubleadhovere2 = cms.vdouble(   0.087,     0.065,     0.087,      0.05),
             cutsubleadr92 = cms.vdouble(    0.94,      0.34,      0.94,      0.29),
  cutsublead_drtotk_25_992 = cms.vdouble(       1,     0.029,     0.021,     0.028),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.887763        fake=0.042122        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0016     end of iteration 6          SuperTight  ---------------------------
      cutsubleadisosumoet3 = cms.vdouble(     3.8,       2.2,      1.77,      1.29),
   cutsubleadisosumoetbad3 = cms.vdouble(    11.7,       3.4,       3.9,      1.84),
    cutsubleadtrkisooetom3 = cms.vdouble(     3.5,       2.2,       2.3,      1.45),
          cutsubleadsieie3 = cms.vdouble(  0.0106,    0.0097,     0.028,     0.027),
         cutsubleadhovere3 = cms.vdouble(   0.082,     0.062,     0.065,     0.048),
             cutsubleadr93 = cms.vdouble(    0.94,      0.36,      0.94,      0.32),
  cutsublead_drtotk_25_993 = cms.vdouble(       1,     0.062,      0.97,      0.97),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.844291        fake=0.0290732        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0032     end of iteration 6          HyperTight1  ---------------------------
      cutsubleadisosumoet4 = cms.vdouble(     3.2,      1.76,      1.39,      1.18),
   cutsubleadisosumoetbad4 = cms.vdouble(     6.1,       2.7,       2.8,      0.66),
    cutsubleadtrkisooetom4 = cms.vdouble(     3.4,      1.86,      1.67,      1.44),
          cutsubleadsieie4 = cms.vdouble(  0.0104,    0.0094,     0.028,     0.025),
         cutsubleadhovere4 = cms.vdouble(   0.076,      0.03,     0.047,     0.046),
             cutsubleadr94 = cms.vdouble(    0.94,      0.41,      0.94,      0.34),
  cutsublead_drtotk_25_994 = cms.vdouble(       1,      0.97,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.774013        fake=0.0190779        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.00625     end of iteration 6          HyperTight2  ---------------------------
      cutsubleadisosumoet5 = cms.vdouble(     2.6,      1.31,      1.33,      0.82),
   cutsubleadisosumoetbad5 = cms.vdouble(     5.1,      1.62,      1.38,     -0.22),
    cutsubleadtrkisooetom5 = cms.vdouble(     2.9,       1.6,      1.55,      1.44),
          cutsubleadsieie5 = cms.vdouble(  0.0101,    0.0093,     0.027,     0.023),
         cutsubleadhovere5 = cms.vdouble(   0.048,    0.0189,     0.032,    0.0085),
             cutsubleadr95 = cms.vdouble(    0.94,      0.47,      0.94,      0.52),
  cutsublead_drtotk_25_995 = cms.vdouble(       1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.676391        fake=0.0121019        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0125     end of iteration 6          HyperTight3  ---------------------------
      cutsubleadisosumoet6 = cms.vdouble(    1.85,      0.96,      1.21,    -0.029),
   cutsubleadisosumoetbad6 = cms.vdouble(     3.7,      0.97,      1.38,     -0.88),
    cutsubleadtrkisooetom6 = cms.vdouble(    1.93,       1.4,      1.48,     0.056),
          cutsubleadsieie6 = cms.vdouble(  0.0099,    0.0092,     0.027,     0.023),
         cutsubleadhovere6 = cms.vdouble(   0.042,    0.0173,     0.023,    0.0085),
             cutsubleadr96 = cms.vdouble(    0.94,      0.69,      0.97,      0.52),
  cutsublead_drtotk_25_996 = cms.vdouble(       1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.524271        fake=0.00631764        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.025     end of iteration 6          HyperTight4  ---------------------------
      cutsubleadisosumoet7 = cms.vdouble(    1.31,       0.3,      1.15,    -0.029),
   cutsubleadisosumoetbad7 = cms.vdouble(    1.72,      0.69,      1.14,     -0.88),
    cutsubleadtrkisooetom7 = cms.vdouble(    1.42,      0.76,      1.48,     0.056),
          cutsubleadsieie7 = cms.vdouble(  0.0098,     0.009,     0.026,     0.023),
         cutsubleadhovere7 = cms.vdouble(   0.037,   0.00049,    0.0198,   0.00024),
             cutsubleadr97 = cms.vdouble(    0.94,      0.69,      0.97,      0.73),
  cutsublead_drtotk_25_997 = cms.vdouble(       1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.383546        fake=0.00362626        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.045     end of iteration 6          HyperTight5  ---------------------------
      cutsubleadisosumoet8 = cms.vdouble(    0.94,     0.112,     0.039,    -0.029),
   cutsubleadisosumoetbad8 = cms.vdouble(    0.86,      0.45,     -0.37,     -0.88),
    cutsubleadtrkisooetom8 = cms.vdouble(    1.21,      0.51,      0.27,     0.056),
          cutsubleadsieie8 = cms.vdouble(  0.0097,     0.009,     0.026,     0.023),
         cutsubleadhovere8 = cms.vdouble(   0.028,   1.4e-05,    0.0198,     7e-06),
             cutsubleadr98 = cms.vdouble(    0.94,      0.69,      0.97,      0.73),
  cutsublead_drtotk_25_998 = cms.vdouble(       1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.265878        fake=0.00218518        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.08     end of iteration 6          HyperTight6  ---------------------------
      cutsubleadisosumoet9 = cms.vdouble(     0.3,     0.072,     -0.04,     -0.97),
   cutsubleadisosumoetbad9 = cms.vdouble(    0.59,      0.42,     -0.84,     -1.77),
    cutsubleadtrkisooetom9 = cms.vdouble(    0.43,       0.4,    0.0075,     0.056),
          cutsubleadsieie9 = cms.vdouble(  0.0094,     0.009,     0.024,     0.023),
         cutsubleadhovere9 = cms.vdouble(  0.0071,         0,   0.00055,         0),
             cutsubleadr99 = cms.vdouble(    0.95,      0.69,      0.97,      0.84),
  cutsublead_drtotk_25_999 = cms.vdouble(       1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.109842        fake=0.000923024        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.15     end of iteration 6          HyperTight7  ---------------------------
      cutsubleadisosumoet10 = cms.vdouble(   0.146,     0.058,     -0.04,     -1.16),
   cutsubleadisosumoetbad10 = cms.vdouble(    0.42,      0.41,     -0.84,     -1.77),
    cutsubleadtrkisooetom10 = cms.vdouble(    0.43,    0.0111,    0.0075,    0.0073),
          cutsubleadsieie10 = cms.vdouble(  0.0094,     0.009,     0.024,     0.023),
         cutsubleadhovere10 = cms.vdouble(  0.0002,         0,   1.6e-05,         0),
             cutsubleadr910 = cms.vdouble(    0.95,      0.69,      0.97,      0.85),
  cutsublead_drtotk_25_9910 = cms.vdouble(       1,         1,         1,         1)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.0898582        fake=0.00079285        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
)
