npr <- c(
    "Adcyap1r1",
    "Avpr1a",
    "Calcr",
    "Calcrl",
    "Cckar",
    "Cckbr",
    "Cntfr",
    "Crhr1",
    "Crhr2",
    "Esr1",
    "Galr1",
    "Galr2",
    "Galr3",
    "Ghr",
    "Ghrhr",
    "Ghsr",
    "Glp1r",
    "Gpr55",
    "Gpr83",
    "Gpr149",
    "Grpr",
    "Hcrtr1",
    "Hcrtr2",
    "Igf1r",
    "Insr",
    "Insrr",
    "Kiss1r",
    "Lepr",
    "Mc1r",
    "Mc3r",
    "Mc4r",
    "Mchr1",
    "Nmbr",
    "Nmur1",
    "Nmur2",
    "Npffr1",
    "Npffr2",
    "Npr1",
    "Npr2",
    "Npr3",
    "Npsr1",
    "Npsr2",
    "Npy1r",
    "Npy2r",
    "Npy5r",
    "Ntrk2",
    "Ntsr1",
    "Ntsr2",
    "Oprd1",
    "Oprk1",
    "Oprl1",
    "Oprm1",
    "Oxtr",
    "Prlhr",
    "Prlr",
    "Prokr2",
    "Qrfpr",
    "Rxfp1",
    "Rxfp2",
    "Sstr1",
    "Sstr2",
    "Sstr3",
    "Tacr1",
    "Tacr3",
    "Trhr",
    "Trhr2",
    "Tshr",
    "Vipr1",
    "Vipr2"
)

np <- c(
    "Adcyap1",
    "Agrp",
    "Avp",
    "Bdnf",
    "Cartpt",
    "Cck",
    "Cntf",
    "Crh",
    "Gal",
    "Ghrh",
    "Ghrl",
    "Grp",
    "Hcrt",
    "Kiss1",
    "Lep",
    "Nmb",
    "Nms",
    "Nmu",
    "Npvf",
    "Npw",
    "Npy",
    "Nts",
    "Oxt",
    "Pdyn",
    "Penk",
    "Pmch",
    "Pnoc",
    "Pomc",
    "Qrfp",
    "Reln",
    "Rln1",
    "Rln3",
    "Sst",
    "Tac1",
    "Tac2",
    "Trh"
)

irs_genes <- c(
    "Alk", "Insr", "Ltk", "Igf1r", "Irs1",
    "Ptn", "Mdk", "Fam150a", "Fam150b",
    "Mc4r", "Lepr", "Sim1", "Lmo4",
    "Slc2a1", "Slc2a3"
)

neurotrans <- c(
    "Slc17a6", "Slc17a7", "Slc17a8", "Slc1a1", "Slc1a2", "Slc1a6",
    "Gad1", "Slc32a1", "Slc6a1"
)
glut <- c("Slc17a6", "Slc17a7", "Slc17a8", "Slc1a1", "Slc1a2", "Slc1a6")
glutr <- c(
    "Gria1", "Gria2", "Gria3", "Gria4", # iGlu AMPA receptors
    "Grid1", "Grid2", # iGlu delta receptors
    "Grik1", "Grik2", "Grik3", "Grik4", "Grik5", # iGlu kainate receptors
    "Grin1", "Grin2a", "Grin2b", "Grin2c", "Grin2d", "Grin3a", "Grin3b", # iGlu NMDA receptors
    "Grm1", "Grm5", # mGluRs 1
    "Grm2", "Grm3", # mGluRs 2
    "Grm4", "Grm6", "Grm7", "Grm8" # mGluRs 3
)
gaba <- c("Gad1", "Gad2", "Slc32a1", "Slc6a1")
gabar <- c(
    "Gabra1", "Gabra2", "Gabra3", "Gabra4", "Gabra5", "Gabra6",
    "Gabrb1", "Gabrb2", "Gabrb3",
    "Gabrg1", "Gabrg2", "Gabrg3",
    "Gabrd", "Gabre", "Gabrp", "Gabrq",
    "Gabrr1", "Gabrr2", "Gabrr3",
    "Gabbr1", "Gabbr2"
)
nmr <- c(
    "Adra1a",
    "Adra1b",
    "Adra1d",
    "Adra2a",
    "Adra2b",
    "Adra2c",
    "Adrb1",
    "Adrb2",
    "Adrb3",
    "Adrbk1",
    "Adrbk2",
    "Adrm1", # adrenergic receptors
    "Adora1",
    "Adora2a",
    "Adora2b",
    "Adora3", # adenosine receptors
    "Chrm1",
    "Chrm2",
    "Chrm3",
    "Chrm4",
    "Chrm5",
    "Chrna1",
    "Chrna2",
    "Chrna3",
    "Chrna4",
    "Chrna5",
    "Chrna6",
    "Chrna7",
    "Chrna9",
    "Chrna10",
    "Chrnb1",
    "Chrnb2",
    "Chrnb3",
    "Chrnd",
    "Chrng", # cholinergic receptors
    "Gria1",
    "Gria2",
    "Gria3",
    "Gria4", # iGlu AMPA receptors
    "Grid1",
    "Grid2", # iGlu delta receptors
    "Grik1",
    "Grik2",
    "Grik3",
    "Grik4",
    "Grik5", # iGlu kainate receptors
    "Grin1",
    "Grin2a",
    "Grin2b",
    "Grin2c",
    "Grin2d",
    "Grin3a",
    "Grin3b", # iGlu NMDA receptors
    "Grm1",
    "Grm5", # mGluRs 1
    "Grm2",
    "Grm3", # mGluRs 2
    "Grm4",
    "Grm6",
    "Grm7",
    "Grm8", # mGluRs 3
    "Gabra1",
    "Gabra2",
    "Gabra3",
    "Gabra4",
    "Gabra5",
    "Gabra6",
    "Gabrb1",
    "Gabrb2",
    "Gabrb3",
    "Gabrg1",
    "Gabrg2",
    "Gabrg3",
    "Gabrd",
    "Gabre",
    "Gabrp",
    "Gabrq",
    "Gabrr1",
    "Gabrr2",
    "Gabrr3",
    "Gabbr1",
    "Gabbr2", # GABA receptors
    "Drd1",
    "Drd2",
    "Drd3",
    "Drd4",
    "Drd5", # dopamine receptors
    "Htr1a",
    "Htr1b",
    "Htr1d",
    "Htr1f",
    "Htr2a",
    "Htr2b",
    "Htr2c",
    "Htr3a",
    "Htr3b",
    "Htr4",
    "Htr5a",
    "Htr5b",
    "Htr6",
    "Htr7", # serotonin receptors
    "Gnas",
    "Gnai1",
    "Gnai2",
    "Gnai3",
    "Gnao1",
    "Gnao2",
    "Gnaq",
    "Gna11",
    "Gna12",
    "Gna13",
    "Gnal",
    "Gnasxl", # G protein alpha subunit
    "Gnb1",
    "Gnb2",
    "Gnb3",
    "Gnb4",
    "Gnb5", # G protein beta subunit
    "Gng2",
    "Gng3",
    "Gng4",
    "Gng5",
    "Gng7",
    "Gng8",
    "Gng10",
    "Gng11",
    "Gng12",
    "Gng13",
    "Gngt1",
    "Gngt2", # G protein gamma subunit
    "P2rx1",
    "P2rx2",
    "P2rx3",
    "P2rx4",
    "P2rx5",
    "P2rx6",
    "P2rx7",
    "P2ry1",
    "P2ry2",
    "P2ry4",
    "P2ry6",
    "P2ry12",
    "P2ry13",
    "P2ry14", # purinergic receptors
    "Ryr1",
    "Ryr2",
    "Ryr3" # ryanodine receptors
)

dopam <- c("Th", "Slc6a3", "Slc18a2", "Ddc", "Slc18a3", "Drd1", "Drd2", "Drd3", "Drd4", "Drd5")
sert <- c("Htr1a", "Htr1b", "Htr1d", "Htr1f", "Htr2a", "Htr2b", "Htr2c", "Htr3a", "Htr3b", "Htr4", "Htr5a", "Htr5b", "Htr6", "Htr7", "Gnai1", "Gnai3", "Gnao1", "Gnaz")
ach <- c("Chat", "Slc18a3", "Ache", "Slc5a7")
mcr_genes <- c("Mc1r", "Mc2r", "Mc3r", "Mc4r")
pvn_genes <- c("Trh", "Crh", "Mbnl3", "Pgf", "Irs4", "Gpr101", "Nr3c2", "Agtr1")
dmh_trh_g <- c("Onecut3", "Cartpt")

mitochondrial <- c(
    "Mfn1", # Mitofusin 1: A GTPase located on the outer mitochondrial membrane that plays a crucial role in mitochondrial fusion. Mfn1 mediates tethering and fusion of the outer mitochondrial membranes, helping to maintain the integrity and function of the mitochondrial network.
    "Mfn2", # Mitofusin 2: Another GTPase on the outer mitochondrial membrane that is involved in mitochondrial fusion. Mfn2 works in concert with Mfn1 to mediate the fusion of outer mitochondrial membranes, maintaining proper mitochondrial function and dynamics.
    "Opa1", # Optic atrophy 1: A dynamin-related GTPase located in the inner mitochondrial membrane that is crucial for the regulation of mitochondrial fusion. Opa1 mediates the fusion of inner mitochondrial membranes, contributing to the maintenance of mitochondrial morphology, cristae structure, and overall mitochondrial function. Additionally, Opa1 plays a role in the regulation of apoptosis and is involved in the proper functioning of the electron transport chain.
    "Dnm1", # Dynamin 1: A GTPase involved in the regulation of membrane remodeling and trafficking, including vesicle formation and endocytosis. It has also been implicated in the regulation of mitochondrial morphology and function.
    "Drp1", # Dynamin-related protein 1: Another GTPase, which plays a central role in mitochondrial fission. It is recruited to the outer mitochondrial membrane during fission and is crucial for the proper division and maintenance of mitochondrial networks.
    "Dnm1l", # Dynamin 1-like, involved in mitochondrial fission
    "Fis1", # Fission 1, involved in mitochondrial fission
    "Mff", # Mitochondrial fission factor, involved in mitochondrial fission
    "Mid49", # Mitochondrial dynamics protein of 49 kDa, involved in mitochondrial fission
    "Mid51", # Mitochondrial dynamics protein of 51 kDa, involved in mitochondrial fission
    "Ppargc1a", # Peroxisome proliferator-activated receptor gamma coactivator 1-alpha, involved in mitochondrial biogenesis
    "Ppargc1b", # Peroxisome proliferator-activated receptor gamma coactivator 1-beta, involved in mitochondrial biogenesis
    "Nrf1", # Nuclear respiratory factor 1, involved in mitochondrial biogenesis
    "Nrf2", # Nuclear respiratory factor 2, involved in mitochondrial biogenesis
    "Tfam", # Mitochondrial transcription factor A, involved in mitochondrial biogenesis
    "Sod2", # Superoxide dismutase 2, involved in mitochondrial antioxidant defense
    "Prdx3", # Peroxiredoxin 3, involved in mitochondrial antioxidant defense
    "Slc25a1", # Solute carrier family 25 member 1, involved in mitochondrial metabolism
    "Slc25a4", # Solute carrier family 25 member 4, involved in mitochondrial metabolism
    "Slc25a5", # Solute carrier family 25 member 5, involved in mitochondrial metabolism
    "Uqcrc1", # Ubiquinol-cytochrome c reductase core protein 1, involved in mitochondrial electron transport chain
    "Uqcrc2", # Ubiquinol-cytochrome c reductase core protein 2, involved in mitochondrial electron transport chain
    "Cox4i1", # Cytochrome c oxidase subunit 4 isoform 1, involved in mitochondrial electron transport chain
    "Cox4i2", # Cytochrome c oxidase subunit 4 isoform 2, involved in mitochondrial electron transport chain
    "Atp5a1", # ATP synthase, H+ transporting, mitochondrial F1 complex, alpha subunit 1, involved in mitochondrial ATP production
    "Atp5b", # ATP synthase, H+ transporting, mitochondrial F1 complex, beta subunit, involved in mitochondrial ATP production
    "Bak1", # BCL2 antagonist/killer 1, involved in mitochondrial-mediated apoptosis
    "Bax", # BCL2 associated X, involved in mitochondrial-mediated apoptosis
    "Bcl2", # B-cell CLL/lymphoma 2, involved in mitochondrial-mediated apoptosis
    "Bcl2l1", # BCL2-like 1, involved in mitochondrial-mediated apoptosis
    "Bnip3", # BCL2/adenovirus E1B 19kDa interacting protein 3, involved in mitoptosis
    "Bnip3l", # BCL2/adenovirus E1B 19kDa interacting protein 3-like, involved in mitoptosis
    "Casp9", # Caspase 9, involved in mitochondrial-mediated apoptosis
    "Cybb", # Cytochrome b-245, beta polypeptide, involved in oxidative stress
    "Gpx1", # Glutathione peroxidase 1, involved in oxidative stress
    "Gpx4", # Glutathione peroxidase 4, involved in oxidative stress
    "Mapk14", # Mitogen-activated protein kinase 14, involved in ischemia and oxidative stress
    "Nfe2l2", # Nuclear factor, erythroid 2 like 2, involved in oxidative stress response
    "Nox4", # NADPH oxidase 4, involved in oxidative stress
    "Polg", # DNA polymerase gamma, involved in mitochondrial DNA replication
    "Polg2", # DNA polymerase gamma 2, involved in mitochondrial DNA replication
    "Pink1", # PTEN-induced putative kinase 1, involved in mitophagy
    "Park2", # Parkin RBR E3 ubiquitin protein ligase, involved in mitophagy
    "Aifm1", # Apoptosis-inducing factor, mitochondrion-associated, 1, involved in DNA damage and apoptosis
    "Ogg1", # 8-oxoguanine DNA glycosylase 1, involved in DNA damage repair
    "Mutyh", # MutY DNA glycosylase, involved in DNA damage repair
    "Sod1", # Superoxide dismutase 1, involved in oxidative stress
    "Sod3", # Superoxide dismutase 3, involved in oxidative stress
    "Ucp1", # Uncoupling protein 1, involved in thermogenesis
    "Ucp2", # Uncoupling protein 2, involved in thermogenesis
    "Ucp3", # Uncoupling protein 3, involved in thermogenesis
    "Ucp4", # Uncoupling protein 4, involved in thermogenesis
    "Ucp5", # Uncoupling protein 5, involved in thermogenesis
    "Fgf21", # Fibroblast growth factor 21, involved in metabolic regulation
    "Klb", # Klotho beta (beta-klotho), co-receptor for Fgf21
    # Additional genes related to brain metabolism switch and different catabolic/anabolic modes
    "Bdh1", # 3-hydroxybutyrate dehydrogenase 1, involved in ketone body metabolism
    "Bdh2", # 3-hydroxybutyrate dehydrogenase 2, involved in ketone body metabolism
    "Hmgcs2", # 3-hydroxy-3-methylglutaryl-CoA synthase 2, involved in ketone body synthesis
    "Slc16a1", # Solute carrier family 16 member 1 (monocarboxylate transporter 1), involved in ketone body transport
    "Slc16a7", # Solute carrier family 16 member 7 (monocarboxylate transporter 2), involved in ketone body transport
    "Slc2a1", # Solute carrier family 2 member 1 (glucose transporter type 1), involved in glucose transport
    "Slc2a3", # Solute carrier family 2 member 3 (glucose transporter type 3), involved in glucose transport
    "Gck", # Glucokinase, involved in glycolysis
    "Pfkp", # Phosphofructokinase, platelet, involved in glycolysis
    "Pkm", # Pyruvate kinase, muscle, involved in glycolysis
    "Pdha1", # Pyruvate dehydrogenase E1 alpha 1, involved in the conversion of pyruvate to acetyl-CoA
    "Pdha2", # Pyruvate dehydrogenase E1 alpha 2, involved in the conversion of pyruvate to acetyl-CoA
    "Pdk1", # Pyruvate dehydrogenase kinase 1, involved in the regulation of pyruvate dehydrogenase complex
    "Pdk2", # Pyruvate dehydrogenase kinase 2, involved in the regulation of pyruvate dehydrogenase complex
    "Pdk3", # Pyruvate dehydrogenase kinase 3, involved in the regulation of pyruvate dehydrogenase complex
    "Pdk4", # Pyruvate dehydrogenase kinase 4, involved in the regulation of pyruvate dehydrogenase complex
    "Acox1", # Acyl-CoA oxidase 1, involved in fatty acid oxidation
    "Cpt1a", # Carnitine palmitoyltransferase 1A, involved in fatty acid oxidation
    "Cpt1b", # Carnitine palmitoyltransferase 1B, involved in fatty acid oxidation
    "Cpt1c", # Carnitine palmitoyltransferase 1C, involved in fatty acid oxidation
    "Cpt2", # Carnitine palmitoyltransferase 2, involved in fatty acid oxidation
    "Acsl1", # Acyl-CoA synthetase long-chain family member 1, involved in fatty acid activation
    "Acsl3", # Acyl-CoA synthetase long-chain family member 3, involved in fatty acid activation
    "Acsl4", # Acyl-CoA synthetase long-chain family member 4, involved in fatty acid activation
    "Acsl5", # Acyl-CoA synthetase long-chain family member 5, involved in fatty acid activation
    "Acsl6", # Acyl-CoA synthetase long-chain family member 6, involved in fatty acid activation
    "Acly", # ATP citrate lyase, involved in fatty acid synthesis
    "Fasn", # Fatty acid synthase, involved in fatty acid synthesis
    "Scd1", # Stearoyl-CoA desaturase 1, involved in fatty acid desaturation
    "Scd2", # Stearoyl-CoA desaturase 2, involved in fatty acid desaturation
    "Acaa1", # Acetyl-CoA acyltransferase 1: Involved in the beta-oxidation of fatty acids in the mitochondria. Catalyzes the last step in the breakdown of fatty acids, converting 3-ketoacyl-CoA to acetyl-CoA.
    "Acaa2", # Acetyl-CoA acyltransferase 2: Also involved in beta-oxidation of fatty acids in the mitochondria. Similar to Acaa1, it catalyzes the last step of fatty acid breakdown but has different substrate specificity.
    "Aldh1a1", # Aldehyde dehydrogenase 1 family, member A1: An enzyme involved in the detoxification of aldehydes (generated from alcohol metabolism and lipid peroxidation) by converting them to carboxylic acids.
    "Eno1", # Enolase 1: A glycolytic enzyme that catalyzes the conversion of 2-phosphoglycerate to phosphoenolpyruvate in the glycolysis pathway, which generates ATP through the breakdown of glucose.
    "Hadhb", # Hydroxyacyl-CoA dehydrogenase/3-ketoacyl-CoA thiolase/enoyl-CoA hydratase, beta subunit: A mitochondrial enzyme that participates in fatty acid beta-oxidation, playing a role in multiple steps of the process.
    "Pygb" # Glycogen phosphorylase, brain form: An enzyme that breaks down glycogen into glucose-1-phosphate, which is then converted to glucose-6-phosphate for entry into glycolysis. This enzyme plays a critical role in energy production, particularly during periods of increased energy demand.
)

genes.embed <- c(
    "Abcd1",
    "Abcd2",
    "Abcd3",
    "Acaa1",
    "Acaa2",
    "Acox1",
    "Agrn",
    "Agt",
    "Alcam",
    "Aldh1a1",
    "Aldh1l1",
    "Aldoc",
    "Angpt1",
    "Apoe",
    "App",
    "Aqp4",
    "Arf1",
    "Bmp7",
    "Bsg",
    "Caf4",
    "Ccl25",
    "Ckb",
    "Cnr1",
    "Cnr2",
    "Col4a5",
    "Cst3",
    "Dagla",
    "Daglb",
    "Decr2",
    "Dnm1",
    "Drp1",
    "Ech1",
    "Efna5",
    "Egfr",
    "Enho",
    "Eno1",
    "Faah",
    "Fgf1",
    "Fgfr3",
    "Fis1",
    "Fos",
    "Fth1",
    "Ftl1",
    "Gfap",
    "Gja1",
    "Gli1",
    "Glul",
    "Gnai2",
    "Gnas",
    "H2-K1",
    "Hacd2",
    "Hadhb",
    "Hbegf",
    "Hepacam",
    "Hif1",
    "Htra1",
    "Igsf1",
    "Il18",
    "Il1rapl1",
    "Itgav",
    "Jam2",
    "Lama2",
    "Lamb2",
    "Lcat",
    "Lgi1",
    "Lgi4",
    "Lpcat3",
    "Lrpap1",
    "Lrrc4b",
    "Lxn",
    "Mdk",
    "Mdv1",
    "Mfn1",
    "Mfn2",
    "Mgll",
    "Mief1",
    "Napepld",
    "Ncam1",
    "Ncan",
    "Ndrg2",
    "Nfasc",
    "Nfia",
    "Nlgn3",
    "Nrxn1",
    "Nrxn2",
    "Ntn1",
    "Ntrk3",
    "Opa1",
    "Otp",
    "Pex1",
    "Pex10",
    "Pex12",
    "Pex13",
    "Pex14",
    "Pex16",
    "Pex2",
    "Pex26",
    "Pex3",
    "Pex6",
    "Pkm",
    "Pla2g7",
    "Plcb1",
    "Psap",
    "Ptn",
    "Pygb",
    "Rgma",
    "Rtn4",
    "S100a1",
    "S100a6",
    "S100b",
    "Scd2",
    "Sdc2",
    "Sema6a",
    "Sema6d",
    "Sgcd",
    "Sirpa",
    "Slc1a2",
    "Slc1a3",
    "Slc38a1",
    "Slc4a4",
    "Slc6a11",
    "Slc7a10",
    "Slit1",
    "Slit2",
    "Slitrk2",
    "Sorbs1",
    "Sox9",
    "Sparc",
    "Spon1",
    "Tafa1",
    "Timp3",
    "Tkt",
    "Trpv1",
    "Vcam1",
    "Vegfa"
)

cnbn <- c(
    "Cnr1",
    "Cnr2",
    "Gpr55",
    "Dagla",
    "Daglb",
    "Mgll",
    "Faah",
    "Napepld",
    "Trpv1",
    "Gde1",
    "Pparg"
)

genes.nature <- c(
    "Ndrg2", "Aqp4", "Gja1",
    "Slc17a9", "Slc17a7", "Slc17a6",
    "Slc16a1", "Slc16a7", "Slc16a3",
    "Panx1", "P2rx7", "Srr",
    "Dao", "Vamp2", "Gab1",
    "Slc18a1", "Slc18a2", "Slc18a3",
    "Slc17a5", "Osmr", "S100a6",
    "Ogn", "Itih5", "Rdh10",
    "Fst", "1500015O10Rik", "Rnf13",
    "A2m", "Rsph9", "Galnt16",
    "Rad23b", "Tgfbr2", "Ppp1r15a",
    "Mlc1", "Slc6a11", "Slc1a3", "S100b",
    "S100b", "Fabp7", "Gab1", "Igfbp5",
    "Hopx", "Igsf1", "Tgfb2", "2810459M11Rik",
    "Rnf13", "Itih5", "Slc1a2", "Cd59a",
    "Vim", "Slc7a10", "Fos"
)

genes.jj <- c(
    "Acsbg1", "Acsl3", "Actb", "Agt", "Aldh1l1", "Aldoc", "Apoe",
    "Aqp4", "Atp1a2", "Bmpr1b", "Cbs", "Ckb", "Clu", "Cnx43",
    "Cpe", "Cst3", "Cth", "Cyp4f14", "Dbi", "Dbx2", "Dctd",
    "Dio2", "Ednrb", "Fgfr3", "Gabbr1", "Gabbr2", "Gjb6", "Gli3",
    "Gm266", "Gpr37l1", "Grhl1", "Grm3", "Gs", "Gucy2c", "Hapln1",
    "Hes5", "Heyl", "Hgf", "Id2", "Itih3", "Lars2", "Lcat",
    "Ldhb", "Malat1", "Mc3r", "Mfge8", "Mlc1", "Mmd2", "Mt1",
    "Mt2", "Ndrg2", "Nfia", "Npas3", "Nrxn1", "Ntrk2", "Ntsr2",
    "Olx1", "Otx2", "Pax6", "Pbxip1", "Phkg1", "Pla2g3", "Pla2g7",
    "Plcd4", "Plce1", "Plp1", "Ppap2b", "Ppia", "Prodh", "Ptprz1",
    "Rfx4", "Rpl41", "S1pr1", "Scd2", "Serpinb1c", "Slc19a3", "Slc1a2",
    "Slc1a3", "Slc39a12", "Slc4a4", "Slc6a11", "Son", "Sox2", "Sox9",
    "Sparc", "Sparcl1", "Tlr3"
)

genes.anatomy.jj <- c(
    "Gfap", "Meg3", "Snhg11", "Ttc3", "Cst3", "Sparcl1",
    "Ndrg2", "Cnx43", "Agt", "S100b", "Aldh1l1", "Lxn",
    "Aqp4", "Slc1a2", "Fabp7", "Glul", "Olig2"
)

# public resources:
housekeeping_mouse <-
    read_lines(file = here(data_dir, "housekeeping_mouse.tsv"))
transcription_factors <-
    read_lines(file = here(data_dir, "mm_tfs.csv"))
astroenriched_mouse <-
    read_lines(file = here(data_dir, "astrocyte_enriched_genes_shared_by_mice.tsv"))

astroprogenitor_humans <-
    read_lines(file = here(data_dir, "top_astrocytes_progenitor_cells_genes_by_humans.tsv")) |>
    str_to_sentence()
astromature_humans <-
    read_lines(file = here(data_dir, "top_mature_astrocyte_genes_by_humans.tsv")) |>
    str_to_sentence()

# aggregate:
genes.manual <- unique(c(
    genes.embed, genes.nature,
    genes.jj, genes.anatomy.jj,
    astroenriched_mouse,
    astromature_humans,
    astroprogenitor_humans
))


gene_int <-
    c(
        npr, np, irs_genes, neurotrans, mcr_genes,
        genes.manual, pvn_genes, dmh_trh_g
    ) %>% unique()
