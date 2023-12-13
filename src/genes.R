npr <- c('Oxtr',  'Avpr1a', 'Avpr1b', 'Avpr2',
         'Tacr1', 'Mc1r', 'Mc3r', 'Mc4r', 'Oprl1',  'Tacr3',
         'Gpr83', 'Galr1', 'Npy1r', 'Npy2r', 'Npy4r', 'Npy4r2', 'Npy5r',
         'Sstr1', 'Sstr2', 'Sstr3', 'Mchr1',  'Oprd1', 'Oprk1', 'Oprm1',
         'Trhr', 'Hcrtr1', 'Hcrtr2', 'Qrfpr',   'Npffr1', 'Npffr2',
         'Prlhr',  'Ghr', 'Ghrhr',  'Grpr', 'Vipr1', 'Vipr2', 'Prokr2',
         'Nmur2',  'Nmur1', 'Nmbr',  'Kiss1r', 'Crhr1', 'Crhr2',
         'Cntfr', 'Cckar', 'Cckbr', 'Galr3')

np <- c('Adcyap1', 'Oxt', 'Avp', 'Tac1', 'Pomc', 'Pnoc', 'Tac2',
        'Nts', 'Gal', 'Agrp', 'Npy', 'Sst', 'Cartpt', 'Pmch',
        'Reln', 'Rxfp1', 'Penk', 'Pdyn', 'Trh', 'Hcrt', 'Qrfp',
        'Npw', 'Npvf', 'Ghrh', 'Grp', 'Vip', 'Nms', 'Nmu', 'Nmb',
        'Kiss1', 'Crh', 'Bdnf', 'Cntf', 'Cck')

irs_genes <- c('Alk', 'Insr', 'Ltk', 'Igf1r', 'Irs1',
               'Ptn', 'Mdk', 'Fam150a', 'Fam150b',
               'Mc4r', 'Lepr', 'Sim1', 'Lmo4',
               'Slc2a1', 'Slc2a3')

neurotrans <- c('Slc17a6', 'Slc17a7', 'Slc17a8', 'Slc1a1', 'Slc1a2', 'Slc1a6',
                'Gad1', 'Slc32a1', 'Slc6a1')
glut       <- c('Slc17a6', 'Slc17a7', 'Slc17a8', 'Slc1a1', 'Slc1a2', 'Slc1a6')
gaba       <- c('Gad1', 'Gad2', 'Slc32a1', 'Slc6a1')
dopam      <- c('Th', 'Slc6a3', 'Slc18a2', 'Ddc',  'Slc18a3')
ach        <- c('Chat', 'Slc18a3', 'Ache', 'Slc5a7')

gene_int <- c(npr, np, irs_genes, neurotrans) %>% unique()
