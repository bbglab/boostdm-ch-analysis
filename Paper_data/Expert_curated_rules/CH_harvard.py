import pandas as pd


def niroula_CH(data, gene):
    
    """
    Select CH mutations according to Harvard/CNIC definition.
    Based in Niroula et al. 2021
    !!! TP53 definition might be wrong in supplementary table
    """
    
    ## Modify database if necessary
    # Create column with AA position if it does not exist already
    if 'Prot_pos' not in data.columns:
        data['Prot_pos'] = [int(i[1:-1]) for i in data['aachange'].to_list()]    
    # Change column name if necessary
    if 'pos' not in data.columns:
        data = data.rename(columns={'POS': 'pos'})
    
    ## Select CH mutations
    if gene == 'ASXL1':
        # Frameshift/nonsense/splice-site in exons 11-12
        # exon 11 coordinates > 32433284 < 32433917 || exon 12 coordinates > 32434432 < 32439319
        CH_mutations = data[(data['pos'] >= 32433284) & (data['pos'] <= 32433917) |
                    (data['pos'] >= 32434432) & (data['pos'] <= 32439319)]
        CH_mutations = CH_mutations[(CH_mutations['csqn_type_nonsense'] == 1) |
                                    (CH_mutations['csqn_type_splicing'] == 1) ]
            
    elif gene == 'PPM1D':
        # Frameshift/stop gain/splice site AA421-605
        CH_mutations = data[(data['Prot_pos'] >= 421) & (data['Prot_pos'] <= 605)]
        CH_mutations = CH_mutations[(CH_mutations['csqn_type_nonsense'] == 1) |
                                    (CH_mutations['csqn_type_splicing'] == 1) ]
        
    elif gene == 'TET2':
        # Frameshift/nonsense/splice-site
        CH_mutations1 = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1)]
        # Nonsynonymous AA 1104-1481, AA 1843-2002
        CH_mutations2 = data[(data['csqn_type_missense'] == 1)]
        CH_mutations2 = CH_mutations2[((CH_mutations2['Prot_pos'] >= 1104) & (CH_mutations2['Prot_pos'] <= 1481)) |
                                      ((CH_mutations2['Prot_pos'] >= 1843) & (CH_mutations2['Prot_pos'] <= 2002))]
        CH_mutations = pd.concat([CH_mutations1, CH_mutations2])

    elif gene == 'DNMT3A':
        # Frameshift/stop gain/splice site
        CH_mutations1 = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1)]
        # Nonsynonymous AA290-390, AA547-911
        CH_mutations2 = data[(data['csqn_type_missense'] == 1)]
        CH_mutations2 = CH_mutations2[((CH_mutations2['Prot_pos'] >= 290) & (CH_mutations2['Prot_pos'] <= 390)) |
                                      ((CH_mutations2['Prot_pos'] >= 547) & (CH_mutations2['Prot_pos'] <= 911))]
        # AA414, AA494, AA497, AA508, AA527, AA529, AA531, AA532, AA537, AA543
        CH_mutations3 = data[(data['csqn_type_missense'] == 1)]
        CH_mutations3 = CH_mutations3[CH_mutations3['Prot_pos'].isin([414, 494, 497, 508, 527, 
                                                                      529, 531, 532, 537, 543])]
        
        # Concatenate all driver mutations
        CH_mutations = pd.concat([CH_mutations1, CH_mutations2])
        CH_mutations = pd.concat([CH_mutations, CH_mutations3])
    
    elif gene == 'CHEK2':
        # Frameshift, nonsense, splice site
        CH_mutations = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1) ]
    
    elif gene == 'GNAS':
        # Nonsynonymous
        CH_mutations = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1) |
                            (data['csqn_type_missense'] == 1)]
        # !!!! In paper AA201, AA227, AA844: corresponds to another transcript
        # !!!! Transformed to AA844, AA870, AA1017 (checked in Bick-CNIC, Uddin)        
        CH_mutations = CH_mutations[CH_mutations['Prot_pos'].isin([844, 870, 1017])]
    
    elif gene == 'IDH2':
        # Nonsynonymous
        CH_mutations = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1) |
                            (data['csqn_type_missense'] == 1)]
        # AA134-146, 164-180
        CH_mutations = CH_mutations[((CH_mutations['Prot_pos'] >= 134) & (CH_mutations['Prot_pos'] <= 146)) |
                                    ((CH_mutations['Prot_pos'] >= 164) & (CH_mutations['Prot_pos'] <= 180))]

    elif gene == 'SF3B1':
        # Nonsynonymous
        CH_mutations = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1) |
                            (data['csqn_type_missense'] == 1)]
        # AA550-800
        CH_mutations = CH_mutations[(CH_mutations['Prot_pos'] >= 550) & (CH_mutations['Prot_pos'] <= 800)]       

    elif gene == 'SRSF2':
        # Nonsynonymous
        CH_mutations = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1) |
                            (data['csqn_type_missense'] == 1)]
        # AA57, AA95, AA85-100
        CH_mutations = CH_mutations[(CH_mutations['Prot_pos'] == 57) | (CH_mutations['Prot_pos'] == 95) |
                                    ((CH_mutations['Prot_pos'] >= 85) & (CH_mutations['Prot_pos'] <= 100))]       

    elif gene == 'TP53':
        # Frameshift/stop gain/splice site, Nonsynonymous
        CH_mutations = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1) |
                            (data['csqn_type_missense'] == 1)]
        
    elif gene == 'U2AF1':
        # Nonsynonymous
        CH_mutations = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1) |
                            (data['csqn_type_missense'] == 1)]
        # AA22, AA28, AA32, AA34, AA154, AA156, AA157
        CH_mutations = CH_mutations[CH_mutations['Prot_pos'].isin([22, 28, 32, 34, 154, 156, 157])]       
   
    return CH_mutations


def bick_CH(data, gene):
    
    """
    Select CH mutations according to Harvard/CNIC definition.
    Based in Bick et al. 2020
    ! GNAS that comes from Fuster/CNIC: It's the same but with other protein position (other transcript)
    """
    ## Modify database if necessary
    # Create column with AA position if it does not exist already
    if 'Prot_pos' not in data.columns:
        data['Prot_pos'] = [int(i[1:-1]) for i in data['aachange'].to_list()]  
    # Change column name if necessary
    if 'pos' not in data.columns:
        data = data.rename(columns={'POS': 'pos'})
    
    ## Select CH mutations
    if gene == 'ASXL1':
        # Frameshift/nonsense/splice-site in exons 11-12
        # exon 11 coordinates > 32433284 < 32433917 || exon 12 coordinates > 32434432 < 32439319
        CH_mutations = data[(data['pos'] >= 32433284) & (data['pos'] <= 32433917) |
                    (data['pos'] >= 32434432) & (data['pos'] <= 32439319)]
        CH_mutations = CH_mutations[(CH_mutations['csqn_type_nonsense'] == 1) |
                                    (CH_mutations['csqn_type_splicing'] == 1) ]
            
    elif gene == 'PPM1D':
        # Frameshift/nonsense in exons 5-6
        CH_mutations = data[(data['pos'] > 60656599) & (data['pos'] < 60656841) |
                            (data['pos'] > 60662995) & (data['pos'] < 60666280)]
        CH_mutations = CH_mutations[(CH_mutations['csqn_type_nonsense'] == 1)]
        
    elif gene == 'TET2':
        # Frameshift/nonsense/splice-site
        CH_mutations1 = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1)]
        # Nonsynonymous AA 1104-1481, AA 1843-2002
        CH_mutations2 = data[(data['csqn_type_missense'] == 1)]
        CH_mutations2 = CH_mutations2[((CH_mutations2['Prot_pos'] >= 1104) & (CH_mutations2['Prot_pos'] <= 1481)) |
                                      ((CH_mutations2['Prot_pos'] >= 1843) & (CH_mutations2['Prot_pos'] <= 2002))]
        CH_mutations = pd.concat([CH_mutations1, CH_mutations2])

    elif gene == 'DNMT3A':
        # Frameshift/nonsense/splice-site
        CH_mutations1 = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1)]
        # F290I, F290C, R366P, R366H, R366G, A368T, A368V, R379H, R379C, I407T, I407N, I407S, F414L, F414S, F414C, A462V, K468R
        CH_mutations2 = data[data['aachange'].\
            isin(['F290I', 'F290C', 'V296M', 'P307S', 'P307R', 'R326H', 'R326L', 'R326C', 'R326S', 'G332R', 
                  'G332E', 'V339A', 'V339M', 'V339G', 'L344Q', 'L344P', 'R366P', 'R366H', 'R366G', 'A368T', 
                  'A368V', 'R379H', 'R379C', 'I407T', 'I407N', 'I407S', 'F414L', 'F414S', 'F414C', 'A462V', 
                  'K468R', 'C497G', 'C497Y', 'Q527H', 'Q527P', 'Y533C', 'S535F', 'C537G', 'C537R', 'G543A', 
                  'G543S', 'G543C', 'L547H', 'L547P', 'L547F', 'M548I', 'M548K', 'G550R', 'W581R', 'W581G', 
                  'W581C', 'R604Q', 'R604W', 'R635W', 'R635Q', 'S638F', 'G646V', 'G646E', 'L653W', 'L653F', 
                  'I655N', 'V657A', 'V657M', 'R659H', 'Y660C', 'V665G', 'V665L', 'M674V', 'R676W', 'R676Q', 
                  'G685R', 'G685E', 'G685A', 'D686Y', 'D686G', 'R688H', 'G699R', 'G699S', 'G699D', 'P700L', 
                  'P700S', 'P700R', 'P700Q', 'P700T', 'P700A', 'D702N', 'D702Y', 'V704M', 'V704G', 'I705F', 
                  'I705T', 'I705S', 'I705N', 'G707D', 'G707V', 'C710S', 'C710Y', 'S714C', 'V716D', 'V716F', 
                  'V716I', 'N717S', 'N717I', 'P718L', 'R720H', 'R720G', 'K721R', 'K721T', 'Y724C', 'R729Q', 
                  'R729W', 'R729G', 'F731C', 'F731L', 'F731Y', 'F731I', 'F732del', 'F732C', 'F732S', 'F732L', 
                  'E733G', 'E733A', 'F734L', 'F734C', 'Y735C', 'Y735N', 'Y735S', 'R736H', 'R736C', 'R736P', 
                  'L737H', 'L737V', 'L737F', 'L737R', 'A741V', 'P742P', 'P743R', 'P743L', 'R749C', 'R749L', 
                  'R749H', 'R749G', 'F751L', 'F751C', 'F752del', 'F752C', 'F752L', 'F752I', 'F752V', 'W753G', 
                  'W753C', 'W753R', 'L754P', 'L754R', 'L754H', 'F755S', 'F755I', 'F755L', 'M761I', 'M761V', 
                  'G762C', 'V763I', 'S770L', 'S770W', 'S770P', 'R771Q', 'F772I', 'F772V', 'L773R', 'L773V', 
                  'E774K', 'E774D', 'E774G', 'I780T', 'D781G', 'R792H', 'W795C', 'W795L', 'G796D', 'G796V', 
                  'N797Y', 'N797H', 'N797S', 'P799S', 'P799R', 'P799H', 'R803S', 'R803W', 'P804L', 'P804S', 
                  'K826R', 'S828N', 'K829R', 'T835M', 'N838D', 'K841Q', 'Q842E', 'P849L', 'D857N', 'W860R', 
                  'E863D', 'F868S', 'G869S', 'G869V', 'M880V', 'S881R', 'S881I', 'R882H', 'R882P', 'R882C', 
                  'R882G', 'A884P', 'A884V', 'Q886R', 'L889P', 'L889R', 'G890D', 'G890R', 'G890S', 'V895M', 
                  'P896L', 'V897G', 'V897D', 'R899L', 'R899H', 'R899C', 'L901R', 'L901H', 'P904L', 'F909C', 
                  'P904Q', 'A910P', 'C911R', 'C911Y'])]
        
        # Concatenate all driver mutations
        CH_mutations = pd.concat([CH_mutations1, CH_mutations2])
    
    elif gene == 'GNAS':
        # Fuster et al!!! It's the same but with other protein position (other transcript)
        # Original definition Bick: R201S,  R201C,  R201H,  R201L,  Q227K,  Q227R,  Q227L,  Q227H,  R374C
        CH_mutations = data[data['aachange']\
            .isin(['R844S',  'R844C',  'R844H',  'R844L',  'Q870K',  'Q870R',  'Q870L',  'Q870H',  'R1017C'])]
    
    elif gene == 'IDH2':
        # R140W,  R140Q,  R140L,  R140G,  R172W,  R172G,  R172K,  R172T,  R172M,  R172N,  R172S
        CH_mutations = data[data['aachange']\
            .isin(['R140W',  'R140Q',  'R140L',  'R140G',  'R172W',  'R172G', 
                   'R172K',  'R172T',  'R172M',  'R172N',  'R172S'])]
         

    elif gene == 'SF3B1':
        # G347V,  R387W,  R387Q,  E592K,  E622D,  Y623C,  R625L,  R625C, R625G,  H662Q,  H662D, T663I, K666N,  K666T,  K666E,  K666R,  K700E,  V701F,  A708T,  G740R,  G740E,  A744P,  D781G,  E783K,  R831Q,  L833F,  E862K,  R957Q
        CH_mutations = data[data['aachange']\
            .isin(['G347V',  'R387W',  'R387Q', 'E592K', 'E622D', 'Y623C', 'R625L', 'R625C',
                   'R625G',  'H662Q',  'H662D', 'T663I', 'K666N',  'K666T',  'K666E',  'K666R',  
                   'K700E',  'V701F',  'A708T',  'G740R',  'G740E',  'A744P',  'D781G',  'E783K', 
                   'R831Q',  'L833F',  'E862K',  'R957Q'])]       

    elif gene == 'SRSF2':
        # Y44H,  P95H,  P95L,  P95T,  P95R,  P95A,  P107H,  P95fs
        CH_mutations = data[data['aachange']\
            .isin(['Y44H',  'P95H',  'P95L',  'P95T',  'P95R',  'P95A',  'P107H',  'P95fs'])]    

    elif gene == 'TP53':
        # Frameshift/nonsense/splice-site
        CH_mutations1 = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1)]
        # Nonsynonymous AA xxx
        CH_mutations2 = data[data['aachange']\
            .isin(['S46F', 'G105C', 'G105R', 'G105D', 'G108S', 'G108C', 'R110L', 'R110C', 'T118A', 'T118R', 
                   'T118I', 'S127F', 'S127Y', 'L130V', 'L130F', 'K132Q', 'K132E', 'K132W', 'K132R', 'K132M', 
                   'K132N', 'F134V', 'F134L', 'F134S', 'C135W', 'C135S', 'C135F', 'C135G', 'C135Y', 'Q136K', 
                   'Q136E', 'Q136P', 'Q136R', 'Q136L', 'Q136H', 'A138P', 'A138V', 'A138A', 'A138T', 'T140I', 
                   'C141R', 'C141G', 'C141A', 'C141Y', 'C141S', 'C141F', 'C141W', 'V143M', 'V143A', 'V143E', 
                   'L145Q', 'W146C', 'W146L', 'L145R', 'V147G', 'P151T', 'P151A', 'P151S', 'P151H', 'P151R', 
                   'P152S', 'P152R', 'P152L', 'T155P', 'T155A', 'V157F', 'R158H', 'R158L', 'A159V', 'A159P', 
                   'A159S', 'A159D', 'A161T', 'A161D', 'Y163N', 'Y163H', 'Y163D', 'Y163S', 'Y163C', 'K164E', 
                   'K164M', 'K164N', 'K164P', 'H168Y', 'H168P', 'H168R', 'H168L', 'H168Q', 'M169I', 'M169T', 
                   'M169V', 'E171K', 'E171Q', 'E171G', 'E171A', 'E171V', 'E171D', 'V172D', 'V173M', 'V173L', 
                   'V173G', 'R174W', 'R175G', 'R175C', 'R175H', 'C176R', 'C176G', 'C176Y', 'C176F', 'C176S', 
                   'P177R', 'P177R', 'P177L', 'H178D', 'H178P', 'H178Q', 'H179Y', 'H179R', 'H179Q', 'R181C', 
                   'R181Y', 'D186G', 'G187S', 'P190L', 'P190T', 'H193N', 'H193P', 'H193L', 'H193R', 'L194F', 
                   'L194R', 'I195F', 'I195N', 'I195T', 'R196P', 'V197L', 'G199V', 'Y205N', 'Y205C', 'Y205H', 
                   'D208V', 'R213Q', 'R213P', 'R213L', 'R213Q', 'H214D', 'H214R', 'S215G', 'S215I', 'S215R', 
                   'V216M', 'V217G', 'Y220N', 'Y220H', 'Y220S', 'Y220C', 'E224D', 'I232F', 'I232N', 'I232T', 
                   'I232S', 'Y234N', 'Y234H', 'Y234S', 'Y234C', 'Y236N', 'Y236H', 'Y236C', 'M237V', 'M237K', 
                   'M237I', 'C238R', 'C238G', 'C238Y', 'C238W', 'N239T', 'N239S', 'S241Y', 'S241C', 'S241F', 
                   'C242G', 'C242Y', 'C242S', 'C242F', 'G244S', 'G244C', 'G244D', 'G245S', 'G245R', 'G245C', 
                   'G245D', 'G245A', 'G245V', 'G245S', 'M246V', 'M246K', 'M246R', 'M246I', 'N247I', 'R248W', 
                   'R248G', 'R248Q', 'R249G', 'R249W', 'R249T', 'R249M', 'P250L', 'I251N', 'L252P', 'I254S', 
                   'I255F', 'I255N', 'I255S', 'L257Q', 'L257P', 'E258K', 'E258Q', 'D259Y', 'S261T', 'G262D', 
                   'G262V', 'L265P', 'G266R', 'G266E', 'G266V', 'R267W', 'R267Q', 'R267P', 'E271K', 'V272M', 
                   'V272L', 'R273S', 'R273G', 'R273C', 'R273H', 'R273P', 'R273L', 'V274F', 'V274D', 'V274A', 
                   'V274G', 'V274L', 'C275Y', 'C275S', 'C275F', 'A276P', 'C277F', 'C277Y', 'P278T', 'P278A', 
                   'P278S', 'P278H', 'P278R', 'P278L', 'G279E', 'R280G', 'R280K', 'R280T', 'R280I', 'R280S', 
                   'D281N', 'D281H', 'D281Y', 'D281G', 'D281E', 'R282G', 'R282W', 'R282Q', 'R282P', 'E285K', 
                   'E285V', 'E286G', 'E286V', 'E286K', 'K320N', 'L330R', 'G334V', 'R337C', 'R337L', 'A347T', 
                   'L348F', 'T377P'])]
        CH_mutations = pd.concat([CH_mutations1, CH_mutations2])
                             
    elif gene == 'U2AF1':
        # D14G,  S34F,  S34Y,  R35L,  R156H,  R156Q,  Q157R,  Q157P
        CH_mutations = data[data['aachange']\
            .isin(['D14G',  'S34F',  'S34Y',  'R35L',  'R156H',  'R156Q',  'Q157R',  'Q157P'])]      
   
    return CH_mutations


def WHO_CH(data, gene):
    
    """
    CHECKED 03/03/2023 by joan enric
    
    Select CH mutations according to OMS definition.
    Based in table sent by Vassiliou
    """
    
    ## Modify database if necessary
    # Create column with AA position if it does not exist already
    if 'Prot_pos' not in data.columns:
        data['Prot_pos'] = [int(i[1:-1]) for i in data['aachange'].to_list()]    
    # Change column name if necessary
    if 'pos' not in data.columns:
        data = data.rename(columns={'POS': 'pos'})
    
    ## Select CH mutations
    if gene == 'ASXL1':
        # Frameshift/nonsense/splice-site in exons 11-12 (NM_015338)
        # exon 11 coordinates > 32433284 < 32433917 || exon 12 coordinates > 32434432 < 32439319
        CH_mutations = data[(data['pos'] >= 32433284) & (data['pos'] <= 32433917) |
                    (data['pos'] >= 32434432) & (data['pos'] <= 32439319)]
        CH_mutations = CH_mutations[(CH_mutations['csqn_type_nonsense'] == 1) |
                                    (CH_mutations['csqn_type_splicing'] == 1) ]
            
    elif gene == 'PPM1D':
        # Frameshift/nonsense/splice-site in exon 5/6 (NM_003620)
        CH_mutations = data[(data['pos'] > 60656599) & (data['pos'] < 60656841) |
                            (data['pos'] > 60662995) & (data['pos'] < 60666280)]
        CH_mutations = CH_mutations[(CH_mutations['csqn_type_nonsense'] == 1) |
                                    (CH_mutations['csqn_type_splicing'] == 1) ]
        
    elif gene == 'TET2':
        # Frameshift/nonsense/splice-site
        CH_mutations1 = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1)]
        # Missense in aa range: p.1104-1481 and p.1843-2002 (NM_001127208)
        CH_mutations2 = data[(data['csqn_type_missense'] == 1)]
        CH_mutations2 = CH_mutations2[((CH_mutations2['Prot_pos'] >= 1104) & (CH_mutations2['Prot_pos'] <= 1481)) |
                                      ((CH_mutations2['Prot_pos'] >= 1843) & (CH_mutations2['Prot_pos'] <= 2002))]
        CH_mutations = pd.concat([CH_mutations1, CH_mutations2])

    elif gene == 'DNMT3A':
        # Frameshift/stop gain/splice site
        CH_mutations1 = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1)]
        # Missense in aa range: p.292-350, p.482-614 and p.634-912 (NM_022552)
        CH_mutations2 = data[(data['csqn_type_missense'] == 1)]
        CH_mutations2 = CH_mutations2[((CH_mutations2['Prot_pos'] >= 292) & (CH_mutations2['Prot_pos'] <= 350)) |
                                      ((CH_mutations2['Prot_pos'] >= 482) & (CH_mutations2['Prot_pos'] <= 614)) |
                                      ((CH_mutations2['Prot_pos'] >= 634) & (CH_mutations2['Prot_pos'] <= 912))]
        # Concatenate all driver mutations
        CH_mutations = pd.concat([CH_mutations1, CH_mutations2])
    
    elif gene == 'GNAS':
        # Missense at p.R201 (NM_016592)
        # !!!! In paper AA201 corresponds to another transcript transformed to AA844 (checked in Bick-CNIC, Uddin) 
        CH_mutations = data[(data['csqn_type_missense'] == 1)]       
        CH_mutations = CH_mutations[CH_mutations['Prot_pos'].isin([844])]
    
    elif gene == 'IDH2':
        # Missense at p.R140 or p.R172 (NM_002168)
        CH_mutations = data[(data['csqn_type_missense'] == 1)]
        CH_mutations = CH_mutations[CH_mutations['Prot_pos'].isin([140, 172])]

    elif gene == 'SF3B1':
        # Missense in terminal HEAT domains (p.529-1201) (NM_012433)
        CH_mutations = data[(data['csqn_type_missense'] == 1)]       
        CH_mutations = CH_mutations[(CH_mutations['Prot_pos'] >= 529) & (CH_mutations['Prot_pos'] <= 1201)]       

    elif gene == 'SRSF2':
        # Missense/in-frame deletion involving P95 (NM_003016)
        CH_mutations = data[(data['csqn_type_missense'] == 1)]       
        CH_mutations = CH_mutations[CH_mutations['Prot_pos'].isin([95])]       

    elif gene == 'TP53':
        # Frameshift/nonsense/splice-site
        CH_mutations1 = data[(data['csqn_type_nonsense'] == 1) |
                            (data['csqn_type_splicing'] == 1)]
        # Missense in aa range: p.72, p.95-288 and p.337 (NM_001126112)
        CH_mutations2 = data[(data['csqn_type_missense'] == 1)]
        CH_mutations2 = CH_mutations2[(CH_mutations2['Prot_pos'] == 72) |
                                      ((CH_mutations2['Prot_pos'] >= 95) & (CH_mutations2['Prot_pos'] <= 288)) |
                                      (CH_mutations2['Prot_pos'] == 337)]
        CH_mutations = pd.concat([CH_mutations1, CH_mutations2])
        
    elif gene == 'U2AF1':
        # Missense at p.S34 / p.R156 / p.Q157 (NM_006758)
        CH_mutations = data[(data['csqn_type_missense'] == 1)]
        CH_mutations = CH_mutations[CH_mutations['Prot_pos'].isin([34, 156, 157])]       
   
    return CH_mutations
