# Retrieve scaffold proteins in organelles of interest
# Based on DrLLPS, LLPSDB and STRINGdb

getScaffolds <- function(organelle=NA, trsh=600, ensembl_host="https://jul2023.archive.ensembl.org") {
 	# get list of LLPS scaffolds from DrLLPS and LLPSDB
  llpsData <- getLLPSdata(host=ensembl_host)

	if(organelle == "Heterochromatin focus") {
		# get list of Cbx5/Mecp2 interactors
	  cbx5_interactors <- getInteractors("cbx5", scaffolds=llpsData, host=ensembl_host, score=trsh, toAdd=c("Hmga2", "Kmt5c", "Ube2i"), toRemove=c("Smarca4", "H1f4", "Sgo1"))
		mecp2_interactors <- getInteractors("mecp2", scaffolds=llpsData, host=ensembl_host, score=trsh, toRemove=c("Smarca4", "Aurkb", "Shank3", "Rbfox1", "Tbl1xr1"))
		## Histone H1 removed because we consider it as part of the chromatin portion of the condensate
		## Smarca4 removed because its specific link with heterochromatin is unclear
		## Tbl1xr1 removed because it is only active upon hormone stimulation
		## Aurkb and Sgo1 removed because they function in mitosis
		## Shank3 and Rbfox1 removed because they are linked to  neuronal functions
		## Hmga2 and Kmt5c added because they localize to heterochromatin in various cell types
		##
		## Found interactors: Cbx1, Suv39h1, Ube2i, Trim28
		
		# get list of scaffolds that interact with Cbx5/Mecp2
		cbx5_scaffolds <- intersect(llpsData, cbx5_interactors)
		mecp2_scaffolds <- intersect(llpsData, mecp2_interactors)

		# return list of scaffolds in heterochromatin foci
		res <- rbind(cbx5_scaffolds, mecp2_scaffolds)
		res <- res[!duplicated(res[,'uniprotswissprot']),]
	}

	if(organelle == "Polycomb body") {
	  # get list of Cbx2 interactors
	  cbx2_interactors <- getInteractors("cbx2", scaffolds=llpsData, host=ensembl_host, score=trsh)
	  ## Found interactors: Phc1, Ube2i, Sumo3
	  
		# get list of scaffolds that interact with Cbx2
		cbx2_scaffolds <- intersect(llpsData, cbx2_interactors) 
		
		# return list of scaffolds in Polycomb bodies
		res <- cbx2_scaffolds
		res <- res[!duplicated(res[,'uniprotswissprot']),]
	}
		
	if(organelle == "Txn condensate") {
	  # get list of Med1/Brd4 interactors
	  med1_interactors <- getInteractors("med1", scaffolds=llpsData, host=ensembl_host, score=trsh, toAdd=c("Ccnt1"), toRemove=c("Ar", "Rxrg", "Tbl1xr1", "Srsf2", "Rara", "Esr1"))
	  brd4_interactors <- getInteractors("brd4", scaffolds=llpsData, host=ensembl_host, score=trsh, toRemove=c("Cdk1", "Ar", "Srsf2", "Vhl"))
	  ## Cdk1 removed because it is not expected to be a core component of transcriptional condensates
	  ## Ar, Rxrg and Tbl1xr1 removed because they are only active upon hormone stimulation
	  ## Srsf2 removed because it is a nuclear speckle marker, Spector and Lamond, 2011, Cold Spring Harb Perspect Biol 3: a000646, PMID 20926517, https://doi.org/10.1101/cshperspect.a000646
	  ## Rara removed because it localizes in the cytoplasm of mESCs, Sart et al, 2014, Tissue Eng Part A 20: 54–66, PMID 23848515, https://doi.org/10.1089/ten.tea.2012.0690
	  ## Esr1 removed because it only forms foci upon stimulation with estradiol, Nair et al, 2020, Nat Struct Mol Biol 26: 193–203, PMID 30833784, https://doi.org/10.1038/s41594-019-0190-5
	  ##
	  ## Found interactors: Cdk9, Pou5f1, Brd3, Smarca4, Polr2a
	  
		# get list of scaffolds that interact with Med1/Brd4
		med1_scaffolds <- intersect(llpsData, med1_interactors)
		brd4_scaffolds <- intersect(llpsData, brd4_interactors)

		# return list of scaffolds in transcriptional condensates
		res <- rbind(med1_scaffolds, brd4_scaffolds)
		res <- res[!duplicated(res[,'uniprotswissprot']),]
	}

	if(organelle == "Nucleolus, GC") {
	  # get list of Npm1 interactors
		npm1_interactors <- getInteractors("npm1", scaffolds=llpsData, host=ensembl_host, score=trsh, toAdd=c("Eloa", "Lin28a", "Rpl23a", "Gnl2"), toRemove=c("Fbl", "Srsf2"))
		## Srsf2 removed because it is a nuclear speckle marker, Spector and Lamond, 2011, Cold Spring Harb Perspect Biol 3: a000646, PMID 20926517, https://doi.org/10.1101/cshperspect.a000646
		## Rpl23a and Gnl2 added, Mitrea et al, 2016, eLife 5: e13571, PMID 26836305, https://doi.org/10.7554/eLife.13571
		## Lin28 added, Sun, Yu, Zhao et al, 2022, Protein Cell 13: 490-512, PMID 34331666, https://doi.org/10.1007/s13238-021-00864-5
		## Eloa added, Ardehali et al, 2021, JBC 296: 100202, PMID 33334895, https://doi.org/10.1074/jbc.RA120.015877

		# remove proteins that are not detected in all three replicates of nucleolar extract
		# Szerlong et al, 2015, J Mol Biol 427: 2056-2071, PMID 25584861, https://doi.org/10.1016/j.jmb.2015.01.001, supplement-3.xlsx
		npm1_interactors_filtered <- intersectNucleolarExtract(npm1_interactors, toAdd=c("Rpl23a", "Gnl2", "Eloa", "Rpl5"))
		
		# get list of scaffolds that interact with Npm1
		npm1_scaffolds <- intersect(llpsData, npm1_interactors_filtered) 
		
		# return list of scaffolds in the GC of the nucleolus
		res <- npm1_scaffolds
		res <- res[!duplicated(res[,'uniprotswissprot']),]
	}

	if(organelle == "Nucleolus, DFC") {
	  # get list of Fbl interactors
	  fbl_interactors <- getInteractors("fbl", scaffolds=llpsData, host=ensembl_host, score=trsh, toAdd=c("Eloa", "Lin28a"), toRemove=c("Ncl", "Npm1", "Srsf2", "Surf6", "Eif4a1"))
	  ## Lin28 added, Sun, Yu, Zhao et al, 2022, Protein Cell 13: 490-512, PMID 34331666, https://doi.org/10.1007/s13238-021-00864-5
	  ## Eloa added, Ardehali et al, 2021, JBC 296: 100202, PMID 33334895, https://doi.org/10.1074/jbc.RA120.015877
	  ## Srsf2 removed because it is a nuclear speckle marker, Spector and Lamond, 2011, Cold Spring Harb Perspect Biol 3: a000646, PMID 20926517, https://doi.org/10.1101/cshperspect.a000646
	  ## Eif4a1 removed because it is mostly cytoplasmic
	  
	  # remove proteins that are not detected in all three replicates of nucleolar extract
	  # Szerlong et al, 2015, J Mol Biol 427: 2056-2071, PMID 25584861, https://doi.org/10.1016/j.jmb.2015.01.001, supplement-3.xlsx
	  fbl_interactors_filtered <- intersectNucleolarExtract(fbl_interactors, toAdd=c("Polr1b", "Eloa"))
	  
	  # get list of scaffolds that interact with Npm1
	  fbl_scaffolds <- intersect(llpsData, fbl_interactors_filtered) 
	  
	  # return list of scaffolds in the GC of the nucleolus
	  res <- fbl_scaffolds
	  res <- res[!duplicated(res[,'uniprotswissprot']),]
	}
  
  # return results
  return(res)
}


## ###################
## accessory functions
## ###################

# query LLPS databases to obtain a list of scaffold proteins underoing LLPS
getLLPSdata <- function(host) {
  # load libraries
  library("biomaRt")
  library("RCurl")
  library("readr")
  library("xlsx")
  
  # get scaffolds from DrLLPS database (Ning et al, 2020, Nucleic Acids Res 48: D288-D295, PMID: 31691822)
  drllps_data <- RCurl::getURLContent("http://llps.biocuckoo.cn/download/LLPS.txt")
  drllps_table <- read.table(textConnection(drllps_data), sep = "\t", header = T, stringsAsFactors = F)
  drllps_scaffolds <- drllps_table[drllps_table$LLPS.Type=="Scaffold",]
  drllps_scaffolds <- drllps_scaffolds[drllps_scaffolds$Species %in% c("Homo sapiens", "Mus musculus", "Rattus norvegicus"),]
  
  # add scaffolds that have been identified after the last release of DrLLPS
  gata3 <- data.frame("-", "P23771", "GATA3", "ENSG00000107485.15", "Homo sapiens", "Droplet", "Scaffold", "30833784", "MEVTADQPRWVSHHHPAVLNGQHPDTHHPGLSHSYMDAAQYPLPEEVDVLFNIDGQGNHVPPYYGNSVRATVQRYPPTHHGSQVCRPPLLHGSLPWLDGGKALGSHHTASPWNLSPFSKTSIHHGSPGPLSVYPPASSSSLSGGHASPHLFTFPPTPPKDVSPDPSLSTPGSAGSARQDEKECLKYQVPLPDSMKLESSHSRGSMTALGGASSSTHHPITTYPPYVPEYSSGLFPPSSLLGGSPTGFGCKSRPKARSSTEGRECVNCGATSTPLWRRDGTGHYLCNACGLYHKMNGQNRPLIKPKRRLSAARRAGTSCANCQTTTTTLWRRNANGDPVCNACGLYYKLHNINRPLTMKKEGIQTRNRKMSSKSKKCKKVHDSLEDFPKNSSFNPAALSRHMSSLSHISPFSHSSHMLTTPTPMHPPSSLSFGPHHPSSMVTAMG")
  names(gata3) <- names(drllps_scaffolds)
  
  fbl <- data.frame("-", "P35550", "Fbl", "ENSMUSP00000037613", "Mus musculus", "Droplet", "Scaffold","", "MKPGFSPRGGGFGGRGGFGDRGGRGGGRGGRGGFGGGRGGFGGGGRGRGGGGGGFRGRGGGGGRGGGFQSGGNRGRGGGRGGKRGNQSGKNVMVEPHRHEGVFICRGKEDALVTKNLVPGESVYGEKRVSISEGDDKIEYRAWNPFRSKLAAAILGGVDQIHIKPGAKVLYLGAASGTTVSHVSDIVGPDGLVYAVEFSHRSGRDLINLAKKRTNIIPVIEDARHPHKYRMLIAMVDVIFADVAQPDQTRIVALNAHTFLRNGGHFVISIKANCIDSTASAEAVFASEVKKMQQENMKPQEQLTLEPYERDHAVVVGVYRPPPKVKN")
  names(fbl) <- names(drllps_scaffolds)
  
  mecp2 <- data.frame("-", "P51608", "MECP2", "ENSG00000169057.19", "Homo sapiens", "Droplet", "Scaffold", "32698189", "MAAAAAAAPSGGGGGGEEERLEEKSEDQDLQGLKDKPLKFKKVKKDKKEEKEGKHEPVQPSAHHSAEPAEAGKAETSEGSGSAPAVPEASASPKQRRSIIRDRGPMYDDPTLPEGWTRKLKQRKSGRSAGKYDVYLINPQGKAFRSKVELIAYFEKVGDTSLDPNDFDFTVTGRGSPSRREQKPPKKPKSPKAPGTGRGRGRPKGSGTTRPKAATSEGVQVKRVLEKSPGKLLVKMPFQTSPGGKAEGGGATTSTQVMVIKRPGRKRKAEADPQAIPKKRGRKPGSVVAAAAAEAKKKAVKESSIRSVQETVLPIKKRKTRETVSIEVKEVVKPLLVSTLGEKSGKGLKTCKSPGRKSKESSPKGRSSSASSPPKKEHHHHHHHSESPKAPVPLLPPLPPPPPEPESSEDPTSPPEPQDLSSSVCKEEKMPRGGSLESDGCPKEPAKTQPAVATAATAAEKYKHRGEGERKDIVSSSMPRPNREEPVDSRTPVTERVS")
  names(mecp2) <- names(drllps_scaffolds)
  
  kmt5c <- data.frame("-", "Q86Y97", "KMT5C", "ENSG00000133247.13", "Homo sapiens", "Droplet", "Scaffold", "33326747", "MGPDRVTARELCENDDLATSLVLDPYLGFRTHKMNVSPVPPLRRQQHLRSALETFLRQRDLEAAYRALTLGGWTARYFQSRGPRQEAALKTHVYRYLRAFLPESGFTILPCTRYSMETNGAKIVSTRAWKKNEKLELLVGCIAELREADEGLLRAGENDFSIMYSTRKRSAQLWLGPAAFINHDCKPNCKFVPADGNAACVKVLRDIEPGDEVTCFYGEGFFGEKNEHCECHTCERKGEGAFRTRPREPALPPRPLDKYQLRETKRRLQQGLDSGSRQGLLGPRACVHPSPLRRDPFCAACQPLRLPACSARPDTSPLWLQWLPQPQPRVRPRKRRRPRPRRAPVLSTHHAARVSLHRWGGCGPHCRLRGEALVALGQPPHARWAPQQDWHWARRYGLPYVVRVDLRRLAPAPPATPAPAGTPGPILIPKQALAFAPFSPPKRLRLVVSHGSIDLDVGGEEL")
  names(kmt5c) <- names(drllps_scaffolds)
  
  cbx4 <- data.frame("-", "O00257", "CBX4", "ENSG00000141582.14", "Homo sapiens", "Droplet", "Scaffold", "31965994", "MELPAVGEHVFAVESIEKKRIRKGRVEYLVKWRGWSPKYNTWEPEENILDPRLLIAFQNRERQEQLMGYRKRGPKPKPLVVQVPTFARRSNVLTGLQDSSTDNRAKLDLGAQGKGQGHQYELNSKKHHQYQPHSKERAGKPPPPGKSGKYYYQLNSKKHHPYQPDPKMYDLQYQGGHKEAPSPTCPDLGAKSHPPDKWAQGAGAKGYLGAVKPLAGAAGAPGKGSEKGPPNGMMPAPKEAVTGNGIGGKMKIVKNKNKNGRIVIVMSKYMENGMQAVKIKSGEVAEGEARSPSHKKRAADERHPPADRTFKKAAGAEEKKVEAPPKRREEEVSGVSDPQPQDAGSRKLSPTKEAFGEQPLQLTTKPDLLAWDPARNTHPPSHHPHPHPHHHHHHHHHHHHAVGLNLSHVRKRCLSETHGEREPCKKRLTARSISTPTCLGGSPAAERPADLPPAAALPQPEVILLDSDLDEPIDLRCVKTRSEAGEPPSSLQVKPETPASAAVAVAAAAAPTTTAEKPPAEAQDEPAESLSEFKPFFGNIIITDVTANCLTVTFKEYVTV")
  names(cbx4) <- names(drllps_scaffolds)
  
  phc1 <- data.frame("-", "Q64028", "Phc1", "ENSMUSG00000040669", "Mus musculus", "", "","", "METESEQNSSSTNGSSSSGASSRPQIAQMSLYERQAVQALQALQRQPNAAQYFHQFMLQQQLSNAQLHSLAAVQQATIAASRQASSPNSSTAQQQTATTQASMNLATTSAAQLISRSQSVSSPSATTLTQSVLLGNTTSPPLNQSQAQMYLRPQLGNLLQVNRTLGRNVPLASQLILMPNGAVAAVQQEVPPAQSPGVHADADQVQNLAVRNQQASAQGPQMPGSTQKAIPPGASPVSGLSQTSSQALAVAQASSGASGQSLNLSQAGGGSGNSLPGSMGPGGGGQAPGGLGQLPSSGLTGGSCPRKGTGVVQPLPAAQTVTVSQGSQTEAESAAAKKAEADGSGQQSVGMNLTRTATPAPSQTLISSATYTQIQPHSLIQQQQQIHLQQKQVVIQQQIAIHHQQQFQHRQSQLLHTATHLQLAQQQQQQQQQQQQQQQQQQQQQQGTTLTAPQPPQVPPTQQVPPSQSQQQAQTLVVQPMLQSSPLTLPPEPTSKPPIPIQSKPPVAPIKPPQLGAAKMSATQQPPPHIPVQVVGTRQPGSAQAQALGLAQLAAAVPTPRGITGAVQPGQAHLASSPPSSQAAPGALQECPPALAAGMTLAPVQGTAHVVKGGPTASSPVVAQVPAAFYMQSVHLPGKAQTLAVKRKAESEEERDDLSALASVLPTKASPAAESPKVIEEKNSLGEKAEPVASLNANPPNSDLVALAPTPSAPPPTLALVSRQMGDSKPPQAIVKPQILTHIIEGFVIQEGAEPFPVGCSQFLKETKKPLQAGLPTGLNESQPSGPLGGDSPSVELEKKANLLKCEYCGKYAPAEQFRGSKRFCSMTCAKRYNVSCSHQFRLKRKKMKEFQEASYARVRRRGPRRSSSDIARAKIQGKRHRGQEDSSRGSDNSSYDEALSPTSPGPLSVRAGHGERDLGNTITTPSTPELQGINPVFLSSNPSQWSVEEVYEFIASLQGCQEIAEEFRSQEIDGQALLLLKEEHLMSAMNIKLGPALKICAKINVLKET" )
  names(phc1) <- names(drllps_scaffolds)
  
  hmga1 <- data.frame("-", "P17096", "HMGA1", "ENSG00000137309", "Homo sapiens", "Droplet", "Scaffold", "", "MSESSSKSSQPLASKQEKDGTEKRGRGRPRKQPPVSPGTALVGSQKEPSEVPTPKRPRGRPKGSKNKGAAKTRKTTTTPGRKPRGRPKKLEKEEEEGISQESSEEEQ")
  names(hmga1) <- names(drllps_scaffolds)
  
  hmga2 <- data.frame("-", "P52927", "Hmga2", "ENSMUSG00000056758", "Mus musculus", "Droplet", "Scaffold", "", "MSARGEGAGQPSTSAQGQPAAPVPQKRGRGRPRKQQQEPTCEPSPKRPRGRPKGSKNKSPSKAAQKKAETIGEKRPRGRPRKWPQQVVQKKPAQETEETSSQESAEED")
  names(hmga2) <- names(drllps_scaffolds)
  
  polr2a <- data.frame("-", "P08775", "Polr2a", "ENSMUSG00000005198", "Mus musculus", "Droplet", "Scaffold", "", "MHGGGPPSGDSACPLRTIKRVQFGVLSPDELKRMSVTEGGIKYPETTEGGRPKLGGLMDPRQGVIERTGRCQTCAGNMTECPGHFGHIELAKPVFHVGFLVKTMKVLRCVCFFCSKLLVDSNNPKIKDILAKSKGQPKKRLTHVYDLCKGKNICEGGEEMDNKFGVEQPEGDEDLTKEKGHGGCGRYQPRIRRSGLELYAEWKHVNEDSQEKKILLSPERVHEIFKRISDEECFVLGMEPRYARPEWMIVTVLPVPPLSVRPAVVMQGSARNQDDLTHKLADIVKINNQLRRNEQNGAAAHVIAEDVKLLQFHVATMVDNELPGLPRAMQKSGRPLKSLKQRLKGKEGRVRGNLMGKRVDFSARTVITPDPNLSIDQVGVPRSIAANMTFAEIVTPFNIDRLQELVRRGNSQYPGAKYIIRDNGDRIDLRFHPKPSDLHLQTGYKVERHMCDGDIVIFNRQPTLHKMSMMGHRVRILPWSTFRLNLSVTTPYNADFDGDEMNLHLPQSLETRAEIQELAMVPRMIVTPQSNRPVMGIVQDTLTAVRKFTKRDVFLERGEVMNLLMFLSTWDGKVPQPAILKPRPLWTGKQIFSLIIPGHINCIRTHSTHPDDEDSGPYKHISPGDTKVVVENGELIMGILCKKSLGTSAGSLVHISYLEMGHDITRLFYSNIQTVINNWLLIEGHTIGIGDSIADSKTYQDIQNTIKKAKQDVIEVIEKAHNNELEPTPGNTLRQTFENQVNRILNDARDKTGSSAQKSLSEYNNFKSMVVSGAKGSKINISQVIAVVGQQNVEGKRIPFGFKHRTLPHFIKDDYGPESRGFVENSYLAGLTPTEFFFHAMGGREGLIDTAVKTAETGYIQRRLIKSMESVMVKYDATVRNSINQVVQLRYGEDGLAGESVEFQNLATLKPSNKAFEKKFRFDYTNERALRRTLQEDLVKDVLSNAHIQNELEREFERMREDREVLRVIFPTGDSKVVLPCNLLRMIWNAQKIFHINPRLPSDLHPIKVVEGVKELSKKLVIVNGDDPLSRQAQENATLLFNIHLRSTLCSRRMAEEFRLSGEAFDWLLGEIESKFNQAIAHPGEMVGALAAQSLGEPATQMTLNTFHYAGVSAKNVTLGVPRLKELINISKKPKTPSLTVFLLGQSARDAERAKDILCRLEHTTLRKVTANTAIYYDPNPQSTVVAEDQEWVNVYYEMPDFDVARISPWLLRVELDRKHMTDRKLTMEQIAEKINAGFGDDLNCIFNDDNAEKLVLRIRIMNSDENKMQEEEEVVDKMDDDVFLRCIESNMLTDMTLQGIEQISKVYMHLPQTDNKKKIIITEDGEFKALQEWILETDGVSLMRVLSEKDVDPVRTTSNDIVEIFTVLGIEAVRKALERELYHVISFDGSYVNYRHLALLCDTMTCRGHLMAITRHGVNRQDTGPLMKCSFEETVDVLMEAAAHGESDPMKGVSENIMLGQLAPAGTGCFDLLLDAEKCKYGMEIPTNIPGLGAAGPTGMFFGSAPSPMGGISPAMTPWNQGATPAYGAWSPSVGSGMTPGAAGFSPSAASDASGFSPGYSPAWSPTPGSPGSPGPSSPYIPSPGGAMSPSYSPTSPAYEPRSPGGYTPQSPSYSPTSPSYSPTSPSYSPTSPNYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPSYSPTSPNYSPTSPNYTPTSPSYSPTSPSYSPTSPNYTPTSPNYSPTSPSYSPTSPSYSPTSPSYSPSSPRYTPQSPTYTPSSPSYSPSSPSYSPTSPKYTPTSPSYSPSSPEYTPASPKYSPTSPKYSPTSPKYSPTSPTYSPTTPKYSPTSPTYSPTSPVYTPTSPKYSPTSPTYSPTSPKYSPTSPTYSPTSPKGSTYSPTSPGYSPTSPTYSLTSPAISPDDSDEEN")
  names(polr2a) <- names(drllps_scaffolds)
  
  polr1b <- data.frame("-", "P70700", "Polr1b", "ENSMUSG00000027395", "Mus musculus", "Droplet", "Scaffold", "", "MDVDGRWRNLPSGPSLKHLTDPSYGIPPEQQKAALQDLTRAHVDSFNYAALEGLSHAVQAIPPFEFAFKDERISLTIVDAVISPPSVPKGTICKDLNVYPAECRGRKSTYRGRLTADISWAVNGVPKGIIKQFLGYVPIMVKSKLCNLYNLPPRVLIEHHEEAEEMGGYFIINGIEKVIRMLIVPRRNFPIAMVRPKWKSRGLGYTQFGVSMRCVREEHSAVNMNLHYVENGTVMLNFIYRKELFFLPLGFALKALVSFSDYQIFQELIKGKEEDSFFRNSVSQMLRIVIEEGCHSQKQVLNYLGECFRVKLSLPDWYPNVEAAEFLLNQCICIHLQSNTDKFYLLCLMTRKLFALARGECMDDNPDSLVNQEVLSPGQLFLMFLKEKMENWLVSIKIVLDKRAQKANVSINNENLMKIFSMGTELTRPFEYLLATGNLRSKTGLGFLQDSGLCVVADKLNFLRYLSHFRCVHRGAAFAKMRTTTVRRLLPESWGFLCPVHTPDGAPCGLLNHLTAVCEVVTKFVYTASIPALLCGLGVTPVDTAPCRPYSDCYPVLLDGVMVGWVDKDLAPEVADTLRRFKVLREKRIPPWMEVALIPMTGKPSLYPGLFLFTTPCRLVRPVQNLELGREELIGTMEQLFMNVAIFEDEVFGGISTHQELFPHSLLSVIANFIPFSDHNQSPRNMYQCQMGKQTMGFPLLTYQNRSDNKLYRLQTPQSPLVRPCMYDFYDMDNYPIGTNAIVAVISYTGYDMEDAMIVNKASWERGFAHGSVYKSEFIDLSEKFKQGEDNLVFGVKPGDPRVMQKLDDDGLPFIGAKLEYGDPYYSYLNLNTGEGFVVYYKSKENCVVDNIKVCSNDMGSGKFKCICITVRIPRNPTIGDKFASRHGQKGILSRLWPAEDMPFTESGMMPDILFNPHGFPSRMTIGMLIESMAGKSAALHGLCHDATPFIFSEENSALEYFGEMLKAAGYNFYGTERLYSGISGMELEADIFIGVVYYQRLRHMVSDKFQVRTTGARDKVTNQPLGGRNVQGGIRFGEMERDALLAHGTSFLLHDRLFNCSDRSVAHMCVECGSLLSPLLEKPPPSWSAMRNRKYNCTVCGRSDTIDTVSVPYVFRYFVAELAAMNIKVKLDVI")
  names(polr1b) <- names(drllps_scaffolds)
  
  ncl <- data.frame("-", "P09405", "ncl", "ENSMUSG00000026234", "Mus musculus", "Droplet", "Scaffold", "", "MVKLAKAGKTHGEAKKMAPPPKEVEEDSEDEEMSEDEDDSSGEEEVVIPQKKGKKATTTPAKKVVVSQTKKAAVPTPAKKAAVTPGKKAVATPAKKNITPAKVIPTPGKKGAAQAKALVPTPGKKGAATPAKGAKNGKNAKKEDSDEDEDEEDEDDSDEDEDDEEEDEFEPPIVKGVKPAKAAPAAPASEDEEDDEDEDDEEDDDEEEEDDSEEEVMEITTAKGKKTPAKVVPMKAKSVAEEEDDEEEDEDDEDEDDEEEDDEDDDEEEEEEEPVKAAPGKRKKEMTKQKEAPEAKKQKVEGSEPTTPFNLFIGNLNPNKSVNELKFAISELFAKNDLAVVDVRTGTNRKFGYVDFESAEDLEKALELTGLKVFGNEIKLEKPKGRDSKKVRAARTLLAKNLSFNITEDELKEVFEDAMEIRLVSQDGKSKGIAYIEFKSEADAEKNLEEKQGAEIDGRSVSLYYTGEKGQRQERTGKTSTWSGESKTLVLSNLSYSATKETLEEVFEKATFIKVPQNPHGKPKGYAFIEFASFEDAKEALNSCNKMEIEGRTIRLELQGSNSRSQPSKTLFVKGLSEDTTEETLKESFEGSVRARIVTDRETGSSKGFGFVDFNSEEDAKAAKEAMEDGEIDGNKVTLDWAKPKGEGGFGGRGGGRGGFGGRGGGRGGRGGFGGRGRGGFGGRGGFRGGRGGGGDFKPQGKKTKFE")
  names(ncl) <- names(drllps_scaffolds)
  
  lin28 <- data.frame("-", "Q8K3Y3", "Lin28a", "ENSMUSG00000050966", "Mus musculus", "Droplet", "Scaffold", "", "MGSVSNQQFAGGCAKAAEKAPEEAPPDAARAADEPQLLHGAGICKWFNVRMGFGFLSMTARAGVALDPPVDVFVHQSKLHMEGFRSLKEGEAVEFTFKKSAKGLESIRVTGPGGVFCIGSERRPKGKNMQKRRSKGDRCYNCGGLDHHAKECKLPPQPKKCHFCQSINHMVASCPLKAQQGPSSQGKPAYFREEEEEIHSPALLPEAQN")
  names(lin28) <- names(drllps_scaffolds)
  
  eloa <- data.frame("-", "Q8CB77", "Eloa", "ENSMUSG00000028668", "Mus musculus", "Droplet", "Scaffold", "", "MAAESALQVVEKLQARLAANPDPKKLLKYLKKLSILPITVDILVETGVGKTVNSFRKHEQVGNFARDLVAQWKKLVPVERNSEAEDQDFEKNNSRKRPRDALQREEELEGNYQESWKPSGSRSYSPEHRQKKHKKLSEPERPHKVAHSHEKRDERKRCHKVSPPYSSDPESSDYGHVQSPPPSSPHQMYTDLSRSPEEDQEPIISHQKPGKVHSNTFQDRLGVSHLGEQGKGAVSHHKQHRSSHKEKHPADAREDEKISAVSREKSHKASSKEESRRLLSGDSAKEKLPSSVVKKDKDREGSSLKKKFSPALDVASDNHFKKPKHKDSEKAKSDKNKQSVDGVDSGRGTGDPLPKAKEKVPNHLKAQEGKVRTNADGKSAGPLHPKAEETDVDDEFERPTMSFESYLSYDQPRKKKKKVVKTSSTALGEKGLKKKDSKSTSKNLNSAQKLPKVNENKSEKLQPAGAEPTRPRKVPTDVLPALPDIPLPAIHANYRPLPSLELIPSFQPKRKAFSSPQEEEEAGFTGRRMNSKMQVYSGSKCAYLPKMMTLHQQCIRVLKNNIDSIFEVGGVPYSVLEPVLERCTPDQLYRIEECNHVLIEETDQLWKVHCHRDFKEERPEEYESWREMYLRLQDAREQRLRLLTNNIRSAHANKPKGRQAKMAFVNSVAKPPRDVRRRQEKFGTGGAAVPEKVRIKPAPYTTGSSHVPASNSSSNFHSSPEELAYDGPSTSSAHLAPVASSSVSYDPRKPAVKKIAPMMAKTIKAFKNRFSRR")
  names(eloa) <- names(drllps_scaffolds)
  
  # update list of scaffolds
  scaffolds_drllps <- rbind(drllps_scaffolds, gata3, mecp2, kmt5c, cbx4, phc1, hmga1, hmga2, fbl, polr2a, polr1b, ncl, lin28, eloa)
  scaffolds_drllps_genenames <- unique(c(tolower(scaffolds_drllps$Gene.name)))
  
  
  # get scaffolds from LLPSDB v2.0 database (Wang et al, 2022, Bioinformatics 38: 2010-2014, PMID: 35025997)
  llpsdb_ambiguous_zip <- RCurl::getURLContent("http://bio-comp.org.cn/llpsdbv2/download/Phase_separation_ambiguous.zip", binary = TRUE)
  llpsdb_unambiguous_zip <- RCurl::getURLContent("http://bio-comp.org.cn/llpsdbv2/download/Phase_separation_unambiguous.zip", binary = TRUE)
  tf <- tempfile()
  readr::write_file(llpsdb_ambiguous_zip, tf)
  uz <- unzip(tf, exdir=dirname(tf), list=TRUE)
  unzip(tf, exdir=dirname(tf))
  llpsdb_ambiguous <- xlsx::read.xlsx(file=paste0(dirname(tf),"/",uz[1,1],"protein.xls"), 1, as.data.frame = TRUE)
  tf <- tempfile()
  readr::write_file(llpsdb_unambiguous_zip, tf)
  uz <- unzip(tf, exdir=dirname(tf), list=TRUE)
  unzip(tf, exdir=dirname(tf))
  llpsdb_unambiguous <- xlsx::read.xlsx(file=paste0(dirname(tf),"/",uz[1,1],"protein.xls"), 1, as.data.frame = TRUE)
  
  llpsdb_unambiguous_data <- llpsdb_unambiguous[which(llpsdb_unambiguous$Uniprot.ID != "" & llpsdb_unambiguous$Species %in% c("Homo sapiens", "Mus musculus", "Rattus norvegicus")), c(2,4,6)]	
  llpsdb_ambiguous_data <- llpsdb_ambiguous[which(llpsdb_ambiguous$UniprotID != "" & llpsdb_ambiguous$Species %in% c("Homo sapiens", "Mus musculus", "Rattus norvegicus")), c(2,3,6)]
  colnames(llpsdb_ambiguous_data) <- colnames(llpsdb_unambiguous_data)
  
  scaffold_llpsdb_data <- rbind(llpsdb_unambiguous_data,llpsdb_ambiguous_data)
  llpsdb_genenames <- as.character(tolower(scaffold_llpsdb_data$Gene.name)); llpsdb_genenames= unlist(strsplit(llpsdb_genenames,",")); llpsdb_genenames= unlist(strsplit(llpsdb_genenames,";"))
  
  # combine scaffold lists from DrLLPS and LLPSDB
  drllps_llpsdb_genenames <- unique(c(tolower(scaffolds_drllps_genenames),tolower(llpsdb_genenames))) 
  
  # convert gene names to UniProt IDs
  mart <- biomaRt::useDataset("mmusculus_gene_ensembl", useMart("ensembl", host=host))
  drllps_llpsdb_IDs <- biomaRt::getBM(filters="external_gene_name", attributes=c("uniprotswissprot", "external_gene_name"), values=drllps_llpsdb_genenames, mart=mart, useCache=FALSE)
  drllps_llpsdb_IDs <- rbind(drllps_llpsdb_IDs, c("Q64028", "Phc1"))
  
  # remove duplicates lacking an UniProt ID
  drllps_llpsdb_IDs_with_ID <- unique(drllps_llpsdb_IDs[drllps_llpsdb_IDs[,1]!="",])
  drllps_llpsdb_IDs_unique <- merge(unique(subset(drllps_llpsdb_IDs, select="external_gene_name")), drllps_llpsdb_IDs_with_ID, all.x = TRUE)
  drllps_llpsdb_IDs_unique[is.na(drllps_llpsdb_IDs_unique$uniprotswissprot),]$uniprotswissprot <- ""
  
  # return list of LLPS scaffolds
  return(drllps_llpsdb_IDs_unique[c("uniprotswissprot","external_gene_name")])
}
        

# query STRINGdb to obtain interactors
getInteractors <- function(marker, scaffolds, host, score=600, toAdd=NA, toRemove=NA) {
  # load libraries
  library("biomaRt")
  library("STRINGdb")
  
  # load biomaRt for gene symbol conversion
  mart <- biomaRt::useDataset("mmusculus_gene_ensembl", useMart("ensembl", host=host))
  
  # load STRINGdb
  string_db <- STRINGdb$new(version="11.5", species=10090, protocol="http", score_threshold=score) # default threshold
  
  # check if all scaffolds have at least one ENSEMBL ID
  scaffolds_ensembl_protein_ids <- biomaRt::getBM(filters="uniprotswissprot", attributes= c("ensembl_peptide_id", "uniprotswissprot"), values=scaffolds$uniprotswissprot, mart=mart, useCache=FALSE)
  scaffolds_ensembl_protein_ids <- sapply(scaffolds[scaffolds$Species=="Homo sapiens",]$UniProt.ID, function(x) {scaffolds_ensembl_protein_ids[scaffolds_ensembl_protein_ids[,2]==x,][1,1]} )
  if(sum(is.na(scaffolds_ensembl_protein_ids)) > 0) {
    print("The following scaffold(s) lack ENSEMBL ID(s):")
    print(names(which(is.na(scaffolds_ensembl_protein_ids))))
  }
 
  # get ENSEMBL IDs of interactors of "marker protein"
  interactors <- string_db$get_neighbors(string_db$mp(marker))
  interactors <- c(interactors, string_db$mp(marker))
  
  # get UniProt IDs of interactors of "marker protein"
  interactors_uniprot <- biomaRt::getBM(filters="ensembl_peptide_id", attributes=c( "uniprotswissprot", "external_gene_name"), values=sub(".*ENSMUSP", "ENSMUSP", interactors), mart=mart, useCache=FALSE) 
  
  # remove proteins without UniProt ID
  interactors_uniprot <- interactors_uniprot[which(interactors_uniprot$uniprotswissprot != ""),]
  nbr_interactors <- nrow(interactors_uniprot)
  
  # manually add proteins to the list (for example because relevant interactions have been documented in the cell type of interest)
  if(length(!is.na(toAdd))>0) {
    toAdd_uniprot <- biomaRt::getBM(filters="external_gene_name", attributes=c("uniprotswissprot", "external_gene_name"), values=toAdd, mart=mart, useCache=FALSE)
    interactors_uniprot <- rbind(interactors_uniprot, toAdd_uniprot)
  }
  
  # manually remove proteins from the list (for example because they not localize in the organelle of interest in the cell type of interest)
  if((length(!is.na(toRemove))>0) & (sum(toRemove %in% interactors_uniprot$external_gene_name)>0)) {
    ind <- which(interactors_uniprot$external_gene_name %in% toRemove)
    interactors_uniprot <- interactors_uniprot[-ind,]
  }
  
  
  # remove lines without UniProt IDs
  if(length(which(interactors_uniprot$uniprotswissprot == ""))>0) {
    interactors_uniprot <- interactors_uniprot[-which(interactors_uniprot$uniprotswissprot == ""),]
  }
  interactors_uniprot <- interactors_uniprot[which(!is.na(interactors_uniprot$uniprotswissprot)),]
  
  # return list of interactors
  return(interactors_uniprot)
}



intersect <- function(llpsData, interactors, getLinkStringNetwork=FALSE){
  # load library to retrieve protein sequences
  library("protr")
  
  # merge lists of interactors and scaffolds
  res <- interactors[which(interactors$uniprotswissprot %in% llpsData$uniprotswissprot),]
  
  # retrieve protein sequences
  uniID <- res$uniprotswissprot
  res$protein_sequence <- unlist(protr::getUniProt(uniID))
  
  # return result
  return(res)
}



intersectNucleolarExtract <- function(input_list, toAdd=NA) {
  # load library
  library("rio")
  
  # get list of proteins in nucleolar extract of mESCs
  # Szerlong et al, 2015, J Mol Biol 427: 2056-2071, PMID 25584861, https://doi.org/10.1016/j.jmb.2015.01.001, supplement-3.xlsx
  nucleolar_extract <- rio::import("https://europepmc.org/articles/PMC4417401/bin/NIHMS655033-supplement-3.xlsx")
  colnames(nucleolar_extract) <- nucleolar_extract[3,]
  nucleolar_extract <- nucleolar_extract[-c(1,2,3), c(1,2,4,5,6)] # remove header (first 3 lines), keep relevant columns (protein, accession number, three replicates in wildtype ESCs)
  
  # convert columns containing protein abundances to numbers
  nucleolar_extract[,3] <- as.numeric(nucleolar_extract[,3])
  nucleolar_extract[,4] <- as.numeric(nucleolar_extract[,4])
  nucleolar_extract[,5] <- as.numeric(nucleolar_extract[,5])
  
  # keep proteins that are present in all three replicates
  nucleolar_extract_filtered <- nucleolar_extract[which(nucleolar_extract[,3]>0 & nucleolar_extract[,4]>0 & nucleolar_extract[,5]>0),]
  
  # retrieve gene names
  nucleolar_extract_filtered_gene_names <- strsplit(as.character(nucleolar_extract_filtered[,1]),"GN")
  nucleolar_extract_filtered_gene_names <- unlist(lapply(nucleolar_extract_filtered_gene_names, function(x) x[[2]]))
  nucleolar_extract_filtered_gene_names <- strsplit(nucleolar_extract_filtered_gene_names," ")
  nucleolar_extract_filtered_gene_names <- unlist(lapply(nucleolar_extract_filtered_gene_names, function(x) x[[1]]))
  nucleolar_extract_filtered_gene_names <- strsplit(nucleolar_extract_filtered_gene_names,"=")
  nucleolar_extract_filtered_gene_names <- unlist(lapply(nucleolar_extract_filtered_gene_names, function(x) x[[2]]))

  # manually add proteins to the list (for example because relevant interactions have been documented in the cell type of interest)
  if(length(!is.na(toAdd))>0) {
    nucleolar_extract_filtered_gene_names <- c(nucleolar_extract_filtered_gene_names, toAdd)
  }
  
  # intersect list of proteins in nucleolar extract and input list
  return(input_list[input_list$external_gene_name %in% nucleolar_extract_filtered_gene_names,])
}


