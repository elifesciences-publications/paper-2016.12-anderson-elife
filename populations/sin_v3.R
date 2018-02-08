#### sin analysis

### sin is at pos 16410774 in r6

### libraries
	library(data.table)
	library(foreach)
	library(doMC)
	registerDoMC(20)
	library(ggplot2)
	
### is sin's geographic distribution expected by chance?

	### first, load in pooled allele frequency data	
		dat <- foreach(i=system("ls ~/sin/pooled/*.delim", intern=T))%do%fread(i, na.strings=c("-", "NA"))
		dat <- rbindlist(dat)
		#dat[nAlt=="-", nAlt:=NA]
		dat[,nAlt := as.numeric(as.character(nAlt))]
		dat[nAlt>rd, nAlt:=rd]
	
	### clean up and identify comparable sites
		dat[,insertion:=grepl("\\+", alt)]
		dat[,length:=nchar(alt)-1]
		dat[,freq := as.numeric(as.character(nAlt))/rd]
		dat[,H := 2*freq*(1-freq)*(rd/(rd-1))]

	### get position information for indels identified in pooled allele frequency data set 
		
		dat[,use := length<=2 & rd < 120 & freq<1 & freq>0 & !is.na(freq) & chr!="X"]

		setkey(dat, chr, pos)
	
		datSites <- na.omit(dat[use==T & !duplicated(dat), c("chr", "pos"), with=F])
		datSites[,focalSite := paste(chr, pos, sep=":")]	
		datSites[,chr := paste("chr", chr, sep="")]
		datSites[,stop := pos+1 ]
		datSites[,pos:= pos -1 ]
		
		setcolorder(datSites, c("chr", "pos", "stop", "focalSite"))
		
		write.table(datSites, file="~/sin/pooled_indels.r6.bed", col.names=F, quote=F, row.names=F, sep="\t")
				
	### calculated nEff and down-sampled frequencies
		### first, pull in pooled meta-data
			load("/mnt/icy_3/2014-mel-seasonality/analysis/analyses_data/00d-sample_metadata.RData")
			
			mapping <- fread("~/sin/rtec_sample_mapping.csv")
			setnames(mapping, names(mapping), c("pop", "sample_name"))
			
			setkey(sample_data, sample_name)
			setkey(mapping, sample_name)
			
			sample_data <- merge(sample_data, mapping)
			setkey(sample_data, pop)
			
		### extract out useable data
			datFreq <- na.omit(dat[use==T])
			setkey(datFreq, pop)
		
		### tack in nChr & pop_name
			datFreq[,nChr := 2*sample_data[datFreq]$individuals]
			datFreq[,pop_name := sample_data[datFreq]$pop_name]
			datFreq[,sample_name := sample_data[datFreq]$sample_name]
		
		### finally, generate nEff & recalculate frequencies
			datFreq[, rd.down := round((rd*nChr-1)/(rd+nChr))]
			datFreq[, freq.down := round(freq*rd.down)/rd.down]
		
		### make r6 position column
			datFreq[,r6_pos:=paste(chr, pos, sep=":")]
		
		### export
			write.table(datFreq, file="~/sin/pooled_indels.r6.down.freq", quote=F, row.names=F, sep="\t")
			
				
	#### write position file and use lift-over to convert to R6 coordinates
		### get lift-over file here: wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/liftOver/dm6ToDm3.over.chain.gz		
		### get liftover program here: wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v287/liftOver
		
		### run lift-over
			system("~/sin/liftOver ~/sin/pooled_indels.r6.bed ~/sin/dm6ToDm3.over.chain ~/sin/pooled_indels.r3.out ~/sin/unmapped.out")
			system("sort -k1,1 -k2,2n  ~/sin/pooled_indels.r3.out > ~/sin/pooled_indels.r3.out.sort")
			system("cat ~/sin/pooled_indels.r3.out.sort | sed 's/^chr2L/3/g' | sed 's/^chr2R/7/g' | sed 's/^chr3L/5/g' | sed 's/^chr3R/8/g'  > ~/sin/pooled_indels.r3.out.sort.chrSwap")
		
		### iterate through each DGN tabix'd vcf file; these pull out files in r3
			fl <- system("ls ~/sin/dpgp/*.vcf.gz", intern=T)
			
			foreach(i=1:length(fl))%do%{
				print(paste(i, length(fl), sep=" / "))
				system(paste("/home/bergland/sin/get_indel_DGN.sh ", fl[i], " > ", fl[i], ".geno", sep=""))
			}
			
	### now, I need to import the data and calculate per-population allele frequencies; remember these are in r3
		### meta-data	
			popInfo <- fread("/home/bergland/sin/dpgp_pops.csv")
			popInfo <- popInfo[nchar(Country)>0]
			setnames(popInfo, names(popInfo)[1], "pop")
			setnames(popInfo, names(popInfo)[6], "lat")
			setnames(popInfo, names(popInfo)[7], "long")
			setnames(popInfo, names(popInfo)[8], "elev")
			setnames(popInfo, names(popInfo)[12], "data")

		### read in each individual's DPGP datafile
			fl <- system("ls ~/sin/dpgp/*.vcf.gz", intern=T)

			dgn <- foreach(i=1:length(fl), .errorhandling="remove")%do%{
				print(paste(i, length(fl), sep=" / "))

				dat <- fread(paste(fl[i], "geno", sep="."))
				setnames(dat, names(dat), c("chr", "pos", "ref", "alt", "indiv"))
				
				dat[,pop:=popInfo$pop[sapply(paste("^", popInfo$pop, sep=""), function(x) grepl(x, dat$indiv[1]))]]
				dat[,region:=popInfo$region[sapply(paste("^", popInfo$pop, sep=""), function(x) grepl(x, dat$indiv[1]))]]
				
			}
			dgn <- rbindlist(dgn)
		
		### pad with zeros when strain/snp combo is absent
			dgn[,allele:=1]
			
			setkey(dgn, chr, pos)			
			allSites <- dgn[!duplicated(dgn)][,c("chr", "pos", "ref", "alt"), with=F]
			
			setkey(dgn, indiv)
			allIndiv.plusMeta <- dgn[!duplicated(dgn)][,c("indiv", "pop", "region"), with=F]
			
			
			all.all <- cbind(data.table(chr=rep(allSites$chr, each=dim(allIndiv.plusMeta)[1]),
										pos=rep(allSites$pos, each=dim(allIndiv.plusMeta)[1]),
										ref=rep(allSites$ref, each=dim(allIndiv.plusMeta)[1]),
										alt=rep(allSites$alt, each=dim(allIndiv.plusMeta)[1])),
							 allIndiv.plusMeta[rep(seq_len(nrow(allIndiv.plusMeta)), dim(allSites)[1]), ])

								  
			setkey(dgn, chr, pos, indiv)
			setkey(all.all, chr, pos, indiv)
			
			all.all.missing <- all.all[!dgn]
			all.all.missing[,allele:=0]
			
			dgn <- rbind(dgn, all.all.missing)
			
			
		### add in admixture data sub-Saharan populations
			admix <- fread("~/sin/TableS5_admixture.csv")
			setnames(admix, "Chr. Arm", "chr")
			setnames(admix, names(admix)[1], "indiv")
			setnames(admix, names(admix)[3], "start")
			setnames(admix, names(admix)[4], "stop")
			admix[,start := as.numeric(start)]
			admix[,stop := as.numeric(stop)]
			
			setkey(dgn, chr, pos, indiv)
			
			registerDoMC(10)
			admixed_regions <- foreach(i=1:dim(admix)[1], .errorhandling="remove")%dopar%{
				print(paste(i, dim(admix)[1], sep=" / "))
				dgn[J(data.table(chr=sub("Chr", "", admix[i]$chr),
							pos=admix[i]$start:admix[i]$stop,
							indiv=admix[i]$indiv,
							key="chr,pos,indiv")), nomatch=0][,c("chr", "pos", "indiv"), with=F]
			}
			
			admixed_regions <- rbindlist(admixed_regions)
			setkey(admixed_regions, chr, pos, indiv)
			admixed_regions[,admix:=T]
			
			dgn <- admixed_regions[J(dgn)]
			dgn[is.na(admix), admix:=F]
				
		
		### calcualte regional allele frequencies
			dgn.freq <- rbind(dgn[pop!="RAL",
								list(nStrains=length(allele),
									 ad=sum(allele),
									 dgrp=F,
									 admixExclude=F),
								list(region, chr, pos)],
			
							dgn[pop%in%c(unique(dgn[region=="Europe"]$pop), "RAL"),
								list(nStrains=length(allele),
									 ad=sum(allele),
									 dgrp=T,
									 admixExclude=T),
								list(region, chr, pos)],
			
							dgn[grepl("Africa", region) & admix==F,
								list(nStrains=length(allele),
									 ad=sum(allele),
									 dgrp=F,
									 admixExclude=T),
								list(region, chr, pos)])
	
			dgn.freq[,freq:=ad/nStrains]

		### save 
			dgn.freq[,r3_pos:=paste(chr, pos, sep=":")]
			save(dgn.freq, file="~/sin/dgn_freq.r3.Rdata")
	
############################	
### analysis starts here ###
############################	
	
	### sin is at 12236496 in r3.dgn
	### sin is at 3R:16410774 in r6 in pooled dataset
	### 16410774 (r6) converts to 
	
	### librarie
		library(data.table)
		library(foreach)
		library(ggplot2)
		library(gridExtra)
	
	### data
		### load pre-made data sets
			load(file="~/sin/dgn_freq.r3.Rdata")
			
			datFreq <- fread("~/sin/pooled_indels.r6.down.freq")
		
		### load in map; [l]ift[O]ver.[r]esutls
			lor <- fread("~/sin/pooled_indels.r3.out.sort")
			setnames(lor, names(lor), c("chr", "start", "stop", "r6_pos"))
			lor[,chr:=gsub("chr", "", chr)]
		
			lor[,pos := start+1 ]
			lor[,r3_pos := paste(chr, start+1, sep=":") ]
		
		### tack r6 into dgn.freq
			setkey(lor, r3_pos)
			setkey(dgn.freq, r3_pos)
		
			dgn.freq[,r6_pos:=lor[dgn.freq]$r6_pos]
		
		
		### identify comprable indels
			### we base this on allele frequency in the DGRP
			compSites <- dgn.freq[dgrp==T][!is.na(r6_pos)][freq>.2 & freq<.3][,c("chr", "pos", "r3_pos", "r6_pos"), with=F]	
	
		### pull out sub-dataframes & eventually merge them together
			### DGN data
				setkey(dgn.freq, r3_pos)
				dgn.sub <- dgn.freq[J(compSites$r3_pos)]
				dgn.sub[,set:="dgn"]
			
			### nescent data
				setkey(datFreq, r6_pos)
				datFreq.temp.sub <- datFreq[J(compSites$r6_pos)]
	
				datFreq.temp.sub[pop_name%in%c("BA_12", "VI_12"), region:= "Europe"]
				datFreq.temp.sub[!pop_name%in%c("BA_12", "VI_12"), region:= pop_name]
				datFreq.temp.sub[grepl("PA", pop_name), region:= "PA"]
	
				datFreq.sub <- na.omit(datFreq.temp.sub[,list(nStrains=sum(rd.down, na.rm=T),
													  freq=sum(freq.down*rd.down, na.rm=T)/sum(rd.down, na.rm=T),
													  dgrp=F,
													  admixExclude=T,
													  set="6d"),
												 list(region, chr, pos, r6_pos)])
							  
			### merge them
				m <- rbind(dgn.sub[,names(datFreq.sub),with=F], datFreq.sub)
				m[set=="dgn" & continent!="Africa", admixExclude := T]
				m[region=="SE_NorthAmerica" & dgrp==T, region:="SE_NorthAmerica_DGRP"]
	
			### define continent column
				NorthAmerica <- c("SE_NorthAmerica", "Carribbean", 
								  "FL_10", "MA_12", "NY_12", "PA", "SC_10", "VA_12",
								  "GA_08", "ME_09", "WI_12")
			
				Africa <- c("W_Africa", "E_Africa", "C_Africa", "S_Africa")
			
			
			
				m[region%in%NorthAmerica, continent:= "NorthAmerica"]
				m[!region%in%NorthAmerica, continent:= region]
				m[region%in%Africa, continent:= "Africa"]

	### Fst	
		### functions
		
			fstFun <- function(p1, p2, d1, d2) {
				phat <- (p1*d1 + p2*d2)/(d1+d2)
				
				Htot <- 2*phat*(1-phat)
				
				Hwith <- 2*p1*(1-p1)*d1/(d1+d2) + 2*p2*(1-p2)*d2/(d1+d2)
				
				(Htot - Hwith)/Htot
			}

			hetFun <- function(p1, p2, d1, d2) {
				phat <- (p1*d1 + p2*d2)/(d1+d2)
				
				Htot <- 2*phat*(1-phat)
				
				Htot
			}

	### within & between continents
		regionPlaces <- unique(m$region)
		continentsPlaces <- unique(m$continent)
		places <- unique(regionPlaces, continentsPlaces)

		
		registerDoMC(20)
		fst.o <- foreach(i=1:(length(places)-1), .errorhandling="remove")%do%{
			if(places[i]%in%regionPlaces) setkey(m, "region")
			if(places[i]%in%continentsPlaces) setkey(m, "continent")
			
			m.i <- m[places[i]]
			m.i[,pair_set:=places[i]]
			m.i.ag <- m.i[admixExclude==T,list(freq=sum(freq*nStrains)/sum(nStrains),
								n=sum(nStrains)),
							list(pair_set, r6_pos, continent)]
							
			o <- foreach(j=(i+1):length(places), .errorhandling="remove")%dopar%{
				print(paste(places[i], places[j], sep=" | "))
				
				if(places[j]%in%regionPlaces) setkey(m, "region")
				if(places[j]%in%continentsPlaces) setkey(m, "continent")
			
				m.j <- m[places[j]]
				m.j[,pair_set:=places[j]]
				m.j.ag <- m.j[admixExclude==T,list(freq=sum(freq*nStrains)/sum(nStrains),
								n=sum(nStrains)),
							list(pair_set, r6_pos, continent)]
			
				
				
				m.ij.ag <- rbind(m.i.ag, m.j.ag)
	

				fst.ij.ag <- na.omit(m.ij.ag[,list(fst=fstFun(p1=freq[pair_set==places[i]],
													p2=freq[pair_set==places[j]],
													d1=n[pair_set==places[i]],
													d2=n[pair_set==places[j]]),
										hetTot=hetFun(p1=freq[pair_set==places[i]],
													p2=freq[pair_set==places[j]],
													d1=n[pair_set==places[i]],
													d2=n[pair_set==places[j]]),
										pop.i=places[i],
										pop.j=places[j],
										continent.i=continent[pair_set==places[i]],
										continent.j=continent[pair_set==places[j]],
										freq.i=(freq[pair_set==places[i]]),
										freq.j=(freq[pair_set==places[j]]),
										nStrains.i=(n[pair_set==places[i]]),
										nStrains.j=(n[pair_set==places[j]])),
									list(r6_pos)])
	
	
				fst.ij.ag[]
			}
			rbindlist(o)
		}
		fst.o <- rbindlist(fst.o)
	
	### save
		save(fst.o, file="~/sin/fstOut.Rdata")
		load(file="~/sin/fstOut.Rdata")
	
	### formal test
		fst.o[,minFreq:=ifelse(freq.i<freq.j, freq.i, freq.j)]
	
	
		fst.quan <- fst.o[,list(z=qnorm(rank(fst)/(length(fst)+1), 0, 1),
								quan=(rank(fst)/(length(fst)+1)),
								logit.fst=qlogis(fst),
								fst=fst,
								Htot=hetTot,
								r6_pos=r6_pos,
								continent.i=continent.i,
								continent.j=continent.j,
								pop.i=pop.i,
								pop.j=pop.j,
								minFreq=minFreq),
							list(pair=paste(pop.i, pop.j, sep=":"))]
		
		fst.quan[,continent.i:=factor(continent.i, levels=c("NorthAmerica", "Europe", "Africa", "Australia"))]
		fst.quan[,continent.j:=factor(continent.j, levels=c("NorthAmerica", "Europe", "Africa", "Australia"))]
		
		fst.quan.ag <- fst.quan[logit.fst!=-Inf,list(z.mu=mean(z),
													  logit.fst.mu=mean(logit.fst)),
												list(r6_pos, 
													continent.i=continent.i, 
													continent.j=continent.j)]
		
		fst.quan.ag[,c1:=ifelse(as.numeric(continent.i) < as.numeric(continent.j), 
								c("NorthAmerica", "Europe", "Africa", "Australia")[continent.i], 
								c("NorthAmerica", "Europe", "Africa", "Australia")[continent.j])]		
		
		fst.quan.ag[,c2:=ifelse(as.numeric(continent.i) < as.numeric(continent.j), 
								c("NorthAmerica", "Europe", "Africa", "Australia")[continent.j], 
								c("NorthAmerica", "Europe", "Africa", "Australia")[continent.i])]
		
		fst.quan.ag.ag <- fst.quan.ag[continent.i!="Australia"][continent.j!="Australia"][,
										list(z.mu=mean(z.mu),
											logit.fst.mu=mean(logit.fst.mu, na.rm=T)),
										list(c1, c2, r6_pos,
											 within=(c1==c2))]
		### plot
			### plot average normalized rank of sin relative to controls
				betweenContinent <- ggplot(fst.quan.ag.ag[within==F], 
											aes(x=z.mu)) + 
									geom_histogram() + 
									facet_grid(c1~c2) +
									geom_vline(data=fst.quan.ag.ag[within==F][r6_pos=="3R:16410774"],
												aes(xintercept=z.mu), color="red") +
									ggtitle("Between")
			
				withinNorthAmerica <- ggplot(fst.quan.ag.ag[within==T][c1=="NorthAmerica"], 
											aes(x=z.mu)) + 
									geom_histogram() + 
									facet_grid(c1~c2) +
									geom_vline(data=fst.quan.ag.ag[within==T][c1=="NorthAmerica"][r6_pos=="3R:16410774"],
												aes(xintercept=z.mu), color="red") +
									ggtitle("Within")
		
			### plot average logit
				betweenContinent <- ggplot(fst.quan.ag.ag[within==F], 
											aes(x=logit.fst.mu)) + 
									geom_histogram() + 
									facet_grid(c1~c2) +
									geom_vline(data=fst.quan.ag.ag[within==F][r6_pos=="3R:16410774"],
												aes(xintercept=logit.fst.mu), color="red") +
									ggtitle("Between")
			
				withinNorthAmerica <- ggplot(fst.quan.ag.ag[within==T][c1=="NorthAmerica"], 
											aes(x=z.mu)) + 
									geom_histogram() + 
									facet_grid(c1~c2) +
									geom_vline(data=fst.quan.ag.ag[within==T][c1=="NorthAmerica"][r6_pos=="3R:16410774"],
												aes(xintercept=z.mu), color="red") +
									ggtitle("Within")
		
			

			ggsave("~/fst.pdf", grid.arrange(betweenContinent, withinNorthAmerica, nrow=1))
			
			
	
	
#### REHH analysis
	library(data.table)
	library(rehh)
	library(doMC)
	registerDoMC(10)
	library(foreach)
	
### data	
	
	### data
		### load pre-made data sets
			load(file="~/sin/dgn_freq.r3.Rdata")
			
			datFreq <- fread("~/sin/pooled_indels.r6.down.freq")
		
		### load in map; [l]ift[O]ver.[r]esutls
			lor <- fread("~/sin/pooled_indels.r3.out.sort")
			setnames(lor, names(lor), c("chr", "start", "stop", "r6_pos"))
			lor[,chr:=gsub("chr", "", chr)]
		
			lor[,pos := start+1 ]
			lor[,r3_pos := paste(chr, start+1, sep=":") ]
		
		### tack r6 into dgn.freq
			setkey(lor, r3_pos)
			setkey(dgn.freq, r3_pos)
		
			dgn.freq[,r6_pos:=lor[dgn.freq]$r6_pos]
		
		
		### identify comprable indels
			### we base this on allele frequency in the DGRP
			compSites_NorthAmerica_30 <- dgn.freq[dgrp==T][!is.na(r6_pos)][freq>.2 & freq<.3][,c("chr", "pos", "r3_pos", "r6_pos"), with=F]	
			compSites_Africa_0 <- dgn.freq[grepl("Africa", region)][admixExclude==T][,list(freq=sum(freq*nStrains, na.rm=T)/sum(nStrains, na.rm=T)),
																							list(r3_pos)][freq<.05]
	
			setkey(compSites_NorthAmerica_30, r3_pos)
			setkey(compSites_Africa_0, r3_pos)
	
			compSites <- merge(compSites_NorthAmerica_30, compSites_Africa_0)

	### load in data for EHH test
		dat <- list()		
		dat$chr2L <- data2haplohh(hap_file="~/sin/dgrp.2L.thap",
						      map_file="~/sin/dgrp.inp",
						      chr.name="2L",
						      haplotype.in.columns=T,
						      recode.allele=TRUE,
						      min_perc_geno.hap=75,
     						  min_perc_geno.snp=75)	

		dat$chr2R <- data2haplohh(hap_file="~/sin/dgrp.2R.thap",
						      map_file="~/sin/dgrp.inp",
						      chr.name="2R",
						      haplotype.in.columns=T,
						      recode.allele=TRUE,
						      min_perc_geno.hap=75,
     						  min_perc_geno.snp=75)	

		dat$chr3L <- data2haplohh(hap_file="~/sin/dgrp.3L.thap",
						      map_file="~/sin/dgrp.inp",
						      chr.name="3L",
						      haplotype.in.columns=T,
						      recode.allele=TRUE,
						      min_perc_geno.hap=75,
     						  min_perc_geno.snp=75)	
	
		dat$chr3R <- data2haplohh(hap_file="~/sin/dgrp.3R.thap",
						      map_file="~/sin/dgrp.inp",
						      chr.name="3R",
						      haplotype.in.columns=T,
						      recode.allele=TRUE,
						      min_perc_geno.hap=75,
     						  min_perc_geno.snp=75)	
	
	
	
	### pull out relevant sites from annoying REHH object	
		chrDt <- data.table(chr=c("2L", "2R", "3L", "3R"), index=c(1:4))
	
		registerDoMC(10)	
		o <- foreach(i=c(1:dim(compSites)[1]))%dopar%{
			print(paste(i, dim(compSites)[1], sep=" / "))
		
			sp <- strsplit(compSites[i]$r3_pos, ":")[[1]]
		
			tmp <- which(dat[[chrDt[chr==sp[1]]$index]]@position>=(as.numeric(sp[2])-10) & 
						 dat[[chrDt[chr==sp[1]]$index]]@position<=(as.numeric(sp[2])+10))
		
			snp.name <- dat[[chrDt[chr==sp[1]]$index]]@snp.name[tmp]
	
			if(any(grepl("INS|DEL", snp.name))) {
		
				ambig <-length( tmp[grepl("INS|DEL", snp.name)] ) > 1
			
				tmp <- tmp[grepl("INS|DEL", snp.name)][1]
			
				return(data.table(r3_pos=compSites[i]$r3_pos,
						chrIndex=chrDt[chr==sp[1]]$index,
						mrk=tmp,
						snpName=dat[[chrDt[chr==sp[1]]$index]]@snp.name[tmp],
						ambig=ambig))
					
			} else {
				return(data.table(r3_pos=compSites[i]$r3_pos,
					chrIndex=NA,
					mrk=NA,
					snpName=NA,
					ambig=NA))

			}
		}
		o <- rbindlist(o)
		o <- na.omit(o)

	### save
		save(o, dat, file="~/sin/rehh_inputData.Rdata")
		
		load(file="~/sin/rehh_inputData.Rdata")
		
	### calculate EHH stats
		registerDoMC(10)
		
		reso <- foreach(i=c(1:dim(o)[1]), .errorhandling="remove")%dopar%{
		
			print(paste(i, dim(o)[1], sep=" / "))
	
			res <- calc_ehh(haplohh=dat[[o$chrIndex[i]]],
							mrk=o$mrk[i], 
							plotehh=F)
	
			ehh <- as.data.table(t(res$ehh))
			setnames(ehh, names(ehh), c("ehh_anc", "ehh_der"))
			ehh$snpName <- colnames(res$ehh)
			ehh$focalSNP <- o$snpName[i]
			ehh$ihs_anc <- res$ihh[1]
			ehh$ihs_der <- res$ihh[2]
			ehh$freq_all1 <- res$freq_all1
	
			### return
			ehh[ehh_anc>0 | ehh_der>0]
	
		}
		reso <- rbindlist(reso)
		
	### save
		save(reso, file="~/sin/reso.Rdata")
	
	
	### load
		load(file="~/sin/reso.Rdata")
		load(file="~/sin/rehh_inputData.Rdata")
	
		reso[, dist := as.numeric(do.call("rbind", strsplit(snpName, "_"))[,2]) - as.numeric(do.call("rbind", strsplit(focalSNP, "_"))[,2])]
	
	### aggregate
		reso.ag <- reso[freq_all1>.7 & freq_all1<.8, list(ihs=log2(mean(ihs_der/ihs_anc))), focalSNP]
		
		
			
	### plot it
		ehh_plot <- ggplot(data=reso[focalSNP=="3R_12236497_INS"]) +
					geom_line(aes(x=dist, y=ehh_anc), color="darkgrey", lwd=1) + 
					geom_line(aes(x=dist, y=ehh_der), color="blue", lwd=1) +
					theme_classic()
							
		
		ihs_dist <- ggplot(data=reso.ag, aes(x=ihs)) + geom_histogram() + 
					geom_vline(xintercept=reso.ag[focalSNP=="3R_12236497_INS"]$ihs, color="red") +
					theme_classic()
	
		ggsave("ehh.pdf", grid.arrange(ehh_plot, ihs_dist, nrow=1))
	
	
	

########################################
### world-wide allele frequency data ###
########################################
	### libraries
		library(data.table)
		library(foreach)
		library(ggplot2)
		library(ggmap)
		library(maps)
		library(mapdata)
		library(gridExtra)
		library(rgdal)
		library(scatterpie)
		
		setwd("~/sin")
		
	### pull out and parse pooled data
		
			### dump and modify in Excel to match sample names
			#	pooled <- na.omit(dat[pos==16410774])
			#	setnames(pooled, "pop", "sample")
			#	pooled[,pop:=sapply(pooled$sample, function(x) strsplit(x, "_")[[1]][1])]
			#	pooled[,season:=sapply(pooled$sample, function(x) strsplit(x, "_")[[1]][1])]
			#	write.csv(pooled, "~/pooled.csv")

		### load meta-data
			load("/mnt/icy_3/2014-mel-seasonality/analysis/analyses_data/00d-sample_metadata.RData")

		### reload
			pooled <- fread("~/sin/pooled_UPDATED.csv")
	
		### merge with sample metadata
			setkey(pooled, sample_name)
			setkey(sample_data, sample_name)
			pooled <- merge(pooled, sample_data)
		
		### aggregate
			pooled.ag <- pooled[,list(A=mean(freq), B=1-mean(freq),
									  freq=mean(freq), n=sum(rd),
									  elev=NA,
									  seqType="pooled",
									  set="all"),
								 list(pop=pop, locality=pop,
								 	  lat=latitude,
								 	  long=longitude)]
		
			pooled.season <- pooled[pop=="PA",list(A=mean(freq), B=1-mean(freq),
									  freq=mean(freq), n=sum(rd),
									  seqType=NA,
									  type="pooled"),
								 list(season, year)]
		
	### DPGP data
		### meta-data	
			popInfo <- fread("/home/bergland/sin/dpgp_pops.csv")
			popInfo <- popInfo[nchar(Country)>0]
			setnames(popInfo, names(popInfo)[1], "pop")
			setnames(popInfo, names(popInfo)[7], "lat")
			setnames(popInfo, names(popInfo)[8], "long")
			setnames(popInfo, names(popInfo)[9], "elev")
			setnames(popInfo, names(popInfo)[12], "data")

		### first, parse bigger file
			popAf <- fread("/home/bergland/sin/dpgp_freqs.csv")
	
			popAf <- foreach(i=1:dim(popInfo)[1], .combine="rbind")%do%{
				tempAf <- popAf[grepl(paste("^", popInfo[i]$pop, sep=""), V1)]
				tempAf[,pop:=popInfo[i]$pop]
				tempAf[,locality:=popInfo[i]$Locality]
				tempAf[,lat:=popInfo[i]$lat]
				tempAf[,long:=popInfo[i]$long]
				tempAf[,elev:=popInfo[i]$elev]
				tempAf
			}
	
		### now, Bergman's data
			popAf_B <- fread("/home/bergland/sin/Bergman_dpgp_freqs.csv")
			popInfo_B <- popInfo[grepl("Bergman", data)]
		
			popAf_B <- foreach(i=1:dim(popInfo_B)[1], .combine="rbind")%do%{
				tempAf <- popAf_B[grepl(paste("^", gsub("_", "", popInfo_B[i]$pop), sep=""), V1)]
				tempAf[,pop:=popInfo_B[i]$pop]
				tempAf[,locality:=popInfo_B[i]$Locality]
				tempAf[,lat:=popInfo_B[i]$lat]
				tempAf[,long:=popInfo_B[i]$long]
				tempAf[,elev:=popInfo_B[i]$elev]
				tempAf
			}
		
		### now, Pool's separate data
			popAf_P <- fread("/home/bergland/sin/Pool_dpgp_freqs.csv")
			popInfo_P <- popInfo[pop%in%unique(substr(popAf_P$V1, 0, 2))]
		
			popAf_P <- foreach(i=1:dim(popInfo_P)[1], .combine="rbind")%do%{
				tempAf <- popAf_P[grepl(paste("^", gsub("_", "", popInfo_P[i]$pop), sep=""), V1)]
				tempAf[,pop:=popInfo_P[i]$pop]
				tempAf[,locality:=popInfo_P[i]$Locality]
				tempAf[,lat:=popInfo_P[i]$lat]
				tempAf[,long:=popInfo_P[i]$long]
				tempAf[,elev:=popInfo_P[i]$elev]
				tempAf
			}

		### And, Clark's separate data
			popAf_C <- fread("/home/bergland/sin/Clark_dpgp_freqs.csv")
			popInfo_C <- popInfo[data=="Grenier et al. 2015"]
		
			popAf_C <- foreach(i=1:dim(popInfo_C)[1], .combine="rbind")%do%{
				tempAf <- popAf_C[grepl(paste("^", gsub("_", "", popInfo_C[i]$pop), sep=""), V1)]
				tempAf[,pop:=popInfo_C[i]$pop]
				tempAf[,locality:=popInfo_C[i]$Locality]
				tempAf[,lat:=popInfo_C[i]$lat]
				tempAf[,long:=popInfo_C[i]$long]
				tempAf[,elev:=popInfo_C[i]$elev]
				tempAf
			}
		
		### merge
			popAf <- rbind(popAf, popAf_B, popAf_P, popAf_C)
	
		### tack in information about admixture in Subsaharan Africa.
			admix <- fread("~/sin/TableS5_admixture.csv")
			setnames(admix, "Chr. Arm", "chr")
			setnames(admix, names(admix)[1], "V1")
			setnames(admix, names(admix)[3], "start")
			setnames(admix, names(admix)[4], "stop")
			admix[,start := as.numeric(start)]
			admix[,stop := as.numeric(stop)]
			
			flagStrains <- admix[chr=="Chr3R"][start<=12236496 & stop>=12236496]
			setkey(flagStrains, V1)
			setkey(popAf, V1)
			
			
			
			#popAf <- flagStrains[popAf][is.na(chr)][,c("V1", "pop", "locality", "lat", "long", "elev", "V2"), with=F]
			
			popAf[,admixed:=ifelse(V1%in%flagStrains$V1, T, F)]
			
		### 
			popAf.ag <- rbind(popAf[admixed==F,list(A=mean(V2), 
													B=1-mean(V2),
													freq=mean(V2), 
													n=length(V2),
													set="admix"), 
												list(pop, locality, lat, long, elev)],
							popAf[,list(A=mean(V2), 
													B=1-mean(V2),
													freq=mean(V2), 
													n=length(V2),
													set="full"), 
												list(pop, locality, lat, long, elev)])

			popAf.ag <- popAf.ag[order(n, decreasing=T)]
			popAf.ag[,seqType:="haplotype"]
	
		
	### mega-merge
		sinAF <- rbind(popAf.ag, pooled.ag)
		sinAF$latOrig <- sinAF$lat
		sinAF$longOrig <- sinAF$long		
		setnames(sinAF, "type", "seqType")
	
	### try again
		world <- map_data("world")
		### Africa
			### africa w/o admixture			
				admixRemove.africa <- 
						ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
						coord_quickmap(xlim=c(-25, 55), ylim=c(-35, 35)) +
						geom_scatterpie(data=sinAF[seqType=="haplotype"][set=="admix"][order(n, decreasing=T)][n>4], 
										aes(x=long, y=lat, group=pop, r=log2(n+1), fill=c("grey", "blue")), 
										cols=LETTERS[1:2], color="black", alpha=.8, lwd=.25)	+
						geom_scatterpie(data=sinAF[seqType=="haplotype"][set=="admix"][order(n, decreasing=T)][n<=4], 
										aes(x=long, y=lat, group=pop, r=log2(n+1), fill=c("grey", "blue")), 
										cols=LETTERS[1:2], color="black", alpha=.8, lwd=.25)	+			
						geom_scatterpie_legend(radius=c(1, 2, 4, 6, 8), x=-20, y=-20, n=8, labeller=function(x) 2^(x))	+
						ggtitle("African haplotypes - admixed haplotypes removed")			
			
			### africa w/ admixture			
				admixInclude.africa <- 
						ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
						coord_quickmap(xlim=c(-25, 55), ylim=c(-35, 35)) +
						geom_scatterpie(data=sinAF[seqType=="haplotype"][set=="full"][order(n, decreasing=T)][n>4], 
										aes(x=long, y=lat, group=pop, r=log2(n+1), fill=c("grey", "blue")), 
										cols=LETTERS[1:2], color="black", alpha=.8, lwd=.25)	+
						geom_scatterpie(data=sinAF[seqType=="haplotype"][set=="full"][order(n, decreasing=T)][n<=4], 
										aes(x=long, y=lat, group=pop, r=log2(n+1), fill=c("grey", "blue")), 
										cols=LETTERS[1:2], color="black", alpha=.8, lwd=.25)	+			
									geom_scatterpie_legend(radius=c(1, 2, 4, 6, 8), x=-20, y=-20, n=8, labeller=function(x) 2^(x))	+
						ggtitle("African haplotypes - admixed haplotypes present")			
					
			###  make composite African map
				ggsave("~/africa_sin.pdf", grid.arrange(admixRemove.africa, admixInclude.africa, nrow=2))
		
		### Europe
			### haplos
				europeHaplos <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
						coord_quickmap(xlim=c(-10, 20), ylim=c(35, 55)) +
						geom_scatterpie(data=sinAF[seqType=="haplotype"][set=="admix"][order(n, decreasing=T)][n>4], 
										aes(x=long, y=lat, group=pop, r=log2(n+1), fill=c("grey", "blue")), 
										cols=LETTERS[1:2], color="black", alpha=.8, lwd=.25)	+
						geom_scatterpie(data=sinAF[seqType=="haplotype"][set=="admix"][order(n, decreasing=T)][n<=4], 
										aes(x=long, y=lat, group=pop, r=log2(n+1), fill=c("grey", "blue")), 
										cols=LETTERS[1:2], color="black", alpha=.8, lwd=.25)	+			
						geom_scatterpie_legend(radius=c(1, 2, 4), x=15, y=40, n=8, labeller=function(x) 2^(x))	+
						ggtitle("European haplotypes")			
			
			### pooled
				europePooled <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
						coord_quickmap(xlim=c(-10, 20), ylim=c(35, 55)) +
						geom_scatterpie(data=sinAF[seqType=="pooled"][set=="all"][order(n, decreasing=T)][n>4], 
										aes(x=long, y=lat, group=pop, r=log(n+1, 5), fill=c("grey", "blue")), 
										cols=LETTERS[1:2], color="black", alpha=.8, lwd=.25)	+			
						geom_scatterpie_legend(radius=c(0:3), x=15, y=40, n=8, labeller=function(x) 5^(x))	+
						ggtitle("European pooled samples")			
		
			### composite European
				ggsave("~/europe_sin.pdf", grid.arrange(europeHaplos, europePooled))
		
		### North America
			### haplos
				namHaplos <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
						coord_quickmap(xlim=c(-100, -50), ylim=c(15, 50)) +
						geom_scatterpie(data=sinAF[seqType=="haplotype"][set=="admix"][order(n, decreasing=T)][n>4], 
										aes(x=long, y=lat, group=pop, r=log(n+1, 5), fill=c("grey", "blue")), 
										cols=LETTERS[1:2], color="black", alpha=.8, lwd=.25)	+
						geom_scatterpie(data=sinAF[seqType=="haplotype"][set=="admix"][order(n, decreasing=T)][n<=4], 
										aes(x=long, y=lat, group=pop, r=log(n+1, 5), fill=c("grey", "blue")), 
										cols=LETTERS[1:2], color="black", alpha=.8, lwd=.25)	+			
						geom_scatterpie_legend(radius=c(0:3), x=-60, y=30, n=8, labeller=function(x) 5^(x))	+
						ggtitle("European haplotypes")			
			
			### pooled
				namPooled <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
						coord_quickmap(xlim=c(-100, -50), ylim=c(15, 50)) +
						geom_scatterpie(data=sinAF[seqType=="pooled"][set=="all"][order(n, decreasing=T)][n>4], 
										aes(x=long, y=lat, group=pop, r=log(n+1, 5), fill=c("grey", "blue")), 
										cols=LETTERS[1:2], color="black", alpha=.8, lwd=.25)	+			
						geom_scatterpie_legend(radius=c(0:3), x=-60, y=30, n=8, labeller=function(x) 5^(x))	+
						ggtitle("European pooled samples")			
		
			### composite European
				ggsave("~/nam_sin.pdf", grid.arrange(namHaplos, namPooled))

		
		
		
		
		
				
		### without african admixture.
			world <- map_data("world")
			world <- world[world$region!="Antarctica",]
			
			admix.haplos <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
					coord_cartesian(ylim = c(-35, 60), xlim=c(-130, 50))+

					geom_scatterpie(data=sinAF[n>=15][seqType=="haplotype"][set=="admix"], aes(x=long, y=lat, group=pop, r=sqrt(4*log(n))), 
									cols=LETTERS[1:2], color="black", alpha=.8, lwd=.5) +
					geom_scatterpie(data=sinAF[n<15][seqType=="haplotype"][set=="admix"], aes(x=long, y=lat, group=pop, r=sqrt(4*log(n))), 
									cols=LETTERS[1:2], color="black", alpha=.8, lwd=.5)	
					
	
	
		### US only
			world <- map_data("world")

			g2 <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
					coord_fixed(1.3) +
					geom_scatterpie(data=sinAF[n>=15][seqType=="pooled"][set=="all"], aes(x=long, y=lat, group=pop, r=sqrt(3*log(n))), 
									cols=LETTERS[1:2], color="black", alpha=.8, lwd=.5) +
					geom_scatterpie(data=sinAF[n<15][seqType=="pooled"][set=="all"], aes(x=long, y=lat, group=pop, r=sqrt(3*log(n))), 
									cols=LETTERS[1:2], color="black", alpha=.8, lwd=.5) +
					coord_cartesian(ylim = c(20, 60), xlim=c(-100, 25))

		### with NO african admixture; africa only.
			world <- map_data("world")
			
			all.haplos.africa <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
				coord_map("cylindrical")+					
					geom_scatterpie(data=sinAF[n>=15][seqType=="haplotype"][set=="full"], aes(x=long, y=lat, group=pop, r=sqrt(3*log(n))), 
									cols=LETTERS[1:2], color="black", alpha=.8, lwd=.5) +
					geom_scatterpie(data=sinAF[n<15][seqType=="haplotype"][set=="full"], aes(x=long, y=lat, group=pop, r=sqrt(3*log(n))), 
									cols=LETTERS[1:2], color="black", alpha=.8, lwd=.5) +				
					coord_cartesian(ylim = c(-30, 20), xlim=c(-25, 50)) +
					ggtitle("No admixture filtering")

		### with filtering; africa only
			filter.haplos.africa <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
				coord_map("cylindrical")+					
					geom_scatterpie(data=sinAF[n>=15][seqType=="haplotype"][set=="admix"], aes(x=long, y=lat, group=pop, r=sqrt(3*log(n))), 
									cols=LETTERS[1:2], color="black", alpha=.8, lwd=.5) +
					geom_scatterpie(data=sinAF[n<15][seqType=="haplotype"][set=="admix"], aes(x=long, y=lat, group=pop, r=sqrt(3*log(n))), 
									cols=LETTERS[1:2], color="black", alpha=.8, lwd=.5) +				
					coord_cartesian(ylim = c(-30, 20), xlim=c(-25, 50)) +
					ggtitle("Admixture filtering")

	
	    ### joint plots
	    	ggsave("~/sin/map_africa.pdf", grid.arrange(all.haplos.africa, filter.haplos.africa, nrow=2), device="pdf")
	    	ggsave("~/sin/map_world.pdf", grid.arrange(admix.haplos, g2, nrow=2), device="pdf")
                       		
		
		
		
		
		
		
		
		
