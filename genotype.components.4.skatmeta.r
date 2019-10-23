
########


the.top<-1:dim(meta.results.burden)[1]




to.unwind<-c(meta.results.burden[the.top,"gene"])# ,meta.results.skatO[1:the.top,"gene"])






################################################################################
################################################################################
################################################################################
################################################################################


genotypes<-a.indel[pass,the.samples.use] ## ordered correctly for phenotypes
snp.names<-key[pass] ## GEFOS ony name with start

#### snpinfo now A different size than a.indel since added pathways!!!  snpinfo[snpinfo[,"gene"]=="KNSTRN",]
snpinfo<-snpinfo.ori[snpinfo.ori[,"Name"] %in% snp.names,]

if(!exists("gene.weights")){
  gene.weights.subset<-1
}else{
gene.weights.subset<-gene.weights[snpinfo.ori[,"Name"] %in% snp.names] # weight in same order as snpinfo.ori
}
snpinfo<-cbind(snpinfo,gene.weights.subset)
#snpinfo[1:5,]
sum(is.na(as.numeric(snpinfo[,"gene.weights.subset"])))

###################################################/media/scratch/software/matlab/network.lic

if( sum(!(snp.names %in% snpinfo.ori[,"Name"]))>0){print("WARINING snp.names not in snpinfo- unusual!")}
dim(snpinfo)
length(snp.names)
dim(genotypes)
dim(genotypes)
print("start QC")
genotypes[genotypes=="NA"]<-NA
genotypes[genotypes=="0/0"]<-0
genotypes[genotypes=="0/1"]<-1
genotypes[genotypes=="1/1"]<-2

########### prevent any averaging
dim(genotypes)
genotypes[is.na(genotypes)]<-0
dim(genotypes)
########### prevent any averaging

num.col<-dim(genotypes)[2]
num.row<-dim(genotypes)[1]

genotypes<-as.numeric(as.matrix(genotypes))
dim(genotypes)<-c(num.row,num.col)
genotypes<-t(genotypes) # samples x SNPS

colnames(genotypes)<-snp.names
rownames(genotypes)<-gsub(".GT$","",the.samples.use)

################################################################################
################################################################################
################################################################################
################################################################################





#to.unwind.name<-to.unwind[1]
#to.unwind.name<-"EVERYTHING"
# to.unwind.name<-"TOP500"
#match(net,meta.results.burden[,"gene"])
# to.unwind.name<-"SYNON_test"
# to.unwind.name<-"Pathways"
# to.unwind.name<-"ALL_significant"
# to.unwind.name<-"ALL_significant"

snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
loci<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,"Name"] # this is IDH1 not IDH1 in cluster # are the snp.names
loci<-unique(loci)
the.genes<-unique(snpinfo.ex[,"cluster"])
the.genes<-unique(snpinfo.ex[,"gene"])
the.genes<-the.genes[!(the.genes %in% clusters.wanted)]

sort(the.genes) #245 ### if used a cluster name need to do back up to (**) the.genes<-c(the.genes,"STAG2")

############repest to clean out cluster names 

## the.genes.burden<-meta.results.burden[meta.results.burden[,"gene"] %in% the.genes,]

## the.genes.burden
## write.table(the.genes.burden,file=paste(to.unwind.name,"conponents:","Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

## the.genes.burden<-meta.results.skatO[meta.results.skatO[,"gene"] %in% the.genes,]
## #the.genes.burden
## write.table(the.genes.burden,file=paste(paste(to.unwind.name,collapse="."),"conponents:","SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



########### single point only
## length(loci)
## meta.results.burden[1:5,]
## loci<-meta.results.burden[1:550,"Name"]
###############


dim(genotypes)
#genotypes[1:5,1:5]
genotypes.ex<-genotypes[,loci]


dim(genotypes.ex)
genotypes.ex[is.na(genotypes.ex)]<-0
dim(genotypes.ex)


snpinfo.ex<-snpinfo[snpinfo[,"cluster"] %in% to.unwind,]
dim(snpinfo.ex)
dim(genotypes.ex)
dim(pheno)
#snpinfo.ex[1:5,]


########### single point p values

cohort.seq.ex <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "Name",verbose=FALSE)
## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
meta.results.burden.ex<-burdenMeta(cohort.seq.ex,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "Name")
#meta.results.burden.ex[1:5,]
#pheno[1:5,]
dim(meta.results.burden.ex)
dim(genotypes.ex)

cohort.seq.test <- skatCohort(genotypes.ex, formula, SNPInfo = snpinfo.ex, data=pheno,aggregateBy = "cluster",verbose=FALSE)

meta.results.burden.test<-burdenMeta(cohort.seq.test,wts=1,mafRange = c(0,1),SNPInfo = snpinfo.ex,aggregateBy = "cluster")
#meta.results.burden.test

## meta.results.skat.ex<-skatMeta(cohort.seq,SNPInfo = snpinfo)
#meta.results.skatO.test<-skatOMeta(cohort.seq.test,burden.wts =1,SNPInfo = snpinfo.ex,aggregateBy="cluster")
#meta.results.skatO.test
figure<- match(loci,key)

#genotypes.PD<-a.indel[figure, c("LPH-001-27_PD.GT",paste(pheno.ori[pheno.ori[,"PD"],"SAMPLE"],".GT",sep="")) ]


if("PD" %in% colnames(pheno.ori)){
PDs.get<-rep(FALSE,times=dim(pheno.ori)[1])
for(ipd in 1:length(PD.group)){
PDs.get<-PDs.get | pheno.ori[,PD.group[ipd]]
}

PDs<-pheno.ori[PDs.get,"SAMPLE"]
PDs<-pheno.ori[pheno.ori[,"AML-Child"] | pheno.ori[,"Asian-AML-Child"] | pheno.ori[,"Asian-AML"]  | pheno.ori[,"AML-NotDiagnosis-Child"] | pheno.ori[, "Asian-AML-NotDiagnosis-Child"],"SAMPLE"]

genotypes.PD<-a.indel[figure, c(paste(PDs,".GT",sep="")) ]
genotypes.PD<-t(genotypes.PD)
genotypes.PD[genotypes.PD=="NA"]<-NA
genotypes.PD[genotypes.PD=="0/0"]<-0
genotypes.PD[genotypes.PD=="0/1"]<-1
genotypes.PD[genotypes.PD=="1/1"]<-2
rownames(genotypes.PD)<-gsub(".GT","",rownames(genotypes.PD))
dim(genotypes.ex)
dim(genotypes.PD)
options(max.print=200)

muts.in.PD<-apply(genotypes.PD,2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
}else{ ## no PD class
    muts.in.PD<-rep("",times=sum(pass))
}


muts.in.cases<-apply(genotypes.ex[pheno[,cancer.group],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
muts.in.controls<-apply(genotypes.ex[pheno[,control.group],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})


 controls<- paste(pheno[pheno[,control.group],"SAMPLE"],".GT",sep="")





########################################################
check<-16

quality.cases<-rep("",times=length(loci))
quality.controls<-rep("",times=length(loci))
quality.PD<-rep("",times=length(loci))

depth.cases<-rep("",times=length(loci))
depth.fad.cases<-rep("",times=length(loci))
dup.cases<-rep("",times=length(loci))

depth.controls<-rep("",times=length(loci))
depth.fad.controls<-rep("",times=length(loci))
dup.controls<-rep("",times=length(loci))

depth.PD<-rep("",times=length(loci))
depth.fad.PD<-rep("",times=length(loci))
dup.PD<-rep("",times=length(loci))

a.indel.sub<-a.indel[figure,]
a.indel.stats.sub<-a.indel.stats[figure,]

########## for somatic guessing
## somatic.matrix.desc.full.sub<-somatic.matrix.desc.full[figure,]
## somatic.matrix.p.full.sub<-somatic.matrix.p.full[figure,]

somatic.cases<-rep("",times=length(loci))
somatic.PD<-rep("",times=length(loci))

somatic.p.cases<-rep("",times=length(loci))
somatic.p.PD<-rep("",times=length(loci))



# a.indel.stats.sub[1:5,1:20]
# dim(a.indel.sub)

check<-1
for(check in 1:length(loci)){

posn<-check


if(muts.in.PD[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"GQ",sep=".")
quality.PD[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

the.ad<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"AD",sep=".")
depth.PD[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

the.fad<-paste(unlist(strsplit(muts.in.PD[check],split=",")),FAD.col,sep=".")
depth.fad.PD[check]<-paste(a.indel.stats.sub[posn,the.fad],collapse=";")

the.dup<-paste(unlist(strsplit(muts.in.PD[check],split=",")),DUP.col,sep=".")
dup.PD[check]<-paste(a.indel.stats.sub[posn,the.dup],collapse=";")

## the.ad.soma<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"GT",sep=".")
## somatic.PD[check]<-paste(somatic.matrix.desc.full.sub[posn,the.ad.soma],collapse=";")

## the.ad.soma<-paste(unlist(strsplit(muts.in.PD[check],split=",")),"GT",sep=".")
## somatic.p.PD[check]<-paste(signif(somatic.matrix.p.full.sub[posn,the.ad.soma],digits=4),collapse=";")


a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}




if(muts.in.cases[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GQ",sep=".")
quality.cases[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

the.ad<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"AD",sep=".")
depth.cases[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

the.fad<-paste(unlist(strsplit(muts.in.cases[check],split=",")),FAD.col,sep=".")
depth.fad.cases[check]<-paste(a.indel.stats.sub[posn,the.fad],collapse=";")

the.dup<-paste(unlist(strsplit(muts.in.cases[check],split=",")),DUP.col,sep=".")
dup.cases[check]<-paste(a.indel.stats.sub[posn,the.dup],collapse=";")


## the.ad.soma<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
## somatic.cases[check]<-paste(somatic.matrix.desc.full.sub[posn,the.ad.soma],collapse=";")

## the.ad.soma<-paste(unlist(strsplit(muts.in.cases[check],split=",")),"GT",sep=".")
## somatic.p.cases[check]<-paste(signif(somatic.matrix.p.full.sub[posn,the.ad.soma],digits=4),collapse=";")



a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}

if(muts.in.controls[check]!=""){
#the.gt<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GT",sep=".")
the.gq<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"GQ",sep=".")
quality.controls[check]<-paste(a.indel.sub[posn,the.gq],collapse=",")

the.ad<-paste(unlist(strsplit(muts.in.controls[check],split=",")),"AD",sep=".")
depth.controls[check]<-paste(a.indel.sub[posn,the.ad],collapse=";")

the.fad<-paste(unlist(strsplit(muts.in.controls[check],split=",")),FAD.col,sep=".")
depth.fad.controls[check]<-paste(a.indel.stats.sub[posn,the.fad],collapse=";")

the.dup<-paste(unlist(strsplit(muts.in.controls[check],split=",")),DUP.col,sep=".")
dup.controls[check]<-paste(a.indel.stats.sub[posn,the.dup],collapse=";")

a.indel[posn,the.gq]
## a.indel[posn,the.gt]
## a.indel[posn,the.dp]
}

} # end check
##########################################################################

                             
#figure
length(figure)
dim(meta.results.burden.ex)
length(muts.in.cases)
length(muts.in.controls)
#pass[figure]
#help[figure,]
## muts.in.cases[1:10]
## x<-muts.in.cases[6]


if("capture" %in% colnames(pheno.ori)){
    
    capture.counts.cases<-apply(as.matrix(muts.in.cases),1,function(x){
        if(x!=""){
            x<-unlist(strsplit(x,split=","))
            x<-pheno.ori[pheno.ori[,"SAMPLE"] %in% x,"capture"]
            x<-table(x)
            x<-x[x!=0]
            x<-paste(paste(names(x),x,sep="="),collapse=" ;")
        }else{x<-""}
        x
    }
                                )


    if("PD" %in% colnames(pheno.ori)){
        capture.counts.PD<-apply(as.matrix(muts.in.PD),1,function(x){
            if(x!=""){
                x<-unlist(strsplit(x,split=","))
                x<-pheno.ori[pheno.ori[,"SAMPLE"] %in% x,"capture"]
                x<-table(x)
                x<-x[x!=0]
                x<-paste(paste(names(x),x,sep="="),collapse=" ;")
            }else{x<-""}
            x
        }
                                 )
    }else{
        capture.counts.PD<-rep("",times=sum(pass))
    }
    
}else{  # no capture
capture.counts.PD<-rep("",times=sum(pass))
capture.counts.cases<-rep("",times=sum(pass))    
}




if("Aligner" %in% colnames(pheno.ori)){
    
    Aligner.counts.cases<-apply(as.matrix(muts.in.cases),1,function(x){
        if(x!=""){
            x<-unlist(strsplit(x,split=","))
            x<-pheno.ori[pheno.ori[,"SAMPLE"] %in% x,"Aligner"]
            x<-table(x)
            x<-x[x!=0]
            x<-paste(paste(names(x),x,sep="="),collapse=" ;")
        }else{x<-""}
        x
    }
                                )

    if("PD" %in% colnames(pheno.ori)){
        Aligner.counts.PD<-apply(as.matrix(muts.in.PD),1,function(x){
            if(x!=""){
                x<-unlist(strsplit(x,split=","))
                x<-pheno.ori[pheno.ori[,"SAMPLE"] %in% x,"Aligner"]
                x<-table(x)
                x<-x[x!=0]
                x<-paste(paste(names(x),x,sep="="),collapse=" ;")
            }else{x<-""}
            x
        }
                                 )
    }else{
        Aligner.counts.PD<-rep("",times=sum(pass))
    }

}else{  # no sligner
Aligner.counts.PD<-rep("",times=sum(pass))
Aligner.counts.cases<-rep("",times=sum(pass))    
} # aligner



if("sample.Source" %in% colnames(pheno.ori)){    
    source.counts.cases<-apply(as.matrix(muts.in.cases),1,function(x){
        if(x!=""){
            x<-unlist(strsplit(x,split=","))
            x<-pheno.ori[pheno.ori[,"SAMPLE"] %in% x,"sample.Source"]
            x<-table(x)
            x<-x[x!=0]
            x<-paste(paste(names(x),x,sep="="),collapse=" ;")
        }else{x<-""}
        x
    }
                               )

 if("PD" %in% colnames(pheno.ori)){
    source.counts.PD<-apply(as.matrix(muts.in.PD),1,function(x){
        if(x!=""){
            x<-unlist(strsplit(x,split=","))
            x<-pheno.ori[pheno.ori[,"SAMPLE"] %in% x,"sample.Source"]
            x<-table(x)
            x<-x[x!=0]
            x<-paste(paste(names(x),x,sep="="),collapse=" ;")
        }else{x<-""}
        x
    }
                            )
      }else{
         source.counts.PD<-rep("",times=sum(pass))
    }
    
}else{  # no sligner
source.counts.PD<-rep("",times=sum(pass))
source.counts.cases<-rep("",times=sum(pass))    
} # aligner





#if(sum( c("capture","Aligner","sample.Source") %in%  colnames(pheno.ori))==3){
extra.lib.info.cases<-cbind(capture.counts.cases,Aligner.counts.cases,source.counts.cases)
extra.lib.info.PD<-cbind(capture.counts.PD,Aligner.counts.PD,source.counts.PD)
#}else{
#    source.counts.PD<-rep("",times=sum(pass))
#source.counts.cases<-rep("",times=sum(pass))
#}
################################ end counting in 



#toString(colnames(a.indel)[c(1:6,8,11,16,28,7,30,34,35,36,37:42,43,14,32,33)])
#colnames(a.indel)[1:60]

 ann.cols<-c("chr","start","end","REF","ALT","TYPE","refGene::type","knownGene::type","ensGene::type","Gene.Names","Genes.mentioned.at.ASH","refGene::location","knownGene::location","ensGene::location","OMIM (Gene::Status::OMIM::description::disease)","Consequence.Embl","Uploaded_variation.Embl","Gene.Embl","Feature.Embl", "Protein_position.Embl", "Amino_acids.Embl" , "ensGene::type","ID::maf","FILTER")# ,"rs.id")
ann.cols[!(ann.cols %in% colnames(a.indel))]
ann.cols<-ann.cols[ann.cols %in% colnames(a.indel)]

annotations<-a.indel[,ann.cols]
dim(annotations)
dim(help)
dim(summary.geno.extra)
dim(a.indel)
dim(poss.model)
length(quality.cases)
length(figure)
dim(meta.results.burden.ex)
gerp.scores<-a.indel[,"gerp.scores"]
                             
#sum(meta.results.burden.ex[,"gene"]!=loci)
## colnames(a.indel)[1:50]

## key[grep("chr17",key)[1:100]]
## grep("chr17:41197708",key)
## key[grep("10088407",key)]
#out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,28,7,30,34,37:42,43)],summary.geno.extra[figure,],high.missing[figure,],help[figure,])
## out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,28,7,30,34,37:42,43,14,32,33)],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],high.missing[figure,])
## summary.geno.extra[figure,]
## annotations[figure,]
## help[figure,]

dim(meta.results.burden.ex)
if(!exists("summary.geno.extra.ori")){summary.geno.extra.ori<-summary.geno.extra}
if(is.null(dim(summary.geno.extra.ori))){summary.geno.extra.ori<-summary.geno.extra}

## if(!exists("pass.old")){pass.old<-pass}
## if(!exists("pass.new")){pass.new<-pass}
#out<-cbind(meta.results.burden.ex,a.indel[figure,c(1:6,16,43,28,7,30,34,37:42)],summary.geno.extra[figure,c("GENO.AML","GENO.Control","GENO.AML.filt","GENO.Control.filt")],help[figure,],muts.in.cases,muts.in.controls)
a.functions<-a.indel[,c("PolyPhen.scores","SIFT.scores","PolyPhen.desc","SIFT.desc")]


posns<-match(key,filt.key)
missing<-is.na(posns)
sum(missing)
filt.sub<-filt[posns,]




if(types[itypes]=="sliding.window"){ ### att the cluster name:
  
posns<-match(meta.results.burden.ex[,"gene"],snpinfo.sliding[,"Name"])
missing<-is.na(posns)
sum(missing)
the.window<-snpinfo.sliding[posns,"cluster"]
meta.results.burden.ex<-cbind(the.window,meta.results.burden.ex)

}


all.GQ<-cbind(the.QG.cancer,the.QG.PD,the.QG.Controls)
colnames(all.GQ)<-paste(colnames(all.GQ),"GQ",sep=".")

#abundance<-cbind(test.nextera,test.trueSeq,test.bwa,test.novoalign,test.PD.nextera,test.PD.trueSeq,test.PD.bwa,test.PD.novoalign)
abundance<-cbind(bias.cancer.pval,bias.control.pval,bias.PD.pval)
abundance[1:5,]
colnames(abundance)<-paste(colnames(abundance),"Pval",sep=".")
#colnames(abundance)<-gsub("^test","EnRiched",colnames(abundance))
colnames(abundance)<-paste("EnRiched",colnames(abundance),sep=".")

truth.table<-cbind(GQ.cancer.pass,GQ.Control.pass,fail.bias.cancer,fail.bias.control,fail.bias.PD)


# filt.sub[figure,c("FILTER_SUMMARY","SUMMARY_CALLED","SUMMARY_NOT_CALLED")]
                             
## out<-cbind(meta.results.burden.ex,a.functions[figure,],gerp.scores[figure],annotations[figure,],maf.lt.all[figure,],is.benign.missense[figure],annotations[figure,],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]], filt.sub[figure,c("FILTER_SUMMARY","SUMMARY_CALLED","SUMMARY_NOT_CALLED")],pass.lose.filters[figure],pass.0.01.use[figure],pass.0.001.use[figure],pass.0.01.bad.loc[figure],pass.0.001.bad.loc[figure],pass.0.01.old[figure],pass.0.001.old[figure],help[figure,],high.missing.table[figure,],validated.posn[figure],poss.model[figure,],poss.model.lib[figure,],muts.in.cases,somatic.cases,somatic.p.cases,quality.cases,depth.fad.cases,depth.cases,dup.cases,muts.in.PD,somatic.PD,somatic.p.PD,quality.PD,depth.fad.PD,depth.PD,muts.in.controls,quality.controls,depth.fad.controls,depth.controls,dup.controls,summary.geno.extra.ori[figure,colnames(summary.geno.extra.ori)[grep("^GENO",colnames(summary.geno.extra.ori))]])

## out<-cbind(meta.results.burden.ex,a.functions[figure,],gerp.scores[figure],annotations[figure,],maf.lt.all[figure,],is.benign.missense[figure],annotations[figure,],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]], filt.sub[figure,c("FILTER_SUMMARY","SUMMARY_CALLED","SUMMARY_NOT_CALLED")],pass[figure],pass.0.01.use[figure],pass.0.005.use[figure],pass.0.001.use[figure],help[figure,],alt.counts.thresh.4.rare.in.Controls,high.missing.table[figure,],poss.model[figure,],poss.model.lib[figure,],muts.in.cases,somatic.cases,somatic.p.cases,quality.cases,depth.fad.cases,depth.cases,dup.cases,muts.in.PD,somatic.PD,somatic.p.PD,quality.PD,depth.fad.PD,depth.PD,muts.in.controls,quality.controls,depth.fad.controls,depth.controls,dup.controls,summary.geno.extra.ori[figure,colnames(summary.geno.extra.ori)[grep("^GENO",colnames(summary.geno.extra.ori))]]) ### use for AMP with one
enum<-1:dim(meta.results.burden.ex)[1]

out<-cbind(enum,meta.results.burden.ex,a.functions[figure,],gerp.scores[figure],annotations[figure,],maf.lt.all[figure,],is.benign.missense[figure],annotations[figure,],summary.geno.extra[figure,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]], filt.sub[figure,c("FILTER_SUMMARY","SUMMARY_CALLED","SUMMARY_NOT_CALLED")],pass[figure],pass.all.cohorts[figure],alt.counts.thresh.4.rare.in.Controls[figure],all.GQ[figure,],abundance[figure,],truth.table[figure,],help[figure,],high.missing.table[figure,],poss.model[figure,],poss.model.lib[figure,],muts.in.cases,somatic.cases,somatic.p.cases,quality.cases,depth.fad.cases,depth.cases,dup.cases,extra.lib.info.cases,muts.in.PD,somatic.PD,somatic.p.PD,quality.PD,depth.fad.PD,depth.PD,extra.lib.info.PD,muts.in.controls,quality.controls,depth.fad.controls,depth.controls,dup.controls,summary.geno.extra.ori[figure,colnames(summary.geno.extra.ori)[grep("^GENO",colnames(summary.geno.extra.ori))]]) ### use for out



#all.data[figure,]
#out<-cbind(meta.results.burden.ex,annotations[figure,],muts.in.cases,muts.in.controls)
dim(out)
#out[,1:13]
# help["chr7:150700484:150700484:G:A:snp",]


## table(out[,"refGene::location"])
## table(out[,"Consequence.Embl"]) # to.unwind.name<-"IDH"
getwd()
setwd(analysis.dir)
paste(paste(to.unwind,collapse="."))
paste(to.unwind.name,collapse=".")
paste(paste(to.unwind.name,collapse="."),p,"GENOTYPE.conponents.","SkatO","clusters",snap.file,"txt",sep=".")

order.by<-order(out[,"p"],decreasing=FALSE)
#enum<-1:dim(meta.results.burden.ex)[1]
out[order.by,][1:10,1:10]
setwd(analysis.dir)

##AN -In
#output.dir<- "/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/burden_test_code_new_by_paul"
#setwd(output.dir)
##AN-out
getwd()
#write.table(out[order.by,],file=paste(paste(paste0("/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/burden_test_code_new_by_paul/",to.unwind.name,collapse=".")),p,"GENOTYPE.conponents.",snap.file,sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(out[order.by,],file=paste(paste(paste0(path.save,type.project,to.unwind.name,collapse=".")),p,"GENOTYPE.conponents.",snap.file,sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

getwd()
