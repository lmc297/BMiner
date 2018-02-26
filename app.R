# BMiner
# Written by Laura Carroll, Cornell University

# load required packages
library(shiny)
library(ggplot2)
library(readr)
library(stringr)
library(vegan)
library(plyr)
library(dplyr)
library(cluster)
library(magrittr)
library(ggrepel)
library(ggtree)

# thanks to daattali for the code for the spinner
# https://github.com/daattali/advanced-shiny/blob/master/plot-spinner/app.R

# begin ui
ui <- fluidPage(
  
  titlePanel(h1("BMiner")),
  
  sidebarLayout(
    sidebarPanel(
      
      conditionalPanel(condition = "input.tabs=='Home'",
                       fileInput("file1","Choose BTyper final results files to analyze",
                                 multiple=TRUE,accept=c("_final_results.txt")),
                       helpText("Final results files are text files that BTyper
                                creates. They have the extension '_final_results.txt'."),
                       fileInput("file1m","Upload a metadata file to accompany your BTyper final results files",
                                 multiple=FALSE,accept=c(".txt")),
                       helpText("A metadata file is a  2-column tab-separated (.txt) file with no header; the first column should contain the name of each final results file,
                                while the second should contain the corresponding metadata.")),# end conditional panel
      conditionalPanel(condition = "input.tabs=='Virulence'",
                       selectInput("vplot",
                                   label="Virulence Typing Analyses",
                                   choices=c(
                                             "Virulence Gene Distribution (Bar Plot)"=1, 
                                             "Non-Metric Multidimensional Scaling (NMDS)"=2,
                                             "Principal Component Analysis (PCA)"=3,
                                             "Virulene Gene Presence/Absence Table"=4),selected=1)),
      conditionalPanel(condition = "input.vplot=='2' && input.tabs=='Virulence'",
                       checkboxInput("nmdsMeta",
                                     label="Overlay Metadata",value = FALSE),
                       radioButtons("dmetric",
                                    label="Dissimilarity Index",
                                    choices=list("Jaccard"="jaccard",
                                                 "Raup-Crick"="raup"),
                                    selected="jaccard")),
      conditionalPanel(condition = "input.vplot=='3' && input.tabs=='Virulence'",# && output.listPC",
                       checkboxInput("pcaMeta",
                                     label="Overlay Metadata", value=FALSE),
                       uiOutput("listPC1"),
                       uiOutput("listPC2"),
                       uiOutput("listPC3")),
      conditionalPanel(condition = "input.tabs=='Antimicrobial Resistance'",
                       selectInput("aplot",
                                   label="AMR Typing Analyses",
                                   choices=c(
                                             "AMR Gene Distribution (Bar Plot)"=1, 
                                             "Non-Metric Multidimensional Scaling (NMDS)"=2,
                                             "Principal Component Analysis (PCA)"=3,
                                             "AMR Gene Presence/Absence Table"=4),selected=1)),
      conditionalPanel(condition = "input.aplot=='2' && input.tabs=='Antimicrobial Resistance'",
                       checkboxInput("amrnmdsMeta",
                                     label="Overlay Metadata",value = FALSE),
                       radioButtons("amrdmetric",
                                    label="Dissimilarity Index",
                                    choices=list("Jaccard"="jaccard",
                                                 "Raup-Crick"="raup"),
                                    selected="jaccard")),
      conditionalPanel(condition = "input.aplot=='3' && input.tabs=='Antimicrobial Resistance'",# && output.listPC",
                       checkboxInput("amrpcaMeta",
                                     label="Overlay Metadata", value=FALSE),
                       uiOutput("amrlistPC1"),
                       uiOutput("amrlistPC2"),
                       uiOutput("amrlistPC3")),
      conditionalPanel(condition = "input.tabs=='MLST'",
                       selectInput("mplot",
                                   label="Multi-Locus Sequence Typing (MLST) Analyses",
                                   choices = c("Sequence Type (ST) Distribution (Bar Plot)"=1,
                                               "Minimum Spanning Tree"=2),
                                   selected = 1)),
      conditionalPanel(condition = "input.mplot=='2' && input.tabs=='MLST'",
                       checkboxInput("mstreemeta",
                                     label="Overlay Metadata", value=FALSE)),
      conditionalPanel(condition = "input.tabs=='rpoB'",
                       selectInput("rplot",
                                   label="rpoB Allelic Typing (AT) Analyses",
                                   choices = c("rpoB Allelic Type (AT) Distribution (Bar Plot)"=1),
                                   selected=1)),
      conditionalPanel(condition = "input.tabs=='panC'",
                       selectInput("pplot",
                                   label="panC Clade Typing Analyses",
                                   choices=c(
                                             "panC Clade Distribution (Bar Plot)"=1),
                                   selected=1)),
      conditionalPanel(condition = "input.tabs=='16S'",
                       selectInput("splot",
                                   label="16S rDNA Typing Analyses",
                                   choices=c(
                                             "Closest 16S Sequence Distribution (Bar Plot)"=1),
                                   selected=1)),
      conditionalPanel(condition = "input.vplot=='1' && input.tabs=='Virulence'", 
                       radioButtons("bcTypev",label = "Y-axis Display",
                                    choices = list("Counts","Percent (%)"),
                                    selected="Counts")),
      conditionalPanel(condition = "input.aplot=='1' && input.tabs=='Antimicrobial Resistance'", 
                       radioButtons("bcTypea",label = "Y-axis Display",
                                    choices = list("Counts","Percent (%)"),
                                    selected="Counts")),
      conditionalPanel(condition = "input.mplot=='1' && input.tabs=='MLST'", 
                       radioButtons("bcTypem",label = "Y-axis Display",
                                    choices = list("Counts","Percent (%)"),
                                    selected="Counts")),
      conditionalPanel(condition = "input.rplot=='1' && input.tabs=='rpoB'", 
                       radioButtons("bcTyper",label = "Y-axis Display",
                                    choices = list("Counts","Percent (%)"),
                                    selected="Counts")),
      conditionalPanel(condition = "input.pplot=='1' && input.tabs=='panC'", 
                       radioButtons("bcTypep",label = "Y-axis Display",
                                    choices = list("Counts","Percent (%)"),
                                    selected="Counts")),
      conditionalPanel(condition = "input.splot=='1' && input.tabs=='16S'", 
                       radioButtons("bcTypes",label = "Y-axis Display",
                                    choices = list("Counts","Percent (%)"),
                                    selected="Counts")),
      
      conditionalPanel(condition = "input.tabs!='Home' && input.tabs != 'Metadata' && input.vplot!='4' && input.aplot!='4'",
                       downloadButton(outputId = "downloadPlot",label="Download this plot")),
      conditionalPanel(condition = "input.tabs=='Virulence' && input.vplot=='4'",
                       downloadButton(outputId = "downloadTable", label="Download this table")),
      conditionalPanel(condition = "input.tabs=='Antimicrobial Resistance' && input.aplot=='4'",
                       downloadButton(outputId = "amrdownloadTable", label="Download this table"))
                       ),#end sidebarpanel
    
    mainPanel(tabsetPanel(id = "tabs", 
                          tabPanel("Home",htmlOutput("hometext1")),
                          tabPanel("Virulence",
                                   plotOutput("virbarchart",click = "plotclick",dblclick="plotdb",brush=brushOpts(id="plotbrush",resetOnNew = TRUE)),
                                   conditionalPanel(condition = "input.vplot!='1' && input.tabs=='Virulence' && input.vplot!='4'",
                                                    fluidRow(column(width=12,h3("Selected genome(s)"),
                                                                    helpText("Click on a point to view genome name(s) and the associated coordinates."),
                                                                    helpText("The 'dist_' column refers to the distance of a point from coordinates of a mouse click"),
                                                                    helpText("Drag mouse and double-click on plot to zoom in; double-click to reset plot."))),
                                                    tableOutput("clickinfo")),
                                   conditionalPanel(condition= "input.tabs=='Virulence' && input.vplot=='4'",
                                                    tableOutput(outputId = "vtable"))
                          ),
                          tabPanel("Antimicrobial Resistance",
                                   plotOutput("amrbarchart",click = "amrplotclick",dblclick="plotdb",brush=brushOpts(id="plotbrush",resetOnNew = TRUE)),
                                   conditionalPanel(condition = "input.aplot!='1' && input.tabs=='Antimicrobial Resistance' && input.aplot!='4'",
                                                    fluidRow(column(width=12,h3("Selected genome(s)"),
                                                                    helpText("Click on a point to view genome name(s) and the associated coordinates."),
                                                                    helpText("The 'dist_' column refers to the distance of a point from coordinates of a mouse click"),
                                                                    helpText("Drag mouse and double-click on plot to zoom in; double-click to reset plot."))),
                                                    tableOutput("amrclickinfo")),
                                   conditionalPanel(condition= "input.tabs=='Antimicrobial Resistance' && input.aplot=='4'",
                                                    tableOutput(outputId = "atable"))
                          ),
                          tabPanel("MLST",plotOutput("mlstbarchart",click="stclick",dblclick="plotdb",brush=brushOpts(id="plotbrush",resetOnNew = TRUE)),
                                   conditionalPanel(condition="input.mplot!='1' && input.tabs=='MLST'",
                                                    fluidRow(column(width=12,h3("Selected genome(s)"),
                                                                    helpText("Click on a point to view genome name(s) and the associated coordinates."),
                                                                    helpText("The 'dist_' column refers to the distance of a point from coordinates of a mouse click"),
                                                                    helpText("Drag mouse and double-click on plot to zoom in; double-click to reset plot."))),
                                                    tableOutput("mlstclick"))),
                          tabPanel("rpoB",plotOutput("rpoBbarchart")),
                          tabPanel("panC",plotOutput("panCbarchart")),
                          tabPanel("16S",plotOutput("S16barchart")),
                          tabPanel("Metadata",tableOutput("metadata"))
    )) # end mainPanel
)) # end sidebar layout+ui


# server
server <- function(input, output) {
  
  # raise maximum file upload size
  options(shiny.maxRequestSize=10000*1024^2)
  
  # set xy ranges for plot clicking
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  # function to create a bar chart
  barchart<-function(vec,xlab,ylab,text.angle,vmat){
    myvec<-as.character(vec)
    myst<-as.data.frame(table(myvec))
    myst<-myst[order(myst$Freq,decreasing = TRUE),]
    if (input$bcTypev=="Counts" && input$tabs=="Virulence"){
      bc<-ggplot(myst,aes(x=reorder(myvec,-Freq),y=Freq,fill=reorder(myvec,-Freq)))+geom_bar(stat="identity",show.legend = FALSE)+xlab(label = xlab)+ylab(label = paste("Total Number of Genomes",ylab))+theme(axis.text.x = element_text(angle=text.angle))}
    else if (input$bcTypea=="Counts" && input$tabs=="Antimicrobial Resistance"){
      bc<-ggplot(myst,aes(x=reorder(myvec,-Freq),y=Freq,fill=reorder(myvec,-Freq)))+geom_bar(stat="identity",show.legend = FALSE)+xlab(label = xlab)+ylab(label = paste("Total Number of Genomes",ylab))+theme(axis.text.x = element_text(angle=text.angle))}
    else if (input$bcTypem=="Counts" && input$tabs=="MLST"){
      bc<-ggplot(myst,aes(x=reorder(myvec,-Freq),y=Freq,fill=reorder(myvec,-Freq)))+geom_bar(stat="identity",show.legend = FALSE)+xlab(label = xlab)+ylab(label = paste("Total Number of Genomes",ylab))+theme(axis.text.x = element_text(angle=text.angle))}
    else if (input$bcTyper=="Counts" && input$tabs=="rpoB"){
      bc<-ggplot(myst,aes(x=reorder(myvec,-Freq),y=Freq,fill=reorder(myvec,-Freq)))+geom_bar(stat="identity",show.legend = FALSE)+xlab(label = xlab)+ylab(label = paste("Total Number of Genomes",ylab))+theme(axis.text.x = element_text(angle=text.angle))}
    else if (input$bcTypep=="Counts" && input$tabs=="panC"){
      bc<-ggplot(myst,aes(x=reorder(myvec,-Freq),y=Freq,fill=reorder(myvec,-Freq)))+geom_bar(stat="identity",show.legend = FALSE)+xlab(label = xlab)+ylab(label = paste("Total Number of Genomes",ylab))+theme(axis.text.x = element_text(angle=text.angle))}
    else if (input$bcTypes=="Counts" && input$tabs=="16S"){
      bc<-ggplot(myst,aes(x=reorder(myvec,-Freq),y=Freq,fill=reorder(myvec,-Freq)))+geom_bar(stat="identity",show.legend = FALSE)+xlab(label = xlab)+ylab(label = paste("Total Number of Genomes",ylab))+theme(axis.text.x = element_text(angle=text.angle))}
    else if (input$tabs!="Virulence"&&input$tabs!="Antimicrobial Resistance"){
      bc<-ggplot(myst,aes(x=reorder(myvec,-Freq),y=(Freq/sum(Freq))*100,fill=reorder(myvec,-Freq)))+geom_bar(stat="identity",show.legend = FALSE)+xlab(label = xlab)+ylab(label = paste("Percentage (%) of Genomes",ylab))+theme(axis.text.x = element_text(angle=text.angle))}
    else{
      bc<-ggplot(myst,aes(x=reorder(myvec,-Freq),y=(Freq/nrow(vmat))*100,fill=reorder(myvec,-Freq)))+geom_bar(stat="identity",show.legend = FALSE)+xlab(label = xlab)+ylab(label = "Percentage (%) of Genomes with Gene")+theme(axis.text.x = element_text(angle=text.angle))}
    return(bc)
  }
  
  
  # f'n to create empty plot if data does not exist
  empty.plot<-function(text){
    ggplot()+
      annotate("text", x = 4, y = 25, size=12, label = text) + 
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())}
  
  # f'n to run nmds
  run.nmds<-function(vmat,metadata,dmetric){
    vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
    vmat<-vmat[rowSums(vmat[, -1])>0,]
    # if metadata is input
    if(!(is.null(metadata))){
      metadata<-data.frame(gsub(pattern = "_final_results.txt",replacement = "",x = metadata$V1),as.factor(metadata$V2))
      ordered_meta<-metadata[match(rownames(vmat),metadata[,1]),]
      mymat<-data.frame(ordered_meta[,2],vmat)
      mymds<-metaMDS(comm=mymat[,2:ncol(mymat)],distance=dmetric,k=2,trymax = 10000)
      mymds.scores <- as.data.frame(scores(mymds)) 
      mymds.scores$site <- rownames(mymds.scores)  
      mymds.scores$metadata <- mymat[,1]
      species.scores <- as.data.frame(scores(mymds, "species"))  
      species.scores$species <- rownames(species.scores) 
      hull.data<-NULL
      for(i in 1:length(levels(mymds.scores$metadata))){
        mylevel<-levels(mymds.scores$metadata)[i]
        print("mymds.scores")
        print(mymds.scores)
        mygroup<-mymds.scores[mymds.scores$metadata == mylevel, ][chull(mymds.scores[mymds.scores$metadata==mylevel,c("NMDS1","NMDS2")]),]
        hull.data<-rbind(hull.data,mygroup)}
      nmdsplot<-ggplot()
      nmdsplot<- nmdsplot+ 
        geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,group=metadata,colour=metadata,fill=metadata),alpha=0.30)+
        geom_text_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5,size=3) +  
        geom_point(data=mymds.scores,aes(x=NMDS1,y=NMDS2,colour=metadata,shape="a"),size=2) +
        guides(shape=FALSE) +
        coord_cartesian(xlim = ranges$x, ylim = ranges$y)
      nmds.final<-mymds.scores[,1:2]
      return(list(nmdsplot,nmds.final))}
    # if there is no metadata
    else{
      mymds<-metaMDS(comm=vmat,distance=dmetric,k=2,trymax = 10000)
      mymds.scores <- as.data.frame(scores(mymds)) 
      mymds.scores$site <- rownames(mymds.scores)  
      species.scores <- as.data.frame(scores(mymds, "species"))  
      species.scores$species <- rownames(species.scores) 
      nmdsplot<-ggplot()
      nmdsplot<- nmdsplot+ 
        geom_text_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5,size=3) +  
        geom_point(data=mymds.scores,aes(x=NMDS1,y=NMDS2,colour=factor(mymds.scores$site),shape="a"),size=2) + 
        coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
        theme(legend.position = "none")
      nmds.final<-mymds.scores[,1:2]
      return(list(nmdsplot,nmds.final))
    }}
  
  
  # create gene matrix
  make.vmat<-function(lines,task){
    # virulence
    isolate.names<-lapply(lines,function(x) strsplit(x[1],split="BTyper Results for ")[[1]][2])
    if(task=="virulence"){
      vgenes<-sapply(lines,function(x) which(vapply(x,function(r) grepl("\\|",r)&&!(grepl("rpoB\\|",r))&&!(grepl("\\(",r)),FUN.VALUE = 1)==1), simplify=FALSE)}
    else if (task=="AMR"){
      vgenes<-sapply(lines,function(x) which(vapply(x,function(r) grepl("\\|",r)&&!(grepl("rpoB\\|",r))&&(grepl("\\(",r)),FUN.VALUE = 1)==1), simplify=FALSE)}
    print("vgenes")
    print(input$tabs)
    if(length(unlist(vgenes))==0){
      return(NULL)
    }
    names(vgenes)<-isolate.names
    finalv<-lapply(1:length(vgenes),function(i) strsplit(names(vgenes[[i]]),"\\|"))
    finalv<-lapply(finalv,lapply,function(i) i[1])
    finalv<-lapply(finalv,lapply,function(i) gsub(pattern = ".*\t",replacement = "",i))
    names(finalv)<-isolate.names
    unique.genes<-unique(unlist(finalv))
    vmat<-matrix(nrow = length(finalv),ncol=length(unique.genes))
    rownames(vmat)<-names(finalv)
    colnames(vmat)<-unique.genes
    for (i in 1:nrow(vmat)){
      for (g in 1:ncol(vmat)){
        isolate<-rownames(vmat)[i]
        gene<-colnames(vmat)[g]
        vgenes<-unlist(finalv[[isolate]])
        if (gene%in%vgenes){
          testgene<-1
        }
        else{
          testgene<-0
        }
        vmat[i,g]<-testgene
      }
    }
    return(vmat)
  }
  
  
  # run pca to get total # of PCs
  test.pca<-function(vmat){
    vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
    pca<-prcomp(vmat,scale = TRUE,center = TRUE)
    nPC<-ncol(pca$x)
    return(nPC)}
  
  # run PCA for plotting
  run.pca<-function(vmat,metadata,pc1,pc2,pc3){
    pc1<-as.integer(pc1)
    pc2<-as.integer(pc2)
    pc3<-as.integer(pc3)
    vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
    pca<-prcomp(vmat,scale = TRUE,center = TRUE)
    pca3col <- as.data.frame(cbind(pca$x[,pc1],pca$x[,pc2],pca$x[,pc3]))
    colnames(pca3col) <- c(paste("PC",pc1,sep=""),paste("PC",pc2,sep=""),paste("PC",pc3,sep=""))
    pc1name<-paste("PC",pc1,sep="")
    pc2name<-paste("PC",pc2,sep="")
    pc3name<-paste("PC",pc3,sep="")
    g <- ggplot(pca3col, aes_string(x=pc1name, y=pc2name, size=pc3name)) 
    if (!is.null(metadata)){
      covar<-metadata[match(rownames(pca3col),metadata[,1]),]
      metadata<-covar[,2]
      # take this line out if you ever want to make it continuous
      metadata<-as.factor(metadata)
      g <- g + geom_point(aes(color=metadata)) + xlab(label = paste("PC",pc1,sep = "")) + ylab(label = paste("PC",pc2,sep="")) + scale_size_continuous(name=paste("PC",pc3,sep="")) + coord_cartesian(xlim = ranges$x, ylim = ranges$y)
    }
    else{
      g <- g + geom_point(aes(color=rownames(pca3col))) + guides(colour=FALSE) + xlab(label = paste("PC",pc1,sep = "")) + ylab(label = paste("PC",pc2,sep="")) + scale_size_continuous(name=paste("PC",pc3,sep="")) + coord_cartesian(xlim = ranges$x, ylim = ranges$y)
    }
    return(list(g,pca3col))}
  
  # home text
  output$hometext1 <- renderText({
    hometext1<-"<h1>Thanks for using BMiner!</h1><br><br><h3>Use the panel to your left to upload BTyper final results files and/or any associated categorical metadata.
    After uploading your data, click on the tabs to view, analyze, and interact with your data. Happy mining!</h3>"
    return(hometext1)
  })
  
  # metadata
  output$metadata <- renderTable({
    infile<-input$file1m
    if (is.null(infile)){
      return(NULL)
    } 
    if (is.null(input$file1)){
      return(NULL)
    }else{
      mtable<-read.table(infile$datapath,header=FALSE,sep = "\t")
      validate(need(ncol(mtable)==2,"Your metadata file needs to have exactly 2 columns (names of BTyper final results files in column 1, corresponding metadata in column 2)."))
      finalfiles<-input$file1
      target<-finalfiles$name
      validate(need(nrow(mtable)==length(target),"It looks like the number of final results files listed in your metadata doesn't match the number of final results files that you've uploaded. Please correct your metadata table and try uploading it again."))
      validate(need(any(is.na(match(target,mtable$V1)))==FALSE,"It looks like the names of the final results files that you supplied don't match those provided in your metadata. Please correct your metadata table and try uploading it again."))
      ordered_table<-mtable[match(target,mtable$V1),]
      return(ordered_table)
    }
  })# end metadata
  
  
  # MLST barchart/MST
  mlst<-function(){
    #output$mlstbarchart <- renderPlot({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      # create mlst barchart
      if(input$mplot=='1'){
        mgenes<-lapply(lines,function(x) which(x%in%"ST\tglp\tgmk\tilv\tpta\tpur\tpyc\ttpi"))
        st<-c()
        glp<-c()
        gmk<-c()
        ilv<-c()
        pta<-c()
        pur<-c()
        pyc<-c()
        tpi<-c()
        for (i in 1:length(lines)){
          testline<-lines[[i]][mgenes[[i]][1]+1]
          print("testline")
          print(testline)
          if(!grepl("Predicted",testline)){
            st<-append(st,strsplit(lines[[i]][mgenes[[i]][1]+1],split="\t")[[1]][1])
            glp<-append(glp,strsplit(lines[[i]][mgenes[[i]][1]+1],split="\t")[[1]][2])
            gmk<-append(gmk,strsplit(lines[[i]][mgenes[[i]][1]+1],split="\t")[[1]][3])
            ilv<-append(ilv,strsplit(lines[[i]][mgenes[[i]][1]+1],split="\t")[[1]][4])
            pta<-append(pta,strsplit(lines[[i]][mgenes[[i]][1]+1],split="\t")[[1]][5])
            pur<-append(pur,strsplit(lines[[i]][mgenes[[i]][1]+1],split="\t")[[1]][6])
            pyc<-append(pyc,strsplit(lines[[i]][mgenes[[i]][1]+1],split="\t")[[1]][7])
            tpi<-append(tpi,strsplit(lines[[i]][mgenes[[i]][1]+1],split="\t")[[1]][8])
          }}
        st<-gsub("\\*","",st)
        st <- st[!is.na(st)]
        if (length(st)>0){
          mbar<-barchart(vec=st,xlab="Sequence Type (ST)", ylab="with ST",text.angle=45,vmat=NULL)}
        else{
          text<-paste("No MLST results were detected","in the selected files.",sep="\n")
          mbar<-empty.plot(text = text)}
        return(mbar)}
      
      # create mlst minimum spanning tree
      else if(input$mplot=='2'){
        element1 <- lines[1]
        print("element1")
        print(element1)
        mgenes<-lapply(lines,function(x) which(x%in%"ST\tglp\tgmk\tilv\tpta\tpur\tpyc\ttpi"))
        isolate.names<-lapply(lines,function(x) strsplit(x[1],split="BTyper Results for ")[[1]][2])
        allst<-c()
        for (i in 1:length(lines)){
          allst<-append(allst,lines[[i]][mgenes[[i]][1]+1])
          names(allst)[i]<-isolate.names[i]}
        stmat<-matrix(nrow = length(allst),ncol = 8)
        stvec<-c()
        counter<-0
        for (a in 1:length(allst)){
          splits<-strsplit(allst[a],split = "\t")[[1]]
          splits<-gsub("\\*","",splits)
          testst<-splits[1]
          if (testst=="?"){
            counter=counter+1
            testst<-paste("?",counter,sep="")}
          stvec<-append(stvec,testst)
          names(stvec)[a]<-names(allst)[a]
          for (t in 1:length(splits)){
            if (t==1 && splits[t]=="?"){
              finalat<-paste(splits[t],counter,sep="")}
            else{
              finalat<-splits[t]}
            stmat[a,t]<-finalat
          }}
        colnames(stmat)<-c("ST","glp","gmk","ilv","pta","pur","pyc","tpi")
        stdf<-stmat[!duplicated(stmat),]
        stdf<-stdf[!duplicated(stdf[,2:ncol(stdf)]),]
        rownames(stdf)<-stdf[,1]
        stdf<-stdf[,2:ncol(stdf)]
        newdf<-as.data.frame(stdf)
        stdist<-daisy(newdf,metric="gower")
        stdist<-as.matrix(stdist,labels=TRUE)
        colnames(stdist)<-rownames(stdist)
        sttree<-spantree(stdist)
        stplot<-plot(sttree,type="t",labels=sttree$labels)
        sites<-as.data.frame(stplot$sites)
        rownames(sites)<-sttree$labels
        xdim<-sites[,1]
        ydim<-sites[,2]
        zz<-as.data.frame(table(stmat[,1]))
        zdim<-zz$Freq
        names(zdim)<-zz$Var1
        print(zdim)
        zdim<-zdim[match(rownames(sites),names(zdim))]
        alldim<-cbind(sites,zdim)
        n<-sttree$n
        children<-sttree$kid
        print(alldim)
        colnames(alldim)<-c("Dimension1","Dimension2", "Size")
        g <- ggplot(data=alldim,aes(x=Dimension1,y=Dimension2))
        mymat<-cbind(sites[-1,1],sites[-1,2],sites[children,1],sites[children,2])
        if(!input$mstreemeta|is.null(input$file1m)){
          for(i in 1:nrow(mymat)){
            xstart<-mymat[i,1]
            ystart<-mymat[i,2]
            xend<-mymat[i,3]
            yend<-mymat[i,4]
            g <- g + geom_segment(x=xstart,y=ystart,xend=xend,yend=yend,color="slategray",alpha=0.2)}
          g <- g + geom_point(aes(color=rownames(alldim),size=Size))
          g <- g + coord_cartesian(xlim = ranges$x, ylim = ranges$y) + guides(color=FALSE)}
        else{
          mimport<-input$file1m
          print(mimport)
          metadata<-read.table(mimport$datapath,header=FALSE,sep="\t")
          validate(need(ncol(metadata)==2,"Your metadata file needs to have exactly 2 columns (names of BTyper final results files in column 1, corresponding metadata in column 2)."))
          target<-names(stvec)
          validate(need(nrow(metadata)==length(target),"It looks like the number of final results files listed in your metadata doesn't match the number of final results files that you've uploaded. Please correct your metadata table and try uploading it again."))
          metadata<-data.frame(gsub(pattern = "_final_results.txt",replacement = "",x = metadata$V1),metadata$V2)
          colnames(metadata)<-c("V1","V2")
          validate(need(any(is.na(match(target,metadata$V1)))==FALSE,"It looks like the names of the final results files that you supplied don't match those provided in your metadata. Please correct your metadata table and try uploading it again."))
          ordered_table<-metadata[match(target,metadata$V1),]
          ordered_table<-cbind(ordered_table,stvec)
          print(ordered_table)
          gg_color_hue <- function(n) {
            hues = seq(15, 375, length = n + 1)
            hcl(h = hues, l = 65, c = 100)[1:n]} # courtesy of user John Colby: http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
          colorvec<-gg_color_hue(length(unique(ordered_table[,2])))
          print(colorvec)
          names(colorvec)<-unique(ordered_table[,2])
          print("colorvec")
          print(colorvec)
          for(i in 1:nrow(mymat)){
            xstart<-mymat[i,1]
            ystart<-mymat[i,2]
            xend<-mymat[i,3]
            yend<-mymat[i,4]
            g <- g + geom_segment(x=xstart,y=ystart,xend=xend,yend=yend,color="slategray",alpha=0.2)
            g <- g + geom_point(aes(color=rep_len(ordered_table[,2],nrow(alldim))))+scale_color_manual("Metadata",values=colorvec)
            g <- g + coord_cartesian(xlim = ranges$x, ylim = ranges$y)}
          g
          for (r in 1:length(rownames(alldim))){
            st<-rownames(alldim)[r]
            metarows<-ordered_table[which(ordered_table[,3]==st),]
            props<-as.data.frame(metarows)
            print(st)
            print(props)
            pie<-ggplot(props,aes(x=factor(1),fill=V2))+geom_bar(width=1)+scale_fill_manual(values=colorvec)
            pie<-pie+coord_polar(theta="y") + theme_tree() + xlab(NULL) + ylab(NULL) + theme_transparent() 
            h<-((0.2-0.1)*(alldim[r,3]-min(alldim[,3])))/(max(alldim[,3])-min(alldim[,3]))+0.1
            g%<>%subview(pie,x=alldim[r,1],y=alldim[r,2],width = h,height=h)
          }
        }
        print(g)
        
        return(g)}
    }
  }# end mlst.barchart
  
  
  # rpoB barchart
  rpoB<-function(){
    #output$rpoBbarchart <- renderPlot({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      # rpoB
      rgenes<-lapply(lines,function(x) which(x%in%"Predicted rpoB Allelic Type\tPercent (%) Identity\tPercent (%) Coverage"))
      rts<-c()
      for (i in 1:length(lines)){
        rt<-strsplit(lines[[i]][rgenes[[i]][1]+1],split="\t")[[1]][1]
        if (!grepl("Predicted",rt)){
          rts<-append(rts,strsplit(rt,split = "\\|")[[1]][2])}
      }
      rts<-gsub("\\*","",rts)
      rts<-gsub("rpoB\\|","",rts)
      ats<-strsplit(rts,"\\|")
      ats<-sapply(ats,"[",-2)
      rts<-sapply(ats,paste,collapse="_")
      rts<-gsub("_","\n",rts)
      testrts<-rts[which(rts!="NA")]
      if(length(testrts)>0){
        rbar<-barchart(vec=rts,xlab="rpoB Allelic Type (AT)", ylab="with AT",text.angle=45,vmat=NULL)}
      else{
        text<-paste("No rpoB allelic typing results were detected","in the selected files.",sep="\n")
        rbar<-empty.plot(text = text)}
      return(rbar)}
    
  } # end rpoB.barchart
  
  # panC barchart
  panC<-function(){
    #output$panCbarchart <- renderPlot({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      # panC
      pgenes<-lapply(lines,function(x) which(x%in%"panC Clade Name\tClosest Strain\tPercent (%) Identity\tPercent (%) Coverage"))
      pts<-c()
      for (i in 1:length(lines)){
        pgroup<-strsplit(lines[[i]][pgenes[[i]][1]+1],split="\t")[[1]][1:2]
        punit<-paste(pgroup[1],pgroup[2],sep="_")
        if (!grepl("Predicted",punit)){
          pts<-append(pts,punit)}
      }
      pts<-gsub("_.*","\n",pts)
      pts<-gsub("\\*","",pts)
      testpts<-pts[which(as.character(pts)!="NA\nNA")]
      print(testpts)
      if (length(testpts)>0){
        pbar<-barchart(vec=pts,xlab="panC Clade", ylab="belonging to panC Clade",text.angle=45,vmat=NULL)}
      else{
        text<-paste("No panC clade typing results were detected","in the selected files.",sep="\n")
        pbar<-empty.plot(text = text)}
      return(pbar)}
    
  } # end panC.barchart
  
  # S16 barchart
  S16<-function(){
    #output$S16barchart <- renderPlot({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      # 16S
      sgenes<-lapply(lines,function(x) which(x%in%"Predicted 16s Type\tPercent (%) Identity\tPercent (%) Coverage"))
      sts<-c()
      for (i in 1:length(lines)){
        sts<-append(sts,strsplit(lines[[i]][sgenes[[i]][1]+1],split="\t")[[1]][1])
      }
      sts<-gsub("16s_","",sts)
      sts<-gsub("_","\n",sts)
      sts<-gsub("\\*","",sts)
      teststs<-sts[which(as.character(sts)!="NA")]
      print(teststs)
      if(length(teststs)>0){
        sbar<-barchart(vec=sts,xlab="Closest 16S Sequence", ylab=NULL,text.angle=45,vmat=NULL)}
      else{
        text<-paste("No 16S  results were detected","in the selected files.",sep="\n")
        sbar<-empty.plot(text = text)}
      return(sbar)}
    
  }# end S16.barchart
  
  # create mlst plot
  output$mlstbarchart <- renderPlot({
    
    print(mlst())
    #return(myplot)
  })
  
  # create rpoB plot
  output$rpoBbarchart <- renderPlot({
    
    print(rpoB())
    #return(myplot)
  })
  
  # create panC plot
  output$panCbarchart <- renderPlot({
    
    print(panC())
    #return(myplot)
  })
  
  # create 16S plot
  output$S16barchart <- renderPlot({
    
    print(S16())
    #return(myplot)
  })
  
  output$yay<- renderText({
    print(input$tabs)
    print(input$vplot)
    print(input$aplot)
    print(input$pplot)
    print(input$rplot)
    print(input$splot)
  })
  
  
  # get x-axis PC for PCA (virulence)
  output$listPC1 <- renderUI({
    infile <- input$file1
    vselect <- input$vplot
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      if (input$tabs=="Virulence"){
        task="virulence"}
      else if (input$tabs=="Antimicrobial Resistance"){
        task="AMR"}
      vmat<-make.vmat(lines=lines,task=task)
      if(vselect=='3'){
        totalPC<-test.pca(vmat = vmat)
        nPC<-c(1:totalPC)
        print(nPC)
        selectInput("pc1",
                    label="X-axis Principal Component (PC)",
                    choices=nPC,
                    selected=1)
        
      }
    }})# end renderTable
  
  # get y axis PC for PCA (virulence)
  output$listPC2 <- renderUI({
    infile <- input$file1
    vselect <- input$vplot
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      if (input$tabs=="Virulence"){
        task="virulence"}
      else if (input$tabs=="Antimicrobial Resistance"){
        task="AMR"}
      vmat<-make.vmat(lines=lines,task=task)
      if(vselect=='3'){
        totalPC<-test.pca(vmat = vmat)
        nPC<-c(1:totalPC)
        print(nPC)
        selectInput("pc2",
                    label="Y-axis Principal Component (PC)",
                    choices=nPC,
                    selected=2)
        
      }
    }})# end renderTable
  
  
  # get z-axis (size) PC for PCA (virulence)
  output$listPC3 <- renderUI({
    infile <- input$file1
    vselect <- input$vplot
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      if (input$tabs=="Virulence"){
        task="virulence"}
      else if (input$tabs=="Antimicrobial Resistance"){
        task="AMR"}
      vmat<-make.vmat(lines=lines,task=task)
      if(vselect=='3'){
        totalPC<-test.pca(vmat = vmat)
        nPC<-c(1:totalPC)
        print(nPC)
        selectInput("pc3",
                    label="Z-axis Principal Component (PC)",
                    choices=nPC,
                    selected=3)
        
      }
    }})# end renderTable
  
  
  # get x-axis PC for PCA (AMR)
  output$amrlistPC1 <- renderUI({
    infile <- input$file1
    vselect <- input$aplot
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      if (input$tabs=="Virulence"){
        task="virulence"}
      else if (input$tabs=="Antimicrobial Resistance"){
        task="AMR"}
      vmat<-make.vmat(lines=lines,task=task)
      if(vselect=='3'){
        totalPC<-test.pca(vmat = vmat)
        nPC<-c(1:totalPC)
        print(nPC)
        selectInput("amrpc1",
                    label="X-axis Principal Component (PC)",
                    choices=nPC,
                    selected=1)
        
      }
    }})# end renderTable
  
  
  # get y-axis PC for PCA (AMR)
  output$amrlistPC2 <- renderUI({
    infile <- input$file1
    vselect <- input$aplot
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      if (input$tabs=="Virulence"){
        task="virulence"}
      else if (input$tabs=="Antimicrobial Resistance"){
        task="AMR"}
      vmat<-make.vmat(lines=lines,task=task)
      if(vselect=='3'){
        totalPC<-test.pca(vmat = vmat)
        nPC<-c(1:totalPC)
        print(nPC)
        selectInput("amrpc2",
                    label="Y-axis Principal Component (PC)",
                    choices=nPC,
                    selected=2)
        
      }
    }})# end renderTable
  
  # get z-axis (size) PC for PCA (AMR)
  output$amrlistPC3 <- renderUI({
    infile <- input$file1
    vselect <- input$aplot
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      if (input$tabs=="Virulence"){
        task="virulence"}
      else if (input$tabs=="Antimicrobial Resistance"){
        task="AMR"}
      vmat<-make.vmat(lines=lines,task=task)
      if(vselect=='3'){
        totalPC<-test.pca(vmat = vmat)
        nPC<-c(1:totalPC)
        print(nPC)
        selectInput("amrpc3",
                    label="Z-axis Principal Component (PC)",
                    choices=nPC,
                    selected=3)
        
      }
    }})# end renderTable
  
  # virulence/AMR plot
  virulence<-function(task){
    #output$virbarchart <- renderPlot({
    if(task=="virulence"){
      infile <- input$file1
      vselect <- input$vplot}
    else if (task=="AMR"){
      infile <- input$file1
      vselect <- input$aplot}
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      vmat<-make.vmat(lines=lines,task=task)
      if (is.null(vmat)){
        text = paste("\n   There don't appear to be any genes\n",
                     "       in your selected files.\n")
        g<-empty.plot(text=text)
        return(g)}
      # if virulence matrix selected
      else if (vselect=="4"){
        return(vmat)}
      # if virulence barchart selected
      else if (vselect=="1"){ 
        vs<-colSums(vmat)
        vts<-strsplit(names(vs),"\\|")
        vts<-lapply(vts, tail, n = 1L)
        vts<-unlist(vts)
        names(vs)<-vts
        finalvvec<-c()
        for (v in 1:length(vs)){
          genename<-names(vs)[v]
          times<-vs[v]
          finalvvec<-append(finalvvec,rep(genename,times))
        }
        vbar<-barchart(vec=finalvvec,xlab="Genes", ylab="with Gene",text.angle=45,vmat=vmat)
        return(vbar)}
      # if NMDS selected
      else if (vselect=="2"){
        meta.in<-input$file1m
        if (task=="virulence"){
          dmetric<-input$dmetric}
        else if (task=="AMR"){
          dmetric<-input$amrdmetric}
        if (is.null(meta.in)){
          vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
          vnmds<-run.nmds(vmat=vmat,metadata = NULL,dmetric = dmetric)
          return(vnmds[[1]])}
        else{
          if (task=="virulence"){
            overlay<-input$nmdsMeta}
          else if (task=="AMR"){
            overlay<-input$amrnmdsMeta}
          if (overlay==FALSE){
            vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
            vnmds<-run.nmds(vmat=vmat,metadata = NULL,dmetric = dmetric)
            return(vnmds[[1]])}
          if (overlay==TRUE){
            mtable<-read.table(meta.in$datapath,header=FALSE,sep="\t")
            validate(need(ncol(mtable)==2,"Your metadata file needs to have exactly 2 columns (names of BTyper final results files in column 1, corresponding metadata in column 2)."))
            finalfiles<-input$file1
            target<-finalfiles$name
            validate(need(nrow(mtable)==length(target),"It looks like the number of final results files listed in your metadata doesn't match the number of final results files that you've uploaded. Please correct your metadata table and try uploading it again."))
            validate(need(any(is.na(match(target,mtable$V1)))==FALSE,"It looks like the names of the final results files that you supplied don't match those provided in your metadata. Please correct your metadata table and try uploading it again."))
            ordered_table<-mtable[match(target,mtable$V1),]
            vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
            vnmds<-run.nmds(vmat=vmat,metadata=ordered_table,dmetric = dmetric)
            return(vnmds[[1]])}
        }}
      # if PCA
      else if (vselect=="3"){
        meta.in<-input$file1m
        if(task=="virulence"){
          if(is.null(input$pc1)){
            pc1<-1
          }
          else{
            pc1<-input$pc1
          }
          if(is.null(input$pc2)){
            pc2<-2
          }
          else{
            pc2<-input$pc2
          }
          if(is.null(input$pc3)){
            pc3<-3
          }else{
            pc3<-input$pc3
          }}
        else if (task=="AMR"){
          if(is.null(input$amrpc1)){
            pc1<-1
          }
          else{
            pc1<-input$amrpc1
          }
          if(is.null(input$amrpc2)){
            pc2<-2
          }
          else{
            pc2<-input$amrpc2
          }
          if(is.null(input$amrpc3)){
            pc3<-3
          }else{
            pc3<-input$amrpc3
          }}
        if (is.null(meta.in)){
          print(c(pc1,pc2,pc3))
          pcaplot<-run.pca(vmat = vmat,metadata = NULL,pc1 = pc1,pc2 = pc2,pc3 = pc3)[[1]]
        }
        else{
          if (input$pcaMeta==FALSE && input$tabs=="Virulence"){
            pcaplot<-run.pca(vmat = vmat,metadata = NULL,pc1 = pc1,pc2 = pc2,pc3 = pc3)[[1]]}
          else if (input$amrpcaMeta==FALSE && input$tabs=="Antimicrobial Resistance"){
            pcaplot<-run.pca(vmat = vmat,metadata = NULL,pc1 = pc1,pc2 = pc2,pc3 = pc3)[[1]]}
          else{
            mtable<-read.table(meta.in$datapath,header=FALSE,sep="\t")
            validate(need(ncol(mtable)==2,"Your metadata file needs to have exactly 2 columns (names of BTyper final results files in column 1, corresponding metadata in column 2)."))
            finalfiles<-input$file1
            target<-finalfiles$name
            validate(need(nrow(mtable)==length(target),"It looks like the number of final results files listed in your metadata doesn't match the number of final results files that you've uploaded. Please correct your metadata table and try uploading it again."))
            validate(need(any(is.na(match(target,mtable$V1)))==FALSE,"It looks like the names of the final results files that you supplied don't match those provided in your metadata. Please correct your metadata table and try uploading it again."))
            ordered_table<-mtable[match(target,mtable$V1),]
            ordered_table<-data.frame(gsub(pattern = "_final_results.txt",replacement = "",x = ordered_table$V1),ordered_table$V2)
            pcaplot<-run.pca(vmat = vmat,metadata = ordered_table,pc1 = pc1,pc2 = pc2,pc3 = pc3)[[1]]
          }}
        
        return(pcaplot)}
    }}# end virbarchart
  
  
  output$virbarchart <- renderPlot({
    
    myplot<-virulence(task="virulence")
    print(myplot)
    myplot
    #return(myplot)
  })
  
  
  output$vtable <- renderTable({
    mytable<-virulence(task="virulence")
    return(mytable)
  }) # end vtable
  
  output$amrbarchart <- renderPlot({
    
    myplot<-virulence(task="AMR")
    print(myplot)
    myplot
    return(myplot)
  })
  
  output$atable <- renderTable({
    mytable<-virulence(task="AMR")
    return(mytable)
  }) # end atable
  
  # get plot coordinates on click
  getclicks <- function(){
    infile <- input$file1
    if(input$tabs=="Virulence"){
      vselect <- input$vplot}
    else if (input$tabs=="Antimicrobial Resistance"){
      vselect <- input$aplot}
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      if (input$tabs=="Virulence"){
        task="virulence"}
      else if (input$tabs=="Antimicrobial Resistance"){
        task="AMR"}
      vmat<-make.vmat(lines=lines,task = task)
      if (input$tabs=="Virulence"){
        dmetric<-input$dmetric}
      else if (input$tabs=="Antimicrobial Resistance"){
        dmetric<-input$amrdmetric}
      vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
      if(vselect=='2'){
        vnmds<-run.nmds(vmat=vmat,metadata = NULL,dmetric = dmetric)
        mypoints<-vnmds[[2]]}
      else if(vselect=='3'){
        if(task=="virulence"){
          if(is.null(input$pc1)){
            pc1<-1
          }
          else{
            pc1<-input$pc1
          }
          if(is.null(input$pc2)){
            pc2<-2
          }
          else{
            pc2<-input$pc2
          }
          if(is.null(input$pc3)){
            pc3<-3
          }else{
            pc3<-input$pc3
          }}
        if(task=="AMR"){
          if(is.null(input$amrpc1)){
            pc1<-1
          }
          else{
            pc1<-input$amrpc1
          }
          if(is.null(input$amrpc2)){
            pc2<-2
          }
          else{
            pc2<-input$amrpc2
          }
          if(is.null(input$amrpc3)){
            pc3<-3
          }else{
            pc3<-input$amrpc3
          }}
        mypoints<-run.pca(vmat = vmat,metadata = NULL,pc1 = pc1,pc2 = pc2,pc3 = pc3)[[2]]}
      return(mypoints)}}# end getclick
  
  
  # output plot clicks
  output$clickinfo <- renderTable({
    mypoints<-getclicks()
    mypoints$Sample<-rownames(mypoints)
    print(head(mypoints))
    print(class(mypoints))
    print(input$plotclick)
    clickinfo<-nearPoints(mypoints,input$plotclick, addDist=TRUE)
    print("clickinfo")
    print(clickinfo)
    return(clickinfo)
  })
  
  # output amr plot clicks
  output$amrclickinfo <- renderTable({
    mypoints<-getclicks()
    mypoints$Sample<-rownames(mypoints)
    print(head(mypoints))
    print(class(mypoints))
    print(input$amrplotclick)
    clickinfo<-nearPoints(mypoints,input$amrplotclick, addDist=TRUE)
    print("clickinfo")
    print(clickinfo)
    return(clickinfo)
  })
  
  # get mlst mst clicks
  getclicksMLST <- function(){
    infile <- input$file1
    mselect <- input$mplot
    if (is.null(infile)){
      return(NULL)
    } else{
      # if there is a file
      allfiles <- c()
      for (i in 1:nrow(infile)){
        rfile <- read_file(infile$datapath[i])
        allfiles <- append(allfiles,rfile)
      }
      if (length(allfiles)<2){
        text<-paste("It looks like you've only selected one file.","Please select 2 or more files to use BMiner.",sep="\n")
        g<-empty.plot(text=text)
        return(g)}
      # split file by new line character
      lines<-strsplit(allfiles,"\n") 
      if(input$mplot=='2'){
        mgenes<-lapply(lines,function(x) which(x%in%"ST\tglp\tgmk\tilv\tpta\tpur\tpyc\ttpi"))
        isolate.names<-lapply(lines,function(x) strsplit(x[1],split="BTyper Results for ")[[1]][2])
        names(mgenes)<-isolate.names
        allst<-c()
        for (i in 1:length(lines)){
          allst<-append(allst,lines[[i]][mgenes[[i]][1]+1])
          names(allst)[i]<-isolate.names[i]}
        stmat<-matrix(nrow = length(allst),ncol = 8)
        stvec<-c()
        counter<-0
        for (a in 1:length(allst)){
          splits<-strsplit(allst[a],split = "\t")[[1]]
          splits<-gsub("\\*","",splits)
          testst<-splits[1]
          if (testst=="?"){
            counter=counter+1
            testst<-paste("?",counter,sep="")}
          stvec<-append(stvec,testst)
          names(stvec)[a]<-names(allst)[a]
          for (t in 1:length(splits)){
            if (t==1 && splits[t]=="?"){
              finalat<-paste(splits[t],counter,sep="")}
            else{
              finalat<-splits[t]}
            stmat[a,t]<-finalat
          }}
        colnames(stmat)<-c("ST","glp","gmk","ilv","pta","pur","pyc","tpi")
        stdf<-stmat[!duplicated(stmat),]
        stdf<-stdf[!duplicated(stdf[,2:ncol(stdf)]),]
        rownames(stdf)<-stdf[,1]
        stdf<-stdf[,2:ncol(stdf)]
        newdf<-as.data.frame(stdf)
        print("newdf")
        print(class(newdf[1,1]))
        print(newdf)
        stdist<-daisy(newdf,metric="gower")
        stdist<-as.matrix(stdist,labels=TRUE)
        colnames(stdist)<-rownames(stdist)
        sttree<-spantree(stdist)
        stplot<-plot(sttree,type="t",labels=sttree$labels)
        sites<-as.data.frame(stplot$sites)
        rownames(sites)<-sttree$labels
        print(sites)
        xdim<-sites[,1]
        ydim<-sites[,2]
        zz<-as.data.frame(table(stmat[,1]))
        zdim<-zz$Freq
        names(zdim)<-zz$Var1
        zdim<-zdim[match(rownames(sites),names(zdim))]
        alldim<-cbind(sites,zdim)
        colnames(alldim)<-c("Dimension1","Dimension2","Size")
        return(alldim)}}}# end getclickMLST
  
  
  # output mlst mst click
  output$mlstclick <- renderTable({
    mypoints<-getclicksMLST()
    mypoints$ST<-rownames(mypoints)
    print(mypoints)
    print(input$stclick)
    mlstclick<-nearPoints(df = mypoints,coordinfo = input$stclick, xvar  = "Dimension1", yvar = "Dimension2", addDist=TRUE, threshold = 100)#, xvar="Dimension1", yvar="Dimension2")
    print("mlstclick")
    print(mlstclick)
    return(mlstclick)
    #return(NULL)
  })
  
  
  # observe double click and zoom
  observeEvent(input$plotdb, {
    brush <- input$plotbrush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })# end double-click
  
  # download plot
  output$downloadPlot<-downloadHandler(
    filename<-function(){
      timext<-gsub(" ","_",Sys.time())
      paste("bminer_plot_",timext,".pdf",sep="")},
    content<-function(file){
      pdf(file,onefile = TRUE,paper = "a4r",width = "10",height = "7")
      if(input$tabs=="Virulence"){
        print(virulence(task="virulence"))}
      else if(input$tabs=="Antimicrobial Resistance"){
        print(virulence(task="AMR"))}
      else if(input$tabs=="MLST"){
        print(mlst())}
      else if(input$tabs=="rpoB"){
        print(rpoB())}
      else if(input$tabs=="panC"){
        print(panC())}
      else if(input$tabs=="16S"){
        print(S16())}
      dev.off()}
  )
  
  # download virulence table
  output$downloadTable<-downloadHandler(
    filename = function() { 
      timext<-gsub(" ","_",Sys.time())
      paste("bminer_table_",timext,".txt",sep="")},
    content = function(file) {
      mytable<-virulence(task="virulence")
      print("mytable")
      print(mytable)
      write.table(mytable, file)
    }
  )
  
  # download AMR table
  output$amrdownloadTable<-downloadHandler(
    filename = function() { 
      timext<-gsub(" ","_",Sys.time())
      paste("bminer_table_",timext,".txt",sep="")},
    content = function(file) {
      mytable<-virulence(task="AMR")
      print("mytable")
      print(mytable)
      write.table(mytable, file)
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)

