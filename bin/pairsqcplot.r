#!/usr/bin/env Rscript
######################################################
## modified from https://github.com/4dn-dcic/pairsqc #
## source+commit:https://github.com/4dn-dcic/pairsqc/blob/33a332bd49589b37b5046a037add48425e38ca5b/plot.r
## download date: 11/18/2021, commit: 33a332b
## This source code is licensed under the MIT license
## modified by Jianhong to combine multiple files into one R file and fix some bugs.
######################################################

args = commandArgs(TRUE)
RE_len = as.numeric(args[1]) # either 4 or 6
if(length(args)>1) { report_dir = args[2] } else { report_dir = './report' }

rainbow_w_offset <- function(L, offset = NA){
    if(is.na(offset)) offset = floor(L/2)
    color = rainbow(L+ offset)[1:L]
    return( color )
}

# set libPath to pwd
pwd <- getwd()
pwd <- file.path(pwd, "lib")
dir.create(pwd)
.libPaths(c(pwd, .libPaths()))

# (this maybe a issue)
# both plotosaurus and interactive_multiline_d3prep is not test on multiple platforms
# library(plotosaurus)
# difficult to load the plotosaurus in docker, copied the source code from https://github.com/SooLee/plotosaurus
COLOR3B <-function() { c('grey80','grey60',1) }
COLOR4 <-function() { c('dark violet','dodgerblue2','mediumseagreen','lightcoral') }
exp_axis<-function(x, axis_ind, n=5){
    at=pretty(x, n)
    at_exp = sapply(at, function(xx)as.expression(bquote(10^ .(xx))))
    axis(axis_ind, at=at, labels=at_exp)
}
get_preset <- function( n=1, layout=''){
    # 1 plot in one page
    if(n==1){
        preset = list(
            plts = list(
                c(0.2,0.95,0.2,0.95)
            ),
            width=7,
            height=7
        )

    # 2 plots in one page
    } else if(n==2){

        if(layout == 'h50' ){
            preset = list(
                plts=list(
                    c(0.1,0.47,0.2,0.95),
                    c(0.58,0.95,0.2,0.95)
                ),
                width=12,
                height=6
            )
        } else if(layout == 'v50' ){
            preset = list(
                plts=list(
                    c(0.2,0.95,0.2,0.5),
                    c(0.2,0.95,0.65,0.95)
                ),
                width=7,
                height=7
            )
        }

    # 3 plots in one page
    } else if(n==3){

        preset = list(
            plts=list(
                c(0.1,0.47,0.2,0.95),
                c(0.58,0.95,0.45,0.95),
                c(0.58,0.95,0.2,0.35)
            ),
            width=12,
            height=6
        )
    }

    return(preset)
}
pngpdf<-function(plot.func,fprefix,width=7,height=7,pointsize=12,add.date=TRUE,skip.png=FALSE,nplots=1,png.res=600, display.only=FALSE){

    if(add.date==TRUE){ date.tag=paste(".",gsub("-","",Sys.Date()),sep="")  ## eg. ".20141209"
    } else { date.tag="" }

    if(display.only) {
        return.value = plot.func()

    } else {

        if(skip.png==FALSE){
            png(paste(fprefix,date.tag,".png",sep=""),width=width,height=height*nplots,pointsize=pointsize*(nplots^0.3),units="in",res=png.res)
            par(mfrow=c(nplots,1))
            plot.func()
            dev.off()
        }

        pdf.options(useDingbats=FALSE)
        pdf(paste(fprefix,date.tag,".pdf",sep=""),width=width,height=height,pointsize=pointsize)
        return.value = plot.func()
        dev.off()

    }

    invisible(return.value)
}
pngpdf_preset<-function (plotfuncs, name, stylefunc=stylefunc0, preset=get_preset(length(plotfuncs)), ...)
{
    pngpdf(function() {
        retn=list()
        stylefunc()
        par(plt = preset$plts[[1]])
        retn[[1]] = plotfuncs[[1]]()
        if (length(plotfuncs) > 1) {
            for (i in 2:length(plotfuncs)) {
                par(new = T, plt = preset$plts[[i]])
                retn[[i]] = plotfuncs[[i]]()
            }
        }
        return(retn);
    }, name, width = preset$width, height = preset$height, ...)
}
stylefunc0<-function(){
    par(las=1, lend=2)
    par(bg='white',fg=1,col.lab=1, col.axis=1)
}
#### end of copy ######

library(stringr)
# source("https://raw.githubusercontent.com/SooLee/d3forNozzleR/master/interactive_multiline_d3prep.r") # for d3 integration
## Author: Soo Lee (duplexa@gmail.com), Department of Biomedical Informatics, Harvard Medical School.
create_d3_js_for_interactive_multiline_plot <- function(
    sample_name = "mysample", # Some kind of plot identification (in case you want to add multiple plots into a single html)
    xmin= 0,
    xmax= 1,
    ymin= 0,
    ymax= 1,
    div_id = paste("d3div_myplot_", sample_name, "___", sep=""), # Your d3 plot will be placed in a div section with this div_id.
    xlab = 'x',
    ylab = 'y',
    report_dir = "report",  ## output directory that will contain html, js, tsv etc.
    js_file_prefix = 'interactive_multiline',

    # tsvfile relative path ( to be used inside js file)
    tsvfile = "sample_tsvfile.tsv", # relative to report dir

    # data frame for column pairs, colors (on and off) and labels.
    tsvcoldata = read.table("sample/sample_tsvcolfile.tsv", comment.char="", sep="\t", header=T, stringsAsFactors=FALSE)
){

    # The function assumes a directory structure of a report_dir which contains all necessary html, js and tsv.
    # the data tsv file must already be in the report_dir.

    # tsvcolfile relative path (to be used inside js file)
    tsvcolfile =paste("./tsvcol.", div_id, ".tsv", sep="") # relative to report dir

    # create tsvcolfilecurr from data frame tsvcoldata
    tsvcolfilecurr =paste(report_dir, "/tsvcol.", div_id, ".tsv", sep="") # relative to current dir
    write.table(tsvcoldata, tsvcolfilecurr, quote=FALSE, row.names=FALSE, sep="\t")

    # create js file
    js_content = c(
        paste("d3.tsv('", tsvcolfile, "', function(tsvcolumns) {", sep=""),
        paste( "interactive_multiline_plot('", tsvfile, "', tsvcolumns, ", xmin, ",", xmax, ",", ymin, ",", ymax ,", '", xlab, "', '", ylab, "', '", div_id , "');", sep=""),
        "});"
    )
    js_file = paste("./", js_file_prefix, ".", sample_name, ".js", sep="") # relative to report dir  (to be used in html)
    js_file_curr = paste(report_dir, "/distanceplot_perchr.", sample_name, ".js", sep="") # relative to current dir
    write.table(js_content, js_file_curr,quote=FALSE, row.names=FALSE, col.names=FALSE)

    # tags for html to call js file
    html_script_tag = paste("<script src=\"", js_file, "\"></script>",sep="")   # relative to report dir
    html_common_script_tag = "<script src=\"https://d3js.org/d3.v3.js\"></script><script src=\"https://rawgit.com/SooLee/d3forNozzleR/master/interactive_multiline.js\"></script>"
    html_div_tag = paste("<div id=\"", div_id, "\"></div>",sep="")

    minihtml = paste("<html><body>",html_common_script_tag,html_div_tag,html_script_tag,"</body></html>",sep="")

    # return necessary things to create an html component
    return(list(div_id=div_id, tsvcolfile=tsvcolfile, js_file=js_file, html_script_tag = html_script_tag, html_common_script_tag=html_common_script_tag, html_div_tag=html_div_tag, minihtml=minihtml))
}

plot_orientation_proportion_vs_distance <- function(x, RE_len, xlim=c(2,5), no_xlabel=FALSE){
    matplot(x$distance,x[,c('proportion.Inner','proportion.Outer','proportion.Right','proportion.Left')],pch=19,type='o',xlab="", ylab="Proportion",lwd=1,lty=1, col=COLOR4(), axes=F, xlim=xlim)
    exp_axis(xlim,1)
    axis(2)
    box()
    if(no_xlabel==FALSE) mtext("distance",side=1,line=2.5)
    legend("topright",c('Inner','Outer','Right','Left'),col=COLOR4(),bty="n",pch=19,lwd=1,lty=1)
    if(RE_len==4) { abline(v=3.5, lty=2, col="grey") # 4-cutters at 3kb
    } else { abline(v=4.5, lty=2, col="grey") } # 6-cutters at 30kb
}


plot_orientation_log10count_vs_distance <- function(x, RE_len, xlim=c(2,5), no_xlabel=FALSE){
    ind = which(x$distance>=xlim[1] & x$distance<=xlim[2] & x$log10count.Inner!=-100 & x$log10count.Outer!=-100 & x$log10count.Right!=-100 & x$log10count.Left!=-100)
    if(length(ind)>0){
        ylim=range(x[ind,c('log10count.Inner','log10count.Outer','log10count.Right','log10count.Left')])
    }else{
        ylim=range(x[,c('log10count.Inner','log10count.Outer','log10count.Right','log10count.Left')])
    }
    matplot(x$distance,x[,c('log10count.Inner','log10count.Outer','log10count.Right','log10count.Left')],pch=19,type='o',xlab="", ylab="Contact frequency",lwd=1,lty=1,ylim=ylim, col=COLOR4(), axes=F, xlim=xlim)
    exp_axis(xlim,1)
    exp_axis(ylim,2)
    box()
    if(no_xlabel==FALSE) mtext("distance",side=1,line=2.5)
    legend("bottomright", c('Inner','Outer','Right','Left'), col=COLOR4(), bty="n", pch=19, lwd=1, lty=1)
    if(RE_len==4) { abline(v=3.5, lty=2, col="grey") # 4-cutters at 3kb
    } else { abline(v=4.5, lty=2, col="grey") } # 6-cutters at 30kb
}


plot_contact_probability_vs_distance <- function(x, xlim=c(3.2,8), no_xlabel=FALSE) {
    ylim=range(x$log10prob)
    plot(x$distance,x$log10prob,pch=19,type='o',xlab="", ylab="Contact probability",lwd=1,lty=1,ylim=ylim,axes=F,xlim=xlim)
    exp_axis(xlim,1)
    exp_axis(ylim,2)
    box()
    if(no_xlabel==FALSE) mtext("distance",side=1,line=2.5)

    abline(v=4, lty=2, col="grey80")
    abline(v=5.5, lty=2, col="grey80")

    # slope
    tad=which(x$distance>=4 & x$distance<=5.5)
    xx=x$distance[tad]
    yy=x$log10prob[tad]
    tad_slope=glm(yy~xx)$coefficients[2]
    tad_intercept=glm(yy~xx)$coefficients[1]
    spacing=abs(yy[1]*0.07)
    text(xx[1], yy[1]+spacing, paste("slope=",round(tad_slope,2),sep=' '), pos=4)
    abline(a=tad_intercept + spacing, b=tad_slope, lty=2)

    return(round(tad_slope,2));
}


# entropy
entropy<-function(xx)-sum(xx*log2(xx))
plot_entropy <- function(x, xlim=c(2,5), plt=c(0.2,0.95,0.2,0.95)){
    par(las=1,lend=2)
    par(plt=plt)
    entropies=apply(x[,c('proportion.Inner','proportion.Outer','proportion.Right','proportion.Left')],1,entropy)
    plot(x$distance, entropies,type='o',pch=19,ylab="entropy",xlab="distance",xlim=xlim,axes=F)
    exp_axis(xlim,1)
    axis(2)
}


plot_sd_with_cutoff <- function(x, xlim=c(2,5)){
    x=x[which(x$distance>=xlim[1] & x$distance<=xlim[2]),];
    sds=apply(x[,c('proportion.Inner','proportion.Outer','proportion.Right','proportion.Left')],1,sd)
    plot(x$distance, log10(sds), type='o',pch=19,ylab="sd",xlab="distance",xlim=xlim,axes=F, col=COLOR3B()[1])
    exp_axis(xlim,1)
    exp_axis(log10(sds), 2, 3)
    box()
    almost_converged = which(sds<0.02)
    well_converged = which(sds<0.005)
    points(x$distance[almost_converged], log10(sds[almost_converged]), col=COLOR3B()[2], pch=19)
    points(x$distance[well_converged], log10(sds[well_converged]), col=COLOR3B()[3], pch=19)
    abline(v=x[almost_converged[1] ,"distance"],col=COLOR3B()[2],lty=2)
    abline(v=x[well_converged[1] ,"distance"],col=COLOR3B()[3],lty=2)
    legend("topright",c("sd<0.02","sd<0.005"),col=COLOR3B()[2:3],pch=19,bty='n')

    ## determination for convergence
    if(RE_len==4) { conv_ind = which(x$distance>=3.5) } else { conv_ind = which(x$distance>=4.5) }
    if(length(setdiff(conv_ind, well_converged))==0) { return('Very Good');
    } else if(length(setdiff(conv_ind, almost_converged))==0) { return('Good');
    } else { return('Not converged'); }

}


plot_contact_frequency_vs_genomic_separation_per_chr <- function(x) {
    y=data.frame(x[,grep('log10prob_per_chr',colnames(x))])
    valid = which(x$distance > 3.5 & x$distance < 7)
    y=y[valid,]
    ylim = range(apply(y,2,function(xx){ ind =which(xx>-90); if(length(ind)==0) return(NA) else return( range(xx[ind])) }), na.rm=T)
    L = ncol(y)
    chrnames = sub("log10prob_per_chr.","",colnames(y))

    plt0 = par('plt', no.readonly=FALSE)
    plt1 = plt0; plt2= plt0
    pltx.plotend = (plt0[2] - plt0[1]) * 7/8 + plt0[1]
    pltx.legendstart = (plt0[2] - plt0[1]) * 29/32 + plt0[1]
    plt1[2] = pltx.plotend
    plt2[1] = pltx.legendstart
    par(plt=plt1)
    matplot(x$distance[valid],y,ylim=ylim, pch=19, type='l', lty=1, lwd=1, col=rainbow_w_offset(L), xlab="genomic separation", ylab="Contact frequency", axes=F)
    exp_axis(x$distance[valid], 1)
    exp_axis(ylim, 2)
    box()
    abline(v=4, lty=2, col="grey80")
    abline(v=5.5, lty=2, col="grey80")
    par(new=T, plt=plt2)
    plot.new(); plot.window(c(0,1),c(0,1), xaxs='i', yaxs='i')
    legend("topleft", chrnames, col = rainbow_w_offset(L), lty=1, lwd=1, bty='n' ,xpd=T, cex=0.8)
}

prep_plot_contact_frequency_vs_genomic_separation_per_chr_d3 <- function(x, sample_name, report_dir, xmin=3.5, xmax=7) {
    ycolumns = colnames(x)[grep('log10prob_per_chr',colnames(x))]
    y = data.frame(x[,ycolumns])
    valid = which(x$distance > xmin & x$distance < xmax)
    y=y[valid,,drop=FALSE]
    ylim = range(apply(y,2,function(xx){ ind =which(xx>-90); if(length(ind)==0) return(NA) else return( range(xx[ind])) }), na.rm=T)
    ymin= ylim[1]
    ymax = ylim[2]
    L = ncol(y)
    chrnames = sub("log10prob_per_chr.","",colnames(y))
    colors = substr(str_to_lower(rainbow_w_offset(L)),1,7)
    xlab = "log10 genomic separation"
    ylab = "log10 contact frequency"
    div_id = paste("d3div_contact_frequency_vs_genomic_separation_per_chr_", sample_name, "___", sep="")
    js_file_prefix = "distanceplot_perchr"
    tsvfile = paste("./", sample_name, ".plot_table.out", sep="")
    tsvcoldata = data.frame(columnx="distance", columny=ycolumns, color=colors, offcolor="#999999", name=chrnames, stringsAsFactors=FALSE)
    return( create_d3_js_for_interactive_multiline_plot(sample_name, xmin, xmax, ymin, ymax, div_id, xlab, ylab, report_dir, js_file_prefix, tsvfile, tsvcoldata) )
}


########################


require( "Nozzle.R1" )

generate_pairsqc_report <- function ( sample_name = NA, d3_prep_res = NULL) {

    if(is.na(sample_name)) {
        report_title = "PairsQC Report"
    } else {
        report_title = paste("PairsQC Report for sample", sample_name, sep=" ")
    }

    # Phase 1: create report elements
    r <- newCustomReport( report_title );
    d3script <- newHtml (d3_prep_res[[1]]$html_common_script_tag);
    s1 <- newSection( "Summary table" );
    s2 <- newSection( "Proportion of read orientation versus genomic separation");
    s3 <- newSection( "Number of reads versus genomic separation, stratified by read orientation" );
    s4 <- newSection( "Contact probability versus genomic separation" );
    s5 <- newSection( "Contact probability versus genomic separation, per chromosome" );
    references <- newSection( "References" );

    # References
    refRao = newCitation( authors= 'Rao et al.', title= 'A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping', publication= 'Cell', issue= '159:7', pages='1665-1680', year='2014' )
    refSanborn = newCitation( authors= 'Sanborn et al.', title= 'Chromatin extrusion explains key features of loop and domain formation in wild-type and engineered genome', publication= 'Proc. Natl. Acad. Sci. USA', issue= '112:47', pages='E6456-E6465', year='2015' )
    refImakaev = newCitation( authors= 'Imakaev et al.', title= 'Iterative correction of Hi-C data reveals hallmarks of chromosome organization', publication= 'Nature Methods', issue= '9:10', pages='999-1003', year='2012' )

    # section 1
    table_files = list.files('.', pattern='*summary.out$')
    sample_names = sub('.?summary.out','',table_files, perl=T, fixed=F)
    tableData1 = do.call(data.frame, sapply(table_files, function(table_file){
        tableData=read.table(table_file,stringsAsFactors=F, header=F, sep="\t")
        colnames(tableData) = c('QC field','value')
        return(tableData)
    }, simplify=F))[,c(1,seq(2,2*length(table_files),2))]
    colnames(tableData1) = c('QC field', sample_names)

    t1 <- newTable( tableData1, significantDigits=0, exportId="TABLE_1", "Summary Table" );
    p1 <- newParagraph( "Cis-to-trans ratio is computed as the ratio of long-range cis reads (>20kb) to trans reads plus long-range cis reads. Typically, cis-to-trans ratio higher than 40% is required. Percentage of long-range cis reads is the ratio of long-range cis reads to total number of reads. Minimum 15% is required and 40% or higher suggests a good library", asReference( refRao ), ". Convergence is determined as standard deviation of proportions of four read orientations to be <0.002 (Very good) or <0.05 (Good) (See below section ", asStrong("Proportion of read orientation versus genomic separation"), "). The slope of log10 contact probability vs distance between 10kb ~ 300kb representing TAD is also provided as well. (See below section ", asStrong("Contact probability versus genomic separation"), ".)");

    # section 2
    # using sample_names from section 1 (same ordering)
    png2_files = paste('plots/', sample_names, '.proportion.png',sep='')
    pdf2_files = gsub('png$','pdf',png2_files)
    f2list <- mapply(function(png2, pdf2, sample_name) {
                        newFigure( png2, fileHighRes=pdf2, exportId=paste("FIGURE_2", sample_name, sep="."),
                                    paste("Proportion of read orientation versus genomic separation : ", asStrong(sample_name), sep=""));
                    }, png2_files, pdf2_files, sample_names, SIMPLIFY=FALSE)
    eq2a = '10^3.5'
    eq2b = '10^4.5'
    p2 <- newParagraph("Contact frequency (number of reads, left) and proportion of reads (right) are shown, stratified by read orientation. Good four-cutter and six-cutter samples would converge at ~3kb (", eq2a, ") and ~30kb (", eq2b, "), respectively", asReference( refRao ), ". Convergence is determined by the standard deviation of the proportions being less than 0.005.")

    # section 3
    #png3= 'plots/log10counts.png'
    #pdf3= 'plots/log10counts.pdf'
    #f3 <- newFigure( png3, fileHighRes=pdf3, exportId="FIGURE_3", "Number of reads versus genomic separation, stratified by read orientation");

    # section 4
    # using sample_names from section 1 (same ordering)
    png4_files = paste('plots/', sample_names, '.log10prob.png',sep='')
    pdf4_files = gsub('png$','pdf',png4_files)
    f4list <- mapply(function(png4, pdf4, sample_name) {
                        newFigure( png4, fileHighRes=pdf4,
                                    paste("Contact probability versus genomic separation : ", asStrong(sample_name), sep=""));
                    }, png4_files, pdf4_files, sample_names, SIMPLIFY=FALSE)
    p4 <- newParagraph( "Contact probability (number of reads, normalized by number of bins and bin size) is shown with respect to genomic separation between mates", asReference( refImakaev), asReference( refSanborn ), ". The slope between distance 10kb ~ 300kb (10^4 ~ 10^5.5) representing a TAD is calculated. A good mitotic sample would have a slope close to ~ -0.76", asReference( refSanborn ), "." );

    # section 5
    # this section has d3 elements
    ss5divlist <- mapply(function(sn, k) newHtml(paste("<div class='figure'><div class='caption'><p><strong>Interactive Figure&nbsp;",k,".&nbsp;</strong>Contact probability versus genomic separation, per chromosome : ", asStrong(sn), "</div><div id='", d3_prep_res[[sn]]$div_id, "'></div></div>",d3_prep_res[[sn]]$html_script_tag, sep="")), sample_names, 1:length(sample_names), SIMPLIFY=FALSE);

    # Phase 2: assemble report structure bottom-up
    s1 <- addTo( s1, p1, t1);
    s2 <- addTo( s2, p2);
    for(f2 in f2list) { s2 <- addTo( s2, f2); }
    #s3 <- addTo( s3, f3);
    s4 <- addTo( s4, p4 );
    for(f4 in f4list) { s4 <- addTo( s4, f4); }
    # s5 <- addTo( s5, f5);
    for(i in 1:length(sample_names)) { s5 <- addTo( s5, ss5divlist[[i]] ); }
    references <- addTo( references, refRao, refImakaev, refSanborn )
    #r <- addTo( r, s1, s2, s3, s4, references );
    #r <- addTo( r, s1, s2, s4, s5, references );
    r <- addTo( r, d3script, s1, s2, s4, s5, references );

    # Phase 3: render report to file
    writeReport( r, filename="pairsqc_report" ); # w/o extension

}

##################

cwd = getwd()
sample_names = gsub('.?plot_table.out$', '', list.files(report_dir, pattern='*.plot_table.out$'), perl=T, fixed=F)

plot_for_sample<-function(sample_name, report_dir) {
    setwd(cwd)
    plot_dir = paste(report_dir,"plots",sep="/")
    plot_table_file = paste(report_dir, "/", sample_name, ".plot_table.out", sep="")

    x=read.table(plot_table_file,sep="\t",stringsAsFactors=F,header=T)
    write.csv(x[, grepl("distance|proportion", colnames(x))], file.path(report_dir, paste0(sample_name, ".distance.vs.proportion.csv")), row.names=FALSE)
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    setwd(plot_dir)
    res = pngpdf_preset( list(
        function() plot_orientation_log10count_vs_distance(x, RE_len),
        function() plot_orientation_proportion_vs_distance(x, RE_len, no_xlabel=T),
        function() plot_sd_with_cutoff(x)
        ), paste(sample_name, "proportion", sep="."), stylefunc0, get_preset(3), add.date=FALSE
    )
    convergence_determinant = res[[3]];

    # pngpdf_preset( list(
    #    function() plot_contact_probability_vs_distance(x),
    #    function() plot_contact_frequency_vs_genomic_separation_per_chr(x)
    #    ), paste(sample_name, "log10prob", sep="."), stylefunc0, get_preset(2,'h50'), add.date=FALSE
    #)

    res = pngpdf_preset( list(
            function() plot_contact_probability_vs_distance(x)
        ), paste(sample_name, "log10prob", sep="."), stylefunc0, get_preset(1), add.date=FALSE
    )
    slope = res[[1]]

    #pngpdf.nodate( function()plot_contact_probability_vs_distance(x) ,"log10prob")
    #pngpdf.nodate( function()plot_contact_frequency_vs_genomic_separation_per_chr(x), "log10prob_chr", height=5.5 )
    #pngpdf.nodate( function()plot_entropy(x), "entropy", height=4)

    # d3 plot
    setwd(cwd)
    d3plot_prep_res = prep_plot_contact_frequency_vs_genomic_separation_per_chr_d3(x, sample_name, report_dir)

    ## printing cis-trans stats and post-plot stats to summary.out
    res_all = as.data.frame(rbind(c("convergence",convergence_determinant), c("slope", slope)), stringsAsFactors=FALSE)
    colnames(res_all)=c('V1','V2')
    cis_file = paste(report_dir, '/', sample_name, ".cis_to_trans.out", sep="");
    res_file = paste(report_dir, '/', sample_name, ".summary.out", sep="");
    original_res = read.table(cis_file, stringsAsFactors=FALSE, header=FALSE, sep="\t");
    colnames(original_res)=c('V1','V2')
    res_all = rbind(original_res, res_all)
    write.table(res_all, res_file, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
    setwd(plot_dir)

    return(d3plot_prep_res)
}

d3_prep_res = sapply(sample_names, plot_for_sample, report_dir = report_dir, simplify=FALSE)
print(d3_prep_res)

setwd(cwd)
setwd(report_dir)
generate_pairsqc_report(d3_prep_res = d3_prep_res)
