//
// This file holds several functions specific to the workflow/hicar.nf in the nf-core/hicar pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowHicar {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        genomeExistsError(params, log)

        if (!params.gtf && !params.gff) {
            Nextflow.error("No GTF or GFF3 annotation specified! The pipeline requires at least one of these files.")
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def citation_text = [
                "Tools used in the workflow included:",
                "BWA (Li H et al. 2009)",
                "BEDTools (Quinlan et al. 2010)",
                "ChIPpeakAnno (Zhu et al. 2010)",
                "circos (Krzywinski et al. 2009)",
                "cooler (Abdennur et al. 2020)",
                "cooltools (Nezar et al.)",
                "diffhic (Lun et al. 2015)",
                "edgeR (Robinson et al. 2010)",
                "FastQC (Andrews 2010),",
                "GenMap (Pockrandt et al. 2020)",
                "HiC-DC+ (Sahin et al. 2021)",
                "HiCExplorer (Ramírez et al. 2018)",
                "Homer (Heinz et al. 2010)",
                "igv.js (Thorvaldsdóttir et al. 2013)",
                "Juicer_tools (Durand et al. 2016)",
                "Kraken 2 (Wood et al. 2019)",
                "MACS2 (Zhang et al. 2008)",
                "MAPS (Juric et al. 2019)",
                "MultiQC (Ewels et al. 2016)",
                "pairsqc (Soo Lee)",
                "pairtools (Nezar et al. 2023)",
                "peakachu (Salameh et al. 2020)",
                "SAMtools (Li et al. 2009)",
                "trackViewer (Ou et al. 2019)",
                "UCSC tools (Kent et al. 2010)",
                "."
            ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def reference_text = [
                "<li>Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics. 2009 Jul 15;25(14):1754-60. doi: 10.1093/bioinformatics/btp324. Epub 2009 May 18. PubMed PMID: 19451168; PubMed Central PMCID: PMC2705234.</li>",
                "<li>Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2. doi: 10.1093/bioinformatics/btq033. Epub 2010 Jan 28. PubMed PMID: 20110278; PubMed Central PMCID: PMC2832824.</li>",
                "<li>Zhu LJ, Gazin C, Lawson ND, Pagès H, Lin SM, Lapointe DS, Green MR. ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics. 2010 May 11;11:237. doi: 10.1186/1471-2105-11-237. PMID: 20459804; PMCID: PMC3098059.</li>",
                "<li>Krzywinski M, Schein J, Birol I, Connors J, Gascoyne R, Horsman D, Jones SJ, Marra MA. Circos: an information aesthetic for comparative genomics. Genome Res. 2009 Sep;19(9):1639-45. doi: 10.1101/gr.092759.109. Epub 2009 Jun 18. PMID: 19541911; PMCID: PMC2752132.</li>",
                "<li>Abdennur N, Mirny LA. Cooler: scalable storage for Hi-C data and other genomically labeled arrays. Bioinformatics. 2020 Jan 1;36(1):311-316. doi: 10.1093/bioinformatics/btz540. PMID: 31290943.</li>",
                "<li>Lun AT, Smyth GK. diffHic: a Bioconductor package to detect differential genomic interactions in Hi-C data. BMC Bioinformatics. 2015 Aug 19;16:258. doi: 10.1186/s12859-015-0683-0. PMID: 26283514; PMCID: PMC4539688.</li>",
                "<li>Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010 Jan 1;26(1):139-40. doi: 10.1093/bioinformatics/btp616. Epub 2009 Nov 11. PMID: 19910308; PMCID: PMC2796818.</li>",
                "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
                "<li>Pockrandt C, Alzamel M, Iliopoulos CS, Reinert K. GenMap: ultra-fast computation of genome mappability. Bioinformatics. 2020 Jun 1;36(12):3687-3692. doi: 10.1093/bioinformatics/btaa222. PMID: 32246826; PMCID: PMC7320602.</li>",
                "<li>Sahin M, Wong W, Zhan Y, Van Deynze K, Koche R, Leslie CS. HiC-DC+ enables systematic 3D interaction calls and differential analysis for Hi-C and HiChIP. Nat Commun. 2021 Jun 7;12(1):3366. doi: 10.1038/s41467-021-23749-x. PMID: 34099725; PMCID: PMC8184932.</li>",
                "<li>Ramírez F, Bhardwaj V, Arrigoni L, Lam KC, Grüning BA, Villaveces J, Habermann B, Akhtar A, Manke T. High-resolution TADs reveal DNA sequences underlying genome organization in flies. Nat Commun. 2018 Jan 15;9(1):189. doi: 10.1038/s41467-017-02525-w. PMID: 29335486; PMCID: PMC5768762.</li>",
                "<li>Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432</li>",
                "<li>Thorvaldsdóttir H, Robinson JT, Mesirov JP. Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration. Brief Bioinform. 2013 Mar;14(2):178-92. doi: 10.1093/bib/bbs017. Epub 2012 Apr 19. PMID: 22517427; PMCID: PMC3603213.</li>",
                "<li>Durand NC, Shamim MS, Machol I, Rao SS, Huntley MH, Lander ES, Aiden EL. Juicer Provides a One-Click System for Analyzing Loop-Resolution Hi-C Experiments. Cell Syst. 2016 Jul;3(1):95-8. doi: 10.1016/j.cels.2016.07.002. PMID: 27467249; PMCID: PMC5846465.</li>",
                "<li>Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2. Genome Biol. 2019 Nov 28;20(1):257. doi: 10.1186/s13059-019-1891-0. PMID: 31779668; PMCID: PMC6883579.</li>",
                "<li>Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS. Model-based analysis of ChIP-Seq (MACS). Genome Biol. 2008;9(9):R137. doi: 10.1186/gb-2008-9-9-r137. Epub 2008 Sep 17. PubMed PMID: 18798982; PubMed Central PMCID: PMC2592715.</li>",
                "<li>Juric I, Yu M, Abnousi A, Raviram R, Fang R, Zhao Y, Zhang Y, Qiu Y, Yang Y, Li Y, Ren B, Hu M. MAPS: Model-based analysis of long-range chromatin interactions from PLAC-seq and HiChIP experiments. PLoS Comput Biol. 2019 Apr 15;15(4):e1006982. doi: 10.1371/journal.pcbi.1006982. PMID: 30986246; PMCID: PMC6483256.</li>",
                "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
                "<li>Salameh TJ, Wang X, Song F, Zhang B, Wright SM, Khunsriraksakul C, Ruan Y, Yue F. A supervised learning framework for chromatin loop detection in genome-wide contact maps. Nat Commun. 2020 Jul 9;11(1):3428. doi: 10.1038/s41467-020-17239-9. PMID: 32647330; PMCID: PMC7347923.</li>",
                "<li>Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PubMed PMID: 19505943; PubMed Central PMCID: PMC2723002.</li>",
                "<li>Ou J, Zhu LJ. trackViewer: a Bioconductor package for interactive and integrative visualization of multi-omics data. Nat Methods. 2019 Jun;16(6):453-454. doi: 10.1038/s41592-019-0430-y. PMID: 31133757.</li>",
                "<li>Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics. 2010 Sep 1;26(17):2204-7. doi: 10.1093/bioinformatics/btq351. Epub 2010 Jul 17. PubMed PMID: 20639541; PubMed Central PMCID: PMC2922891.</li>",
            ].join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        //meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        //meta["tool_bibliography"] = toolBibliographyText(params)


        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }
}
