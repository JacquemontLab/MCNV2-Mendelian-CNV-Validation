#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

params <- getOption("MCNV2.params", default = list())

bedtools_path <- params$bedtools_path %||% "bedtools" 
results_dir   <- params$results_dir   %||% tempdir()

# Define server logic required to draw a histogram
function(input, output, session) {
	#https://stackoverflow.com/questions/60150956/attaching-python-script-while-building-r-package
	#system.file("python", "compute_inheritance.py", package = "MCNV2")
	#options(bedtools.path = "/opt/homebrew/bin/bedtools") -> getOption("bedtools.path")
	
	cnv_check <- reactiveVal(FALSE)
	ped_check <- reactiveVal(FALSE)
	
	observeEvent(input$cnv_tsv, {
		req(input$cnv_tsv)  # Attend que le fichier soit chargé
		
		filepath <- input$cnv_tsv$datapath
		
		# Essaie de valider le fichier
		ret <- check_input_file(filepath, required_cols = c("SampleID", "Chr", "Type"))
		
		cnv_check(ret$status)
		if(cnv_check()){
			output$cnv_tsv_status <- renderText(ret$msg)
		} else {
			output$cnv_tsv_status <- renderText(ret$msg)
		}
	})
	
	observeEvent(input$ped_tsv, {
		req(input$ped_tsv)  # Attend que le fichier soit chargé
		
		filepath <- input$ped_tsv$datapath
		
		# Essaie de valider le fichier
		ret <- check_input_file(filepath, required_cols = NULL)
		
		ped_check(ret$status)
		if(ped_check()){
			output$ped_tsv_status <- renderText(ret$msg)
		} else {
			output$ped_tsv_status <- renderText(ret$msg)
		}
		
	})
	
	observeEvent(cnv_check() | ped_check(), {
		if (cnv_check() & ped_check()) {
			updateActionButton(session, "submit_preprocess", disabled = FALSE)
		} else {
			updateActionButton(session, "submit_preprocess", disabled = TRUE)
		}
	})
	
	observeEvent(input$submit_preprocess, {
		req(input$cnv_tsv)  # Attend que le fichier soit chargé
		
		cnvs_file <- input$cnv_tsv$datapath
		
		gene_annotation_script <- system.file("python",
																					"gene_annotation.py", 
																					package = "MCNV2")
		gene_resource_file <- system.file("resources",
																			"gene_resources.tsv", 
																			package = "MCNV2")
		prob_regions_file <- system.file("resources",
																		 "problematic_regions_GRCh38.bed", 
																		 package = "MCNV2")
		output_file <- file.path(results_dir, "cnvs_annotated_by_genes.tsv")
		
		cmd = paste("python3", gene_annotation_script, "--cnv", cnvs_file,
		"--gene_resource", gene_resource_file, "--prob_regions", prob_regions_file,
		"--out", output_file, 
		"--genome_version", input$build, "--bedtools_path", bedtools_path)

		ret <- system(command = cmd, intern = FALSE)
		
		if(ret == 0){
			output$annot_tsv_status <- renderText(output_file)
			output$preview_preproc_tbl <- renderDataTable({
				dat <- readr::read_tsv(output_file, show_col_types = FALSE)
				datatable(dat, options = list(dom = 't', scrollX = TRUE))
			})
		} else {
			output$annot_tsv_status <- renderText("❌ File annotation failed. Check logs")
		}

		
	})
	
	# observeEvent(input$submit_preprocess, {
	# 	output$preview_preproc_tbl <- renderDataTable({
	# 		txt <- readr::read_tsv("~/workspace/projects/Support/BIOINF-219/data/merged_WGSCNV_trios_30.tsv")
	# 		datatable(data.frame(message = txt), options = list(dom = 't', scrollX = TRUE))
	# 	})
	# })
}
