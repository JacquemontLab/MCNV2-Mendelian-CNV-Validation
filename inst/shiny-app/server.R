#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

# Define server logic required to draw a histogram
function(input, output, session) {
	#https://stackoverflow.com/questions/60150956/attaching-python-script-while-building-r-package
	#system.file("python", "compute_inheritance.py", package = "MCNV2")
	output$distPlot <- renderPlot({
		
		# generate bins based on input$bins from ui.R
		x    <- faithful[, 2]
		bins <- seq(min(x), max(x), length.out = input$bins + 1)
		
		# draw the histogram with the specified number of bins
		hist(x, breaks = bins, col = 'darkgray', border = 'white',
				 xlab = 'Waiting time to next eruption (in mins)',
				 main = 'Histogram of waiting times')
		
	})
	
	
	observeEvent(input$submit_preprocess, {
		output$preview_preproc_tbl <- renderDataTable({
			txt <- readr::read_tsv("~/workspace/projects/Support/BIOINF-219/data/merged_WGSCNV_trios_30.tsv")
			datatable(data.frame(message = txt), options = list(dom = 't', scrollX = TRUE))
		})
	})
}
