#' @export
check_input_file <- function(filepath, file_type = c("cnv", "ped", "annot"), sep = "\t") {
	file_type <- match.arg(file_type)
	
	if (!file.exists(filepath)) {
		return(list(status = FALSE, msg = paste0("❌ File not found: ", filepath)))
	}
	
	# Read header
	header <- names(read.table(filepath, header = TRUE, sep = sep, nrows = 1, check.names = FALSE))
	header_lower <- tolower(header)
	
	# Define accepted column variants per file type
	schema_list <- list(
		cnv = list(
			ordered = list(
				chr   = c("chr", "chrom", "seqnames"),
				start = c("start", "begin", "pos_start"),
				end   = c("end", "stop", "pos_end")
			),
			unordered = list(
				sample = c("sampleid", "sample_id", "id_sample"),
				type   = c("type", "svtype", "variant_type")
			)
		),
		ped = list(
			ordered = list(
				iid = c("iid", "individual", "sampleid", "childid", "id")
			),
			unordered = list(
				fid = c("fid", "family", "familyid", "fam_id"),
				father = c("pid", "father", "dad", "fatherid"),
				mother = c("mid", "mother", "mom", "motherid")
			)
		),
		other = list(
			ordered = list(),
			unordered = list()
		)
	)
	
	schema <- schema_list[[file_type]]
	
	# --- Check ordered columns ---
	if (length(schema$ordered) > 0) {
		required_n <- length(schema$ordered)
		if (length(header_lower) < required_n) {
			return(list(status = FALSE, msg = paste0("❌ File has fewer than ", required_n, " columns.")))
		}
		
		ordered_names <- names(schema$ordered)
		for (i in seq_along(ordered_names)) {
			colname <- ordered_names[i]
			variants <- schema$ordered[[colname]]
			if (!header_lower[i] %in% variants) {
				return(list(
					status = FALSE,
					msg = paste0("❌ Column ", i, " must be one of: ", paste(variants, collapse = ", "))
				))
			}
		}
	}
	
	# --- Check unordered columns ---
	if (length(schema$unordered) > 0) {
		missing_cols <- c()
		for (colname in names(schema$unordered)) {
			variants <- schema$unordered[[colname]]
			if (!any(header_lower %in% variants)) {
				missing_cols <- c(missing_cols, colname)
			}
		}
		if (length(missing_cols) > 0) {
			return(list(status = FALSE, msg = paste0("❌ Missing required columns: ", paste(missing_cols, collapse = ", "))))
		}
	}
	
	return(list(status = TRUE, msg = "✅ File check passed!"))
}

