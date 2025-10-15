library(tidyverse)

dir_name <- "high_seed_low_everything"

db_names <- c("Dry_seed_n_all_Stages", "Dry_seed", "Seeds_Stage_3_10_w", "Seeds_Stage_3_5_w", "Seeds_Stage_6_8_w", "Seeds_Stage_9_10_w")

setwd(paste0("C:/Users/YonatanY/Migal/Rachel Amir Team - General/Arabidopsis_db/seed_specific_expression_genes/", dir_name))

for (i_db in db_names) {
    fn <- paste0("ExpressionAnglerOutput/ExpressionAnglerOutput_", i_db, ".txt")

    # Read as tab-delimited; keep original column names
    x <- read.delim(fn, check.names = FALSE, stringsAsFactors = FALSE)

    # Drop metadata rows (e.g., "#tissue")
    x <- x[!grepl("^#", x[[1]]), ]
    names(x)[1:2] <- c("ID", "Info")

    # Parse the Info column: "r_value|GeneChip_ID|Annotation|..."
    split <- strsplit(x$Info, "\\|")
    x$r_value <- as.numeric(trimws(vapply(split, `[`, "", 1)))
    x$GeneChip_ID <- trimws(ifelse(lengths(split) >= 2, vapply(split, `[`, "", 2), NA))
    x$Annotation <- trimws(ifelse(lengths(split) >= 3, vapply(split, function(z) paste(z[-(1:2)], collapse = "|"), ""), NA))

    # Pull TAIR (AGI) locus from ID, GeneChip_ID, or Annotation (first match wins)
    tair_pat <- "AT[1-5MC]G\\d{5}(?:\\.\\d+)?"
    grab_agi <- function(id, chip, ann) {
        for (s in c(id, chip, ann)) {
            if (!is.na(s) && grepl(tair_pat, s, ignore.case = TRUE)) {
                return(toupper(regmatches(s, regexpr(tair_pat, s, ignore.case = TRUE))))
            }
        }
        NA_character_
    }
    tair_id <- mapply(grab_agi, x$ID, x$GeneChip_ID, x$Annotation)

    # Keep only TAIR ID + r_value, drop NAs
    out <- data.frame(tair_id = tair_id, r_value = x$r_value, stringsAsFactors = FALSE)
    out <- out[!is.na(out$tair_id) & !is.na(out$r_value), ]

    # Save + preview
    write.csv(out, paste0(i_db, "_tair_rValues.csv"), row.names = FALSE)
    head(out)
}
