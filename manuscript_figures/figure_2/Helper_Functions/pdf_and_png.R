library(ggplot2)
library(Cairo)

# Saves a PDF and PNG file for a provided ggplot object, with visually similar scaling

pdf_and_png <- function(p, outputDir, filename, pdfWidth = 7, pdfHeight = 7, splitDirectories = F, SAVE_RDS = F, SKIP_IF_EXISTS = F, scale = 1) {
    if (splitDirectories) {
        outputDir_png <- file.path(outputDir, "png")
        outputDir_pdf <- file.path(outputDir, "pdf")
        outputDir_rds <- file.path(outputDir, "rds")

        if (!dir.exists(outputDir_png)) {
            dir.create(outputDir_png, showWarnings = F)
        }
        if (!dir.exists(outputDir_pdf)) {
            dir.create(outputDir_pdf, showWarnings = F)
        }
        if (!dir.exists(outputDir_rds)) {
            dir.create(outputDir_rds, showWarnings = F)
        }

        pdf_path <- file.path(outputDir_pdf, paste0(filename, ".pdf"))
        png_path <- file.path(outputDir_png, paste0(filename, ".png"))
        rds_path <- file.path(outputDir_rds, paste0(filename, ".rds"))
    } else {
        pdf_path <- file.path(outputDir, paste0(filename, ".pdf"))
        png_path <- file.path(outputDir, paste0(filename, ".png"))
        rds_path <- file.path(outputDir, paste0(filename, ".rds"))
    }

    if (!SKIP_IF_EXISTS || !file.exists(pdf_path)) {
        ggsave(pdf_path, plot = p, width = pdfWidth, height = pdfHeight, dpi = 300, limitsize = F, scale = scale, device = cairo_pdf)
    } else {
        print(paste0("SKIPPING: ", pdf_path))
    }

    if (!SKIP_IF_EXISTS || !file.exists(png_path)) {
        ggsave(png_path, plot = p, width = pdfWidth, height = pdfHeight, dpi = 300, limitsize = F, scale = scale)
    } else {
        print(paste0("SKIPPING: ", png_path))
    }

    if (SAVE_RDS) {
        path_output <- list("pdf" = pdf_path, "png" = png_path, "rds" = rds_path)
        if (!SKIP_IF_EXISTS || !file.exists(rds_path)) {
            saveRDS(p, rds_path)
        }
    } else {
        path_output <- list("pdf" = pdf_path, "png" = png_path)
    }

    return(path_output)
}
