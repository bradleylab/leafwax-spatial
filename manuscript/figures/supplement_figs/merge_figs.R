library(pdftools)
library(magick)

# Set the input and output
input_dir <- "."  # update this
output_pdf <- "Figures_all.pdf"

# Get list of PDF files
pdf_files <- list.files(input_dir, pattern = "\\.pdf$", full.names = TRUE)
pdf_files <- sort(pdf_files)

# Temporary folder to store labeled PDFs
temp_dir <- file.path(tempdir(), "labeled_pages")
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

# Initialize vector to store labeled PDF paths
labeled_files <- c()

for (pdf_file in pdf_files) {
  fname <- basename(pdf_file)
  message("Processing: ", fname)
  
  # Convert first (and only) page to image
  img <- tryCatch({
    image_read_pdf(pdf_file, density = 150)
  }, error = function(e) {
    message("Failed to read: ", fname)
    return(NULL)
  })
  
  if (!is.null(img)) {
    # Annotate with filename
    img_annotated <- image_annotate(
      img, text = fname, gravity = "southwest",
      size = 15, color = "black", location = "+10+10"
    )
    
    # Save to temporary PDF file
    output_path <- file.path(temp_dir, fname)
    tryCatch({
      image_write(img_annotated, output_path, format = "pdf")
      labeled_files <- c(labeled_files, output_path)
    }, error = function(e) {
      message("Failed to write: ", fname)
    })
  }
}

# Final check
if (length(labeled_files) == 0) {
  stop("No labeled PDFs were created. Check input files and ImageMagick config.")
}

# Method 1: Try using pdf_combine from pdftools (safer than system calls)
tryCatch({
  message("Attempting to merge PDFs using pdftools...")
  pdf_combine(labeled_files, output = output_pdf)
  cat("✅ Merged PDF saved as:", output_pdf, "\n")
}, error = function(e) {
  message("pdftools merge failed, trying alternative method...")
  
  # Method 2: Use qpdf with proper quoting
  tryCatch({
    # Build command with proper quoting
    cmd <- sprintf("qpdf --empty --pages %s -- \"%s\"", 
                   paste(shQuote(labeled_files), collapse = " "), 
                   output_pdf)
    system(cmd)
    cat("✅ Merged PDF saved as:", output_pdf, "\n")
  }, error = function(e2) {
    message("qpdf also failed, trying one more method...")
    
    # Method 3: Create a response file for qpdf (handles special characters better)
    response_file <- file.path(temp_dir, "file_list.txt")
    writeLines(labeled_files, response_file)
    
    cmd <- sprintf("qpdf --empty --pages @\"%s\" -- \"%s\"", 
                   response_file, output_pdf)
    system(cmd)
    cat("✅ Merged PDF saved as:", output_pdf, "\n")
  })
})