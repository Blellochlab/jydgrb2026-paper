# scRNA-seq Pre-processing and MULTI-seq Demultiplexing
# Author: Deniz Goekbuget

# --- 1. Setup ---
rm(list = ls())
library(Seurat)
library(deMULTIplex2)
library(readxl)
library(ggrastr)

# Define relative paths for the repository structure
data_dir <- "./data"
fastq_dir <- "./fastq"
out_dir <- "./results"

if(!dir.exists(out_dir)) dir.create(out_dir)

# --- 2. Define Functions ---

# Function A: Demultiplex MULTI-seq Barcodes
get_multiseq_bcs <- function(bc_fastq_folder, cellranger_folder, multiseq_bcs, bcfastq_prefix, output_folder) {
  
  # Load reference oligos and filter based on metadata input
  data("multiseq_oligos")
  bcs <- multiseq_oligos[multiseq_oligos %in% multiseq_bcs]
  
  # Get cell barcodes from CellRanger GEX output
  gex <- Read10X(cellranger_folder)
  gex <- CreateSeuratObject(counts = gex, project = bcfastq_prefix, min.cells = 3, min.features = 200)
  cell_bcs <- colnames(gex)
  cell_bcs <- as.character(gsub("-1", "", cell_bcs))
  
  # Load barcode reads
  read_table <- readTags(dir = bc_fastq_folder,
                         name = bcfastq_prefix,
                         barcode.type = "MULTIseq",
                         assay = "RNA",
                         filter.cells = cell_bcs)
  
  # Align and Demultiplex
  tag_mtx <- alignTags(read_table, bcs)
  
  res <- demultiplexTags(tag_mtx,
                         plot.path = output_folder,
                         plot.name = paste0(bcfastq_prefix, "demux2_"),
                         plot.diagnostics = TRUE,
                         min.tag.show = 50)
  
  # Export classifications
  saveRDS(res$final_assign, file.path(output_folder, paste0(bcfastq_prefix, "classifications.rds")))
  return(res$final_assign)
}

# Function B: Initialize Seurat Object and Apply Demultiplexing
initialize_seurat_object <- function(cellranger_folder, samplesheet, name, bcfastq_prefix, output_folder) {
  
  # Load 10x Data
  infile <- Read10X(cellranger_folder)
  dat <- CreateSeuratObject(counts = infile, project = name, min.cells = 3, min.features = 200)
  
  # Remove human cells (isolating murine TILs)
  mouse_cells <- !grepl("GRCh38-", rownames(dat))
  dat <- dat[mouse_cells, ]
  
  # Clean up formatting for matching
  colnames(dat) <- gsub("-1", "", colnames(dat))
  rownames(dat) <- gsub("mm10---", "", rownames(dat))
  
  # Load deMULTIplex2 classifications generated in Step A
  class_file <- file.path(output_folder, paste0(bcfastq_prefix, "classifications.rds"))
  classifications <- readRDS(class_file)
  classifications <- data.frame(classifications)
  classifications$sample <- NA
  
  # Map barcodes to sample names using the provided metadata
  data("multiseq_oligos")
  ind <- match(classifications$classifications, names(multiseq_oligos))
  ind2 <- match(multiseq_oligos[ind], samplesheet$`BC Sequence`)
  classifications$sample <- samplesheet$Sample[ind2]
  
  # Apply classifications to Seurat metadata
  dat[["multiBC"]] <- ifelse(ind, classifications[colnames(dat), ]$classifications, "Empty")
  dat[["sample"]] <- ifelse(ind, classifications[colnames(dat), ]$sample, "Empty")
  
  # Filter out negative, multiplet, and empty droplets
  valid_cells <- dat[["multiBC"]] != "negative" & 
    dat[["multiBC"]] != "multiplet" & 
    dat[["multiBC"]] != "Empty"
  
  dat <- dat[, valid_cells]
  
  # Calculate QC metrics
  dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^mt-")
  
  # Save processed output
  out_path <- file.path(output_folder, paste0(name, "-nodups-seurat.rds"))
  saveRDS(dat, out_path)
  
  cat("\nSuccessfully saved filtered Seurat object to:", out_path, "\n")
  return(dat)
}


# --- 3. Execution Pipeline ---

# Define Library/Project details
library_name <- "jyrb1"
prefix <- "JYRB1-BC_"

# 3.1 Load Metadata
# Assumes 'libraries.xlsx' has columns: `Sample` and `BC Sequence`
meta <- read_excel(file.path(data_dir, "libraries.xlsx"), sheet = "jyrb1")

# Define input paths (Relative to repository)
path_bc_fastq <- file.path(fastq_dir, "bc")
path_cellranger <- file.path(data_dir, "JYRB1-GEX/outs/filtered_feature_bc_matrix/")

# 3.2 Run Demultiplexing
cat("Starting MULTI-seq demultiplexing...\n")
classifications <- get_multiseq_bcs(
  bc_fastq_folder = path_bc_fastq,
  cellranger_folder = path_cellranger,
  multiseq_bcs = meta$`BC Sequence`,
  bcfastq_prefix = prefix,
  output_folder = out_dir
)

# 3.3 Initialize and Filter Seurat Object
cat("Starting Seurat initialization and filtering...\n")
seurat_obj <- initialize_seurat_object(
  cellranger_folder = path_cellranger,
  samplesheet = meta,
  name = library_name,
  bcfastq_prefix = prefix,
  output_folder = out_dir
)