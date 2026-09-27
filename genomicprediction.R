###############################################################################
# GBLUP for two traits, autotetraploid dosage data, rrBLUP
#
# Detects genotype file orientation and ID columns automatically by matching
# against the phenotype IDs, so annotation columns and markers-in-rows layouts
# are handled without editing the script.
###############################################################################

library(rrBLUP)
write.csv(genotypes, file="genotypes.csv")
genotest <- read.csv("genotypes.csv")
## ---- Settings --------------------------------------------------------------
GENO_FILE  <- "genotypes.csv"
PHENO_FILE <- "phenotypes.csv"
PED_FILE   <- "pedigree.csv"

PHENO_ID_COL <- "id"
TRAIT1 <- "DryWeight"
TRAIT2 <- "SpecificGravity"

PED_SIRE_COL <- "father"
PED_DAM_COL  <- "mother"

FIXED_EFFECTS <- NULL     # e.g. c("env", "rep")
MAF_MIN <- 0.05
PLOIDY  <- 4
W1 <- 0.5                 # selection index weights
W2 <- 0.5

## ---- Helpers ---------------------------------------------------------------
norm_id  <- function(x) trimws(as.character(x))
loose_id <- function(x) tolower(gsub("[^A-Za-z0-9]", "", as.character(x)))

## ---- Phenotypes ------------------------------------------------------------
pheno <- as.data.frame(read.csv(PHENO_FILE, check.names = FALSE,
                                stringsAsFactors = FALSE))
names(pheno) <- trimws(names(pheno))

for (nm in c(PHENO_ID_COL, TRAIT1, TRAIT2, FIXED_EFFECTS))
  if (!nm %in% names(pheno))
    stop("Column '", nm, "' not in phenotype file. Found: ",
         paste(names(pheno), collapse = ", "))

pheno$id <- norm_id(pheno[[PHENO_ID_COL]])
pheno[[TRAIT1]] <- as.numeric(pheno[[TRAIT1]])
pheno[[TRAIT2]] <- as.numeric(pheno[[TRAIT2]])
pheno_ids <- unique(pheno$id)

## ---- Genotypes: find the individuals automatically -------------------------
geno <- read.csv(GENO_FILE, check.names = FALSE, stringsAsFactors = FALSE)

# How many phenotype IDs appear in each column's values, and in the header?
col_hits  <- sapply(geno, function(x) sum(norm_id(x) %in% pheno_ids))
name_hits <- sum(norm_id(names(geno)) %in% pheno_ids)

if (name_hits > max(col_hits)) {
  # Markers in rows, individuals in columns. Keep only the sample columns,
  # which drops CHROM/POS/REF/ALT and any other annotation automatically.
  message("Markers in rows: transposing. Sample columns found: ", name_hits)
  keep <- norm_id(names(geno)) %in% pheno_ids
  M <- t(as.matrix(geno[, keep, drop = FALSE]))
  rownames(M) <- norm_id(names(geno))[keep]
} else {
  # Individuals in rows. The ID column is whichever matches the phenotypes.
  id_col <- which.max(col_hits)
  message("Individuals in rows: ID column is '", names(geno)[id_col],
          "' with ", max(col_hits), " matches")
  M <- as.matrix(geno[, -id_col, drop = FALSE])
  rownames(M) <- norm_id(geno[[id_col]])
}

suppressWarnings(storage.mode(M) <- "numeric")
M <- M[, colSums(!is.na(M)) > 0, drop = FALSE]   # drop all-NA (text) columns

# Duplicate sample IDs (repeated checks, re-runs) break the factor levels used
# later, so average replicate genotype calls into a single row per individual.
if (any(duplicated(rownames(M)))) {
  message("Averaging ", sum(duplicated(rownames(M))), " duplicate sample(s): ",
          paste(head(unique(rownames(M)[duplicated(rownames(M))]), 5),
                collapse = ", "))
  n_per <- table(rownames(M))
  M <- rowsum(M, rownames(M), na.rm = TRUE)
  M <- M / as.vector(n_per[rownames(M)])
}

if (nrow(M) == 0)
  stop("No genotype IDs matched the phenotypes.\n",
       "  Genotype row names: ", paste(head(rownames(M), 3), collapse = " | "), "\n",
       "  Genotype col names: ", paste(head(names(geno), 3), collapse = " | "), "\n",
       "  Phenotype IDs:      ", paste(head(pheno_ids, 3), collapse = " | "))

message("Genotype matrix: ", nrow(M), " individuals x ", ncol(M), " markers")
message("Dosage range: ", paste(range(M, na.rm = TRUE), collapse = " to "))

## ---- G matrix --------------------------------------------------------------
M <- apply(M, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
p <- colMeans(M) / PLOIDY
M <- M[, !is.na(p) & p > MAF_MIN & p < (1 - MAF_MIN), drop = FALSE]

p <- colMeans(M) / PLOIDY
W <- sweep(M, 2, PLOIDY * p)
G <- tcrossprod(W) / sum(PLOIDY * p * (1 - p))
G <- G + diag(1e-5, nrow(G))
rownames(G) <- colnames(G) <- rownames(M)

message("G: ", nrow(G), " x ", ncol(G),
        " | mean diagonal ", round(mean(diag(G)), 3), " (expect ~1)")

## ---- Align phenotypes to G -------------------------------------------------
pheno <- pheno[pheno$id %in% rownames(G), ]
if (nrow(pheno) == 0) stop("No phenotype record matches a genotyped individual.")
stopifnot(!any(duplicated(rownames(G))))
message("Records: ", nrow(pheno),
        " | genotyped individuals with no phenotype: ",
        length(setdiff(rownames(G), pheno$id)))

## ---- Fit -------------------------------------------------------------------
fit_one <- function(trait) {
  dat <- pheno[!is.na(pheno[[trait]]), ]
  if (nrow(dat) < 10) stop("Fewer than 10 usable records for ", trait)
  kin.blup(data = dat, geno = "id", pheno = trait,
           K = G, GAUSS = FALSE, fixed = FIXED_EFFECTS, PEV = TRUE)
}

fit1 <- fit_one(TRAIT1)
fit2 <- fit_one(TRAIT2)

## ---- Breeding values -------------------------------------------------------
ids  <- rownames(G)
pick <- function(v) v[match(ids, names(v))]

gebv <- data.frame(
  id     = ids,
  GEBV_1 = pick(fit1$g),
  GEBV_2 = pick(fit2$g),
  rel_1  = pmax(0, pmin(1, 1 - pick(fit1$PEV) / fit1$Vg)),
  rel_2  = pmax(0, pmin(1, 1 - pick(fit2$PEV) / fit2$Vg)),
  n_records = as.integer(table(factor(pheno$id, levels = ids))),
  row.names = NULL
)
gebv$index <- W1 * scale(gebv$GEBV_1) + W2 * scale(gebv$GEBV_2)
names(gebv)[2:5] <- c(paste0("GEBV_", TRAIT1), paste0("GEBV_", TRAIT2),
                      paste0("rel_", TRAIT1), paste0("rel_", TRAIT2))

gebv <- gebv[order(-gebv$index), ]
write.csv(gebv, "GEBV_all.csv", row.names = FALSE)

## ---- Parents ---------------------------------------------------------------
ped <- as.data.frame(read.csv(PED_FILE, check.names = FALSE,
                              stringsAsFactors = FALSE))
names(ped) <- trimws(names(ped))
ped <- ped[, names(ped) != "", drop = FALSE]

for (nm in c(PED_SIRE_COL, PED_DAM_COL))
  if (!nm %in% names(ped))
    stop("Pedigree column '", nm, "' not found. Found: ",
         paste(names(ped), collapse = ", "))

parents <- unique(norm_id(c(ped[[PED_SIRE_COL]], ped[[PED_DAM_COL]])))
parents <- parents[!parents %in% c("", "0", "NA", ".", "-", "*") & !is.na(parents)]
matched <- intersect(parents, rownames(G))

# Fall back to case/punctuation-insensitive matching if nothing lines up
if (length(matched) == 0) {
  lp <- loose_id(parents); lg <- loose_id(rownames(G))
  if (length(intersect(lp, lg)) > 0) {
    message("Parent IDs differ by case or punctuation: matching loosely")
    matched <- rownames(G)[lg %in% lp]
  } else {
    message("No parent IDs matched G.\n",
            "  Pedigree parents: ", paste(head(parents, 3), collapse = " | "), "\n",
            "  Genotype IDs:     ", paste(head(rownames(G), 3), collapse = " | "))
  }
}

write.csv(gebv[gebv$id %in% matched, ], "GEBV_parents.csv", row.names = FALSE)

## ---- Summary ---------------------------------------------------------------
cat("\nh2", TRAIT1, ":", round(fit1$Vg / (fit1$Vg + fit1$Ve), 3), "\n")
cat("h2", TRAIT2, ":", round(fit2$Vg / (fit2$Vg + fit2$Ve), 3), "\n")
cat("GEBV_all.csv    :", nrow(gebv), "individuals\n")
cat("GEBV_parents.csv:", sum(gebv$id %in% matched), "parents\n")


gebv_sg_raw <- data.frame(
  id   = names(fit2$g),
  GEBV = as.vector(fit2$g),
  PEV  = as.vector(fit2$PEV[names(fit2$g)]),
  Pred = fit2$pred,
  row.names = NULL
)
write.csv(gebv_yield_raw, "gebv_yield_raw.csv", row.names = FALSE)
write.csv(gebv_sg_raw, "gebv_sg_raw.csv", row.names = FALSE)


gebv_yield_raw <- gebv_yield_raw[order(-gebv_yield_raw$GEBV), ] 