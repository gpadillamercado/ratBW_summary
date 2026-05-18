# Extracting Animal Body Weight Data
# First download the unformatted Text File:
# Recommendations for & Documentation of Biological Values for Use in Risk Assistance 	1988
library(ggplot2)
full_doc <- readLines("500022JL.txt")
full_doc <- full_doc[!(full_doc %in% "")] # Get rid of empty lines
# Split in pages
pgsplit <- grep("-------\f", full_doc)
splitidx <- cumsum(seq_along(full_doc) %in% pgsplit) + 1L
doc <- split(full_doc, splitidx)
rm(full_doc, pgsplit, splitidx)

# FUNCTIONS----------------------------
get_between <- function(x, ix) {
  bix <- ix[1L:length(ix) - 1L] + c(1L, 1L, 2L, 2L, 2L, 1L, 1L)
  eix <- ix[2L:length(ix)] - 1L
  stopifnot(length(bix) == length(eix))
  lst <- lapply(seq_along(bix), \(j) {
    x[bix[j]:eix[j]]
  })
  names(lst) <- c("Species", "Sex", "N", "Age", "Weight", "Variance", "Reference")
  return(lst)
}

rem_nonum <- function(x) {
  nonum <- grep("^(\\d|\\.)", x, invert = TRUE)
  nonum <- min(nonum[nonum > 35]) - 1L
  x[seq(1, nonum)]
}

join_by_decimal <- function(x) {
  x[grep("^NS$", x)] <- NA
  this_pre <- grep("^\\w$", x)
  this_post <- grep("^\\.", x)
  stopifnot(length(this_pre) == length(this_post))
  x[this_pre] <- paste0(x[this_pre], x[this_post])
  x <- x[-this_post]
  return(x)
}

join_by_space <- function(x, lhs, rhs) {
  this_lhs <- which(x %in% lhs)
  this_rhs <- which(x %in% rhs) # Usually later in the vector
  stopifnot(length(this_lhs) == length(this_rhs))
  x[this_lhs] <- paste(x[this_lhs], x[this_rhs])
  x <- x[-this_rhs]
  return(x)
}

find_table_cols <- function(x) {
  x <- x[!(x %in% "")]
  idx <- grep(cols, x)
  idx <- c(idx, length(x))
  idx[3] <- idx[3] - 1L # "No. of", "Animals"
  # Sometimes... Species is not right next to where this column starts
  if (!grepl("^\\w{3,}", x[idx[1]])) {
    idx[1] <- min(grep("^\\w{3,}", x[idx[1]:idx[2]])) + 1L
  }
  xlist <- get_between(x, idx)
  return(xlist)
}

fill_down <- function(x, val, len) {
  stopifnot(length(x) < len)
  ol <- length(x) + 1L
  x[ol:len] <- val
  return(x)
}

shiftswap_vec <- function(x, start, size, shift) {
  xlen <- length(x)
  stopifnot(start + size + shift - 1 <= xlen)
  rshift <- x[start:(start + size - 1)]
  print(rshift)
  plen <- length(rshift)
  lshift <- x[(start + size):(start + size + shift - 1)]
  x[start:(start + shift - 1)] <- lshift
  x[(start + shift):(start + shift + size - 1)] <- rshift
  x
}

# RATS--------------------------------
# For instance, I know Rat data is on page 110
rat_docs <- doc[110:128]
# Looks like Strain Prefix starts at index 117
# but let's find things programmatically
cols <- "(Species|Sex|Animals|Age|(W|H)eight|Variance|Reference)$"
# Evaluate the potential rows of data
check_data_rows <- function(x) {
  x <- x[!(x %in% "")]
  idx <- grep(cols, x)
  # Sometimes... Species is not right next to where this column starts
  idx <- c(idx, length(x))
  idx[3] <- idx[3] - 1L
  diff(idx)
}


tlst <- lapply(rat_docs, find_table_cols)
lenlst <- lapply(tlst, lengths) # Get actual row lengths here
rowlens <- c(42, 43, 42, 44, 41, 43, 43, 45, 44, 42, 41, 44, 38, 44, 44, 43, 42, 44, 32)

c("Species", "Sex", "N", "Age", "Weight", "Variance", "Reference")

# Manual Changes
tlst[["110"]][["Weight"]] <- join_by_decimal(tlst[["110"]][["Weight"]])
tlst[["110"]][["Variance"]] <- join_by_decimal(rem_nonum(tlst[["110"]][["Variance"]]))
tlst[["112"]][["Species"]] <- join_by_space(tlst[["112"]][["Species"]], "August", "28807/Cr")
tlst[["112"]][["N"]] <- fill_down(tlst[["112"]][["N"]], "1", 42)
tlst[["112"]][["Weight"]] <- join_by_decimal(tlst[["112"]][["Weight"]])
tlst[["112"]][["Variance"]] <- join_by_decimal(tlst[["112"]][["Variance"]]) # Note. Variance gets shifted. Fix in post
tlst[["112"]][["Variance"]] <- shiftswap_vec(tlst[["112"]][["Variance"]], 14, 6, 8)
tlst[["112"]][["Reference"]] <- tlst[["112"]][["Reference"]][1:42]
tlst[["114"]][["Reference"]] <- tlst[["114"]][["Reference"]][1:41]
tlst[["116"]][["Reference"]] <- join_by_space(
  tlst[["116"]][["Reference"]],
  c("Cameron et al."), "1985"
)
tlst[["120"]][["Reference"]] <- tlst[["120"]][["Reference"]][1:41]
tlst[["124"]][["Reference"]] <- tlst[["124"]][["Reference"]][1:44]
tlst[["126"]][["Species"]] <- join_by_space(
  tlst[["126"]][["Species"]],
  c(
    "Hlstar/Furth", "Htstar/Lewts",
    "Ulstar/Lewls", "Hlstar/Leuls",
    "Hlstar/Lewls"
  ),
  "Cr"
)
tlst[["126"]][["Weight"]] <- join_by_decimal(tlst[["126"]][["Weight"]])
tlst[["126"]][["Variance"]] <- join_by_decimal(
  fill_down(
    rem_nonum(tlst[["126"]][["Variance"]]),
    ".E-00X", length(rem_nonum(tlst[["126"]][["Variance"]])) + 4L
  )
)
tlst[["127"]][["Species"]] <- tlst[["127"]][["Species"]][2:45]
tlst[["127"]][["Reference"]] <- join_by_space(
  tlst[["127"]][["Reference"]],
  c("Deyl et al.", "Leong et al"),
  c("1975", ", 1964")
)
tlst[["128"]][["Species"]] <- tlst[["128"]][["Species"]][2:33]
tlst[["128"]][["Reference"]] <- join_by_space(
  tlst[["128"]][["Reference"]][1:39],
  c("Leong et al..", "Wlberg et al.", "Ulberg et al."),
  c("1964", ", 1966")
)

all(sapply(tlst, \(x) length(unique(lengths(x))) == 1L))

rat_dflst <- lapply(tlst, as.data.frame)
rat_df <- do.call(rbind, rat_dflst)
rat_df$Variance <- trimws(gsub(" ", "", rat_df$Variance))
rat_df$Variance[grep("^NS$", rat_df$Variance)] <- NA
rat_df$Species[grep("^NS$", rat_df$Species)] <- NA

## FIX SPECIES NAMES
unique(rat_df$Species)

rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c("9935/Cr")),
  "ACP 9935/Cr"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c("AC! 9935/Cr", "AC1 9935/Cr")),
  "AC1 9935/Cr"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c("August 28B07/Cr")),
  "August 28807/Cr"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c("BH/Cr")),
  "BN/Cr"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c("Cpb-.WU", "Cpb:UU")),
  "Cpb:WU"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c("Cr:MGAPS(OH)", "CrtHGAPS(OH)", "Cr:HGAPS(ON)", "Cr:HGAPS(OH)")),
  "Cr:MGAPS(OM)"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c("Cr.-RAR(SD)")),
  "Cr:RAR(SD)"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c("Fischer F344 .")),
  "Fischer F344"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% paste0("F334/C", c("M Lov", "r 1 Lov", "H Lov"))),
  "F334/Crl Lov"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% "F 334/Cr 1 Lov"),
  "F334/Crl Lov"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c("Osborne-Hendel")),
  "Osborne-Mendel"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% paste0("Sprague-Dawl", c(" ey/HCr", "ey/HCr", "ey/RCr"))),
  "Sprague-Dawley/MCr"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c("S5B/PlCr", "SSB/PICr", "SSB/PlCr")),
  "S5B/P1Cr"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c("Hcrypt/Ztm", "Hcrypt/Ztra")),
  "Wcrypt/Ztm"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c(
    "Hlstar/Furth Cr"
  )),
  "Wistar/Furth Cr"
)
rat_df$Species <- replace(
  rat_df$Species,
  which(rat_df$Species %in% c(
    "Htstar/Lewts Cr",
    "Wlstar/Lewls Cr",
    "Ulstar/Lewls Cr",
    "Hlstar/Leuls Cr",
    "Hlstar/Lewls Cr"
  )),
  "Wistar/Lewis Cr"
)
rat_df$Species <- replace(
  rat_df$Species,
  grep("(W|H|U)lstar$", rat_df$Species),
  "Wistar"
)
rat_df$Species <- replace(
  rat_df$Species,
  grep("Yoshlda", rat_df$Species),
  "Yoshida/Cr"
)
unique(rat_df$Species)

## FIX References
rat_df$Reference <- gsub("\\. (?=\\d{4})", ", ", rat_df$Reference, perl = TRUE)
rat_df$Reference <- gsub(" ,", ",", rat_df$Reference, perl = TRUE)
unique(rat_df$Reference)
rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% "1972"),
  "Poiley, 1972"
)
rat_df$Reference <- replace(
  rat_df$Reference,
  grep("Po.{4,7}1972", rat_df$Reference),
  "Poiley, 1972"
)
rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% "Inazu et al."),
  "Inazu et al., 1984"
)
rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% c("HAS, 1971", "MAS, 1971")),
  "NAS, 1971"
)
rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% "Horrlssey and Norred, 1984"),
  "Morrissey and Norred, 1984"
)
rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% c("Nauderley, 1986", "nauderly, 1986", "Hauderly, 1986", "Nauderly, 1986")),
  "Mauderly, 1986"
)
rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% c("Harrlman, 1969b", "HarMman, 1969b")),
  "Harriman, 1969b"
)
rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% "Holloszy and Smith, 1985"),
  "Holloszy and Smith, 1986"
)
rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% "Holloszy and Smith, 1985"),
  "Holloszy and Smith, 1986"
)
rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% "Rlos et al., 1986a"),
  "Rios et al., 1986a"
)
rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% c("Rlos et al., 1986b", "Rlos ct al., 1986b")),
  "Rios et al., 1986b"
)
rat_df$Reference <- gsub(
  "Hortola",
  "Mortola",
  rat_df$Reference
)
rat_df$Reference <- gsub(
  "Horden",
  "Worden",
  rat_df$Reference
)
rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% "Pond et al."),
  "Pond et al, 1985"
)

rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% "Quasi et al., 1983"),
  "Quast et al., 1983"
)

rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% c("a!., 1985", "al., 1985")),
  "Alt et al., 1985"
)

rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% c("Wlberg et al., 1966", "Ulberg et al., 1966")),
  "Wiberg et al., 1966"
)

rat_df$Reference <- replace(
  rat_df$Reference,
  which(rat_df$Reference %in% "Leong et al, 1964"),
  "Leong et al., 1964"
)

unique(rat_df$Reference)

### Make Sex lowercase and abbreviate
rat_df$Sex[which(rat_df$Sex %in% "NS")] <- "both"
rat_df$Sex <- tolower(rat_df$Sex)

## FIX Variance (letters read instead of numbers)
rat_df$Variance <- gsub("Z", "2", rat_df$Variance)
rat_df$Variance <- gsub("B", "8", rat_df$Variance)
rat_df$Variance <- gsub("O", "0", rat_df$Variance)
rat_df$Variance <- gsub("l", "1", rat_df$Variance)
rat_df$Variance <- gsub("I", "1", rat_df$Variance)
rat_df$Variance <- gsub("D", "0", rat_df$Variance)
rat_df$Variance <- gsub("S", "5", rat_df$Variance)
rat_df$Variance <- gsub("F", "E", rat_df$Variance)
rat_df$Variance <- gsub("r,", "5", rat_df$Variance)
rat_df$Variance <- gsub("a", "8", rat_df$Variance)
rat_df$Variance <- replace(
  rat_df$Variance,
  which(rat_df$Variance %in% "2.E-00X"),
  "2.29E-004"
)
rat_df$Variance <- replace(
  rat_df$Variance,
  which(rat_df$Variance %in% "6.E-00X"),
  "6.44E-004"
)
rat_df$Variance <- replace(
  rat_df$Variance,
  which(rat_df$Variance %in% "1.E-00X"),
  "1.16E-004"
)
rat_df$Variance <- replace(
  rat_df$Variance,
  which(rat_df$Variance %in% "3.E-00X"),
  "3.15E-004"
)

## Fix minor things in other numeric columns
rat_df$Age <- gsub("\\. ", "", rat_df$Age, perl = TRUE)
rat_df$Age <- gsub("B", "8", rat_df$Age)
rat_df$N <- gsub("B", "8", rat_df$N)
rat_df$N <- gsub("\\?", "2", rat_df$N, perl = TRUE)
rat_df$N[grep("^NS$", rat_df$N)] <- NA

## TYPE the columns
rat_df$Age <- as.integer(rat_df$Age)
rat_df$N <- as.integer(rat_df$N)
rat_df$Weight <- as.numeric(rat_df$Weight)
rat_df$Variance <- as.numeric(rat_df$Variance)

# rat_df$Variance[is.na(as.numeric(rat_df$Variance))] # This is to check for non-numeric coercion
rat_df |>
  # subset(grepl("Sprague", Species)) |>
  ggplot(aes(as.numeric(Age), as.numeric(Weight), color = Sex)) +
  geom_point() +
  labs(x = "Age (days)", y = "Mass (kg)") +
  geom_vline(xintercept = c(21, 56)) +
  # facet_wrap(~Species) +
  theme_bw() +
  theme(legend.position = "")

rat_df <- tibble::tibble(rat_df)
readr::write_tsv(rat_df, "ratBW_data1.tsv")
