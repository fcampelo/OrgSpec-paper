library(dplyr)
library(epitopes)
library(seqinr)

# Predictions table: O. volvulus
X     <- readRDS("./output/analysis.rds")
prots <- readRDS("../00_general_datasets/00_proteins_20201007.rds")
  
X$mypreds$Seq <- mapply(
  function(prot, st, en){
    substr(prots$TSeq_sequence[prots$UID == prot], st, en)
  },
  X$mypreds$Info_UID,
  X$mypreds$start_pos,
  X$mypreds$end_pos)

names(X$mypreds) <- c("Protein", "Start pos", "End pos", "Length", "Average Prob.", "Sequence")
X$mypreds[2:4] <- lapply(X$mypreds[2:4], as.integer)
X <- X$mypreds %>%
  dplyr::filter(`Average Prob.` >= 0.75) %>%
  dplyr::mutate(CandidateID = paste0("Pred", sprintf("%02d", seq_along(Protein))))

# Get positive observations (O. volvulus epitopes) from IEDB
epits <- readRDS("../00_general_datasets/00_epitopes_20201006.rds") %>%
  epitopes::filter_epitopes(orgIDs = 6282, 
                            tax_load_file = "../00_general_datasets/00_taxonomy_20201007.rds") %>%
  as.data.frame() %>%
  dplyr::filter(n_Positive > n_Negative)

# prepare and run BLAST
BLAST_path <- "./data/BLAST"
if(!dir.exists(BLAST_path)) dir.create(BLAST_path)
fpr <- paste0(BLAST_path, "/preds.fasta")
ftr <- paste0(BLAST_path, "/epits.fasta")
seqinr::write.fasta(as.list(X$Sequence), 
                    names = X$CandidateID,
                    file.out = fpr, as.string = TRUE, nbchar = 100000)
seqinr::write.fasta(as.list(epits$epit_seq), 
                    names = epits$epitope_id,
                    file.out = ftr, as.string = TRUE, nbchar = 100000)

system(paste0("makeblastdb -in ", ftr, " -dbtype prot"))
system(paste0("blastp -query ", fpr, " -db ", ftr, " -seg no ",
              "-outfmt '6 qseqid sseqid length qlen slen nident pident qcovhsp mismatch gaps qstart qend sstart send evalue score' > ",
              fpr, "-BLAST"))

blast <- read.csv(paste0(fpr, "-BLAST"), sep = "\t",
                  header = FALSE,
                  stringsAsFactors = FALSE)
names(blast) <- c("QueryID", "SubjectID", "Alignment_length",
                  "Query_length", "Subject_length", "Num_matches",
                  "Perc_identity", "Query_coverage", "Num_mismatches",
                  "Num_gaps", "Query_match_start", "Query_match_end",
                  "Subject_match_start", "Subject_match_end",
                  "E_value", "Score")

blast <- blast %>%
  filter(Perc_identity > 80) %>%
  arrange(desc(Perc_identity)) %>%
  mutate(Match = paste0(SubjectID, " (", Perc_identity, "\\%)")) %>%
  group_by(QueryID) %>%
  summarise(Matches = paste(Match, collapse = "; ")) %>%
  ungroup()

breaklines <- function(x, n=20){
  l   <- nchar(x)
  str <- seq(1, l, by = n)
  stp <- pmin(str + n - 1, l)
  br <- mapply(substr, 
               start = str, 
               stop = stp, 
               MoreArgs = list(x = x))
  return(paste(br, collapse = "\\newline "))
}

X <- X %>%
  left_join(blast, by = c("CandidateID" = "QueryID")) %>%
  dplyr::select(-CandidateID) %>%
  dplyr::mutate(Sequence = sapply(Sequence, breaklines, n = 30))

kableExtra::kable(X, format = "latex", longtable = TRUE, digits = 2, escape = FALSE)
