
#' Load and preprocess C. elegans neuron-level data
#'
#' Reads the raw neuron-level CSV, splits the combined `Treatment` field
#' into `Generation`, `Treatment`, and `Retreatment`, constructs a unique
#' worm identifier, converts categorical variables to factors.
#'
#' @return A data frame with columns: `Replicate`, `Treatment`, `NeuronScore`,
#'   `UniqueWormID`, `Retreatment`, `Generation`, `logNeuronScore`, `Worm`.
#'
#' @references
#' [Bergemann et al. (2026) 'Progeny effects of rotenone exposure
#' depend on parental toxicity', Toxicological Sciences, Volume 209, Issue 3, 
#' March 2026, kfag011](https://doi.org/10.1093/toxsci/kfag011)
#'
#' @export
#' 
get_data <- function() {

  # Load data
  neuron <- neuron_level_data
  
  neuron$UniqueWormID <- paste0(neuron$Replicate, "_", neuron$Treatment, "_", neuron$Worm)
  neuron_ <- separate(neuron, col = Treatment, into = c("Generation", "Treatment", "Retreatment"), sep = "_")
  
  # Subset data
  raw_data <- neuron_[, c("Replicate", "Treatment", "NeuronScore", "UniqueWormID", "Retreatment", "Generation")]
  
  # Convert categorical variables to factors
  raw_data$UniqueWormID<- as.factor(raw_data$UniqueWormID)
  raw_data$Treatment <- as.factor(raw_data$Treatment)
  raw_data$Replicate <- as.factor(raw_data$Replicate)
  raw_data$Retreatment <- as.factor(raw_data$Retreatment)
  raw_data$Generation <- as.factor(raw_data$Generation)
  
  # Set reference level for retreatment to be control
  raw_data$Retreatment <- relevel(raw_data$Retreatment, ref = "control")
  
  # Transformed variables
  raw_data$NeuronScore <- raw_data$NeuronScore + 1
  raw_data$logNeuronScore <- log(as.numeric(raw_data$NeuronScore))
  raw_data$Worm <- as.numeric(sub(".*_", "", raw_data$UniqueWormID))
  
  return(raw_data)
}

#' Build model matrices for the Bayesian rank-based isotonic regression
#'
#' Constructs the response vector and design matrices required by
#' `fit_BayesRank` from preprocessed neuron-level data. `X` encodes a
#' cell-means representation of the parental early Treatment
#' ("Low" and "High vs. Low") crossed with Generation and Retreatment.
#' `Q` is the worm-level random-effect incidence matrix, and
#' `W` is the replicate (batch) design matrix.
#'
#' @param raw_data Preprocessed data frame, as returned by `get_data`.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{y}{Ordinal neuron damage score.}
#'     \item{X}{Design matrix of three-way interaction terms
#'      (parental early Treatment crossed with Generation and Retreatment}
#'     \item{W}{Replicate design matrix.}
#'     \item{Q}{Worm-level random-effect incidence matrix.}
#'   }
#'
#' @export
#' 
get_data_BayesRank <- function(raw_data) {

  # Damage score
  y <- raw_data$NeuronScore 
  
  # Create design matrix for treatment effects
  X <- matrix(0, nrow = dim(raw_data)[1], ncol = 2 * 4)
  for (i in 1:dim(raw_data)[1]) {
    if (raw_data$Generation[i] == "F1" & raw_data$Retreatment[i] == "control") {
      if (raw_data$Treatment[i] == "003uM" | raw_data$Treatment[i] == "050uM") X[i, 1] <- 1
      if (raw_data$Treatment[i] == "050uM") X[i, 2] <- 1
    } else if (raw_data$Generation[i] == "P0" & raw_data$Retreatment[i] == "control") {
      if (raw_data$Treatment[i] == "003uM" | raw_data$Treatment[i] == "050uM") X[i, 3] <- 1 
      if (raw_data$Treatment[i] == "050uM") X[i, 4] <- 1
    } else if (raw_data$Generation[i] == "F1" & raw_data$Retreatment[i] == "25uM") {
      if (raw_data$Treatment[i] == "003uM" | raw_data$Treatment[i] == "050uM") X[i, 5] <- 1
      if (raw_data$Treatment[i] == "050uM") X[i, 6] <- 1 
    } else if (raw_data$Generation[i] == "P0" & raw_data$Retreatment[i] == "25uM") {
      if (raw_data$Treatment[i] == "003uM" | raw_data$Treatment[i] == "050uM") X[i, 7] <- 1
      if (raw_data$Treatment[i] == "050uM") X[i, 8] <- 1
    }
  }
  colnames(X) <- c(
    "Low, F1, Un-Rechallenged",  "High vs. Low, F1, Un-Rechallenged",
    "Low, P0, Un-Rechallenged",  "High vs. Low, P0, Un-Rechallenged",
    "Low, F1, Rechallenged",     "High vs. Low, F1, Rechallenged",
    "Low, P0, Rechallenged",     "High vs. Low, P0, Rechallenged")

  # Create design matrix for worms random effects
  Q <- model.matrix(y ~ -1 + raw_data$UniqueWormID)
  
  # Create design matrix for replicates batch effects
  W <- model.matrix(y ~ -1 + raw_data$Replicate)
  colnames(W) <- c("Replicate 1", "Replicate 2", "Replicate 3")
  
  return(list(y = y, X = X, W = W, Q = Q))
}

#' Build model matrix and data frame for the CLMM comparison model
#'
#' Constructs the fixed-effect design matrix for the cumulative link mixed
#' model (`clmm`) benchmark. Response variables and the worm-level
#' grouping factor are appended for use with `ordinal::clmm`.
#'
#' @param raw_data Preprocessed data frame, as returned by `get_data`.
#'
#' @return A data frame with the fixed-effect dummy-coding columns for Replicate 
#'   membership and interactions as in `Treatment:Retreatment:Generation`.
#'   For the three-way interactions, the columns involving the control Treatment group
#'   serve as the implicit within-stratum reference, and are thus dropped.
#'   Additional info include NeuronScore, logNeuronScore, factNeuronScore, and
#'   UniqueWormID.
#'
#' @export
#' 
get_data_clmm <- function(raw_data){
  raw_data_multinom <- model.matrix(~ -1 + Replicate + Treatment:Retreatment:Generation, data=raw_data)
  raw_data_multinom <- data.frame(raw_data_multinom[, -c(6, 9, 12, 15)]) # remove all interactions involving control groups
  raw_data_multinom$NeuronScore     <- raw_data$NeuronScore
  raw_data_multinom$logNeuronScore  <- raw_data$logNeuronScore
  raw_data_multinom$factNeuronScore <- as.factor(raw_data$NeuronScore)
  raw_data_multinom$UniqueWormID    <- raw_data$UniqueWormID
  raw_data_multinom
}