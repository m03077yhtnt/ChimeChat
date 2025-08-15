# silence R CMD check notes from dplyr's non-standard evaluation
utils::globalVariables(c(
  ".data", "score", "condition", "human_symbol", "mouse_symbol"
))
