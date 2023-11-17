# Internals, not doccumented
# TODO:: this is beyond moronic
getLinksBlurb <- function(x) {
  x <- x[which(x == "Links:    "):length(x)]
  x <- x[!(1:length(x) %% 2)]
  unlist(strsplit(x, split = "\\("))[1:(2*length(x)) %% 2]
}
