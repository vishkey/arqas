bTree <- function() {
    obj <- list()
    oldClass(obj) <- "bTree"
    return(obj)
}

"add<-" <- function(x, value) {UseMethod("add<-", x)}
"add<-.bTree" <- function(x, value) {
    if (length(x) == 0)  {
      x <- list(node=value, left=bTree(), right=bTree())
      oldClass(x) <- "bTree"
    }
    else {
      if (value >= x$node)
        add(x$right) <- value
      else
        add(x$left) <- value
    }
    return(x)
}