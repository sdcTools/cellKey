## helpers
# check if input is a Node object
check_for_node <- function(node) {
  stopifnot("nodedim" %in% class(node))
}
# check if a node exists in a tree
check_refnode_exists <- function(node, reference_node) {
  cur_node <- FindNode(node, reference_node)
  if (is.null(cur_node)) {
    stop("The provided reference node does not exist! Check your hierarchy!\n")
  }
  cur_node
}

# checking inputs
check_node_inputs <- function(node, node_labs, reference_node) {
  check_for_node(node)
  stopifnot(is_scalar_character(reference_node))
  stopifnot(is_character(node_labs))

  if (any(duplicated(node_labs))) {
    stop("duplicated values in argument 'node_labs' detected\n")
  }
}


##' @name ck_manage_hierarchies
##' @rdname ck_manage_hierarchies
##'
##' @title Create and modify the structure of hierarchies
##'
##' @description Functions \code{ck_create_nodes()}, \code{ck_add_nodes()}
##' and \code{ck_delete_nodes()} allow to define and modify hierarchical structures
##' represented as trees. These objects can be used in \code{\link{perturbTable}} to
##' define the (hierarchical) structure of tables.
##'
##' @param total_lab the name of the overall total (summation over all contributing)
##' @param node an object as created in \code{ck_create_node()} or returned from
##' \code{ck_add_nodes()} or \code{ck_delete_nodes()}.
##' @param node_labs character name(s) of new elements that should be inserted to or
##' deleted from a hierarchical structure
##' @param reference_node character name of an existing node in the hierarchical
##' structure. When using \code{ck_add_nodes()}, the new elements are created as children
##' of the reference node. In \code{ck_delete_nodes()}, all children of the reference node that
##' match the names with argument \code{node_labs} are deleted from the hierarchy.
##' @examples
##' dim <- ck_create_node(total_lab="Total")
##' dim <- ck_add_nodes(dim, reference_node="Total", node_labs=LETTERS[1:4])
##' print(dim)
##'
##' ## add some levels below "A" and "C"
##' dim <- ck_add_nodes(dim, reference_node="A", node_labs=paste0("a", 1:5))
##' dim <- ck_add_nodes(dim, reference_node="C", node_labs=paste0("c", 1:5))
##' print(dim)
##'
##' ## delete some specific levels
##' dim <- ck_delete_nodes(dim, reference_node="A", node_labs=c("a1", "a4"))
##' print(dim)
##'
##' ## delete entire subtree
##' dim <- ck_delete_nodes(dim, reference_node="Total", node_labs=c("C"))
##' print(dim)
##' # plot(dim)
##' @return a \code{Node} object that can be used to specify a hierarchy used
##' to define table inputs in \code{\link{perturbTable}}.
NULL


##' @rdname ck_manage_hierarchies
##' @export
ck_create_node <- function(total_lab="Total") {
  stopifnot(is_scalar_character(total_lab))
  node <- Node$new(total_lab)
  class(node) = c(class(node), "nodedim")
  node
}

##' @rdname ck_manage_hierarchies
##' @export
ck_add_nodes <- function(node, node_labs, reference_node) {
  check_node_inputs(node, node_labs, reference_node)
  cur_node <- check_refnode_exists(node, reference_node)

  for (lab in node_labs) {
    if (!is.null(FindNode(cur_node, lab))) {
      cat("Node",lab,"already exists under",shQuote(reference_node),"--> skipping\n")
    } else {
      cur_node$AddChild(lab)
    }
  }
  return(node)
}

##' @rdname ck_manage_hierarchies
##' @export
ck_delete_nodes <- function(node, node_labs, reference_node) {
  check_node_inputs(node, node_labs, reference_node)
  cur_node <- check_refnode_exists(node, reference_node)
  for (lab in node_labs) {
    if (is.null(FindNode(cur_node, lab))) {
      cat("Node",lab,"does not exists under",shQuote(reference_node),"--> nothing to delete\n")
    } else {
      Prune(cur_node, function(x) x$name != lab)
    }
  }
  return(node)
}

# internal function used in hierarchies before sdcTable
# --> structure required for sdcTable (already works)
node_to_sdcinput <- function(inp, addNumLevels=FALSE) {
  num_levels <- NULL
  check_for_node(inp)

  xx <- Traverse(inp)

  codes <- sapply(xx, function(x) x$name)
  levels <- sapply(xx, function(x) x$level)
  ats <- sapply(1:length(levels), function(x) {
    paste(rep("@", levels[x]), collapse="")
  })
  dt <- data.table(levels=ats, codes=codes)
  if (addNumLevels==TRUE) {
    dt[,num_levels:=levels]
  }
  dt[]
}
