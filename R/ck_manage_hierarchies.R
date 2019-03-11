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
##' @param node_labs_new a character vector of new node names replacing the existing node names given in
##' argument \code{node_labs} when using \code{ck_rename_nodes}.
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
ck_create_node <- function(total_lab="Total", node_labs=NULL) {
  return(hier_create(root = total_lab, nodes = node_labs))
}

##' @rdname ck_manage_hierarchies
##' @export
ck_add_nodes <- function(node, node_labs, reference_node) {
  return(hier_add(tree = node, root = reference_node, nodes = node_labs))
}

##' @rdname ck_manage_hierarchies
##' @export
ck_delete_nodes <- function(node, node_labs, reference_node = NULL) {
  return(hier_delete(tree = node, nodes = node_labs))
}

##' @rdname ck_manage_hierarchies
##' @export
ck_rename_nodes <- function(node, node_labs, node_labs_new) {
  nn <- node_labs_new
  names(nn) <- node_labs
  return(hier_rename(tree = node, nodes = nn))
}
