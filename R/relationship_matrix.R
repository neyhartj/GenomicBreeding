H_mat <- function(A, G, w) {
  # ID of genotyped
  id <- row.names(G)
  # ID in Amat
  id1 <- row.names(A)
  # Ungenotyped in the Amat
  id2 <- setdiff(id1, id)
  n2 <- length(id2)
  A <- A[c(id, id2), c(id, id2)]
  A11 <- A[id, id]
  if (n2 > 0) {
    A.inv <- solve(A)
    A11.inv <- solve(A11)
  }
  nw <- length(w)
  H <- vector("list", nw)
  Hinv <- vector("list", nw)
  for (i in 1:nw) {
    Gw <- (1 - w[i]) * G + w[i] * A11
    if (n2 == 0) {
      H[[i]] <- Gw
    } else {
      Hinv[[i]] <- as(bdiag(solve(Gw) - A11.inv, Matrix(0, nrow = n2, ncol = n2)) + forceSymmetric(A.inv), "symmetricMatrix")
    }
  }
  id <- c(id, id2)

  H.list <- vector("list", nw)
  for (i in 1:nw) {
    if (!is.null(H[[i]])) {
      eigen.H <- eigen(H[[i]], symmetric = TRUE)
      eigen.H$vectors <- Matrix(eigen.H$vectors, dimnames = list(id, id))
    } else {
      eigen.H <- eigen(Hinv[[i]], symmetric = TRUE)
      eigen.H$values <- 1/eigen.H$values
      eigen.H$vectors <- Matrix(eigen.H$vectors, dimnames = list(id, id))
      H[[i]] <- tcrossprod(eigen.H$vectors %*% Diagonal(n = nrow(Hinv[[i]]),  x = sqrt(eigen.H$values)))
    }
    class(eigen.H) <- "list"
    H.list[[i]] <- list(H = H[[i]], eigen.H = eigen.H)
  }

  # Return H.list
  return(H.list)
}

